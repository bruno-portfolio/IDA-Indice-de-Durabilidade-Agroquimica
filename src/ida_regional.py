import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .config import RESULTS_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class RegionalMutation:
    mutation_id: str
    position: int
    original_aa: str
    mutant_aa: str
    region: str
    first_reported: str
    frequency: str
    cross_resistance_group: str = ""

@dataclass
class MutationPenalty:
    penalty_score: float
    n_mutations: int
    mutation_details: list[dict] = field(default_factory=list)
    interpretation: str = ""

@dataclass
class SelectionPressure:
    region_id: str
    region_name: str
    season: str
    pressure_score: float
    components: dict = field(default_factory=dict)

@dataclass
class TargetSaturation:
    saturation_score: float
    raw_proportion: float = 0.0
    data_quality: str = "no_data"
    interpretation: str = ""

@dataclass
class IDARegionalResult:
    layer: str = "IDA_REGIONAL"
    description: str = "Ajuste do IDA por contexto regional de seleção"
    region_id: str = ""
    region_name: str = ""
    season: str = ""
    date_generated: str = ""
    version: str = "2.0"

    ida_molecular: float = 0.0
    ida_regional: float = 0.0

    mutation_penalty: MutationPenalty = field(default_factory=lambda: MutationPenalty(0, 0))
    selection_pressure: SelectionPressure = field(default_factory=lambda: SelectionPressure("", "", "", 0))
    target_saturation: TargetSaturation = field(default_factory=lambda: TargetSaturation(0))

    adjustment_factors: dict = field(default_factory=dict)

    interpretation: str = ""
    confidence_level: str = ""

    disclaimer: str = (
        "IDA Regional ajusta o IDA Molecular pelo contexto de seleção "
        "específico de cada região. Regiões com mais mutações documentadas, "
        "maior pressão de uso e maior saturação do MoA terão IDA Regional "
        "mais baixo, indicando maior risco de erosão de eficácia."
    )

    def to_dict(self) -> dict:
        return {
            "metadata": {
                "layer": self.layer,
                "description": self.description,
                "region_id": self.region_id,
                "region_name": self.region_name,
                "season": self.season,
                "date_generated": self.date_generated,
                "version": self.version,
            },
            "scores": {
                "ida_molecular": round(self.ida_molecular, 4),
                "ida_regional": round(self.ida_regional, 4),
                "degradation_pct": round(
                    (self.ida_molecular - self.ida_regional) / self.ida_molecular * 100
                    if self.ida_molecular > 0 else 0,
                    1
                ),
            },
            "adjustment_factors": {
                "mutation_penalty": self.mutation_penalty.penalty_score,
                "pressure_score": self.selection_pressure.pressure_score,
                "saturation_score": self.target_saturation.saturation_score,
                "combined_factor": self.adjustment_factors.get("combined_factor", 0),
            },
            "mutation_details": {
                "n_mutations": self.mutation_penalty.n_mutations,
                "details": self.mutation_penalty.mutation_details,
            },
            "interpretation": self.interpretation,
            "confidence_level": self.confidence_level,
            "disclaimer": self.disclaimer,
        }

def load_regional_mutations(mutations_csv: str) -> dict[str, list[RegionalMutation]]:
    df = pd.read_csv(mutations_csv)

    mutations_by_region = {}

    for _, row in df.iterrows():
        mutation = RegionalMutation(
            mutation_id=row.get("mutation_id", f"mut_{row.name}"),
            position=int(row["position"]),
            original_aa=row["original_aa"],
            mutant_aa=row["mutant_aa"],
            region=row["region"],
            first_reported=str(row.get("first_reported", "")),
            frequency=row.get("frequency", "medium"),
            cross_resistance_group=row.get("cross_resistance_group", ""),
        )

        if mutation.region not in mutations_by_region:
            mutations_by_region[mutation.region] = []
        mutations_by_region[mutation.region].append(mutation)

    logger.info(
        f"Mutações carregadas: {len(df)} mutações em {len(mutations_by_region)} regiões"
    )

    return mutations_by_region

def calculate_mutation_penalty(
    regional_mutations: list[RegionalMutation],
    ida_molecular_per_residue: dict[int, float],
) -> MutationPenalty:
    if not regional_mutations:
        return MutationPenalty(
            penalty_score=0.0,
            n_mutations=0,
            interpretation="Sem mutações documentadas nesta região"
        )

    penalties = []

    for mut in regional_mutations:
        freq_penalty = {"low": 0.1, "medium": 0.3, "high": 0.5}.get(mut.frequency, 0.2)

        ida_at_position = ida_molecular_per_residue.get(mut.position, 0.5)
        barrier_breach_penalty = ida_at_position * 0.3

        total_penalty = freq_penalty + barrier_breach_penalty

        penalties.append({
            "mutation": f"{mut.original_aa}{mut.position}{mut.mutant_aa}",
            "frequency": mut.frequency,
            "first_reported": mut.first_reported,
            "penalty": total_penalty
        })

    aggregate_penalty = min(1.0, sum(p["penalty"] for p in penalties) ** 0.7)

    interpretation = _interpret_mutation_penalty(aggregate_penalty)

    return MutationPenalty(
        penalty_score=aggregate_penalty,
        n_mutations=len(penalties),
        mutation_details=penalties,
        interpretation=interpretation,
    )

def _interpret_mutation_penalty(penalty: float) -> str:
    if penalty < 0.1:
        return "Risco baixo - sem mutações significativas documentadas"
    elif penalty < 0.3:
        return "Risco moderado - algumas mutações de baixa frequência"
    elif penalty < 0.5:
        return "Risco elevado - mutações estabelecidas na região"
    else:
        return "Risco alto - múltiplas mutações de alta frequência"

def calculate_selection_pressure(
    region_data: dict,
    weights: Optional[dict] = None,
) -> SelectionPressure:
    if weights is None:
        weights = {
            "area_weight": 0.25,
            "years_use_weight": 0.35,
            "severity_weight": 0.25,
            "intensity_weight": 0.15
        }

    area_norm = min(1.0, region_data.get("soy_area_ha", 0) / 15_000_000)

    years_norm = min(1.0, region_data.get("sdhi_years_use", 0) / 12)

    severity_map = {"low": 0.3, "medium": 0.6, "high": 0.9}
    severity_norm = severity_map.get(region_data.get("disease_severity", "medium"), 0.5)

    intensity_norm = region_data.get("application_intensity", 0.5)

    pressure_score = (
        weights["area_weight"] * area_norm +
        weights["years_use_weight"] * years_norm +
        weights["severity_weight"] * severity_norm +
        weights["intensity_weight"] * intensity_norm
    )

    return SelectionPressure(
        region_id=region_data.get("region_id", ""),
        region_name=region_data.get("region_name", ""),
        season=region_data.get("season", ""),
        pressure_score=pressure_score,
        components={
            "area_contribution": area_norm * weights["area_weight"],
            "years_contribution": years_norm * weights["years_use_weight"],
            "severity_contribution": severity_norm * weights["severity_weight"],
            "intensity_contribution": intensity_norm * weights["intensity_weight"]
        }
    )

def calculate_target_saturation(
    moa_usage_data: dict,
    target_moa: str = "SDHI",
) -> TargetSaturation:
    if not moa_usage_data:
        return TargetSaturation(
            saturation_score=0.5,
            data_quality="no_data",
            interpretation="Dados de market share não disponíveis"
        )

    total_applications = sum(moa_usage_data.values())
    target_applications = moa_usage_data.get(target_moa, 0)

    if total_applications == 0:
        saturation = 0.5
    else:
        saturation = target_applications / total_applications

    saturation_score = min(1.0, saturation / 0.5)

    interpretation = _interpret_saturation(saturation_score)

    return TargetSaturation(
        saturation_score=saturation_score,
        raw_proportion=saturation,
        data_quality="available",
        interpretation=interpretation,
    )

def _interpret_saturation(saturation: float) -> str:
    if saturation < 0.3:
        return "Baixa saturação - alvo não sobrecarregado"
    elif saturation < 0.6:
        return "Saturação moderada - atenção à diversificação"
    else:
        return "Alta saturação - risco sistêmico por concentração de MoA"

def calculate_ida_regional(
    ida_molecular: float,
    mutation_penalty: MutationPenalty,
    selection_pressure: SelectionPressure,
    target_saturation: TargetSaturation,
    factor_weights: Optional[dict] = None,
) -> IDARegionalResult:
    if factor_weights is None:
        factor_weights = {
            "mutation": 0.40,
            "pressure": 0.35,
            "saturation": 0.25
        }

    mut_penalty = mutation_penalty.penalty_score
    pressure_score = selection_pressure.pressure_score
    sat_score = target_saturation.saturation_score

    combined_factor = (
        factor_weights["mutation"] * mut_penalty +
        factor_weights["pressure"] * pressure_score * 0.5 +
        factor_weights["saturation"] * sat_score * 0.3
    )

    ida_regional = ida_molecular * (1 - combined_factor)

    interpretation = _interpret_ida_regional(ida_molecular, ida_regional)

    confidence = _assess_regional_confidence(mutation_penalty, target_saturation)

    logger.info(
        f"IDA Regional calculado: {selection_pressure.region_name}, "
        f"IDA_mol={ida_molecular:.3f} -> IDA_reg={ida_regional:.3f}"
    )

    return IDARegionalResult(
        region_id=selection_pressure.region_id,
        region_name=selection_pressure.region_name,
        season=selection_pressure.season,
        date_generated=datetime.now().isoformat(),
        ida_molecular=ida_molecular,
        ida_regional=ida_regional,
        mutation_penalty=mutation_penalty,
        selection_pressure=selection_pressure,
        target_saturation=target_saturation,
        adjustment_factors={
            "mutation_penalty": mut_penalty,
            "pressure_score": pressure_score,
            "saturation_score": sat_score,
            "combined_factor": combined_factor
        },
        interpretation=interpretation,
        confidence_level=confidence,
    )

def _interpret_ida_regional(ida_mol: float, ida_reg: float) -> str:
    degradation = (ida_mol - ida_reg) / ida_mol * 100 if ida_mol > 0 else 0

    if degradation < 10:
        return f"IDA Regional próximo ao Molecular - região com baixa pressão acumulada"
    elif degradation < 25:
        return f"IDA Regional moderadamente reduzido ({degradation:.0f}%) - monitoramento recomendado"
    elif degradation < 40:
        return f"IDA Regional significativamente reduzido ({degradation:.0f}%) - atenção redobrada ao manejo"
    else:
        return f"IDA Regional muito reduzido ({degradation:.0f}%) - região de alto risco, priorizar rotação de MoA"

def _assess_regional_confidence(
    mutation_penalty: MutationPenalty,
    target_saturation: TargetSaturation,
) -> str:
    if mutation_penalty.n_mutations >= 3 and target_saturation.data_quality == "available":
        return "high"
    elif mutation_penalty.n_mutations >= 1 or target_saturation.data_quality == "available":
        return "medium"
    else:
        return "low"

def compare_regions(regional_results: list[IDARegionalResult]) -> dict:
    if not regional_results:
        return {"error": "Sem resultados para comparar"}

    sorted_results = sorted(regional_results, key=lambda x: x.ida_regional)

    ranking = []
    for i, result in enumerate(sorted_results, 1):
        degradation = (
            (result.ida_molecular - result.ida_regional) / result.ida_molecular * 100
            if result.ida_molecular > 0 else 0
        )
        ranking.append({
            "rank": i,
            "region": result.region_name,
            "ida_molecular": round(result.ida_molecular, 4),
            "ida_regional": round(result.ida_regional, 4),
            "degradation_pct": round(degradation, 1),
            "risk_level": (
                "high" if result.ida_regional < 0.4
                else "medium" if result.ida_regional < 0.6
                else "low"
            )
        })

    return {
        "ranking": ranking,
        "highest_risk_region": sorted_results[0].region_name,
        "lowest_risk_region": sorted_results[-1].region_name,
        "avg_ida_regional": round(
            sum(r.ida_regional for r in regional_results) / len(regional_results),
            4
        ),
        "n_regions": len(regional_results),
    }

def export_ida_regional_results(
    regional_results: list[IDARegionalResult],
    output_path: Optional[Path] = None,
) -> dict:
    if output_path is None:
        output_path = RESULTS_DIR / "ida_regional_results.json"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results = {
        "metadata": {
            "layer": "IDA_REGIONAL",
            "description": "Ajuste do IDA por contexto regional de seleção",
            "date_generated": datetime.now().isoformat(),
            "version": "2.0"
        },
        "summary": compare_regions(regional_results),
        "regions": [r.to_dict() for r in regional_results],
        "disclaimer": (
            "IDA Regional ajusta o IDA Molecular pelo contexto de seleção "
            "específico de cada região. Regiões com mais mutações documentadas, "
            "maior pressão de uso e maior saturação do MoA terão IDA Regional "
            "mais baixo, indicando maior risco de erosão de eficácia."
        )
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    logger.info(f"IDA Regional exportado: {output_path}")
    return results

def print_ida_regional_summary(result: IDARegionalResult) -> None:
    degradation = (
        (result.ida_molecular - result.ida_regional) / result.ida_molecular * 100
        if result.ida_molecular > 0 else 0
    )

    print("\n" + "=" * 70)
    print("IDA REGIONAL - CONTEXTO DE SELEÇÃO")
    print("=" * 70)
    print(f"Região: {result.region_name} ({result.region_id})")
    print(f"Safra: {result.season}")
    print(f"Data: {result.date_generated}")
    print("-" * 70)
    print("SCORES:")
    print(f"  IDA Molecular: {result.ida_molecular:.4f}")
    print(f"  IDA Regional:  {result.ida_regional:.4f}")
    print(f"  Degradação:    {degradation:.1f}%")
    print("-" * 70)
    print("FATORES DE AJUSTE:")
    print(f"  Penalidade mutações:  {result.mutation_penalty.penalty_score:.3f} ({result.mutation_penalty.n_mutations} mutações)")
    print(f"  Pressão de seleção:   {result.selection_pressure.pressure_score:.3f}")
    print(f"  Saturação de alvo:    {result.target_saturation.saturation_score:.3f}")
    print(f"  Fator combinado:      {result.adjustment_factors.get('combined_factor', 0):.3f}")
    print("-" * 70)
    print(f"Interpretação: {result.interpretation}")
    print(f"Confiança: {result.confidence_level}")
    print("-" * 70)
    print("NOTA: " + result.disclaimer[:100] + "...")
    print("=" * 70)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cálculo do IDA Regional")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        ida_molecular = 0.72

        ida_per_residue = {145: 0.85, 230: 0.90, 272: 0.65, 273: 0.70}

        mutations_mt = [
            RegionalMutation(
                mutation_id="FRAC_001",
                position=272,
                original_aa="P",
                mutant_aa="L",
                region="BR-MT",
                first_reported="2019",
                frequency="high",
            ),
            RegionalMutation(
                mutation_id="FRAC_002",
                position=273,
                original_aa="I",
                mutant_aa="V",
                region="BR-MT",
                first_reported="2021",
                frequency="medium",
            ),
        ]

        mut_penalty = calculate_mutation_penalty(mutations_mt, ida_per_residue)

        region_data = {
            "region_id": "BR-MT",
            "region_name": "Mato Grosso",
            "season": "2024/25",
            "soy_area_ha": 12_000_000,
            "sdhi_years_use": 10,
            "disease_severity": "high",
            "application_intensity": 0.7,
        }
        pressure = calculate_selection_pressure(region_data)

        moa_usage = {"SDHI": 45, "QoI": 30, "DMI": 25}
        saturation = calculate_target_saturation(moa_usage)

        result = calculate_ida_regional(
            ida_molecular=ida_molecular,
            mutation_penalty=mut_penalty,
            selection_pressure=pressure,
            target_saturation=saturation,
        )

        print_ida_regional_summary(result)

        mutations_rs = [
            RegionalMutation(
                mutation_id="FRAC_003",
                position=272,
                original_aa="P",
                mutant_aa="L",
                region="BR-RS",
                first_reported="2022",
                frequency="low",
            ),
        ]

        region_data_rs = {
            "region_id": "BR-RS",
            "region_name": "Rio Grande do Sul",
            "season": "2024/25",
            "soy_area_ha": 6_000_000,
            "sdhi_years_use": 8,
            "disease_severity": "medium",
        }

        mut_penalty_rs = calculate_mutation_penalty(mutations_rs, ida_per_residue)
        pressure_rs = calculate_selection_pressure(region_data_rs)
        saturation_rs = calculate_target_saturation({"SDHI": 25, "QoI": 40, "DMI": 35})

        result_rs = calculate_ida_regional(
            ida_molecular=ida_molecular,
            mutation_penalty=mut_penalty_rs,
            selection_pressure=pressure_rs,
            target_saturation=saturation_rs,
        )

        print_ida_regional_summary(result_rs)

        print("\n" + "=" * 70)
        print("COMPARATIVO ENTRE REGIÕES")
        print("=" * 70)
        comparison = compare_regions([result, result_rs])
        print(f"Região de maior risco: {comparison['highest_risk_region']}")
        print(f"Região de menor risco: {comparison['lowest_risk_region']}")
        print(f"IDA Regional médio: {comparison['avg_ida_regional']}")
        print("\nRanking:")
        for r in comparison['ranking']:
            print(f"  {r['rank']}. {r['region']}: IDA_reg={r['ida_regional']:.3f} ({r['risk_level']})")