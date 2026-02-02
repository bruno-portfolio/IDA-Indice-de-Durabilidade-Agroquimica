import json
import logging
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from .config import SDH_TARGET, RESULTS_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class IDAMolecularWeights:
    srs: float = 0.30
    sce: float = 0.35
    robustness: float = 0.35

    def __post_init__(self) -> None:
        total = self.srs + self.sce + self.robustness
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Pesos devem somar 1.0, soma atual: {total}")

@dataclass
class ResidueIDAMolecular:
    res_seq: int
    res_name: str
    chain_id: str

    srs: float
    sce: float
    robustness: float
    fitness_cost_proxy: float
    sce_weight: float
    plddt: float
    functional_role: str

    ida_base: float
    ida_confidence_adjusted: float
    ida_molecular: float

    is_binding_site: bool = True
    frac_documented: bool = False

@dataclass
class IDAMolecularResult:
    layer: str = "IDA_MOLECULAR"
    description: str = "Barreira evolutiva intrínseca do alvo"
    target: str = ""
    organism: str = ""
    defensivo_class: str = ""
    date_generated: str = ""
    version: str = "2.0"

    residue_scores: list[ResidueIDAMolecular] = field(default_factory=list)

    ida_molecular_aggregate: float = 0.0
    ida_molecular_mean: float = 0.0
    ida_molecular_std: float = 0.0
    ida_molecular_min: float = 0.0
    ida_molecular_max: float = 0.0

    avg_srs: float = 0.0
    avg_sce: float = 0.0
    avg_robustness: float = 0.0
    avg_fitness_cost: float = 0.0
    avg_plddt: float = 0.0
    min_plddt: float = 0.0
    n_residues: int = 0

    barrier_level: str = ""
    risk_level: str = ""
    interpretation: str = ""
    recommendation: str = ""
    confidence_level: str = ""

    weights: IDAMolecularWeights = field(default_factory=IDAMolecularWeights)

    disclaimer: str = (
        "IDA Molecular mede a barreira INTRÍNSECA ao escape por mutação. "
        "A durabilidade de campo depende também de fatores regionais "
        "(IDA Regional) e de manejo (IDA Programa). Este score NÃO é "
        "uma previsão de vida útil."
    )

    def to_dict(self) -> dict:
        return {
            "metadata": {
                "layer": self.layer,
                "description": self.description,
                "target": self.target,
                "organism": self.organism,
                "defensivo_class": self.defensivo_class,
                "date_generated": self.date_generated,
                "version": self.version,
            },
            "aggregate_score": {
                "ida_molecular": round(self.ida_molecular_aggregate, 4),
                "interpretation": {
                    "barrier_level": self.barrier_level,
                    "risk_level": self.risk_level,
                    "interpretation": self.interpretation,
                    "recommendation": self.recommendation,
                },
            },
            "component_summary": {
                "avg_srs": round(self.avg_srs, 4),
                "avg_sce": round(self.avg_sce, 4),
                "avg_robustness": round(self.avg_robustness, 4),
                "avg_fitness_cost": round(self.avg_fitness_cost, 4),
            },
            "confidence": {
                "avg_plddt": round(self.avg_plddt, 1),
                "min_plddt": round(self.min_plddt, 1),
                "n_residues": self.n_residues,
                "confidence_level": self.confidence_level,
            },
            "weights": asdict(self.weights),
            "residue_details": [
                {
                    "res_seq": r.res_seq,
                    "res_name": r.res_name,
                    "chain_id": r.chain_id,
                    "srs": round(r.srs, 4),
                    "sce": round(r.sce, 4),
                    "robustness": round(r.robustness, 4),
                    "fitness_cost_proxy": round(r.fitness_cost_proxy, 4),
                    "ida_molecular": round(r.ida_molecular, 4),
                    "plddt": round(r.plddt, 1),
                    "functional_role": r.functional_role,
                    "frac_documented": r.frac_documented,
                }
                for r in self.residue_scores
            ],
            "disclaimer": self.disclaimer,
        }

    def get_per_residue_scores(self) -> dict[int, float]:
        return {r.res_seq: r.ida_molecular for r in self.residue_scores}

def _normalize_scores(scores: dict[int, float]) -> dict[int, float]:
    if not scores:
        return {}

    values = list(scores.values())
    min_val = min(values)
    max_val = max(values)

    if max_val == min_val:
        return {k: 0.5 for k in scores}

    return {
        k: (v - min_val) / (max_val - min_val)
        for k, v in scores.items()
    }

def _interpret_ida_molecular(ida_score: float) -> dict:
    if ida_score >= 0.8:
        return {
            "barrier_level": "high_barrier",
            "risk_level": "low",
            "interpretation": "Barreira evolutiva alta - sítio altamente restrito e conservado",
            "recommendation": "Alvo intrinsecamente robusto contra escape por mutação"
        }
    elif ida_score >= 0.6:
        return {
            "barrier_level": "moderate_high_barrier",
            "risk_level": "low_moderate",
            "interpretation": "Barreira moderada-alta - bom perfil estrutural e evolutivo",
            "recommendation": "Alvo com boa robustez, mas monitoramento é prudente"
        }
    elif ida_score >= 0.4:
        return {
            "barrier_level": "moderate_barrier",
            "risk_level": "moderate",
            "interpretation": "Barreira moderada - algumas vias de escape possíveis",
            "recommendation": "Atenção ao manejo de resistência desde o início"
        }
    elif ida_score >= 0.2:
        return {
            "barrier_level": "low_moderate_barrier",
            "risk_level": "moderate_high",
            "interpretation": "Barreira baixa-moderada - múltiplas vias de escape",
            "recommendation": "Priorizar uso em mistura e rotação com outros MoAs"
        }
    else:
        return {
            "barrier_level": "low_barrier",
            "risk_level": "high",
            "interpretation": "Barreira baixa - escape por mutação provável",
            "recommendation": "Evitar uso isolado, sempre em programa de manejo robusto"
        }

def _get_confidence_level(avg_plddt: float, min_plddt: float) -> str:
    if avg_plddt >= 90 and min_plddt >= 70:
        return "muito alta"
    elif avg_plddt >= 70:
        return "alta"
    elif avg_plddt >= 50:
        return "média"
    else:
        return "baixa"

def calculate_ida_molecular(
    binding_site_residues: list[int],
    residue_names: dict[int, str],
    chain_ids: dict[int, str],
    srs_scores: dict[int, float],
    sce_scores: dict[int, float],
    robustness_scores: dict[int, float],
    fitness_cost_scores: dict[int, float],
    sce_weights: dict[int, float],
    plddt_scores: dict[int, float],
    functional_roles: Optional[dict[int, str]] = None,
    frac_positions: Optional[set[int]] = None,
    weights: Optional[IDAMolecularWeights] = None,
    normalize: bool = True,
    apply_fitness_discount: bool = True,
) -> IDAMolecularResult:
    if weights is None:
        weights = IDAMolecularWeights()

    if functional_roles is None:
        functional_roles = {r: "unknown" for r in binding_site_residues}

    if frac_positions is None:
        frac_positions = set()

    if normalize:
        srs_scores = _normalize_scores(srs_scores)
        sce_scores = _normalize_scores(sce_scores)
        robustness_scores = _normalize_scores(robustness_scores)

    residue_scores = []

    for res_seq in binding_site_residues:
        srs = srs_scores.get(res_seq, 0.5)
        sce = sce_scores.get(res_seq, 0.5)
        robustness = robustness_scores.get(res_seq, 0.5)
        fitness_cost = fitness_cost_scores.get(res_seq, 0.5)
        sce_weight = sce_weights.get(res_seq, 0.8)
        plddt = plddt_scores.get(res_seq, 70.0)
        role = functional_roles.get(res_seq, "unknown")

        ida_base = (
            weights.srs * srs +
            weights.sce * sce +
            weights.robustness * robustness
        )

        ida_confidence_adjusted = ida_base * sce_weight

        if apply_fitness_discount:
            fitness_discount = 1 - (1 - fitness_cost) * 0.3
            ida_molecular = ida_confidence_adjusted * fitness_discount
        else:
            ida_molecular = ida_confidence_adjusted

        residue_scores.append(ResidueIDAMolecular(
            res_seq=res_seq,
            res_name=residue_names.get(res_seq, "UNK"),
            chain_id=chain_ids.get(res_seq, "A"),
            srs=srs,
            sce=sce,
            robustness=robustness,
            fitness_cost_proxy=fitness_cost,
            sce_weight=sce_weight,
            plddt=plddt,
            functional_role=role,
            ida_base=ida_base,
            ida_confidence_adjusted=ida_confidence_adjusted,
            ida_molecular=ida_molecular,
            frac_documented=res_seq in frac_positions,
        ))

    if not residue_scores:
        logger.error("Nenhum resíduo processado")
        return IDAMolecularResult(
            target=SDH_TARGET.enzyme,
            organism=SDH_TARGET.organism,
            defensivo_class=SDH_TARGET.defensivo_class,
            date_generated=datetime.now().isoformat(),
            residue_scores=[],
            interpretation="Erro: nenhum resíduo processado",
            confidence_level="none",
            weights=weights,
        )

    ida_values = [r.ida_molecular for r in residue_scores]
    sce_weight_values = [r.sce_weight for r in residue_scores]

    if sum(sce_weight_values) > 0:
        ida_aggregate = sum(
            r.ida_molecular * r.sce_weight for r in residue_scores
        ) / sum(sce_weight_values)
    else:
        ida_aggregate = np.mean(ida_values)

    ida_mean = np.mean(ida_values)
    ida_std = np.std(ida_values)

    avg_plddt = np.mean([r.plddt for r in residue_scores])
    min_plddt = min(r.plddt for r in residue_scores)
    avg_srs = np.mean([r.srs for r in residue_scores])
    avg_sce = np.mean([r.sce for r in residue_scores])
    avg_robustness = np.mean([r.robustness for r in residue_scores])
    avg_fitness_cost = np.mean([r.fitness_cost_proxy for r in residue_scores])

    interp = _interpret_ida_molecular(ida_aggregate)

    logger.info(
        f"IDA Molecular calculado: {len(residue_scores)} resíduos, "
        f"IDA_mol={ida_aggregate:.3f}, risco={interp['risk_level']}"
    )

    return IDAMolecularResult(
        target=SDH_TARGET.enzyme,
        organism=SDH_TARGET.organism,
        defensivo_class=SDH_TARGET.defensivo_class,
        date_generated=datetime.now().isoformat(),
        residue_scores=residue_scores,
        ida_molecular_aggregate=ida_aggregate,
        ida_molecular_mean=ida_mean,
        ida_molecular_std=ida_std,
        ida_molecular_min=min(ida_values),
        ida_molecular_max=max(ida_values),
        avg_srs=avg_srs,
        avg_sce=avg_sce,
        avg_robustness=avg_robustness,
        avg_fitness_cost=avg_fitness_cost,
        avg_plddt=avg_plddt,
        min_plddt=min_plddt,
        n_residues=len(residue_scores),
        barrier_level=interp["barrier_level"],
        risk_level=interp["risk_level"],
        interpretation=interp["interpretation"],
        recommendation=interp["recommendation"],
        confidence_level=_get_confidence_level(avg_plddt, min_plddt),
        weights=weights,
    )

def export_ida_molecular_json(
    result: IDAMolecularResult,
    output_path: Optional[Path] = None,
) -> Path:
    if output_path is None:
        output_path = RESULTS_DIR / "ida_molecular_results.json"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(result.to_dict(), f, indent=2, ensure_ascii=False)

    logger.info(f"IDA Molecular exportado: {output_path}")
    return output_path

def export_ida_molecular_csv(
    result: IDAMolecularResult,
    output_path: Optional[Path] = None,
) -> Path:
    if output_path is None:
        output_path = RESULTS_DIR / "ida_molecular_residues.csv"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    header = (
        "res_seq,res_name,chain_id,srs,sce,robustness,fitness_cost_proxy,"
        "sce_weight,plddt,functional_role,ida_base,ida_molecular,frac_documented"
    )
    lines = [header]

    for r in result.residue_scores:
        lines.append(
            f"{r.res_seq},{r.res_name},{r.chain_id},"
            f"{r.srs:.4f},{r.sce:.4f},{r.robustness:.4f},{r.fitness_cost_proxy:.4f},"
            f"{r.sce_weight:.4f},{r.plddt:.1f},{r.functional_role},"
            f"{r.ida_base:.4f},{r.ida_molecular:.4f},{r.frac_documented}"
        )

    output_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"CSV exportado: {output_path}")
    return output_path

def print_ida_molecular_summary(result: IDAMolecularResult) -> None:
    print("\n" + "=" * 70)
    print("IDA MOLECULAR - BARREIRA EVOLUTIVA INTRÍNSECA")
    print("=" * 70)
    print(f"Alvo: {result.target}")
    print(f"Organismo: {result.organism}")
    print(f"Classe: {result.defensivo_class}")
    print(f"Data: {result.date_generated}")
    print("-" * 70)
    print(f"IDA Molecular: {result.ida_molecular_aggregate:.4f}")
    print(f"  Média: {result.ida_molecular_mean:.4f} (±{result.ida_molecular_std:.4f})")
    print(f"  Range: [{result.ida_molecular_min:.4f}, {result.ida_molecular_max:.4f}]")
    print("-" * 70)
    print("Componentes médios:")
    print(f"  SRS (Restrição):     {result.avg_srs:.4f}")
    print(f"  SCE (Conservação):   {result.avg_sce:.4f}")
    print(f"  Robustez:            {result.avg_robustness:.4f}")
    print(f"  Fitness Cost Proxy:  {result.avg_fitness_cost:.4f}")
    print("-" * 70)
    print("Confiança:")
    print(f"  pLDDT médio: {result.avg_plddt:.1f}")
    print(f"  pLDDT mínimo: {result.min_plddt:.1f}")
    print(f"  Nível: {result.confidence_level}")
    print("-" * 70)
    print("Interpretação:")
    print(f"  Barreira: {result.barrier_level}")
    print(f"  Risco: {result.risk_level}")
    print(f"  {result.interpretation}")
    print(f"  Recomendação: {result.recommendation}")
    print("-" * 70)
    print(f"Resíduos analisados: {result.n_residues}")
    print("-" * 70)
    print("NOTA: " + result.disclaimer)
    print("=" * 70)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Cálculo do IDA Molecular (Barreira Evolutiva Intrínseca)"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        binding_site = [145, 146, 147, 230, 231, 272, 273]
        residue_names = {
            145: "TYR", 146: "GLY", 147: "ALA",
            230: "HIS", 231: "SER", 272: "PRO", 273: "ILE"
        }
        chain_ids = {r: "B" for r in binding_site}

        srs_scores = {
            145: 0.85, 146: 0.72, 147: 0.68,
            230: 0.91, 231: 0.65, 272: 0.78, 273: 0.82
        }
        sce_scores = {
            145: 0.92, 146: 0.88, 147: 0.75,
            230: 0.95, 231: 0.70, 272: 0.85, 273: 0.80
        }
        robustness_scores = {
            145: 0.75, 146: 0.65, 147: 0.55,
            230: 0.80, 231: 0.50, 272: 0.70, 273: 0.60
        }
        fitness_cost_scores = {
            145: 0.78, 146: 0.65, 147: 0.48,
            230: 0.85, 231: 0.55, 272: 0.72, 273: 0.62
        }
        sce_weights = {r: 0.9 for r in binding_site}
        plddt_scores = {
            145: 92.5, 146: 88.3, 147: 85.1,
            230: 94.2, 231: 82.7, 272: 90.1, 273: 87.5
        }
        functional_roles = {
            145: "direct_contact", 146: "structural", 147: "peripheral",
            230: "catalytic", 231: "peripheral", 272: "direct_contact", 273: "structural"
        }
        frac_positions = {272}

        result = calculate_ida_molecular(
            binding_site_residues=binding_site,
            residue_names=residue_names,
            chain_ids=chain_ids,
            srs_scores=srs_scores,
            sce_scores=sce_scores,
            robustness_scores=robustness_scores,
            fitness_cost_scores=fitness_cost_scores,
            sce_weights=sce_weights,
            plddt_scores=plddt_scores,
            functional_roles=functional_roles,
            frac_positions=frac_positions,
        )

        print_ida_molecular_summary(result)

        json_path = export_ida_molecular_json(result)
        csv_path = export_ida_molecular_csv(result)
        print(f"\nExportado: {json_path}")
        print(f"Exportado: {csv_path}")