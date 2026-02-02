import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from .config import RESULTS_DIR, LOG_FORMAT, LOG_DATE_FORMAT
from .ida_regional import IDARegionalResult

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

MOA_CODES = {
    "SDHI": "7",
    "QoI": "11",
    "DMI": "3",
    "multisite": "M",
    "MBC": "1",
}

def _infer_moa_code(product_type: str) -> str:
    return MOA_CODES.get(product_type, "U")

@dataclass
class Application:
    timing: str
    product_type: str
    active_ingredient: str
    partner: Optional[str]
    moa_code: str = ""

    def __post_init__(self):
        if not self.moa_code:
            self.moa_code = _infer_moa_code(self.product_type)

@dataclass
class ManagementProgram:
    program_id: str
    name: str
    applications: list[Application] = field(default_factory=list)

    total_applications: int = 0
    sdhi_applications: int = 0
    sdhi_proportion: float = 0.0
    multisite_applications: int = 0
    multisite_proportion: float = 0.0
    moa_diversity: int = 0

    def __post_init__(self):
        self._calculate_metrics()

    def _calculate_metrics(self):
        self.total_applications = len(self.applications)

        if self.total_applications > 0:
            self.sdhi_applications = sum(
                1 for a in self.applications if a.product_type == "SDHI"
            )
            self.sdhi_proportion = self.sdhi_applications / self.total_applications

            self.multisite_applications = sum(
                1 for a in self.applications if a.product_type == "multisite"
            )
            self.multisite_proportion = self.multisite_applications / self.total_applications

            self.moa_diversity = len(set(a.moa_code for a in self.applications))

@dataclass
class MixtureProtection:
    protection_score: float
    components: dict = field(default_factory=dict)
    interpretation: str = ""

@dataclass
class RotationScore:
    rotation_score: float
    max_sdhi_consecutive: int = 0
    intercalation_bonus: float = 0.0
    interpretation: str = ""

@dataclass
class IDAProgramaResult:
    layer: str = "IDA_PROGRAMA"
    description: str = "IDA completo incluindo estratégia de manejo"
    program_id: str = ""
    program_name: str = ""
    region_id: str = ""
    region_name: str = ""
    season: str = ""
    date_generated: str = ""
    version: str = "2.0"

    ida_molecular: float = 0.0
    ida_regional: float = 0.0

    mixture_protection: float = 0.0
    rotation_score: float = 0.0
    program_bonus: float = 0.0

    ida_programa: float = 0.0

    mixture_details: MixtureProtection = field(default_factory=lambda: MixtureProtection(0))
    rotation_details: RotationScore = field(default_factory=lambda: RotationScore(0))

    interpretation: str = ""
    recommendations: list[str] = field(default_factory=list)

    disclaimer: str = (
        "IDA Programa representa a avaliação COMPLETA de risco de erosão "
        "de eficácia, integrando: (1) barreira molecular intrínseca, "
        "(2) contexto regional de seleção, e (3) estratégia de manejo. "
        "Este NÃO é uma previsão de vida útil, mas uma estimativa comparativa "
        "de risco que pode orientar decisões de manejo de resistência."
    )

    def to_dict(self) -> dict:
        return {
            "metadata": {
                "layer": self.layer,
                "description": self.description,
                "program_id": self.program_id,
                "program_name": self.program_name,
                "region_id": self.region_id,
                "region_name": self.region_name,
                "season": self.season,
                "date_generated": self.date_generated,
                "version": self.version,
            },
            "scores": {
                "ida_molecular": round(self.ida_molecular, 4),
                "ida_regional": round(self.ida_regional, 4),
                "ida_programa": round(self.ida_programa, 4),
            },
            "program_analysis": {
                "mixture_protection": round(self.mixture_protection, 4),
                "rotation_score": round(self.rotation_score, 4),
                "program_bonus": round(self.program_bonus, 4),
                "effect": "protetor" if self.program_bonus > 0 else "acelerador" if self.program_bonus < 0 else "neutro",
            },
            "interpretation": self.interpretation,
            "recommendations": self.recommendations,
            "disclaimer": self.disclaimer,
        }

def load_programs(programs_json: str) -> list[ManagementProgram]:
    with open(programs_json, encoding="utf-8") as f:
        data = json.load(f)

    programs = []
    for prog_data in data.get("programs", []):
        applications = [
            Application(
                timing=app["timing"],
                product_type=app["product_type"],
                active_ingredient=app.get("active", app.get("active_ingredient", "")),
                partner=app.get("partner"),
                moa_code=app.get("moa_code", _infer_moa_code(app["product_type"]))
            )
            for app in prog_data.get("applications", [])
        ]

        programs.append(ManagementProgram(
            program_id=prog_data.get("program_id", f"prog_{len(programs)}"),
            name=prog_data.get("name", "Programa sem nome"),
            applications=applications
        ))

    logger.info(f"Carregados {len(programs)} programas de manejo")
    return programs

def calculate_mixture_protection_score(program: ManagementProgram) -> MixtureProtection:
    if program.total_applications == 0:
        return MixtureProtection(
            protection_score=0.5,
            interpretation="Programa vazio"
        )

    multisite_score = min(0.3, program.multisite_proportion * 0.6)

    diversity_score = min(0.3, (program.moa_diversity - 1) * 0.1)

    if program.sdhi_proportion <= 0.25:
        dilution_score = 0.4
    elif program.sdhi_proportion <= 0.50:
        dilution_score = 0.3
    elif program.sdhi_proportion <= 0.75:
        dilution_score = 0.15
    else:
        dilution_score = 0.0

    solo_sdhi = sum(
        1 for a in program.applications
        if a.product_type == "SDHI" and not a.partner
    )
    solo_penalty = min(0.3, solo_sdhi * 0.15)

    protection_score = multisite_score + diversity_score + dilution_score - solo_penalty
    protection_score = max(0.0, min(1.0, protection_score))

    interpretation = _interpret_mixture_protection(protection_score)

    return MixtureProtection(
        protection_score=protection_score,
        components={
            "multisite_contribution": multisite_score,
            "diversity_contribution": diversity_score,
            "dilution_contribution": dilution_score,
            "solo_penalty": solo_penalty
        },
        interpretation=interpretation,
    )

def _interpret_mixture_protection(score: float) -> str:
    if score >= 0.7:
        return "Excelente proteção - programa bem balanceado com multissítios e diversidade"
    elif score >= 0.5:
        return "Boa proteção - uso adequado de misturas"
    elif score >= 0.3:
        return "Proteção moderada - considerar aumentar diversidade de MoA"
    else:
        return "Baixa proteção - programa concentrado em poucos MoAs, alto risco de seleção"

def calculate_rotation_score(program: ManagementProgram) -> RotationScore:
    if program.total_applications < 2:
        return RotationScore(
            rotation_score=0.5,
            interpretation="Programa muito curto para avaliar rotação"
        )

    moa_sequence = [a.product_type for a in program.applications]

    sdhi_consecutive = 0
    max_sdhi_consecutive = 0
    for moa in moa_sequence:
        if moa == "SDHI":
            sdhi_consecutive += 1
            max_sdhi_consecutive = max(max_sdhi_consecutive, sdhi_consecutive)
        else:
            sdhi_consecutive = 0

    if max_sdhi_consecutive >= 3:
        rotation_penalty = 0.5
    elif max_sdhi_consecutive == 2:
        rotation_penalty = 0.25
    else:
        rotation_penalty = 0.0

    intercalation_bonus = 0.0
    for i in range(len(moa_sequence) - 1):
        if moa_sequence[i] == "SDHI" and moa_sequence[i + 1] == "multisite":
            intercalation_bonus += 0.1
    intercalation_bonus = min(0.3, intercalation_bonus)

    rotation_score = 0.7 - rotation_penalty + intercalation_bonus
    rotation_score = max(0.0, min(1.0, rotation_score))

    interpretation = _interpret_rotation(rotation_score, max_sdhi_consecutive)

    return RotationScore(
        rotation_score=rotation_score,
        max_sdhi_consecutive=max_sdhi_consecutive,
        intercalation_bonus=intercalation_bonus,
        interpretation=interpretation,
    )

def _interpret_rotation(score: float, max_consecutive: int) -> str:
    if max_consecutive >= 3:
        return f"Rotação inadequada - {max_consecutive} aplicações de SDHI consecutivas"
    elif max_consecutive == 2:
        return "Rotação subótima - evitar SDHI consecutivo"
    elif score >= 0.7:
        return "Boa rotação - MoAs bem alternados"
    else:
        return "Rotação aceitável"

def calculate_ida_programa(
    ida_regional_result: IDARegionalResult,
    program: ManagementProgram,
) -> IDAProgramaResult:
    mixture = calculate_mixture_protection_score(program)
    rotation = calculate_rotation_score(program)

    avg_protection = (mixture.protection_score + rotation.rotation_score) / 2

    program_bonus = (avg_protection - 0.5) * 0.4

    ida_programa = ida_regional_result.ida_regional * (1 + program_bonus)
    ida_programa = max(0.0, min(1.0, ida_programa))

    recommendations = _generate_recommendations(program, mixture, rotation)

    interpretation = _interpret_ida_programa(ida_programa, program_bonus)

    logger.info(
        f"IDA Programa calculado: {program.name}, "
        f"IDA_reg={ida_regional_result.ida_regional:.3f} -> IDA_prog={ida_programa:.3f} "
        f"(bonus={program_bonus:+.3f})"
    )

    return IDAProgramaResult(
        program_id=program.program_id,
        program_name=program.name,
        region_id=ida_regional_result.region_id,
        region_name=ida_regional_result.region_name,
        season=ida_regional_result.season,
        date_generated=datetime.now().isoformat(),
        ida_molecular=ida_regional_result.ida_molecular,
        ida_regional=ida_regional_result.ida_regional,
        mixture_protection=mixture.protection_score,
        rotation_score=rotation.rotation_score,
        program_bonus=program_bonus,
        ida_programa=ida_programa,
        mixture_details=mixture,
        rotation_details=rotation,
        interpretation=interpretation,
        recommendations=recommendations,
    )

def _interpret_ida_programa(ida: float, bonus: float) -> str:
    effect = "protege" if bonus > 0 else "acelera" if bonus < 0 else "neutro"

    if ida >= 0.7:
        return f"Risco baixo - programa {effect} a durabilidade"
    elif ida >= 0.5:
        return f"Risco moderado - programa {effect}, monitoramento recomendado"
    elif ida >= 0.3:
        return f"Risco elevado - programa {effect}, ajustes necessários"
    else:
        return f"Risco alto - programa {effect} significativamente, revisar estratégia"

def _generate_recommendations(
    program: ManagementProgram,
    mixture: MixtureProtection,
    rotation: RotationScore,
) -> list[str]:
    recs = []

    if mixture.protection_score < 0.5:
        if program.multisite_proportion < 0.2:
            recs.append("Incluir pelo menos 1 aplicação de multissítio (mancozeb, clorotalonil)")
        if program.moa_diversity < 3:
            recs.append("Diversificar mecanismos de ação no programa")

    if rotation.max_sdhi_consecutive >= 2:
        recs.append("Evitar aplicações consecutivas de SDHI - intercalar com outro MoA")

    if program.sdhi_proportion > 0.5:
        recs.append("Reduzir proporção de SDHI para menos de 50% do programa")

    solo_sdhi = sum(
        1 for a in program.applications
        if a.product_type == "SDHI" and not a.partner
    )
    if solo_sdhi > 0:
        recs.append(f"{solo_sdhi} aplicação(ões) de SDHI sem parceiro - usar sempre em mistura")

    if not recs:
        recs.append("Programa bem estruturado - manter monitoramento de eficácia")

    return recs

def compare_programs(program_results: list[IDAProgramaResult]) -> dict:
    if not program_results:
        return {"error": "Sem programas para comparar"}

    sorted_results = sorted(
        program_results,
        key=lambda x: x.ida_programa,
        reverse=True
    )

    comparison = {
        "ranking": [
            {
                "rank": i + 1,
                "program": r.program_name,
                "ida_programa": round(r.ida_programa, 4),
                "program_bonus": round(r.program_bonus, 4),
                "effect": "protetor" if r.program_bonus > 0 else "acelerador" if r.program_bonus < 0 else "neutro",
                "key_recommendations": r.recommendations[:2]
            }
            for i, r in enumerate(sorted_results)
        ],
        "best_program": sorted_results[0].program_name,
        "worst_program": sorted_results[-1].program_name,
        "ida_range": {
            "max": round(sorted_results[0].ida_programa, 4),
            "min": round(sorted_results[-1].ida_programa, 4),
            "spread": round(sorted_results[0].ida_programa - sorted_results[-1].ida_programa, 4)
        },
        "n_programs": len(program_results),
    }

    return comparison

def export_ida_programa_results(
    program_results: list[IDAProgramaResult],
    output_path: Optional[Path] = None,
) -> dict:
    if output_path is None:
        output_path = RESULTS_DIR / "ida_programa_results.json"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results = {
        "metadata": {
            "layer": "IDA_PROGRAMA",
            "description": "IDA completo incluindo estratégia de manejo",
            "date_generated": datetime.now().isoformat(),
            "version": "2.0"
        },
        "comparison": compare_programs(program_results),
        "programs": [r.to_dict() for r in program_results],
        "disclaimer": (
            "IDA Programa representa a avaliação COMPLETA de risco de erosão "
            "de eficácia, integrando: (1) barreira molecular intrínseca, "
            "(2) contexto regional de seleção, e (3) estratégia de manejo. "
            "Este NÃO é uma previsão de vida útil, mas uma estimativa comparativa "
            "de risco que pode orientar decisões de manejo de resistência."
        )
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    logger.info(f"IDA Programa exportado: {output_path}")
    return results

def print_ida_programa_summary(result: IDAProgramaResult) -> None:
    print("\n" + "=" * 70)
    print("IDA PROGRAMA - AVALIAÇÃO COMPLETA")
    print("=" * 70)
    print(f"Programa: {result.program_name} ({result.program_id})")
    print(f"Região: {result.region_name} ({result.region_id})")
    print(f"Safra: {result.season}")
    print(f"Data: {result.date_generated}")
    print("-" * 70)
    print("SCORES (TODAS AS CAMADAS):")
    print(f"  IDA Molecular: {result.ida_molecular:.4f}")
    print(f"  IDA Regional:  {result.ida_regional:.4f}")
    print(f"  IDA Programa:  {result.ida_programa:.4f}")
    print("-" * 70)
    print("ANÁLISE DO PROGRAMA:")
    print(f"  Proteção por mistura: {result.mixture_protection:.3f}")
    print(f"  Score de rotação:     {result.rotation_score:.3f}")
    print(f"  Bônus do programa:    {result.program_bonus:+.3f}")
    effect = "protetor" if result.program_bonus > 0 else "acelerador" if result.program_bonus < 0 else "neutro"
    print(f"  Efeito:               {effect}")
    print("-" * 70)
    print(f"Interpretação: {result.interpretation}")
    print("-" * 70)
    print("Recomendações:")
    for i, rec in enumerate(result.recommendations, 1):
        print(f"  {i}. {rec}")
    print("-" * 70)
    print("NOTA: " + result.disclaimer[:80] + "...")
    print("=" * 70)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cálculo do IDA Programa")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        ida_regional = IDARegionalResult(
            region_id="BR-MT",
            region_name="Mato Grosso",
            season="2024/25",
            ida_molecular=0.72,
            ida_regional=0.58,
        )

        prog1 = ManagementProgram(
            program_id="P001",
            name="Rotação SDHI + Triazol + Multissítio",
            applications=[
                Application("V8", "multisite", "Mancozeb", None),
                Application("R1", "SDHI", "Fluxapiroxade", "Piraclostrobina"),
                Application("R3", "DMI", "Protioconazol", "Trifloxistrobina"),
                Application("R5.1", "multisite", "Clorotalonil", None),
            ]
        )

        prog2 = ManagementProgram(
            program_id="P002",
            name="SDHI Intensivo",
            applications=[
                Application("V8", "SDHI", "Bixafen", None),
                Application("R1", "SDHI", "Fluxapiroxade", "Piraclostrobina"),
                Application("R3", "SDHI", "Benzovindiflupir", None),
                Application("R5.1", "DMI", "Protioconazol", None),
            ]
        )

        result1 = calculate_ida_programa(ida_regional, prog1)
        result2 = calculate_ida_programa(ida_regional, prog2)

        print_ida_programa_summary(result1)
        print_ida_programa_summary(result2)

        print("\n" + "=" * 70)
        print("COMPARATIVO ENTRE PROGRAMAS")
        print("=" * 70)
        comparison = compare_programs([result1, result2])
        print(f"Melhor programa: {comparison['best_program']}")
        print(f"Pior programa: {comparison['worst_program']}")
        print(f"Spread de IDA: {comparison['ida_range']['spread']:.4f}")
        print("\nRanking:")
        for r in comparison['ranking']:
            print(f"  {r['rank']}. {r['program']}: IDA_prog={r['ida_programa']:.3f} ({r['effect']})")