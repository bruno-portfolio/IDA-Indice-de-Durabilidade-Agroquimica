import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from .config import RESULTS_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class StructuralConfidence:
    score: float
    issues: list[str] = field(default_factory=list)
    mean_plddt: float = 0.0
    min_plddt: float = 0.0
    n_low_confidence_residues: int = 0

@dataclass
class EvolutionaryConfidence:
    score: float
    issues: list[str] = field(default_factory=list)
    n_sequences: int = 0
    alignment_quality: float = 0.0
    avg_gap_fraction: float = 0.0

@dataclass
class RegionalConfidence:
    score: float
    issues: list[str] = field(default_factory=list)
    n_mutations_documented: int = 0
    years_of_monitoring: int = 0
    n_data_sources: int = 0

@dataclass
class ConfidenceReport:

    structural_confidence: StructuralConfidence
    evolutionary_confidence: EvolutionaryConfidence
    regional_confidence: RegionalConfidence

    overall_confidence: float
    confidence_level: str

    recommendations: list[str] = field(default_factory=list)

    date_generated: str = ""
    version: str = "2.0"

    def to_dict(self) -> dict:
        return {
            "overall": {
                "score": round(self.overall_confidence, 3),
                "level": self.confidence_level,
            },
            "breakdown": {
                "structural": {
                    "score": round(self.structural_confidence.score, 3),
                    "issues": self.structural_confidence.issues,
                    "metrics": {
                        "mean_plddt": self.structural_confidence.mean_plddt,
                        "min_plddt": self.structural_confidence.min_plddt,
                        "n_low_confidence_residues": self.structural_confidence.n_low_confidence_residues,
                    }
                },
                "evolutionary": {
                    "score": round(self.evolutionary_confidence.score, 3),
                    "issues": self.evolutionary_confidence.issues,
                    "metrics": {
                        "n_sequences": self.evolutionary_confidence.n_sequences,
                        "alignment_quality": self.evolutionary_confidence.alignment_quality,
                        "avg_gap_fraction": self.evolutionary_confidence.avg_gap_fraction,
                    }
                },
                "regional": {
                    "score": round(self.regional_confidence.score, 3),
                    "issues": self.regional_confidence.issues,
                    "metrics": {
                        "n_mutations_documented": self.regional_confidence.n_mutations_documented,
                        "years_of_monitoring": self.regional_confidence.years_of_monitoring,
                        "n_data_sources": self.regional_confidence.n_data_sources,
                    }
                },
            },
            "recommendations": self.recommendations,
            "interpretation": self._get_interpretation(),
            "metadata": {
                "date_generated": self.date_generated,
                "version": self.version,
            }
        }

    def _get_interpretation(self) -> str:
        level_text = {
            "alta": "Confiança ALTA: resultados robustos, baixa incerteza",
            "média": "Confiança MÉDIA: resultados indicativos, mas com incerteza moderada",
            "baixa": "Confiança BAIXA: resultados devem ser interpretados com cautela",
            "muito baixa": "Confiança MUITO BAIXA: resultados preliminares, incerteza significativa",
        }

        base = level_text.get(self.confidence_level, "Confiança indeterminada")

        weak_areas = []
        if self.structural_confidence.score < 0.6:
            weak_areas.append("estrutura")
        if self.evolutionary_confidence.score < 0.6:
            weak_areas.append("dados evolutivos")
        if self.regional_confidence.score < 0.6:
            weak_areas.append("dados regionais")

        if weak_areas:
            base += f" - principais fontes de incerteza: {', '.join(weak_areas)}"

        return base

def assess_structural_confidence(
    structural_metrics: dict,
) -> StructuralConfidence:
    mean_plddt = structural_metrics.get("mean_plddt", 70)
    min_plddt = structural_metrics.get("min_plddt", 50)
    n_low_conf = structural_metrics.get("n_low_confidence_residues", 0)

    if mean_plddt >= 90:
        conf_score = 1.0
    elif mean_plddt >= 70:
        conf_score = 0.8
    elif mean_plddt >= 50:
        conf_score = 0.5
    else:
        conf_score = 0.3

    issues = []

    if min_plddt < 70:
        issues.append(f"pLDDT mínimo no sítio = {min_plddt:.1f} (< 70)")
        conf_score *= 0.9

    if n_low_conf > 2:
        issues.append(f"{n_low_conf} resíduos com baixa confiança estrutural")
        conf_score *= 0.95

    if mean_plddt < 70:
        issues.append(f"pLDDT médio = {mean_plddt:.1f} - estrutura pode não ser confiável")

    return StructuralConfidence(
        score=conf_score,
        issues=issues,
        mean_plddt=mean_plddt,
        min_plddt=min_plddt,
        n_low_confidence_residues=n_low_conf,
    )

def assess_evolutionary_confidence(
    evolutionary_metrics: dict,
) -> EvolutionaryConfidence:
    n_sequences = evolutionary_metrics.get("n_sequences", 0)
    alignment_quality = evolutionary_metrics.get("alignment_quality", 0.5)
    gap_fraction = evolutionary_metrics.get("avg_gap_fraction", 0)

    if n_sequences >= 50 and alignment_quality >= 0.8:
        conf_score = 1.0
    elif n_sequences >= 20 and alignment_quality >= 0.6:
        conf_score = 0.7
    elif n_sequences >= 10:
        conf_score = 0.5
    else:
        conf_score = 0.3

    issues = []

    if n_sequences < 20:
        issues.append(f"Apenas {n_sequences} sequências no MSA (recomendado: 50+)")

    if gap_fraction > 0.3:
        issues.append(f"Alta fração de gaps no alinhamento ({gap_fraction:.0%})")
        conf_score *= 0.9

    if alignment_quality < 0.6:
        issues.append(f"Qualidade do alinhamento baixa ({alignment_quality:.2f})")

    return EvolutionaryConfidence(
        score=conf_score,
        issues=issues,
        n_sequences=n_sequences,
        alignment_quality=alignment_quality,
        avg_gap_fraction=gap_fraction,
    )

def assess_regional_confidence(
    regional_metrics: dict,
) -> RegionalConfidence:
    n_mutations = regional_metrics.get("n_mutations_documented", 0)
    years_monitoring = regional_metrics.get("years_of_monitoring", 0)
    data_sources = regional_metrics.get("n_data_sources", 0)

    if n_mutations >= 5 and years_monitoring >= 5:
        conf_score = 1.0
    elif n_mutations >= 2 and years_monitoring >= 3:
        conf_score = 0.7
    elif n_mutations >= 1:
        conf_score = 0.5
    else:
        conf_score = 0.3

    issues = []

    if years_monitoring < 3:
        issues.append(f"Apenas {years_monitoring} anos de monitoramento documentado")

    if data_sources < 2:
        issues.append("Poucas fontes de dados regionais")

    if n_mutations == 0:
        issues.append("Nenhuma mutação documentada - dados podem ser incompletos")

    return RegionalConfidence(
        score=conf_score,
        issues=issues,
        n_mutations_documented=n_mutations,
        years_of_monitoring=years_monitoring,
        n_data_sources=data_sources,
    )

def assess_confidence(
    structural_metrics: dict,
    evolutionary_metrics: dict,
    regional_metrics: dict,
    weights: Optional[dict] = None,
) -> ConfidenceReport:
    if weights is None:
        weights = {
            "structural": 0.40,
            "evolutionary": 0.35,
            "regional": 0.25,
        }

    structural = assess_structural_confidence(structural_metrics)
    evolutionary = assess_evolutionary_confidence(evolutionary_metrics)
    regional = assess_regional_confidence(regional_metrics)

    overall = (
        weights["structural"] * structural.score +
        weights["evolutionary"] * evolutionary.score +
        weights["regional"] * regional.score
    )

    if overall >= 0.8:
        level = "alta"
    elif overall >= 0.6:
        level = "média"
    elif overall >= 0.4:
        level = "baixa"
    else:
        level = "muito baixa"

    recommendations = _generate_confidence_recommendations(
        structural, evolutionary, regional, level
    )

    logger.info(f"Confiança avaliada: {level} ({overall:.2f})")

    return ConfidenceReport(
        structural_confidence=structural,
        evolutionary_confidence=evolutionary,
        regional_confidence=regional,
        overall_confidence=overall,
        confidence_level=level,
        recommendations=recommendations,
        date_generated=datetime.now().isoformat(),
    )

def _generate_confidence_recommendations(
    structural: StructuralConfidence,
    evolutionary: EvolutionaryConfidence,
    regional: RegionalConfidence,
    level: str,
) -> list[str]:
    recommendations = []

    if level == "very_low":
        recommendations.append(
            "ATENÇÃO: Confiança muito baixa - resultados devem ser "
            "interpretados com extrema cautela"
        )

    if structural.min_plddt < 70:
        recommendations.append(
            "Validar estrutura do sítio de ligação com dados experimentais se disponíveis"
        )

    if evolutionary.n_sequences < 20:
        recommendations.append("Expandir busca de homólogos para melhorar MSA")

    if regional.score < 0.5:
        recommendations.append(
            "Dados regionais limitados - IDA Regional pode não refletir situação real"
        )

    if regional.years_of_monitoring < 3:
        recommendations.append(
            "Monitoramento recente - validar com dados de mais safras"
        )

    if not recommendations:
        recommendations.append(
            "Dados de boa qualidade - resultados podem ser considerados confiáveis"
        )

    return recommendations

@dataclass
class IDAResultWithConfidence:

    ida_molecular: float
    ida_regional: float
    ida_programa: float

    interpretation: str

    confidence: ConfidenceReport

    def to_dict(self) -> dict:
        return {
            "ida_result": {
                "ida_molecular": round(self.ida_molecular, 4),
                "ida_regional": round(self.ida_regional, 4),
                "ida_programa": round(self.ida_programa, 4),
                "interpretation": self.interpretation,
            },
            "confidence": self.confidence.to_dict(),
        }

def create_ida_with_confidence(
    ida_molecular: float,
    ida_regional: float,
    ida_programa: float,
    interpretation: str,
    structural_metrics: dict,
    evolutionary_metrics: dict,
    regional_metrics: dict,
) -> IDAResultWithConfidence:
    confidence = assess_confidence(
        structural_metrics,
        evolutionary_metrics,
        regional_metrics,
    )

    return IDAResultWithConfidence(
        ida_molecular=ida_molecular,
        ida_regional=ida_regional,
        ida_programa=ida_programa,
        interpretation=interpretation,
        confidence=confidence,
    )

def export_confidence_report(
    report: ConfidenceReport,
    output_path: Optional[Path] = None,
) -> Path:
    if output_path is None:
        output_path = RESULTS_DIR / "confidence_report.json"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(report.to_dict(), f, indent=2, ensure_ascii=False)

    logger.info(f"Reporte de confiança exportado: {output_path}")
    return output_path

def print_confidence_summary(report: ConfidenceReport) -> None:
    print("\n" + "=" * 70)
    print("REPORTE DE CONFIANÇA")
    print("=" * 70)
    print(f"Confiança Geral: {report.overall_confidence:.2f} ({report.confidence_level.upper()})")
    print("-" * 70)
    print("BREAKDOWN:")
    print(f"  Estrutural:  {report.structural_confidence.score:.2f}")
    if report.structural_confidence.issues:
        for issue in report.structural_confidence.issues:
            print(f"    - {issue}")
    print(f"  Evolutivo:   {report.evolutionary_confidence.score:.2f}")
    if report.evolutionary_confidence.issues:
        for issue in report.evolutionary_confidence.issues:
            print(f"    - {issue}")
    print(f"  Regional:    {report.regional_confidence.score:.2f}")
    if report.regional_confidence.issues:
        for issue in report.regional_confidence.issues:
            print(f"    - {issue}")
    print("-" * 70)
    print("RECOMENDAÇÕES:")
    for i, rec in enumerate(report.recommendations, 1):
        print(f"  {i}. {rec}")
    print("=" * 70)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Sistema de Reporte de Confiança")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        structural_metrics = {
            "mean_plddt": 85.5,
            "min_plddt": 68.2,
            "n_low_confidence_residues": 2,
        }

        evolutionary_metrics = {
            "n_sequences": 32,
            "alignment_quality": 0.72,
            "avg_gap_fraction": 0.15,
        }

        regional_metrics = {
            "n_mutations_documented": 3,
            "years_of_monitoring": 4,
            "n_data_sources": 2,
        }

        report = assess_confidence(
            structural_metrics,
            evolutionary_metrics,
            regional_metrics,
        )

        print_confidence_summary(report)

        result = create_ida_with_confidence(
            ida_molecular=0.72,
            ida_regional=0.58,
            ida_programa=0.64,
            interpretation="Risco moderado - programa protetor, monitoramento recomendado",
            structural_metrics=structural_metrics,
            evolutionary_metrics=evolutionary_metrics,
            regional_metrics=regional_metrics,
        )

        print("\n" + "=" * 70)
        print("RESULTADO IDA COMPLETO COM CONFIANÇA")
        print("=" * 70)
        print(json.dumps(result.to_dict(), indent=2, ensure_ascii=False))