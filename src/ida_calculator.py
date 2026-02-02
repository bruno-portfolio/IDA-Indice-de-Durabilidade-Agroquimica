import json
import logging
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

from .config import SDH_TARGET, RESULTS_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class IDAWeights:
    srs: float = 1/3
    sce: float = 1/3
    sel: float = 1/3

    def __post_init__(self) -> None:
        total = self.srs + self.sce + self.sel
        if abs(total - 1.0) > 0.01:
            logger.warning(f"Pesos não somam 1.0: {total}")

@dataclass
class ResidueIDA:
    res_seq: int
    res_name: str
    chain_id: str

    srs: float
    sce: float
    sel: float
    sce_weight: float
    plddt: float

    ida_raw: float
    ida: float

    is_binding_site: bool = True

@dataclass
class IDAResult:
    target: str
    organism: str
    defensivo_class: str
    date_generated: str

    residue_scores: list[ResidueIDA]

    ida_aggregate: float
    ida_mean: float
    ida_std: float
    ida_min: float
    ida_max: float

    n_residues: int
    avg_plddt: float
    avg_srs: float
    avg_sce: float
    avg_sel: float

    interpretation: str
    confidence_level: str

    weights: IDAWeights = field(default_factory=IDAWeights)

    def to_dict(self) -> dict:
        return {
            "metadata": {
                "target": self.target,
                "organism": self.organism,
                "defensivo_class": self.defensivo_class,
                "date_generated": self.date_generated,
            },
            "aggregate_score": {
                "ida": round(self.ida_aggregate, 4),
                "ida_mean": round(self.ida_mean, 4),
                "ida_std": round(self.ida_std, 4),
                "ida_min": round(self.ida_min, 4),
                "ida_max": round(self.ida_max, 4),
                "interpretation": self.interpretation,
                "confidence_level": self.confidence_level,
            },
            "component_summary": {
                "avg_srs": round(self.avg_srs, 4),
                "avg_sce": round(self.avg_sce, 4),
                "avg_sel": round(self.avg_sel, 4),
                "avg_plddt": round(self.avg_plddt, 1),
                "n_residues": self.n_residues,
            },
            "weights": asdict(self.weights),
            "residue_details": [
                {
                    "res_seq": r.res_seq,
                    "res_name": r.res_name,
                    "chain_id": r.chain_id,
                    "srs": round(r.srs, 4),
                    "sce": round(r.sce, 4),
                    "sel": round(r.sel, 4),
                    "ida": round(r.ida, 4),
                    "plddt": round(r.plddt, 1),
                }
                for r in self.residue_scores
            ],
        }

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

def _interpret_ida(ida_score: float) -> str:
    if ida_score >= 0.8:
        return "Alta durabilidade esperada - sítio altamente restrito e conservado"
    elif ida_score >= 0.6:
        return "Durabilidade moderada-alta - bom perfil estrutural e evolutivo"
    elif ida_score >= 0.4:
        return "Durabilidade moderada - monitoramento de eficácia recomendado"
    elif ida_score >= 0.2:
        return "Durabilidade baixa-moderada - atenção ao manejo e rotação"
    else:
        return "Durabilidade baixa - alto risco, priorizar rotação de mecanismos"

def _get_confidence_level(avg_plddt: float) -> str:
    if avg_plddt >= 90:
        return "very_high"
    elif avg_plddt >= 70:
        return "high"
    elif avg_plddt >= 50:
        return "medium"
    else:
        return "low"

def calculate_ida(
    binding_site_residues: list[int],
    residue_names: dict[int, str],
    chain_ids: dict[int, str],
    srs_scores: dict[int, float],
    sce_scores: dict[int, float],
    sel_scores: dict[int, float],
    sce_weights: dict[int, float],
    plddt_scores: dict[int, float],
    weights: Optional[IDAWeights] = None,
    normalize: bool = True,
) -> IDAResult:
    if weights is None:
        weights = IDAWeights()

    if normalize:
        srs_scores = _normalize_scores(srs_scores)
        sce_scores = _normalize_scores(sce_scores)
        sel_scores = _normalize_scores(sel_scores)

    residue_scores = []

    for res_seq in binding_site_residues:
        srs = srs_scores.get(res_seq, 0.5)
        sce = sce_scores.get(res_seq, 0.5)
        sel = sel_scores.get(res_seq, 0.5)
        sce_weight = sce_weights.get(res_seq, 0.8)
        plddt = plddt_scores.get(res_seq, 70.0)

        ida_raw = (
            weights.srs * srs +
            weights.sce * sce +
            weights.sel * (1 - sel)
        )

        ida = ida_raw * sce_weight

        residue_scores.append(ResidueIDA(
            res_seq=res_seq,
            res_name=residue_names.get(res_seq, "UNK"),
            chain_id=chain_ids.get(res_seq, "A"),
            srs=srs,
            sce=sce,
            sel=sel,
            sce_weight=sce_weight,
            plddt=plddt,
            ida_raw=ida_raw,
            ida=ida,
        ))

    if not residue_scores:
        logger.error("Nenhum resíduo processado")
        return IDAResult(
            target=SDH_TARGET.enzyme,
            organism=SDH_TARGET.organism,
            defensivo_class=SDH_TARGET.defensivo_class,
            date_generated=datetime.now().isoformat(),
            residue_scores=[],
            ida_aggregate=0.0,
            ida_mean=0.0,
            ida_std=0.0,
            ida_min=0.0,
            ida_max=0.0,
            n_residues=0,
            avg_plddt=0.0,
            avg_srs=0.0,
            avg_sce=0.0,
            avg_sel=0.0,
            interpretation="Erro: nenhum resíduo processado",
            confidence_level="none",
            weights=weights,
        )

    ida_values = [r.ida for r in residue_scores]
    sce_weight_values = [r.sce_weight for r in residue_scores]

    if sum(sce_weight_values) > 0:
        ida_aggregate = sum(r.ida * r.sce_weight for r in residue_scores) / sum(sce_weight_values)
    else:
        ida_aggregate = sum(ida_values) / len(ida_values)

    ida_mean = sum(ida_values) / len(ida_values)
    ida_std = (sum((x - ida_mean) ** 2 for x in ida_values) / len(ida_values)) ** 0.5

    avg_plddt = sum(r.plddt for r in residue_scores) / len(residue_scores)
    avg_srs = sum(r.srs for r in residue_scores) / len(residue_scores)
    avg_sce = sum(r.sce for r in residue_scores) / len(residue_scores)
    avg_sel = sum(r.sel for r in residue_scores) / len(residue_scores)

    logger.info(
        f"IDA calculado: {len(residue_scores)} resíduos, "
        f"IDA={ida_aggregate:.3f}, pLDDT médio={avg_plddt:.1f}"
    )

    return IDAResult(
        target=SDH_TARGET.enzyme,
        organism=SDH_TARGET.organism,
        defensivo_class=SDH_TARGET.defensivo_class,
        date_generated=datetime.now().isoformat(),
        residue_scores=residue_scores,
        ida_aggregate=ida_aggregate,
        ida_mean=ida_mean,
        ida_std=ida_std,
        ida_min=min(ida_values),
        ida_max=max(ida_values),
        n_residues=len(residue_scores),
        avg_plddt=avg_plddt,
        avg_srs=avg_srs,
        avg_sce=avg_sce,
        avg_sel=avg_sel,
        interpretation=_interpret_ida(ida_aggregate),
        confidence_level=_get_confidence_level(avg_plddt),
        weights=weights,
    )

def export_ida_json(
    ida_result: IDAResult,
    output_path: Optional[Path] = None,
) -> Path:
    if output_path is None:
        output_path = RESULTS_DIR / "ida_results.json"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(ida_result.to_dict(), f, indent=2, ensure_ascii=False)

    logger.info(f"Resultado exportado: {output_path}")
    return output_path

def export_ida_csv(
    ida_result: IDAResult,
    output_path: Optional[Path] = None,
) -> Path:
    if output_path is None:
        output_path = RESULTS_DIR / "ida_residues.csv"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = ["res_seq,res_name,chain_id,srs,sce,sel,sce_weight,plddt,ida_raw,ida"]

    for r in ida_result.residue_scores:
        lines.append(
            f"{r.res_seq},{r.res_name},{r.chain_id},"
            f"{r.srs:.4f},{r.sce:.4f},{r.sel:.4f},"
            f"{r.sce_weight:.4f},{r.plddt:.1f},{r.ida_raw:.4f},{r.ida:.4f}"
        )

    output_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info(f"CSV exportado: {output_path}")
    return output_path

def print_ida_summary(ida_result: IDAResult) -> None:
    print("\n" + "=" * 60)
    print("ÍNDICE DE DURABILIDADE AGROQUÍMICA (IDA)")
    print("=" * 60)
    print(f"Alvo: {ida_result.target}")
    print(f"Organismo: {ida_result.organism}")
    print(f"Classe: {ida_result.defensivo_class}")
    print(f"Data: {ida_result.date_generated}")
    print("-" * 60)
    print(f"IDA Agregado: {ida_result.ida_aggregate:.4f}")
    print(f"IDA Médio: {ida_result.ida_mean:.4f} (±{ida_result.ida_std:.4f})")
    print(f"IDA Range: [{ida_result.ida_min:.4f}, {ida_result.ida_max:.4f}]")
    print("-" * 60)
    print(f"Componentes médios:")
    print(f"  SRS (Restrição): {ida_result.avg_srs:.4f}")
    print(f"  SCE (Conservação): {ida_result.avg_sce:.4f}")
    print(f"  SEL (Estabilidade): {ida_result.avg_sel:.4f}")
    print(f"  pLDDT médio: {ida_result.avg_plddt:.1f}")
    print("-" * 60)
    print(f"Interpretação: {ida_result.interpretation}")
    print(f"Nível de confiança: {ida_result.confidence_level}")
    print(f"Resíduos analisados: {ida_result.n_residues}")
    print("=" * 60)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cálculo do IDA")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        binding_site = [145, 146, 147, 230, 231, 272, 273]
        residue_names = {145: "TYR", 146: "GLY", 147: "ALA", 230: "HIS", 231: "SER", 272: "PRO", 273: "ILE"}
        chain_ids = {r: "B" for r in binding_site}

        srs_scores = {145: 0.85, 146: 0.72, 147: 0.68, 230: 0.91, 231: 0.65, 272: 0.78, 273: 0.82}
        sce_scores = {145: 0.92, 146: 0.88, 147: 0.75, 230: 0.95, 231: 0.70, 272: 0.85, 273: 0.80}
        sel_scores = {145: 0.25, 146: 0.35, 147: 0.45, 230: 0.20, 231: 0.50, 272: 0.30, 273: 0.40}
        sce_weights = {r: 0.9 for r in binding_site}
        plddt_scores = {145: 92.5, 146: 88.3, 147: 85.1, 230: 94.2, 231: 82.7, 272: 90.1, 273: 87.5}

        result = calculate_ida(
            binding_site_residues=binding_site,
            residue_names=residue_names,
            chain_ids=chain_ids,
            srs_scores=srs_scores,
            sce_scores=sce_scores,
            sel_scores=sel_scores,
            sce_weights=sce_weights,
            plddt_scores=plddt_scores,
        )

        print_ida_summary(result)

        json_path = export_ida_json(result)
        csv_path = export_ida_csv(result)
        print(f"\nExportado: {json_path}")
        print(f"Exportado: {csv_path}")