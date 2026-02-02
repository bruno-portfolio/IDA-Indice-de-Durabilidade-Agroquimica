import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .config import STRUCTURES_DIR, ANALYSIS_PARAMS, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class ConfidenceLevel:
    VERY_HIGH: float = 90.0
    CONFIDENT: float = 70.0
    LOW: float = 50.0

@dataclass
class ResidueConfidence:
    residue_idx: int
    plddt: float
    confidence_level: str
    sce_weight: float

    @property
    def is_reliable(self) -> bool:
        return self.plddt >= ConfidenceLevel.CONFIDENT

@dataclass
class ConfidenceData:
    residue_scores: list[ResidueConfidence]
    pae_matrix: Optional[list[list[float]]] = None

    @property
    def n_residues(self) -> int:
        return len(self.residue_scores)

    @property
    def mean_plddt(self) -> float:
        if not self.residue_scores:
            return 0.0
        return sum(r.plddt for r in self.residue_scores) / len(self.residue_scores)

    @property
    def n_reliable(self) -> int:
        return sum(1 for r in self.residue_scores if r.is_reliable)

    @property
    def pct_reliable(self) -> float:
        if not self.residue_scores:
            return 0.0
        return self.n_reliable / len(self.residue_scores) * 100

    def get_residue(self, idx: int) -> Optional[ResidueConfidence]:
        if 1 <= idx <= len(self.residue_scores):
            return self.residue_scores[idx - 1]
        return None

    def get_sce_weights(self, residue_indices: list[int]) -> dict[int, float]:
        return {
            idx: self.residue_scores[idx - 1].sce_weight
            for idx in residue_indices
            if 1 <= idx <= len(self.residue_scores)
        }

def _classify_confidence(plddt: float) -> str:
    if plddt >= ConfidenceLevel.VERY_HIGH:
        return "muito alta"
    elif plddt >= ConfidenceLevel.CONFIDENT:
        return "alta"
    elif plddt >= ConfidenceLevel.LOW:
        return "baixa"
    else:
        return "muito baixa"

def _calculate_sce_weight(plddt: float) -> float:
    if plddt >= ANALYSIS_PARAMS.plddt_very_high:
        return 1.0
    elif plddt >= ANALYSIS_PARAMS.plddt_confident:
        return 0.8
    elif plddt >= ANALYSIS_PARAMS.plddt_low:
        return 0.5
    else:
        return 0.2

def load_confidence_json(json_path: Path) -> Optional[ConfidenceData]:
    json_path = Path(json_path)

    if not json_path.exists():
        logger.error(f"Arquivo não encontrado: {json_path}")
        return None

    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        logger.error(f"Erro ao decodificar JSON: {e}")
        return None
    except OSError as e:
        logger.error(f"Erro ao ler arquivo: {e}")
        return None

    plddt_values = data.get("confidenceScore") or data.get("plddt", [])

    if not plddt_values:
        logger.warning("Nenhum score pLDDT encontrado no JSON")
        return None

    residue_scores = []
    for idx, plddt in enumerate(plddt_values, start=1):
        residue_scores.append(ResidueConfidence(
            residue_idx=idx,
            plddt=float(plddt),
            confidence_level=_classify_confidence(plddt),
            sce_weight=_calculate_sce_weight(plddt),
        ))

    pae_matrix = data.get("pae") or data.get("predicted_aligned_error")

    logger.info(
        f"Carregados {len(residue_scores)} scores pLDDT, "
        f"média={sum(r.plddt for r in residue_scores)/len(residue_scores):.1f}"
    )

    return ConfidenceData(
        residue_scores=residue_scores,
        pae_matrix=pae_matrix,
    )

def extract_plddt_from_pdb(pdb_path: Path) -> Optional[ConfidenceData]:
    pdb_path = Path(pdb_path)

    if not pdb_path.exists():
        logger.error(f"Arquivo não encontrado: {pdb_path}")
        return None

    residue_plddt: dict[int, list[float]] = {}

    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue

            try:
                res_seq = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())

                if res_seq not in residue_plddt:
                    residue_plddt[res_seq] = []
                residue_plddt[res_seq].append(b_factor)
            except (ValueError, IndexError):
                continue

    if not residue_plddt:
        logger.warning("Nenhum B-factor extraído do PDB")
        return None

    residue_scores = []
    for res_seq in sorted(residue_plddt.keys()):
        values = residue_plddt[res_seq]
        avg_plddt = sum(values) / len(values)

        residue_scores.append(ResidueConfidence(
            residue_idx=res_seq,
            plddt=avg_plddt,
            confidence_level=_classify_confidence(avg_plddt),
            sce_weight=_calculate_sce_weight(avg_plddt),
        ))

    logger.info(
        f"Extraídos {len(residue_scores)} scores pLDDT do PDB, "
        f"média={sum(r.plddt for r in residue_scores)/len(residue_scores):.1f}"
    )

    return ConfidenceData(residue_scores=residue_scores)

def load_confidence(
    pdb_path: Optional[Path] = None,
    json_path: Optional[Path] = None,
) -> Optional[ConfidenceData]:
    if json_path:
        result = load_confidence_json(json_path)
        if result:
            return result

    if pdb_path:
        return extract_plddt_from_pdb(pdb_path)

    logger.error("Nenhum arquivo fornecido")
    return None

def analyze_binding_site_confidence(
    confidence_data: ConfidenceData,
    binding_site_residues: list[int],
) -> dict:
    if not binding_site_residues:
        return {"error": "Nenhum resíduo no sítio de ligação"}

    scores = []
    for res_idx in binding_site_residues:
        residue = confidence_data.get_residue(res_idx)
        if residue:
            scores.append(residue)

    if not scores:
        return {"error": "Nenhum resíduo encontrado nos dados de confiança"}

    plddt_values = [r.plddt for r in scores]
    sce_values = [r.sce_weight for r in scores]

    return {
        "n_residues": len(scores),
        "mean_plddt": sum(plddt_values) / len(plddt_values),
        "min_plddt": min(plddt_values),
        "max_plddt": max(plddt_values),
        "mean_sce_weight": sum(sce_values) / len(sce_values),
        "n_very_high": sum(1 for r in scores if r.confidence_level == "very_high"),
        "n_confident": sum(1 for r in scores if r.confidence_level == "confident"),
        "n_low": sum(1 for r in scores if r.confidence_level == "low"),
        "n_very_low": sum(1 for r in scores if r.confidence_level == "very_low"),
        "pct_reliable": sum(1 for r in scores if r.is_reliable) / len(scores) * 100,
    }

def get_confidence_summary(confidence_data: ConfidenceData) -> dict:
    if not confidence_data.residue_scores:
        return {"error": "Sem dados de confiança"}

    plddt_values = [r.plddt for r in confidence_data.residue_scores]
    levels = [r.confidence_level for r in confidence_data.residue_scores]

    return {
        "n_residues": len(plddt_values),
        "mean_plddt": sum(plddt_values) / len(plddt_values),
        "min_plddt": min(plddt_values),
        "max_plddt": max(plddt_values),
        "std_plddt": (sum((x - sum(plddt_values)/len(plddt_values))**2 for x in plddt_values) / len(plddt_values)) ** 0.5,
        "n_very_high": levels.count("very_high"),
        "n_confident": levels.count("confident"),
        "n_low": levels.count("low"),
        "n_very_low": levels.count("very_low"),
        "pct_reliable": confidence_data.pct_reliable,
        "has_pae": confidence_data.pae_matrix is not None,
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extração de confiança AlphaFold")
    parser.add_argument(
        "-j", "--json",
        type=Path,
        help="Arquivo JSON de confiança",
    )
    parser.add_argument(
        "-p", "--pdb",
        type=Path,
        help="Arquivo PDB (extrai do B-factor)",
    )
    parser.add_argument(
        "-r", "--residues",
        type=int,
        nargs="+",
        help="Resíduos do sítio de ligação para análise",
    )

    args = parser.parse_args()

    if not args.json and not args.pdb:
        parser.error("Forneça --json ou --pdb")

    confidence = load_confidence(pdb_path=args.pdb, json_path=args.json)

    if confidence:
        summary = get_confidence_summary(confidence)
        print("\nResumo de confiança:")
        for key, value in summary.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.2f}")
            else:
                print(f"  {key}: {value}")

        if args.residues:
            print("\nAnálise do sítio de ligação:")
            bs_analysis = analyze_binding_site_confidence(confidence, args.residues)
            for key, value in bs_analysis.items():
                if isinstance(value, float):
                    print(f"  {key}: {value:.2f}")
                else:
                    print(f"  {key}: {value}")
    else:
        print("Erro ao carregar dados de confiança")