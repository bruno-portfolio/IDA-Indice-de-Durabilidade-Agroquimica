import logging
from dataclasses import dataclass
from typing import Optional

from .config import ANALYSIS_PARAMS, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

BLOSUM62: dict[str, dict[str, int]] = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1,
          'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3,
          'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3,
          'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3,
          'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1,
          'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3,
          'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3,
          'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4,
          'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3,
          'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4,
          'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2,
          'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3,
          'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1,
          'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0,
          'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3,
          'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2,
          'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1,
          'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3,
          'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1,
          'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3,
          'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4},
}

HYDROPHOBICITY: dict[str, float] = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

@dataclass
class StabilityScore:
    res_seq: int
    res_name: str
    blosum_avg: float
    tolerance: float
    sel_score: float

    @property
    def is_sensitive(self) -> bool:
        return self.sel_score < 0.3

@dataclass
class StabilityResult:
    scores: dict[int, StabilityScore]

    @property
    def mean_sel(self) -> float:
        if not self.scores:
            return 0.0
        return sum(s.sel_score for s in self.scores.values()) / len(self.scores)

    def get_sel_dict(self) -> dict[int, float]:
        return {pos: s.sel_score for pos, s in self.scores.items()}

def _get_blosum_score(aa1: str, aa2: str) -> int:
    aa1 = aa1.upper()
    aa2 = aa2.upper()

    if aa1 not in BLOSUM62 or aa2 not in BLOSUM62:
        return 0

    return BLOSUM62[aa1].get(aa2, 0)

def _calculate_substitution_tolerance(aa: str) -> float:
    if aa not in BLOSUM62:
        return 0.5

    scores = []
    for other_aa in AMINO_ACIDS:
        if other_aa != aa:
            scores.append(_get_blosum_score(aa, other_aa))

    if not scores:
        return 0.5

    avg_score = sum(scores) / len(scores)

    normalized = (avg_score + 4) / 8
    return max(0.0, min(1.0, normalized))

def calculate_stability_scores(
    sequence: str,
    binding_site_residues: list[int],
    burial_scores: Optional[dict[int, float]] = None,
    burial_weight: float = 0.5,
) -> StabilityResult:
    scores = {}

    for res_pos in binding_site_residues:
        if res_pos < 1 or res_pos > len(sequence):
            logger.warning(f"Posição {res_pos} fora do range da sequência")
            continue

        aa = sequence[res_pos - 1].upper()

        if aa not in AMINO_ACIDS:
            logger.warning(f"Aminoácido inválido na posição {res_pos}: {aa}")
            continue

        tolerance = _calculate_substitution_tolerance(aa)

        blosum_scores = [
            _get_blosum_score(aa, other)
            for other in AMINO_ACIDS
            if other != aa
        ]
        blosum_avg = sum(blosum_scores) / len(blosum_scores) if blosum_scores else 0

        if burial_scores and res_pos in burial_scores:
            burial = burial_scores[res_pos]
            adjusted_tolerance = tolerance * (1 - burial * burial_weight)
        else:
            adjusted_tolerance = tolerance

        sel_score = adjusted_tolerance

        scores[res_pos] = StabilityScore(
            res_seq=res_pos,
            res_name=aa,
            blosum_avg=blosum_avg,
            tolerance=tolerance,
            sel_score=sel_score,
        )

    logger.info(
        f"Estabilidade calculada para {len(scores)} resíduos, "
        f"SEL médio={sum(s.sel_score for s in scores.values())/len(scores):.3f}" if scores else ""
    )

    return StabilityResult(scores=scores)

def calculate_position_sensitivity(
    sequence: str,
    position: int,
) -> dict[str, float]:
    if position < 1 or position > len(sequence):
        return {}

    original_aa = sequence[position - 1].upper()

    if original_aa not in BLOSUM62:
        return {}

    sensitivities = {}
    for new_aa in AMINO_ACIDS:
        if new_aa != original_aa:
            score = _get_blosum_score(original_aa, new_aa)
            sensitivity = (4 - score) / 8
            sensitivities[new_aa] = max(0.0, min(1.0, sensitivity))

    return sensitivities

def get_hydrophobicity_score(residues: list[str]) -> float:
    if not residues:
        return 0.0

    scores = []
    for aa in residues:
        aa = aa.upper()
        if aa in HYDROPHOBICITY:
            scores.append(HYDROPHOBICITY[aa])

    if not scores:
        return 0.0

    avg = sum(scores) / len(scores)

    normalized = (avg + 4.5) / 9.0
    return max(0.0, min(1.0, normalized))

def analyze_stability_distribution(
    stability_result: StabilityResult,
) -> dict:
    if not stability_result.scores:
        return {"error": "Sem scores de estabilidade"}

    sel_values = [s.sel_score for s in stability_result.scores.values()]

    return {
        "n_residues": len(sel_values),
        "mean_sel": sum(sel_values) / len(sel_values),
        "min_sel": min(sel_values),
        "max_sel": max(sel_values),
        "n_sensitive": sum(1 for s in stability_result.scores.values() if s.is_sensitive),
        "pct_sensitive": sum(1 for s in stability_result.scores.values() if s.is_sensitive) / len(sel_values) * 100,
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cálculo de estabilidade")
    parser.add_argument(
        "-s", "--sequence",
        required=True,
        help="Sequência de aminoácidos",
    )
    parser.add_argument(
        "-r", "--residues",
        type=int,
        nargs="+",
        required=True,
        help="Resíduos do sítio de ligação",
    )

    args = parser.parse_args()

    result = calculate_stability_scores(args.sequence, args.residues)

    print("\nScores de estabilidade:")
    for pos, score in sorted(result.scores.items()):
        print(
            f"  {score.res_name}{pos}: "
            f"tolerance={score.tolerance:.3f}, "
            f"SEL={score.sel_score:.3f}"
        )

    analysis = analyze_stability_distribution(result)
    print("\nAnálise:")
    for key, value in analysis.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value}")