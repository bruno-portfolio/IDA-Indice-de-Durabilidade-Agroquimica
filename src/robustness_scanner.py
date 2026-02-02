import logging
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .config import LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

BLOSUM62_SCORES = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4},
}

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

@dataclass
class MutationEffect:
    position: int
    original_aa: str
    mutant_aa: str
    blosum_score: int
    estimated_degradation: float
    is_conservative: bool

@dataclass
class ResidueRobustness:
    position: int
    residue_name: str

    burial_score: float
    conservation_score: float
    mean_blosum_tolerance: float

    n_mutations_tested: int

    robustness_score: float
    fragility_score: float

    interpretation: str

    mutations: list[MutationEffect] = field(default_factory=list)

@dataclass
class RobustnessResult:
    method: str
    n_residues: int
    per_residue: dict[int, ResidueRobustness] = field(default_factory=dict)
    mean_robustness: float = 0.0
    min_robustness: float = 0.0
    max_robustness: float = 0.0
    fragile_positions: list[int] = field(default_factory=list)
    robust_positions: list[int] = field(default_factory=list)

def calculate_robustness_heuristic(
    binding_site_residues: list[int],
    residue_names: dict[int, str],
    conservation_scores: dict[int, float],
    burial_scores: dict[int, float],
) -> RobustnessResult:
    per_residue = {}

    for pos in binding_site_residues:
        res_name = residue_names.get(pos, "UNK")
        aa_1letter = AA_3TO1.get(res_name, "X")

        burial = burial_scores.get(pos, 0.5)
        conservation = conservation_scores.get(pos, 0.5)

        blosum_tolerance, mutations = _calculate_blosum_tolerance(aa_1letter, pos)

        fragility = burial * conservation * (1 - blosum_tolerance)

        robustness = 1 - fragility

        interpretation = _interpret_robustness(robustness, burial, conservation)

        per_residue[pos] = ResidueRobustness(
            position=pos,
            residue_name=res_name,
            burial_score=burial,
            conservation_score=conservation,
            mean_blosum_tolerance=blosum_tolerance,
            n_mutations_tested=len(mutations),
            mutations=mutations,
            robustness_score=robustness,
            fragility_score=fragility,
            interpretation=interpretation,
        )

    robustness_values = [r.robustness_score for r in per_residue.values()]

    fragile = [pos for pos, r in per_residue.items() if r.robustness_score < 0.4]
    robust = [pos for pos, r in per_residue.items() if r.robustness_score >= 0.7]

    logger.info(
        f"Robustez calculada (heurística): {len(per_residue)} resíduos, "
        f"média={np.mean(robustness_values):.3f}"
    )

    return RobustnessResult(
        method="heuristic",
        n_residues=len(per_residue),
        per_residue=per_residue,
        mean_robustness=np.mean(robustness_values) if robustness_values else 0.0,
        min_robustness=min(robustness_values) if robustness_values else 0.0,
        max_robustness=max(robustness_values) if robustness_values else 0.0,
        fragile_positions=fragile,
        robust_positions=robust,
    )

def _calculate_blosum_tolerance(
    original_aa: str,
    position: int,
) -> tuple[float, list[MutationEffect]]:
    if original_aa not in BLOSUM62_SCORES:
        return 0.5, []

    mutations = []
    blosum_scores = []

    for mutant_aa, score in BLOSUM62_SCORES[original_aa].items():
        if mutant_aa == original_aa:
            continue

        is_conservative = score >= 0

        degradation = max(0, min(1, (4 - score) / 8))

        mutations.append(MutationEffect(
            position=position,
            original_aa=original_aa,
            mutant_aa=mutant_aa,
            blosum_score=score,
            estimated_degradation=degradation,
            is_conservative=is_conservative,
        ))
        blosum_scores.append(score)

    if not blosum_scores:
        return 0.5, []

    mean_score = np.mean(blosum_scores)
    tolerance = (mean_score + 4) / 8
    tolerance = max(0, min(1, tolerance))

    return tolerance, mutations

def _interpret_robustness(
    robustness: float,
    burial: float,
    conservation: float,
) -> str:
    parts = []

    if robustness >= 0.7:
        parts.append("Robusto")
    elif robustness >= 0.4:
        parts.append("Moderado")
    else:
        parts.append("Frágil")

    if burial >= 0.7:
        parts.append("enterrado")
    elif burial <= 0.3:
        parts.append("exposto")

    if conservation >= 0.9:
        parts.append("extremamente conservado")
    elif conservation >= 0.7:
        parts.append("conservado")

    return " - ".join(parts)

@dataclass
class DockingVariant:
    mutation: str
    structure_path: str

@dataclass
class DockingResult:
    mutation: str
    score: float
    delta_score: float
    degradation_pct: float

def calculate_robustness_docking(
    wildtype_structure: str,
    variants: list[DockingVariant],
    ligand: str,
    docking_tool: str = "autodock_vina",
) -> RobustnessResult:
    logger.warning(
        "calculate_robustness_docking é um placeholder. "
        "Use calculate_robustness_heuristic para análise sem docking."
    )

    return RobustnessResult(
        method="docking",
        n_residues=0,
        per_residue={},
        mean_robustness=0.5,
        min_robustness=0.5,
        max_robustness=0.5,
        fragile_positions=[],
        robust_positions=[],
    )

def get_robustness_scores(result: RobustnessResult) -> dict[int, float]:
    return {pos: r.robustness_score for pos, r in result.per_residue.items()}

def summarize_robustness(result: RobustnessResult) -> dict:
    return {
        "method": result.method,
        "n_residues": result.n_residues,
        "mean_robustness": result.mean_robustness,
        "min_robustness": result.min_robustness,
        "max_robustness": result.max_robustness,
        "n_fragile": len(result.fragile_positions),
        "n_robust": len(result.robust_positions),
        "fragile_positions": result.fragile_positions,
        "robust_positions": result.robust_positions,
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Scanner de Robustez")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        positions = [145, 230, 272, 273]
        residue_names = {145: "TYR", 230: "HIS", 272: "PRO", 273: "ILE"}
        conservation_scores = {145: 0.92, 230: 0.95, 272: 0.72, 273: 0.80}
        burial_scores = {145: 0.85, 230: 0.90, 272: 0.65, 273: 0.75}

        result = calculate_robustness_heuristic(
            binding_site_residues=positions,
            residue_names=residue_names,
            conservation_scores=conservation_scores,
            burial_scores=burial_scores,
        )

        print("\n" + "=" * 60)
        print("ANÁLISE DE ROBUSTEZ")
        print("=" * 60)
        print(f"Método: {result.method}")
        print(f"Resíduos: {result.n_residues}")
        print("-" * 60)

        for pos, robustness in result.per_residue.items():
            print(f"\nPosição {pos} ({robustness.residue_name}):")
            print(f"  Enterramento:       {robustness.burial_score:.3f}")
            print(f"  Conservação:        {robustness.conservation_score:.3f}")
            print(f"  Tolerância BLOSUM:  {robustness.mean_blosum_tolerance:.3f}")
            print(f"  ROBUSTEZ:           {robustness.robustness_score:.3f}")
            print(f"  Interpretação:      {robustness.interpretation}")

        print("\n" + "-" * 60)
        print("SUMÁRIO:")
        summary = summarize_robustness(result)
        print(f"  Média: {summary['mean_robustness']:.3f}")
        print(f"  Range: [{summary['min_robustness']:.3f}, {summary['max_robustness']:.3f}]")
        print(f"  Posições frágeis ({summary['n_fragile']}): {summary['fragile_positions']}")
        print(f"  Posições robustas ({summary['n_robust']}): {summary['robust_positions']}")
        print("=" * 60)