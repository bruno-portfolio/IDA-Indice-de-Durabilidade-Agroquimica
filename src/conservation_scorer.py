import logging
import math
from collections import Counter
from dataclasses import dataclass
from typing import Optional

from .msa_runner import MSAResult, map_alignment_to_structure
from .config import LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
N_AMINO_ACIDS = 20

@dataclass
class ConservationScore:
    alignment_pos: int
    structure_pos: Optional[int]
    shannon_entropy: float
    conservation: float
    gap_fraction: float
    dominant_aa: str
    dominant_freq: float

    @property
    def is_conserved(self) -> bool:
        return self.conservation > 0.8

    @property
    def is_variable(self) -> bool:
        return self.conservation < 0.5

@dataclass
class ConservationResult:
    scores: list[ConservationScore]
    alignment_length: int
    n_sequences: int

    @property
    def mean_conservation(self) -> float:
        if not self.scores:
            return 0.0
        return sum(s.conservation for s in self.scores) / len(self.scores)

    @property
    def n_conserved(self) -> int:
        return sum(1 for s in self.scores if s.is_conserved)

    @property
    def n_variable(self) -> int:
        return sum(1 for s in self.scores if s.is_variable)

    def get_score(self, structure_pos: int) -> Optional[ConservationScore]:
        for score in self.scores:
            if score.structure_pos == structure_pos:
                return score
        return None

    def get_scores_dict(self) -> dict[int, float]:
        return {
            s.structure_pos: s.conservation
            for s in self.scores
            if s.structure_pos is not None
        }

def _calculate_shannon_entropy(column: list[str], include_gaps: bool = False) -> tuple[float, float]:
    if not column:
        return 0.0, 0.0

    if include_gaps:
        valid_chars = column
    else:
        valid_chars = [aa for aa in column if aa != "-" and aa in AMINO_ACIDS]

    if not valid_chars:
        return 0.0, 1.0

    counts = Counter(valid_chars)
    total = len(valid_chars)
    freqs = [count / total for count in counts.values()]

    entropy = -sum(p * math.log2(p) for p in freqs if p > 0)

    gap_fraction = column.count("-") / len(column) if column else 0.0

    return entropy, gap_fraction

def _entropy_to_conservation(entropy: float, max_entropy: Optional[float] = None) -> float:
    if max_entropy is None:
        max_entropy = math.log2(N_AMINO_ACIDS)

    if max_entropy == 0:
        return 1.0

    conservation = 1.0 - (entropy / max_entropy)
    return max(0.0, min(1.0, conservation))

def calculate_conservation(
    msa_result: MSAResult,
    query_id: Optional[str] = None,
    include_gaps: bool = False,
) -> ConservationResult:
    if not msa_result.success or not msa_result.sequences:
        logger.error("MSA inválido ou vazio")
        return ConservationResult(
            scores=[],
            alignment_length=0,
            n_sequences=0,
        )

    aln_to_struct = {}
    if query_id:
        aln_to_struct = map_alignment_to_structure(msa_result, query_id)

    scores = []
    max_entropy = math.log2(N_AMINO_ACIDS)

    for aln_pos, column in msa_result.iter_columns():
        entropy, gap_fraction = _calculate_shannon_entropy(column, include_gaps)

        conservation = _entropy_to_conservation(entropy, max_entropy)

        valid_aas = [aa for aa in column if aa != "-" and aa in AMINO_ACIDS]
        if valid_aas:
            aa_counts = Counter(valid_aas)
            dominant_aa, dominant_count = aa_counts.most_common(1)[0]
            dominant_freq = dominant_count / len(valid_aas)
        else:
            dominant_aa = "-"
            dominant_freq = 0.0

        struct_pos = aln_to_struct.get(aln_pos)

        scores.append(ConservationScore(
            alignment_pos=aln_pos,
            structure_pos=struct_pos,
            shannon_entropy=entropy,
            conservation=conservation,
            gap_fraction=gap_fraction,
            dominant_aa=dominant_aa,
            dominant_freq=dominant_freq,
        ))

    logger.info(
        f"Conservação calculada: {len(scores)} posições, "
        f"média={sum(s.conservation for s in scores)/len(scores):.3f}"
    )

    return ConservationResult(
        scores=scores,
        alignment_length=msa_result.alignment_length,
        n_sequences=msa_result.n_sequences,
    )

def get_binding_site_conservation(
    conservation_result: ConservationResult,
    binding_site_residues: list[int],
) -> dict[int, float]:
    scores = {}

    for res_pos in binding_site_residues:
        score = conservation_result.get_score(res_pos)
        if score:
            scores[res_pos] = score.conservation

    return scores

def analyze_conservation_distribution(
    conservation_result: ConservationResult,
    binding_site_residues: Optional[list[int]] = None,
) -> dict:
    if not conservation_result.scores:
        return {"error": "Sem scores de conservação"}

    all_scores = [s.conservation for s in conservation_result.scores if s.structure_pos]

    result = {
        "total_positions": len(all_scores),
        "mean_conservation": sum(all_scores) / len(all_scores) if all_scores else 0,
        "n_highly_conserved": sum(1 for s in all_scores if s > 0.9),
        "n_conserved": sum(1 for s in all_scores if s > 0.7),
        "n_variable": sum(1 for s in all_scores if s < 0.5),
    }

    if binding_site_residues:
        bs_scores = get_binding_site_conservation(conservation_result, binding_site_residues)

        if bs_scores:
            bs_values = list(bs_scores.values())
            result["binding_site"] = {
                "n_residues": len(bs_scores),
                "mean_conservation": sum(bs_values) / len(bs_values),
                "min_conservation": min(bs_values),
                "max_conservation": max(bs_values),
                "n_highly_conserved": sum(1 for s in bs_values if s > 0.9),
            }

    return result

def calculate_ecs_score(
    conservation_result: ConservationResult,
    binding_site_residues: list[int],
) -> dict[int, float]:
    bs_conservation = get_binding_site_conservation(
        conservation_result, binding_site_residues
    )

    if not bs_conservation:
        return {}

    values = list(bs_conservation.values())
    min_val = min(values)
    max_val = max(values)

    if max_val > min_val:
        return {
            pos: (val - min_val) / (max_val - min_val)
            for pos, val in bs_conservation.items()
        }

    return bs_conservation

if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from .msa_runner import _parse_fasta_alignment

    parser = argparse.ArgumentParser(description="Cálculo de conservação")
    parser.add_argument("alignment_file", type=Path, help="Arquivo de alinhamento FASTA")
    parser.add_argument(
        "-q", "--query",
        help="ID da sequência query para mapeamento",
    )
    parser.add_argument(
        "-r", "--residues",
        type=int,
        nargs="+",
        help="Resíduos do sítio de ligação",
    )

    args = parser.parse_args()

    sequences = _parse_fasta_alignment(args.alignment_file)
    alignment_length = sequences[0].length if sequences else 0

    msa_result = MSAResult(
        sequences=sequences,
        alignment_length=alignment_length,
        tool_used="file",
        success=True,
        output_path=args.alignment_file,
    )

    result = calculate_conservation(msa_result, args.query)

    analysis = analyze_conservation_distribution(result, args.residues)

    print(f"\nAnálise de conservação:")
    for key, value in analysis.items():
        if isinstance(value, dict):
            print(f"  {key}:")
            for k, v in value.items():
                if isinstance(v, float):
                    print(f"    {k}: {v:.3f}")
                else:
                    print(f"    {k}: {v}")
        elif isinstance(value, float):
            print(f"  {key}: {value:.3f}")
        else:
            print(f"  {key}: {value}")