import logging
import subprocess
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterator

from .config import SEQUENCES_DIR, ALIGNMENTS_DIR, ANALYSIS_PARAMS, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class AlignedSequence:
    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def ungapped_length(self) -> int:
        return len(self.sequence.replace("-", ""))

    @property
    def gap_count(self) -> int:
        return self.sequence.count("-")

    @property
    def gap_percent(self) -> float:
        return self.gap_count / self.length * 100 if self.length > 0 else 0.0

@dataclass
class MSAResult:
    sequences: list[AlignedSequence]
    alignment_length: int
    tool_used: str
    success: bool
    output_path: Optional[Path] = None
    error_message: Optional[str] = None

    @property
    def n_sequences(self) -> int:
        return len(self.sequences)

    def get_column(self, position: int) -> list[str]:
        if not 0 <= position < self.alignment_length:
            return []
        return [seq.sequence[position] for seq in self.sequences]

    def iter_columns(self) -> Iterator[tuple[int, list[str]]]:
        for pos in range(self.alignment_length):
            yield pos, self.get_column(pos)

def _check_tool_installed(tool: str) -> bool:
    return shutil.which(tool) is not None

def _parse_fasta_alignment(fasta_path: Path) -> list[AlignedSequence]:
    sequences = []
    current_id = ""
    current_desc = ""
    current_seq = []

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences.append(AlignedSequence(
                        id=current_id,
                        description=current_desc,
                        sequence="".join(current_seq),
                    ))

                header = line[1:]
                parts = header.split(None, 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line)

    if current_id:
        sequences.append(AlignedSequence(
            id=current_id,
            description=current_desc,
            sequence="".join(current_seq),
        ))

    return sequences

class MAFFTRunner:

    def __init__(self, output_dir: Path = ALIGNMENTS_DIR) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._available = _check_tool_installed("mafft")

        if not self._available:
            logger.warning(
                "MAFFT não encontrado. Instale via: "
                "sudo apt-get install mafft (Linux) ou brew install mafft (macOS)"
            )

    @property
    def is_available(self) -> bool:
        return self._available

    def run(
        self,
        input_fasta: Path,
        output_name: Optional[str] = None,
        strategy: str = "auto",
        max_iterations: int = 1000,
    ) -> MSAResult:
        input_fasta = Path(input_fasta)

        if not input_fasta.exists():
            return MSAResult(
                sequences=[],
                alignment_length=0,
                tool_used="mafft",
                success=False,
                error_message=f"Arquivo não encontrado: {input_fasta}",
            )

        if not self._available:
            return MSAResult(
                sequences=[],
                alignment_length=0,
                tool_used="mafft",
                success=False,
                error_message="MAFFT não instalado",
            )

        output_name = output_name or f"{input_fasta.stem}_aligned.fasta"
        output_path = self.output_dir / output_name

        cmd = ["mafft"]

        if strategy == "auto":
            cmd.append("--auto")
        elif strategy == "accurate":
            cmd.extend(["--maxiterate", str(max_iterations), "--localpair"])
        elif strategy == "fast":
            cmd.append("--retree 1")

        cmd.append(str(input_fasta))

        logger.info(f"Executando MAFFT ({strategy}): {input_fasta.name}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800,
            )

            if result.returncode != 0:
                logger.error(f"MAFFT erro: {result.stderr}")
                return MSAResult(
                    sequences=[],
                    alignment_length=0,
                    tool_used="mafft",
                    success=False,
                    error_message=result.stderr[:500],
                )

            output_path.write_text(result.stdout, encoding="utf-8")

            sequences = _parse_fasta_alignment(output_path)

            if not sequences:
                return MSAResult(
                    sequences=[],
                    alignment_length=0,
                    tool_used="mafft",
                    success=False,
                    error_message="Nenhuma sequência no alinhamento",
                )

            alignment_length = sequences[0].length if sequences else 0

            logger.info(
                f"Alinhamento concluído: {len(sequences)} seqs, "
                f"{alignment_length} posições"
            )

            return MSAResult(
                sequences=sequences,
                alignment_length=alignment_length,
                tool_used="mafft",
                success=True,
                output_path=output_path,
            )

        except subprocess.TimeoutExpired:
            return MSAResult(
                sequences=[],
                alignment_length=0,
                tool_used="mafft",
                success=False,
                error_message="Timeout",
            )

class ClustalOmegaRunner:

    def __init__(self, output_dir: Path = ALIGNMENTS_DIR) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._available = _check_tool_installed("clustalo")

    @property
    def is_available(self) -> bool:
        return self._available

    def run(
        self,
        input_fasta: Path,
        output_name: Optional[str] = None,
        n_threads: int = 4,
    ) -> MSAResult:
        input_fasta = Path(input_fasta)

        if not input_fasta.exists():
            return MSAResult(
                sequences=[],
                alignment_length=0,
                tool_used="clustalo",
                success=False,
                error_message=f"Arquivo não encontrado: {input_fasta}",
            )

        if not self._available:
            return MSAResult(
                sequences=[],
                alignment_length=0,
                tool_used="clustalo",
                success=False,
                error_message="Clustal Omega não instalado",
            )

        output_name = output_name or f"{input_fasta.stem}_aligned.fasta"
        output_path = self.output_dir / output_name

        cmd = [
            "clustalo",
            "-i", str(input_fasta),
            "-o", str(output_path),
            "--outfmt=fasta",
            f"--threads={n_threads}",
            "--force",
        ]

        logger.info(f"Executando Clustal Omega: {input_fasta.name}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800,
            )

            if result.returncode != 0:
                return MSAResult(
                    sequences=[],
                    alignment_length=0,
                    tool_used="clustalo",
                    success=False,
                    error_message=result.stderr[:500],
                )

            sequences = _parse_fasta_alignment(output_path)
            alignment_length = sequences[0].length if sequences else 0

            return MSAResult(
                sequences=sequences,
                alignment_length=alignment_length,
                tool_used="clustalo",
                success=True,
                output_path=output_path,
            )

        except subprocess.TimeoutExpired:
            return MSAResult(
                sequences=[],
                alignment_length=0,
                tool_used="clustalo",
                success=False,
                error_message="Timeout",
            )

def run_msa(
    input_fasta: Path,
    output_dir: Path = ALIGNMENTS_DIR,
    preferred_tool: str = "mafft",
) -> MSAResult:
    if preferred_tool == "mafft":
        runner = MAFFTRunner(output_dir)
        if runner.is_available:
            return runner.run(input_fasta)
        logger.info("MAFFT não disponível, tentando Clustal Omega")

    runner_clustal = ClustalOmegaRunner(output_dir)
    if runner_clustal.is_available:
        return runner_clustal.run(input_fasta)

    runner_mafft = MAFFTRunner(output_dir)
    if runner_mafft.is_available:
        return runner_mafft.run(input_fasta)

    return MSAResult(
        sequences=[],
        alignment_length=0,
        tool_used="none",
        success=False,
        error_message="Nenhuma ferramenta de MSA disponível (MAFFT ou Clustal Omega)",
    )

def validate_alignment(
    msa_result: MSAResult,
    min_sequences: int = ANALYSIS_PARAMS.min_homologs,
    max_gap_percent: float = ANALYSIS_PARAMS.max_gap_percent,
) -> dict:
    if not msa_result.success:
        return {"valid": False, "error": msa_result.error_message}

    issues = []

    if msa_result.n_sequences < min_sequences:
        issues.append(
            f"Poucas sequências: {msa_result.n_sequences} < {min_sequences}"
        )

    high_gap_seqs = []
    for seq in msa_result.sequences:
        if seq.gap_percent > max_gap_percent:
            high_gap_seqs.append(f"{seq.id} ({seq.gap_percent:.1f}%)")

    if high_gap_seqs:
        issues.append(f"{len(high_gap_seqs)} sequências com gaps excessivos")

    gap_percents = [s.gap_percent for s in msa_result.sequences]

    return {
        "valid": len(issues) == 0,
        "issues": issues,
        "n_sequences": msa_result.n_sequences,
        "alignment_length": msa_result.alignment_length,
        "mean_gap_percent": sum(gap_percents) / len(gap_percents) if gap_percents else 0,
        "max_gap_percent": max(gap_percents) if gap_percents else 0,
        "n_high_gap_seqs": len(high_gap_seqs),
    }

def map_alignment_to_structure(
    msa_result: MSAResult,
    query_id: str,
) -> dict[int, int]:
    query_seq = None
    for seq in msa_result.sequences:
        if query_id in seq.id:
            query_seq = seq
            break

    if not query_seq:
        logger.warning(f"Query {query_id} não encontrado no alinhamento")
        return {}

    mapping = {}
    struct_pos = 0

    for aln_pos, aa in enumerate(query_seq.sequence):
        if aa != "-":
            struct_pos += 1
            mapping[aln_pos] = struct_pos

    return mapping

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Execução de MSA")
    parser.add_argument("fasta_file", type=Path, help="Arquivo FASTA de entrada")
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=ALIGNMENTS_DIR,
        help="Diretório de saída",
    )
    parser.add_argument(
        "-t", "--tool",
        choices=["mafft", "clustalo"],
        default="mafft",
        help="Ferramenta de MSA",
    )

    args = parser.parse_args()

    result = run_msa(args.fasta_file, args.output, args.tool)

    if result.success:
        validation = validate_alignment(result)
        print(f"\nAlinhamento: {result.output_path}")
        print(f"Ferramenta: {result.tool_used}")
        print(f"Sequências: {result.n_sequences}")
        print(f"Comprimento: {result.alignment_length}")
        print(f"Válido: {validation['valid']}")

        if validation['issues']:
            print("Problemas:")
            for issue in validation['issues']:
                print(f"  - {issue}")
    else:
        print(f"Erro: {result.error_message}")