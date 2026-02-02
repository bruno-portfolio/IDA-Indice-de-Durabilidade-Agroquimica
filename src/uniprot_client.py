import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Iterator

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import UNIPROT_CONFIG, SEQUENCES_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class FastaSequence:
    accession: str
    entry_name: str
    description: str
    organism: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_fasta(self) -> str:
        header = f">{self.accession}|{self.entry_name} {self.description} OS={self.organism}"
        seq_lines = [self.sequence[i:i+60] for i in range(0, len(self.sequence), 60)]
        return header + "\n" + "\n".join(seq_lines)

@dataclass
class SearchResult:
    query: str
    sequences: list[FastaSequence] = field(default_factory=list)
    total_found: int = 0
    success: bool = True
    error_message: Optional[str] = None

class UniProtClient:

    def __init__(
        self,
        output_dir: Path = SEQUENCES_DIR,
        timeout: int = UNIPROT_CONFIG.timeout_seconds,
        max_results: int = UNIPROT_CONFIG.max_results,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.max_results = max_results
        self.session = self._create_session()

    def _create_session(self) -> requests.Session:
        session = requests.Session()
        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        return session

    def _parse_fasta(self, fasta_text: str) -> Iterator[FastaSequence]:
        if not fasta_text.strip():
            return

        header_pattern = re.compile(
            r">(?:sp|tr)\|([A-Z0-9]+)\|(\S+)\s+(.+?)\s+OS=([^=]+?)(?:\s+\w+=|$)"
        )

        entries = fasta_text.strip().split(">")[1:]

        for entry in entries:
            lines = entry.strip().split("\n")
            if not lines:
                continue

            header_line = ">" + lines[0]
            match = header_pattern.match(header_line)

            if match:
                accession, entry_name, description, organism = match.groups()
                sequence = "".join(lines[1:]).replace(" ", "").replace("\n", "")

                if sequence:
                    yield FastaSequence(
                        accession=accession.strip(),
                        entry_name=entry_name.strip(),
                        description=description.strip(),
                        organism=organism.strip(),
                        sequence=sequence.upper(),
                    )

    def search_homologs(
        self,
        protein_name: str,
        taxonomy_id: int = UNIPROT_CONFIG.fungi_taxon_id,
        reviewed_only: bool = True,
        limit: Optional[int] = None,
    ) -> SearchResult:
        if not protein_name or not protein_name.strip():
            return SearchResult(
                query=protein_name,
                success=False,
                error_message="Nome da proteína vazio",
            )

        limit = limit or self.max_results

        query_parts = [f'({protein_name})']
        query_parts.append(f'(taxonomy_id:{taxonomy_id})')
        if reviewed_only:
            query_parts.append('(reviewed:true)')

        query = " AND ".join(query_parts)

        url = f"{UNIPROT_CONFIG.base_url}/search"
        params = {
            "query": query,
            "format": "fasta",
            "size": limit,
        }

        logger.info(f"Buscando: {query}")

        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()

            sequences = list(self._parse_fasta(response.text))

            logger.info(f"Encontradas {len(sequences)} sequências")

            return SearchResult(
                query=query,
                sequences=sequences,
                total_found=len(sequences),
                success=True,
            )

        except requests.exceptions.HTTPError as e:
            logger.error(f"Erro HTTP: {e}")
            return SearchResult(
                query=query,
                success=False,
                error_message=f"HTTP {e.response.status_code}",
            )
        except requests.exceptions.Timeout:
            logger.error("Timeout na busca")
            return SearchResult(
                query=query,
                success=False,
                error_message="Timeout",
            )
        except requests.exceptions.ConnectionError:
            logger.error("Erro de conexão")
            return SearchResult(
                query=query,
                success=False,
                error_message="Erro de conexão",
            )

    def search_sdh_homologs(self, limit: Optional[int] = None) -> SearchResult:
        return self.search_homologs(
            protein_name="succinate dehydrogenase",
            taxonomy_id=UNIPROT_CONFIG.fungi_taxon_id,
            reviewed_only=True,
            limit=limit,
        )

    def save_fasta(
        self,
        sequences: list[FastaSequence],
        filename: str,
        output_dir: Optional[Path] = None,
    ) -> Path:
        output_dir = output_dir or self.output_dir
        output_path = output_dir / filename

        if not filename.endswith((".fasta", ".fa", ".faa")):
            output_path = output_dir / f"{filename}.fasta"

        fasta_content = "\n".join(seq.to_fasta() for seq in sequences)

        output_path.write_text(fasta_content, encoding="utf-8")
        logger.info(f"Salvo: {output_path} ({len(sequences)} sequências)")

        return output_path

    def get_sequence_by_id(self, uniprot_id: str) -> Optional[FastaSequence]:
        if not uniprot_id or not uniprot_id.strip():
            return None

        uniprot_id = uniprot_id.strip().upper()
        url = f"{UNIPROT_CONFIG.base_url}/{uniprot_id}.fasta"

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()

            sequences = list(self._parse_fasta(response.text))
            return sequences[0] if sequences else None

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                logger.warning(f"UniProt ID não encontrado: {uniprot_id}")
            else:
                logger.error(f"Erro HTTP: {e}")
            return None
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            logger.error(f"Erro de conexão: {e}")
            return None

def calculate_sequence_stats(sequences: list[FastaSequence]) -> dict:
    if not sequences:
        return {"count": 0}

    lengths = [seq.length for seq in sequences]
    organisms = set(seq.organism for seq in sequences)

    return {
        "count": len(sequences),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "avg_length": sum(lengths) / len(lengths),
        "unique_organisms": len(organisms),
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Busca de sequências homólogas no UniProt"
    )
    parser.add_argument(
        "-p", "--protein",
        default="succinate dehydrogenase",
        help="Nome da proteína para buscar",
    )
    parser.add_argument(
        "-t", "--taxonomy",
        type=int,
        default=UNIPROT_CONFIG.fungi_taxon_id,
        help="ID de taxonomia (default: 4751 = Fungi)",
    )
    parser.add_argument(
        "-n", "--limit",
        type=int,
        default=100,
        help="Número máximo de resultados",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=SEQUENCES_DIR,
        help="Diretório de saída",
    )

    args = parser.parse_args()

    client = UniProtClient(output_dir=args.output, max_results=args.limit)
    result = client.search_homologs(
        protein_name=args.protein,
        taxonomy_id=args.taxonomy,
        limit=args.limit,
    )

    if result.success and result.sequences:
        safe_name = args.protein.replace(" ", "_").lower()
        filename = f"{safe_name}_homologs.fasta"
        output_path = client.save_fasta(result.sequences, filename)

        stats = calculate_sequence_stats(result.sequences)
        print(f"\nEstatísticas:")
        for key, value in stats.items():
            print(f"  {key}: {value}")
    else:
        print(f"Erro: {result.error_message}")