import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import ALPHAFOLD_CONFIG, STRUCTURES_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class AlphaFoldStructure:
    uniprot_id: str
    pdb_path: Optional[Path]
    confidence_path: Optional[Path]
    success: bool
    error_message: Optional[str] = None

class AlphaFoldDownloadError(Exception):
    pass

class AlphaFoldClient:

    def __init__(
        self,
        output_dir: Path = STRUCTURES_DIR,
        timeout: int = ALPHAFOLD_CONFIG.timeout_seconds,
        max_retries: int = ALPHAFOLD_CONFIG.max_retries,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.session = self._create_session(max_retries)

    def _create_session(self, max_retries: int) -> requests.Session:
        session = requests.Session()
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        session.mount("http://", adapter)
        return session

    def _build_urls(self, uniprot_id: str) -> tuple[str, str]:
        base = ALPHAFOLD_CONFIG.base_url
        version = ALPHAFOLD_CONFIG.model_version
        pdb_url = f"{base}/AF-{uniprot_id}-F1-model_{version}.pdb"
        json_url = f"{base}/AF-{uniprot_id}-F1-confidence_{version}.json"
        return pdb_url, json_url

    def _download_file(self, url: str, output_path: Path) -> bool:
        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            output_path.write_bytes(response.content)
            logger.info(f"Download concluído: {output_path.name}")
            return True
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                logger.warning(f"Estrutura não encontrada: {url}")
            else:
                logger.error(f"Erro HTTP {e.response.status_code}: {url}")
            return False
        except requests.exceptions.Timeout:
            logger.error(f"Timeout ao baixar: {url}")
            return False
        except requests.exceptions.ConnectionError:
            logger.error(f"Erro de conexão: {url}")
            return False

    def download_structure(self, uniprot_id: str) -> AlphaFoldStructure:
        if not uniprot_id or not uniprot_id.strip():
            return AlphaFoldStructure(
                uniprot_id=uniprot_id,
                pdb_path=None,
                confidence_path=None,
                success=False,
                error_message="UniProt ID vazio ou inválido",
            )

        uniprot_id = uniprot_id.strip().upper()
        pdb_url, json_url = self._build_urls(uniprot_id)

        pdb_path = self.output_dir / f"{uniprot_id}.pdb"
        json_path = self.output_dir / f"{uniprot_id}_confidence.json"

        logger.info(f"Iniciando download para: {uniprot_id}")

        pdb_success = self._download_file(pdb_url, pdb_path)
        json_success = self._download_file(json_url, json_path)

        if pdb_success and json_success:
            return AlphaFoldStructure(
                uniprot_id=uniprot_id,
                pdb_path=pdb_path,
                confidence_path=json_path,
                success=True,
            )
        elif pdb_success:
            return AlphaFoldStructure(
                uniprot_id=uniprot_id,
                pdb_path=pdb_path,
                confidence_path=None,
                success=True,
                error_message="Arquivo de confiança não disponível",
            )
        else:
            return AlphaFoldStructure(
                uniprot_id=uniprot_id,
                pdb_path=None,
                confidence_path=None,
                success=False,
                error_message="Estrutura não encontrada no AlphaFold DB",
            )

    def download_multiple(self, uniprot_ids: list[str]) -> list[AlphaFoldStructure]:
        if not uniprot_ids:
            logger.warning("Lista de UniProt IDs vazia")
            return []

        results = []
        for uid in uniprot_ids:
            result = self.download_structure(uid)
            results.append(result)

        success_count = sum(1 for r in results if r.success)
        logger.info(f"Download concluído: {success_count}/{len(uniprot_ids)} estruturas")

        return results

def load_confidence_scores(json_path: Path) -> Optional[dict]:
    if not json_path.exists():
        logger.error(f"Arquivo não encontrado: {json_path}")
        return None

    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)

        plddt = data.get("confidenceScore") or data.get("plddt", [])

        return {
            "plddt": plddt,
            "pae": data.get("pae"),
        }
    except json.JSONDecodeError as e:
        logger.error(f"Erro ao decodificar JSON: {e}")
        return None
    except OSError as e:
        logger.error(f"Erro ao ler arquivo: {e}")
        return None

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Download de estruturas do AlphaFold Database"
    )
    parser.add_argument(
        "uniprot_ids",
        nargs="+",
        help="UniProt IDs para download (ex: P12345 Q67890)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=STRUCTURES_DIR,
        help="Diretório de saída",
    )

    args = parser.parse_args()

    client = AlphaFoldClient(output_dir=args.output)
    results = client.download_multiple(args.uniprot_ids)

    for r in results:
        status = "OK" if r.success else "FALHA"
        print(f"[{status}] {r.uniprot_id}: {r.pdb_path or r.error_message}")