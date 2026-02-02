from pathlib import Path
from dataclasses import dataclass
from typing import Final

PROJECT_ROOT: Final[Path] = Path(__file__).parent.parent
DATA_DIR: Final[Path] = PROJECT_ROOT / "dados"
RAW_DIR: Final[Path] = DATA_DIR / "brutos"
PROCESSED_DIR: Final[Path] = DATA_DIR / "processados"

STRUCTURES_DIR: Final[Path] = RAW_DIR / "estruturas"
SEQUENCES_DIR: Final[Path] = RAW_DIR / "sequencias"
FIELD_DATA_DIR: Final[Path] = RAW_DIR / "dados_campo"

ALIGNMENTS_DIR: Final[Path] = PROCESSED_DIR / "alinhamentos"
POCKETS_DIR: Final[Path] = PROCESSED_DIR / "cavidades"
SCORES_DIR: Final[Path] = PROCESSED_DIR / "scores"

RESULTS_DIR: Final[Path] = PROJECT_ROOT / "resultados"
FIGURES_DIR: Final[Path] = RESULTS_DIR / "figuras"
TABLES_DIR: Final[Path] = RESULTS_DIR / "tabelas"

EXAMPLE_DATA_DIR: Final[Path] = DATA_DIR / "exemplo"

FRAC_DATA_DIR: Final[Path] = RAW_DIR / "frac"
REGIONAL_DATA_DIR: Final[Path] = RAW_DIR / "regionais"
PROGRAMS_DATA_DIR: Final[Path] = RAW_DIR / "programas"

@dataclass(frozen=True)
class AlphaFoldConfig:
    base_url: str = "https://alphafold.ebi.ac.uk/files"
    model_version: str = "v4"
    timeout_seconds: int = 60
    max_retries: int = 3

@dataclass(frozen=True)
class UniProtConfig:
    base_url: str = "https://rest.uniprot.org/uniprotkb"
    timeout_seconds: int = 60
    max_results: int = 100
    fungi_taxon_id: int = 4751

@dataclass(frozen=True)
class SidraConfig:
    base_url: str = "https://apisidra.ibge.gov.br/values"
    table_pam: int = 1612
    soy_product_code: int = 2713
    yield_variable: int = 112
    timeout_seconds: int = 120

ALPHAFOLD_CONFIG: Final[AlphaFoldConfig] = AlphaFoldConfig()
UNIPROT_CONFIG: Final[UniProtConfig] = UniProtConfig()
SIDRA_CONFIG: Final[SidraConfig] = SidraConfig()

@dataclass(frozen=True)
class SDHTarget:
    organism: str = "Phakopsora pachyrhizi"
    enzyme: str = "Succinate dehydrogenase"
    subunits: tuple = ("SdhB", "SdhC", "SdhD")
    defensivo_class: str = "Carboxamidas (SDHI)"

SDH_TARGET: Final[SDHTarget] = SDHTarget()

@dataclass(frozen=True)
class AnalysisParams:
    binding_site_distance: float = 5.0

    plddt_very_high: float = 90.0
    plddt_confident: float = 70.0
    plddt_low: float = 50.0

    min_homologs: int = 50
    max_gap_percent: float = 50.0

    blosum_matrix: str = "BLOSUM62"

ANALYSIS_PARAMS: Final[AnalysisParams] = AnalysisParams()

LOG_FORMAT: Final[str] = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
LOG_DATE_FORMAT: Final[str] = "%Y-%m-%d %H:%M:%S"