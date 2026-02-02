"""Fixtures compartilhadas para testes."""

import pytest
from pathlib import Path
import sys
import os

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "dados"
RESULTS_DIR = PROJECT_ROOT / "resultados"
EXAMPLE_DIR = DATA_DIR / "exemplo"
STRUCTURES_DIR = DATA_DIR / "brutos" / "estruturas"
SEQUENCES_DIR = DATA_DIR / "brutos" / "sequencias"


@pytest.fixture
def project_root() -> Path:
    """Retorna o diretório raiz do projeto."""
    return PROJECT_ROOT


@pytest.fixture
def example_pdb() -> Path:
    """Retorna caminho do PDB de exemplo."""
    return STRUCTURES_DIR / "2FBW_chainB.pdb"


@pytest.fixture
def example_fasta() -> Path:
    """Retorna caminho do FASTA de exemplo."""
    return SEQUENCES_DIR / "sdh_b_fungi_homologs.fasta"


@pytest.fixture
def example_confidence() -> Path:
    """Retorna caminho do JSON de confiança."""
    return STRUCTURES_DIR / "2FBW_chainB_confidence.json"


@pytest.fixture
def sdh_binding_site() -> list[int]:
    """Retorna resíduos do sítio de ligação do SDH."""
    return [169, 170, 172, 173, 216, 218]


@pytest.fixture
def cyp51_binding_site() -> list[int]:
    """Retorna resíduos do sítio de ligação do CYP51."""
    return [118, 121, 307, 376, 379, 508, 509]


@pytest.fixture
def cytb_binding_site() -> list[int]:
    """Retorna resíduos do sítio de ligação do Cyt b."""
    return [53, 180, 181, 182, 183, 184, 185]
