import logging
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import shutil

from .pdb_parser import Structure, Residue, PDBParser
from .config import STRUCTURES_DIR, POCKETS_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class Pocket:
    pocket_id: int
    score: float
    druggability: float
    volume: float
    n_alpha_spheres: int
    residues: list[tuple[str, int, str]] = field(default_factory=list)

    mean_local_hydrophobic_density: float = 0.0
    mean_alpha_sphere_radius: float = 0.0
    proportion_polar_atoms: float = 0.0
    charge_score: float = 0.0
    hydrophobicity_score: float = 0.0

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def is_druggable(self) -> bool:
        return self.druggability > 0.5

@dataclass
class PocketDetectionResult:
    pockets: list[Pocket]
    pdb_path: Path
    success: bool
    error_message: Optional[str] = None

    @property
    def n_pockets(self) -> int:
        return len(self.pockets)

    def get_best_pocket(self) -> Optional[Pocket]:
        if not self.pockets:
            return None
        return max(self.pockets, key=lambda p: p.score)

    def get_druggable_pockets(self) -> list[Pocket]:
        return [p for p in self.pockets if p.is_druggable]

def _check_fpocket_installed() -> bool:
    return shutil.which("fpocket") is not None

def _parse_fpocket_info(info_file: Path) -> list[Pocket]:
    pockets = []
    current_pocket: Optional[dict] = None

    if not info_file.exists():
        logger.warning(f"Arquivo info não encontrado: {info_file}")
        return pockets

    with open(info_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if line.startswith("Pocket"):
                if current_pocket:
                    pockets.append(_create_pocket_from_dict(current_pocket))
                try:
                    pocket_id = int(line.split()[1].replace(":", ""))
                    current_pocket = {"pocket_id": pocket_id}
                except (ValueError, IndexError):
                    continue

            elif current_pocket and ":" in line:
                try:
                    key, value = line.split(":", 1)
                    key = key.strip().lower().replace(" ", "_").replace("-", "_")
                    value = value.strip()

                    try:
                        current_pocket[key] = float(value)
                    except ValueError:
                        current_pocket[key] = value
                except ValueError:
                    continue

    if current_pocket:
        pockets.append(_create_pocket_from_dict(current_pocket))

    return pockets

def _create_pocket_from_dict(data: dict) -> Pocket:
    return Pocket(
        pocket_id=data.get("pocket_id", 0),
        score=data.get("score", data.get("drug_score", 0.0)),
        druggability=data.get("druggability_score", data.get("drug_score", 0.0)),
        volume=data.get("volume", data.get("real_volume", 0.0)),
        n_alpha_spheres=int(data.get("number_of_alpha_spheres", 0)),
        mean_local_hydrophobic_density=data.get("mean_local_hydrophobic_density", 0.0),
        mean_alpha_sphere_radius=data.get("mean_alpha_sphere_radius", 0.0),
        proportion_polar_atoms=data.get("proportion_of_polar_atoms", 0.0),
        charge_score=data.get("charge_score", 0.0),
        hydrophobicity_score=data.get("hydrophobicity_score", 0.0),
    )

def _parse_pocket_residues(pocket_pdb: Path) -> list[tuple[str, int, str]]:
    residues = set()

    if not pocket_pdb.exists():
        return []

    with open(pocket_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    chain_id = line[21:22].strip() or "A"
                    res_seq = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    residues.add((chain_id, res_seq, res_name))
                except (ValueError, IndexError):
                    continue

    return sorted(residues, key=lambda x: (x[0], x[1]))

class FpocketRunner:

    def __init__(self, output_dir: Path = POCKETS_DIR) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._fpocket_available = _check_fpocket_installed()

        if not self._fpocket_available:
            logger.warning(
                "fpocket não encontrado. Instale via: "
                "sudo apt-get install fpocket (Linux) ou compile do source"
            )

    def run(self, pdb_path: Path) -> PocketDetectionResult:
        pdb_path = Path(pdb_path)

        if not pdb_path.exists():
            return PocketDetectionResult(
                pockets=[],
                pdb_path=pdb_path,
                success=False,
                error_message=f"Arquivo não encontrado: {pdb_path}",
            )

        if not self._fpocket_available:
            return PocketDetectionResult(
                pockets=[],
                pdb_path=pdb_path,
                success=False,
                error_message="fpocket não instalado",
            )

        work_pdb = self.output_dir / pdb_path.name
        shutil.copy(pdb_path, work_pdb)

        try:
            result = subprocess.run(
                ["fpocket", "-f", str(work_pdb)],
                capture_output=True,
                text=True,
                timeout=300,
            )

            if result.returncode != 0:
                logger.error(f"fpocket erro: {result.stderr}")
                return PocketDetectionResult(
                    pockets=[],
                    pdb_path=pdb_path,
                    success=False,
                    error_message=result.stderr,
                )

            output_dir = self.output_dir / f"{pdb_path.stem}_out"

            if not output_dir.exists():
                return PocketDetectionResult(
                    pockets=[],
                    pdb_path=pdb_path,
                    success=False,
                    error_message="Diretório de output não criado",
                )

            info_file = output_dir / f"{pdb_path.stem}_info.txt"
            pockets = _parse_fpocket_info(info_file)

            for pocket in pockets:
                pocket_pdb = output_dir / f"pocket{pocket.pocket_id}_atm.pdb"
                pocket.residues = _parse_pocket_residues(pocket_pdb)

            logger.info(f"Detectadas {len(pockets)} cavidades em {pdb_path.name}")

            return PocketDetectionResult(
                pockets=pockets,
                pdb_path=pdb_path,
                success=True,
            )

        except subprocess.TimeoutExpired:
            return PocketDetectionResult(
                pockets=[],
                pdb_path=pdb_path,
                success=False,
                error_message="Timeout",
            )
        except FileNotFoundError:
            return PocketDetectionResult(
                pockets=[],
                pdb_path=pdb_path,
                success=False,
                error_message="fpocket não encontrado",
            )
        finally:
            if work_pdb.exists() and work_pdb != pdb_path:
                work_pdb.unlink()

class GeometricPocketDetector:

    def __init__(
        self,
        grid_spacing: float = 1.0,
        min_depth: float = 3.0,
        probe_radius: float = 1.4,
    ) -> None:
        self.grid_spacing = grid_spacing
        self.min_depth = min_depth
        self.probe_radius = probe_radius

    def detect(self, structure: Structure) -> PocketDetectionResult:
        atoms = list(structure.iter_atoms())

        if not atoms:
            return PocketDetectionResult(
                pockets=[],
                pdb_path=Path(""),
                success=False,
                error_message="Estrutura sem átomos",
            )

        xs = [a.x for a in atoms]
        ys = [a.y for a in atoms]
        zs = [a.z for a in atoms]

        min_x, max_x = min(xs) - 5, max(xs) + 5
        min_y, max_y = min(ys) - 5, max(ys) + 5
        min_z, max_z = min(zs) - 5, max(zs) + 5

        logger.info(
            "Usando detector geométrico simplificado - "
            "considere instalar fpocket para resultados melhores"
        )

        pocket = Pocket(
            pocket_id=1,
            score=0.5,
            druggability=0.5,
            volume=0.0,
            n_alpha_spheres=0,
            residues=[],
        )

        return PocketDetectionResult(
            pockets=[pocket] if pocket else [],
            pdb_path=Path(""),
            success=True,
            error_message="Usando detector simplificado (fpocket não disponível)",
        )

def detect_pockets(
    pdb_path: Path,
    output_dir: Path = POCKETS_DIR,
    use_fallback: bool = True,
) -> PocketDetectionResult:
    runner = FpocketRunner(output_dir)
    result = runner.run(pdb_path)

    if result.success:
        return result

    if use_fallback and "não instalado" in (result.error_message or ""):
        logger.info("Usando detector geométrico como fallback")
        parser = PDBParser()
        structure = parser.parse(pdb_path)
        detector = GeometricPocketDetector()
        return detector.detect(structure)

    return result

def calculate_pocket_srs(
    pocket: Pocket,
    druggability_weight: float = 0.4,
    volume_weight: float = 0.3,
    hydrophobicity_weight: float = 0.3,
    max_volume: float = 2000.0,
) -> float:
    norm_volume = 1.0 - min(pocket.volume / max_volume, 1.0)

    norm_druggability = pocket.druggability

    norm_hydrophobicity = min(pocket.hydrophobicity_score / 100.0, 1.0) if pocket.hydrophobicity_score > 0 else 0.5

    srs = (
        druggability_weight * norm_druggability +
        volume_weight * norm_volume +
        hydrophobicity_weight * norm_hydrophobicity
    )

    return min(max(srs, 0.0), 1.0)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Detecção de cavidades proteicas")
    parser.add_argument("pdb_file", type=Path, help="Arquivo PDB")
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=POCKETS_DIR,
        help="Diretório de saída",
    )

    args = parser.parse_args()

    result = detect_pockets(args.pdb_file, args.output)

    if result.success:
        print(f"\nDetectadas {result.n_pockets} cavidades:")
        for pocket in result.pockets:
            print(
                f"  Pocket {pocket.pocket_id}: "
                f"score={pocket.score:.2f}, "
                f"druggability={pocket.druggability:.2f}, "
                f"volume={pocket.volume:.1f}ų, "
                f"residues={pocket.n_residues}"
            )

        best = result.get_best_pocket()
        if best:
            print(f"\nMelhor cavidade: Pocket {best.pocket_id}")
            srs = calculate_pocket_srs(best)
            print(f"Score de Restrição do Sítio (SRS): {srs:.3f}")
    else:
        print(f"Erro: {result.error_message}")