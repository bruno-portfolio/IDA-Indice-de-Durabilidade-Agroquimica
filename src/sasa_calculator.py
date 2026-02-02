import logging
import math
from dataclasses import dataclass
from typing import Optional

from .pdb_parser import Structure, Residue, Atom, calculate_distance
from .config import LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

VDW_RADII: dict[str, float] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "H": 1.20,
    "P": 1.80,
    "FE": 1.47,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 1.97,
    "MN": 1.39,
}

MAX_SASA_BY_RESIDUE: dict[str, float] = {
    "ALA": 113.0, "ARG": 241.0, "ASN": 158.0, "ASP": 151.0, "CYS": 140.0,
    "GLN": 189.0, "GLU": 183.0, "GLY": 85.0, "HIS": 194.0, "ILE": 182.0,
    "LEU": 180.0, "LYS": 211.0, "MET": 204.0, "PHE": 218.0, "PRO": 143.0,
    "SER": 122.0, "THR": 146.0, "TRP": 259.0, "TYR": 229.0, "VAL": 160.0,
}

PROBE_RADIUS: float = 1.4

@dataclass
class SASAResult:
    res_seq: int
    res_name: str
    chain_id: str
    sasa: float
    max_sasa: float
    rsa: float
    burial: float

    @property
    def is_buried(self) -> bool:
        return self.rsa < 0.25

    @property
    def is_exposed(self) -> bool:
        return self.rsa > 0.5

def _get_vdw_radius(atom: Atom) -> float:
    element = atom.element.strip().upper()
    if not element:
        element = atom.name.strip()[0].upper()
    return VDW_RADII.get(element, 1.70)

def _generate_sphere_points(n_points: int = 92) -> list[tuple[float, float, float]]:
    points = []
    phi = math.pi * (3.0 - math.sqrt(5.0))

    for i in range(n_points):
        y = 1 - (i / (n_points - 1)) * 2
        radius = math.sqrt(1 - y * y)
        theta = phi * i

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        points.append((x, y, z))

    return points

_SPHERE_POINTS = _generate_sphere_points(92)

def calculate_atom_sasa(
    atom: Atom,
    all_atoms: list[Atom],
    probe_radius: float = PROBE_RADIUS,
    n_points: int = 92,
) -> float:
    vdw_radius = _get_vdw_radius(atom)
    test_radius = vdw_radius + probe_radius

    neighbor_cutoff = test_radius + max(VDW_RADII.values()) + probe_radius
    neighbors = []

    for other in all_atoms:
        if other.serial == atom.serial:
            continue
        dist = calculate_distance(atom.coords, other.coords)
        if dist < neighbor_cutoff:
            neighbors.append((other, _get_vdw_radius(other) + probe_radius))

    accessible_count = 0
    sphere_points = _SPHERE_POINTS[:n_points]

    for px, py, pz in sphere_points:
        test_x = atom.x + px * test_radius
        test_y = atom.y + py * test_radius
        test_z = atom.z + pz * test_radius

        is_accessible = True
        for neighbor, neighbor_radius in neighbors:
            dist = math.sqrt(
                (test_x - neighbor.x) ** 2 +
                (test_y - neighbor.y) ** 2 +
                (test_z - neighbor.z) ** 2
            )
            if dist < neighbor_radius:
                is_accessible = False
                break

        if is_accessible:
            accessible_count += 1

    fraction_accessible = accessible_count / n_points
    sphere_area = 4.0 * math.pi * test_radius ** 2

    return fraction_accessible * sphere_area

def calculate_residue_sasa(
    residue: Residue,
    all_atoms: list[Atom],
    probe_radius: float = PROBE_RADIUS,
) -> SASAResult:
    total_sasa = 0.0

    for atom in residue.atoms:
        if atom.element.strip().upper() == "H":
            continue
        total_sasa += calculate_atom_sasa(atom, all_atoms, probe_radius)

    max_sasa = MAX_SASA_BY_RESIDUE.get(residue.res_name, 200.0)
    rsa = min(total_sasa / max_sasa, 1.0) if max_sasa > 0 else 0.0

    return SASAResult(
        res_seq=residue.res_seq,
        res_name=residue.res_name,
        chain_id=residue.chain_id,
        sasa=total_sasa,
        max_sasa=max_sasa,
        rsa=rsa,
        burial=1.0 - rsa,
    )

def calculate_structure_sasa(
    structure: Structure,
    residue_subset: Optional[list[int]] = None,
) -> dict[int, SASAResult]:
    all_atoms = list(structure.iter_atoms())
    results = {}

    residues = list(structure.iter_residues())
    if residue_subset:
        residues = [r for r in residues if r.res_seq in residue_subset]

    total = len(residues)
    for i, residue in enumerate(residues):
        if (i + 1) % 50 == 0:
            logger.info(f"Calculando SASA: {i + 1}/{total} resíduos")

        result = calculate_residue_sasa(residue, all_atoms)
        results[residue.res_seq] = result

    logger.info(f"SASA calculada para {len(results)} resíduos")
    return results

def calculate_burial_score(
    sasa_results: dict[int, SASAResult],
    binding_site_residues: list[int],
) -> dict[int, float]:
    scores = {}

    for res_seq in binding_site_residues:
        if res_seq in sasa_results:
            scores[res_seq] = sasa_results[res_seq].burial

    if scores:
        max_burial = max(scores.values())
        min_burial = min(scores.values())

        if max_burial > min_burial:
            scores = {
                k: (v - min_burial) / (max_burial - min_burial)
                for k, v in scores.items()
            }

    return scores

def get_sasa_statistics(sasa_results: dict[int, SASAResult]) -> dict:
    if not sasa_results:
        return {"count": 0}

    sasas = [r.sasa for r in sasa_results.values()]
    rsas = [r.rsa for r in sasa_results.values()]
    buried = sum(1 for r in sasa_results.values() if r.is_buried)
    exposed = sum(1 for r in sasa_results.values() if r.is_exposed)

    return {
        "count": len(sasa_results),
        "total_sasa": sum(sasas),
        "avg_sasa": sum(sasas) / len(sasas),
        "avg_rsa": sum(rsas) / len(rsas),
        "n_buried": buried,
        "n_exposed": exposed,
        "pct_buried": buried / len(sasa_results) * 100,
    }

if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from .pdb_parser import PDBParser

    parser = argparse.ArgumentParser(description="Cálculo de SASA")
    parser.add_argument("pdb_file", type=Path, help="Arquivo PDB")
    parser.add_argument(
        "-r", "--residues",
        type=int,
        nargs="+",
        help="Resíduos específicos para calcular",
    )

    args = parser.parse_args()

    pdb_parser = PDBParser()
    structure = pdb_parser.parse(args.pdb_file)

    results = calculate_structure_sasa(structure, args.residues)
    stats = get_sasa_statistics(results)

    print(f"\nEstatísticas SASA:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.2f}")
        else:
            print(f"  {key}: {value}")

    if args.residues:
        print(f"\nDetalhes por resíduo:")
        for res_seq, result in sorted(results.items()):
            print(
                f"  {result.res_name}{res_seq}: "
                f"SASA={result.sasa:.1f}Ų, RSA={result.rsa:.2f}, "
                f"Burial={result.burial:.2f}"
            )