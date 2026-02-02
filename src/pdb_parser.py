import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

from .config import STRUCTURES_DIR, LOG_FORMAT, LOG_DATE_FORMAT, ANALYSIS_PARAMS

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    alt_loc: str
    res_name: str
    chain_id: str
    res_seq: int
    x: float
    y: float
    z: float
    occupancy: float
    temp_factor: float
    element: str

    @property
    def coords(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_backbone(self) -> bool:
        return self.name.strip() in ("N", "CA", "C", "O")

    @property
    def is_sidechain(self) -> bool:
        return not self.is_backbone and self.name.strip() != "H"

@dataclass
class Residue:
    res_seq: int
    res_name: str
    chain_id: str
    atoms: list[Atom] = field(default_factory=list)

    AA_MAP: dict = field(default_factory=lambda: {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }, repr=False)

    @property
    def one_letter(self) -> str:
        return self.AA_MAP.get(self.res_name, "X")

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    def get_ca(self) -> Optional[Atom]:
        for atom in self.atoms:
            if atom.name.strip() == "CA":
                return atom
        return None

    def get_centroid(self) -> Optional[tuple[float, float, float]]:
        if not self.atoms:
            return None
        x = sum(a.x for a in self.atoms) / len(self.atoms)
        y = sum(a.y for a in self.atoms) / len(self.atoms)
        z = sum(a.z for a in self.atoms) / len(self.atoms)
        return (x, y, z)

@dataclass
class Chain:
    chain_id: str
    residues: dict[int, Residue] = field(default_factory=dict)

    @property
    def sequence(self) -> str:
        sorted_residues = sorted(self.residues.values(), key=lambda r: r.res_seq)
        return "".join(r.one_letter for r in sorted_residues)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    def get_residue(self, res_seq: int) -> Optional[Residue]:
        return self.residues.get(res_seq)

@dataclass
class Structure:
    name: str
    chains: dict[str, Chain] = field(default_factory=dict)
    header: dict = field(default_factory=dict)

    @property
    def n_chains(self) -> int:
        return len(self.chains)

    @property
    def n_residues(self) -> int:
        return sum(c.n_residues for c in self.chains.values())

    @property
    def n_atoms(self) -> int:
        return sum(
            r.n_atoms
            for c in self.chains.values()
            for r in c.residues.values()
        )

    def get_chain(self, chain_id: str) -> Optional[Chain]:
        return self.chains.get(chain_id)

    def iter_residues(self) -> Iterator[Residue]:
        for chain in self.chains.values():
            for residue in sorted(chain.residues.values(), key=lambda r: r.res_seq):
                yield residue

    def iter_atoms(self) -> Iterator[Atom]:
        for residue in self.iter_residues():
            yield from residue.atoms

    def get_sequence(self, chain_id: Optional[str] = None) -> str:
        if chain_id:
            chain = self.get_chain(chain_id)
            return chain.sequence if chain else ""
        return "".join(c.sequence for c in self.chains.values())

class PDBParseError(Exception):
    pass

class PDBParser:

    def __init__(self, strict: bool = False) -> None:
        self.strict = strict

    def _parse_atom_line(self, line: str) -> Optional[Atom]:
        try:
            record_type = line[0:6].strip()
            if record_type not in ("ATOM", "HETATM"):
                return None

            serial = int(line[6:11].strip())
            name = line[12:16]
            alt_loc = line[16:17].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22]
            res_seq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            try:
                occupancy = float(line[54:60].strip())
            except (ValueError, IndexError):
                occupancy = 1.0

            try:
                temp_factor = float(line[60:66].strip())
            except (ValueError, IndexError):
                temp_factor = 0.0

            try:
                element = line[76:78].strip()
            except IndexError:
                element = name.strip()[0]

            return Atom(
                serial=serial,
                name=name,
                alt_loc=alt_loc,
                res_name=res_name,
                chain_id=chain_id,
                res_seq=res_seq,
                x=x,
                y=y,
                z=z,
                occupancy=occupancy,
                temp_factor=temp_factor,
                element=element,
            )

        except (ValueError, IndexError) as e:
            if self.strict:
                raise PDBParseError(f"Erro ao parsear linha: {line}") from e
            return None

    def parse(self, pdb_path: Path) -> Structure:
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {pdb_path}")

        structure = Structure(name=pdb_path.stem)

        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("HEADER"):
                    structure.header["classification"] = line[10:50].strip()
                    structure.header["date"] = line[50:59].strip()
                    structure.header["id"] = line[62:66].strip()

                elif line.startswith("TITLE"):
                    title = structure.header.get("title", "")
                    structure.header["title"] = (title + " " + line[10:].strip()).strip()

                elif line.startswith(("ATOM", "HETATM")):
                    atom = self._parse_atom_line(line)
                    if atom is None:
                        continue

                    if atom.res_name in ("HOH", "WAT", "DOD"):
                        continue

                    if atom.chain_id not in structure.chains:
                        structure.chains[atom.chain_id] = Chain(chain_id=atom.chain_id)
                    chain = structure.chains[atom.chain_id]

                    if atom.res_seq not in chain.residues:
                        chain.residues[atom.res_seq] = Residue(
                            res_seq=atom.res_seq,
                            res_name=atom.res_name,
                            chain_id=atom.chain_id,
                        )
                    residue = chain.residues[atom.res_seq]

                    residue.atoms.append(atom)

        logger.info(
            f"Parsed {pdb_path.name}: {structure.n_chains} chains, "
            f"{structure.n_residues} residues, {structure.n_atoms} atoms"
        )

        return structure

    def parse_multiple(self, pdb_dir: Path) -> dict[str, Structure]:
        pdb_dir = Path(pdb_dir)
        structures = {}

        for pdb_file in pdb_dir.glob("*.pdb"):
            try:
                structure = self.parse(pdb_file)
                structures[pdb_file.stem] = structure
            except PDBParseError as e:
                logger.error(f"Erro ao parsear {pdb_file}: {e}")

        return structures

def calculate_distance(
    coord1: tuple[float, float, float],
    coord2: tuple[float, float, float],
) -> float:
    return (
        (coord1[0] - coord2[0]) ** 2 +
        (coord1[1] - coord2[1]) ** 2 +
        (coord1[2] - coord2[2]) ** 2
    ) ** 0.5

def find_nearby_residues(
    structure: Structure,
    center: tuple[float, float, float],
    distance: float = ANALYSIS_PARAMS.binding_site_distance,
) -> list[Residue]:
    nearby = []

    for residue in structure.iter_residues():
        centroid = residue.get_centroid()
        if centroid and calculate_distance(center, centroid) <= distance:
            nearby.append(residue)

    return nearby

def find_binding_site_residues(
    structure: Structure,
    ligand_residues: list[int],
    ligand_chain: str,
    distance: float = ANALYSIS_PARAMS.binding_site_distance,
) -> list[Residue]:
    chain = structure.get_chain(ligand_chain)
    if not chain:
        logger.warning(f"Cadeia {ligand_chain} não encontrada")
        return []

    ligand_atoms = []
    for res_seq in ligand_residues:
        residue = chain.get_residue(res_seq)
        if residue:
            ligand_atoms.extend(residue.atoms)

    if not ligand_atoms:
        logger.warning("Nenhum átomo de ligante encontrado")
        return []

    binding_site = set()

    for residue in structure.iter_residues():
        if residue.res_seq in ligand_residues and residue.chain_id == ligand_chain:
            continue

        for res_atom in residue.atoms:
            for lig_atom in ligand_atoms:
                dist = calculate_distance(res_atom.coords, lig_atom.coords)
                if dist <= distance:
                    binding_site.add((residue.chain_id, residue.res_seq))
                    break

    result = []
    for chain_id, res_seq in sorted(binding_site):
        chain = structure.get_chain(chain_id)
        if chain:
            residue = chain.get_residue(res_seq)
            if residue:
                result.append(residue)

    logger.info(f"Encontrados {len(result)} resíduos no sítio de ligação")
    return result

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Parse de arquivos PDB")
    parser.add_argument("pdb_file", type=Path, help="Arquivo PDB")
    parser.add_argument(
        "-s", "--sequence",
        action="store_true",
        help="Mostrar sequência",
    )

    args = parser.parse_args()

    pdb_parser = PDBParser()
    structure = pdb_parser.parse(args.pdb_file)

    print(f"Nome: {structure.name}")
    print(f"Cadeias: {structure.n_chains}")
    print(f"Resíduos: {structure.n_residues}")
    print(f"Átomos: {structure.n_atoms}")

    if args.sequence:
        for chain_id, chain in structure.chains.items():
            print(f"\nCadeia {chain_id} ({chain.n_residues} res):")
            print(chain.sequence)