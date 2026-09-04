"""Fast, dependency-light PDB loading for molecular context display."""

from collections import defaultdict
from pathlib import Path

import numpy as np

from vorpy.workbench.domain import AnalysisResult, Atom, Bond

DISPLAY_RADII = {
    "H": 0.23, "C": 0.36, "N": 0.34, "O": 0.33,
    "P": 0.42, "S": 0.40, "F": 0.32, "CL": 0.40,
    "BR": 0.43, "I": 0.46, "MG": 0.40, "ZN": 0.40,
}
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
    "P": 1.07, "S": 1.05, "F": 0.57, "CL": 1.02,
    "BR": 1.20, "I": 1.39, "MG": 1.30, "ZN": 1.22,
}


def _element_from_record(line: str, atom_name: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    name = "".join(character for character in atom_name if character.isalpha()).upper()
    if len(name) >= 2 and name[:2] in COVALENT_RADII:
        return name[:2]
    return name[:1] or "C"


def load_pdb(path: Path) -> AnalysisResult:
    """Load the first PDB model and construct explicit or inferred bonds."""
    atoms: list[Atom] = []
    serial_to_index: dict[int, int] = {}
    conect: set[tuple[int, int]] = set()
    in_first_model = True
    saw_model = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    in_first_model = False
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if not in_first_model:
                continue
            if record in {"ATOM", "HETATM"}:
                altloc = line[16:17]
                if altloc not in {"", " ", "A", "1"}:
                    continue
                try:
                    serial = int(line[6:11])
                    atom_name = line[12:16].strip()
                    position = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except (ValueError, IndexError):
                    continue
                element = _element_from_record(line, atom_name)
                index = len(atoms)
                atoms.append(Atom(
                    index=index,
                    serial=serial,
                    name=atom_name,
                    element=element,
                    position=position,
                    residue_name=line[17:20].strip(),
                    residue_sequence=line[22:26].strip(),
                    chain=line[21:22].strip(),
                    radius=DISPLAY_RADII.get(element, 0.36),
                ))
                serial_to_index[serial] = index
            elif record == "CONECT":
                try:
                    serials = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
                except ValueError:
                    continue
                for other in serials[1:]:
                    conect.add(tuple(sorted((serials[0], other))))

    if not atoms:
        raise ValueError(f"No ATOM or HETATM records were found in {path.name}")

    explicit = {
        tuple(sorted((serial_to_index[a], serial_to_index[b])))
        for a, b in conect if a in serial_to_index and b in serial_to_index
    }
    bond_pairs = explicit or _infer_bonds(atoms)
    return AnalysisResult(
        source=path,
        name=path.stem,
        atoms=atoms,
        bonds=[Bond(a, b) for a, b in sorted(bond_pairs)],
        complete_cells=0,
        surface_count=0,
    )


def _infer_bonds(atoms: list[Atom]) -> set[tuple[int, int]]:
    """Infer bonds with a spatial hash, avoiding an O(n²) atom scan."""
    cell_size = 2.4
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    positions = np.asarray([atom.position for atom in atoms], dtype=float)
    for index, position in enumerate(positions):
        cells[tuple(np.floor(position / cell_size).astype(int))].append(index)

    bonds: set[tuple[int, int]] = set()
    offsets = [(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)]
    for index, atom in enumerate(atoms):
        key = tuple(np.floor(positions[index] / cell_size).astype(int))
        radius_a = COVALENT_RADII.get(atom.element, 0.76)
        for offset in offsets:
            neighbor_key = tuple(key[axis] + offset[axis] for axis in range(3))
            for other in cells.get(neighbor_key, ()):
                if other <= index:
                    continue
                radius_b = COVALENT_RADII.get(atoms[other].element, 0.76)
                distance = float(np.linalg.norm(positions[index] - positions[other]))
                if 0.4 <= distance <= radius_a + radius_b + 0.45:
                    bonds.add((index, other))
    return bonds


