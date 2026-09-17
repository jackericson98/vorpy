"""Fast, dependency-light PDB loading for molecular context display."""

from pathlib import Path
from dataclasses import replace

from vorpy.workbench.chemistry import ION_RESIDUES, SOLVENT_RESIDUES

import numpy as np
from scipy.spatial import cKDTree

from vorpy.workbench.domain import AnalysisResult, Atom, Bond
from vorpy.workbench.services.trajectory import index_pdb_frames, read_frame_lines
from vorpy.src.chemistry import element_radii, my_masses

# Use the same default radii as molecular atoms created by the solver.
# Keep the public name for the atomic-radii editor and existing callers.
DISPLAY_RADII = dict(element_radii)
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


def _pdb_charge(text: str) -> float | None:
    text = text.strip()
    if not text:
        return None
    try:
        if text[-1] in "+-":
            return float(text[-1] + text[:-1])
        return float(text)
    except ValueError:
        return None


def load_pdb(path: Path, progress=None, preview=None, *, frame_index=1, frame_ranges=None) -> AnalysisResult:
    """Load one indexed PDB frame and construct explicit or inferred bonds."""
    atoms: list[Atom] = []
    records: list[Atom | str] = []
    serial_to_index: dict[int, int] = {}
    conect: set[tuple[int, int]] = set()
    report = progress or (lambda label, value: None)
    ranges = frame_ranges or index_pdb_frames(path, report)
    if not ranges:
        raise ValueError(f"No ATOM or HETATM records were found in {path.name}")
    if not 1 <= frame_index <= len(ranges):
        raise ValueError(f'Frame {frame_index} is outside 1–{len(ranges)}')
    metadata = dict(frame_index=frame_index, frame_count=len(ranges), frame_ranges=ranges)
    for line in read_frame_lines(path, ranges[frame_index - 1], report):
        record = line[:6].strip().upper()
        if record in {"ATOM", "HETATM"}:
            altloc = line[16:17]
            if altloc not in {"", " ", "A", "1"}:
                continue
            if preview is not None and line[17:20].strip().upper() in SOLVENT_RESIDUES:
                records.append(line)
                continue
            atom = _atom_from_line(line, len(atoms))
            if atom is not None:
                atoms.append(atom)
                records.append(atom)
        elif record == "CONECT":
            try:
                serials = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
            except ValueError:
                continue
            for other in serials[1:]:
                conect.add(tuple(sorted((serials[0], other))))

    if preview is not None and atoms and len(records) > len(atoms):
        # A separate list owns the compact preview indices. The full result
        # below restores file order, including interleaved and wrapped solvent.
        primary_serials = {atom.serial: atom.index for atom in atoms}
        pairs = {tuple(sorted((primary_serials[a], primary_serials[b]))) for a, b in conect
                 if a in primary_serials and b in primary_serials}
        pairs = pairs or _infer_bonds(atoms)
        preview(AnalysisResult(source=path, name=path.stem, atoms=atoms,
                               bonds=[Bond(a, b) for a, b in sorted(pairs)], **metadata))

    if preview is not None:
        atoms = []
        count = len(records)
        for row, record in enumerate(records):
            if row % 4096 == 0:
                report(f"Loading solvent · {row:,}/{count:,} atom records", 70 + row * 2 // max(count, 1))
            atom = (_atom_from_line(record, len(atoms)) if isinstance(record, str)
                    else record if record.index == len(atoms) else replace(record, index=len(atoms)))
            if atom is not None:
                atoms.append(atom)
    if not atoms:
        raise ValueError(f"No ATOM or HETATM records were found in {path.name}")
    serial_to_index = {atom.serial: atom.index for atom in atoms}

    report("Resolving atom connectivity…", 72)
    explicit = {
        tuple(sorted((serial_to_index[a], serial_to_index[b])))
        for a, b in conect if a in serial_to_index and b in serial_to_index
    }
    report("Finding bonds between nearby atoms…", 75)
    bond_pairs = explicit or _infer_bonds(atoms)
    report(f"Preparing {len(atoms):,} atoms and {len(bond_pairs):,} bonds…", 95)
    result = AnalysisResult(
        source=path,
        name=path.stem,
        atoms=atoms,
        bonds=[Bond(a, b) for a, b in sorted(bond_pairs)],
        complete_cells=0,
        surface_count=0,
        **metadata,
    )
    report("Structure read", 100)
    return result


def _atom_from_line(line, index):
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        position = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except (ValueError, IndexError):
        return None
    element = _element_from_record(line, name)
    return Atom(index=index, serial=serial, name=name, element=element, position=position,
                residue_name=line[17:20].strip(), residue_sequence=line[22:26].strip(),
                chain=line[21:22].strip(), radius=DISPLAY_RADII.get(element, DISPLAY_RADII["C"]),
                mass=my_masses.get(element.lower(), 1.0), charge=_pdb_charge(line[78:80]))


def _infer_bonds(atoms: list[Atom]) -> set[tuple[int, int]]:
    """Find nearby pairs in compiled code, then apply covalent cutoffs."""
    if len(atoms) < 2:
        return set()
    positions = np.asarray([atom.position for atom in atoms], dtype=float)
    radii = np.asarray([COVALENT_RADII.get(atom.element, 0.76) for atom in atoms])
    pairs = cKDTree(positions).query_pairs(2 * radii.max() + 0.45, output_type="ndarray")
    if not len(pairs):
        return set()
    delta = positions[pairs[:, 0]] - positions[pairs[:, 1]]
    distance_squared = np.einsum("ij,ij->i", delta, delta)
    cutoffs = radii[pairs[:, 0]] + radii[pairs[:, 1]] + 0.45
    valid = (distance_squared >= 0.4 ** 2) & (distance_squared <= cutoffs ** 2)
    # Exclude coordination contacts; retain bonds within polyatomic ions.
    ionic = np.asarray([atom.residue_name.strip().upper() in ION_RESIDUES for atom in atoms])
    residue_keys = [(atom.chain, atom.residue_sequence, atom.residue_name.strip().upper())
                    for atom in atoms]
    for row in np.flatnonzero(valid & (ionic[pairs[:, 0]] | ionic[pairs[:, 1]])):
        a, b = pairs[row]
        valid[row] = residue_keys[a] == residue_keys[b]
    return {tuple(map(int, pair)) for pair in pairs[valid]}
