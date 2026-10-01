"""
cazals_pdb.py

Adapter between a conventional PDB structure and the Cazals/Chothia
power-alpha-complex input used by VorPy validation analyses.

Repository location:
    vorpy/src/analyze/Part_II_Molecular_Analysis/SuplimentalInformation/cazals_pdb.py

The adapter deliberately does not modify VorPy's core PDB reader. It uses
the same fixed-width PDB fields and converts selected atoms into:

    points          (N, 3) Cartesian coordinates [A]
    base_radii      (N,)   Cazals/Chothia base group radii [A]
    expanded_radii  (N,)   base radius + probe radius [A]
    weights         (N,)   expanded_radii**2 [A^2]
    group_a         zero-based indices into the returned arrays
    group_b         zero-based indices into the returned arrays

For the Cazals construction:
    R_i = r_i + 1.4 A
    w_i = R_i^2

The returned indices are NOT PDB serial numbers. Original PDB serials and
metadata are retained in AtomRecord objects and exports.

Standard use:
    prepared = prepare_cazals_pdb(
        "vorpy/data/2kai.pdb",
        group_a_chains={"A"},
        group_b_chains={"B"},
        include_hydrogens=False,
        strict_heavy_fallbacks=True,
    )
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from vorpy.src.chemistry.Cazals import (
    OTHER_RADIUS,
    PROBE_RADIUS,
    base_radius,
    special_radii,
)


@dataclass(frozen=True)
class AtomRecord:
    index: int
    pdb_serial: int
    record_type: str
    atom_name: str
    altloc: str
    residue_name: str
    chain: str
    residue_number: int
    insertion_code: str
    element: str
    x: float
    y: float
    z: float
    radius_class: str
    explicit_assignment: bool
    base_radius: float
    expanded_radius: float
    weight: float
    group: str


@dataclass
class PreparedCazalsPDB:
    source: Path
    points: np.ndarray
    base_radii: np.ndarray
    expanded_radii: np.ndarray
    weights: np.ndarray
    group_a: np.ndarray
    group_b: np.ndarray
    atoms: list[AtomRecord]
    probe_radius: float
    skipped_altlocs: int
    skipped_hydrogens: int
    skipped_unselected_chains: int
    fallback_atoms: list[AtomRecord]

    @property
    def radii(self) -> np.ndarray:
        """Alias expected by the original array-based GUDHI script."""
        return self.expanded_radii


def _infer_element(atom_name: str, element_field: str) -> str:
    element = element_field.strip().upper()
    if element:
        return element

    name = atom_name.strip().upper().lstrip("0123456789")
    if not name:
        return ""

    # Correct for standard protein names such as CA/CB/CG: carbon, not Ca.
    if name[0] in {"C", "N", "O", "S", "H", "P"}:
        return name[0]
    return name[:2]


def _parse_pdb_line(line: str) -> dict:
    """Parse the fields used by this adapter using VorPy/PDB fixed columns."""
    atom_name = line[12:16].strip()
    element = _infer_element(
        atom_name,
        line[76:78] if len(line) >= 78 else "",
    )
    return {
        "record_type": line[0:6].strip().upper(),
        "pdb_serial": int(line[6:11].strip()),
        "atom_name": atom_name.upper(),
        "altloc": line[16:17].strip().upper(),
        "residue_name": line[17:20].strip().upper(),
        "chain": line[21:22].strip(),
        "residue_number": int(line[22:26].strip()),
        "insertion_code": line[26:27].strip(),
        "x": float(line[30:38].strip()),
        "y": float(line[38:46].strip()),
        "z": float(line[46:54].strip()),
        "element": element,
    }


def _is_hydrogen(atom: dict) -> bool:
    return (
        atom["element"].upper() == "H"
        or atom["atom_name"].upper().lstrip("0123456789").startswith("H")
    )


def _classification(residue: str, atom_name: str, element: str) -> tuple[str, bool]:
    """
    Return a human-readable class and whether the assignment was explicit.

    This mirrors Cazals.py without silently pretending an element fallback is
    a chemically specific C/N classification.
    """
    residue = residue.upper()
    atom_name = atom_name.upper()
    element = element.upper()

    explicit = residue in special_radii and atom_name in special_radii[residue]
    if not explicit:
        if element == "O":
            return "oxygen_element_fallback", False
        if element == "S":
            return "sulfur_element_fallback", False
        if element == "C":
            return "ambiguous_carbon_fallback", False
        if element == "N":
            return "ambiguous_nitrogen_fallback", False
        return "other_fallback", False

    r = special_radii[residue][atom_name]
    if element == "C":
        if np.isclose(r, 1.87):
            return "aliphatic_C", True
        if np.isclose(r, 1.76):
            return "trigonal_C", True
    elif element == "N":
        if np.isclose(r, 1.65):
            return "neutral_N", True
        if np.isclose(r, 1.50):
            return "charged_N", True
    elif element == "O":
        return "oxygen", True
    elif element == "S":
        return "sulfur", True

    return "explicit_residue_atom", True


def _choose_altlocs(raw_atoms: list[dict]) -> tuple[list[dict], int]:
    """
    Keep one conformer per PDB atom identity.

    Preference:
        blank altloc > A > first encountered other altloc.

    The identity includes chain/residue/insertion-code/atom-name.
    """
    selected: dict[tuple, dict] = {}
    priority = {"": 0, "A": 1}

    for atom in raw_atoms:
        key = (
            atom["chain"],
            atom["residue_number"],
            atom["insertion_code"],
            atom["residue_name"],
            atom["atom_name"],
        )
        old = selected.get(key)
        if old is None:
            selected[key] = atom
            continue

        old_p = priority.get(old["altloc"], 2)
        new_p = priority.get(atom["altloc"], 2)
        if new_p < old_p:
            selected[key] = atom

    chosen = list(selected.values())
    chosen.sort(key=lambda a: a["pdb_serial"])
    return chosen, len(raw_atoms) - len(chosen)


def prepare_cazals_pdb(
    pdb_path: str | Path,
    group_a_chains: Iterable[str],
    group_b_chains: Iterable[str],
    *,
    probe_radius: float = PROBE_RADIUS,
    include_hydrogens: bool = False,
    include_hetatm: bool = False,
    strict_heavy_fallbacks: bool = False,
    verbose: bool = True,
) -> PreparedCazalsPDB:
    """
    Convert two PDB chain groups into Cazals power-complex arrays.

    Only atoms belonging to group_a_chains or group_b_chains are returned.
    Groups must be disjoint.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.is_file():
        raise FileNotFoundError(pdb_path)

    a_chains = {str(x).strip() for x in group_a_chains}
    b_chains = {str(x).strip() for x in group_b_chains}
    if not a_chains or not b_chains:
        raise ValueError("Both chain groups must be non-empty.")
    overlap = a_chains & b_chains
    if overlap:
        raise ValueError(f"Chain groups overlap: {sorted(overlap)}")

    raw = []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[0:6].strip().upper()
            if record == "ATOM" or (include_hetatm and record == "HETATM"):
                try:
                    raw.append(_parse_pdb_line(line))
                except (ValueError, IndexError) as exc:
                    raise ValueError(
                        f"Could not parse PDB atom line in {pdb_path}: {line.rstrip()}"
                    ) from exc

    raw, skipped_altlocs = _choose_altlocs(raw)

    records: list[AtomRecord] = []
    skipped_h = 0
    skipped_chains = 0

    for atom in raw:
        chain = atom["chain"]
        if chain not in a_chains and chain not in b_chains:
            skipped_chains += 1
            continue
        if not include_hydrogens and _is_hydrogen(atom):
            skipped_h += 1
            continue

        group = "A" if chain in a_chains else "B"
        cls, explicit = _classification(
            atom["residue_name"], atom["atom_name"], atom["element"]
        )
        r = float(
            base_radius(
                atom["residue_name"],
                atom["atom_name"],
                atom["element"],
            )
        )
        R = r + float(probe_radius)
        w = R * R
        idx = len(records)

        records.append(
            AtomRecord(
                index=idx,
                pdb_serial=atom["pdb_serial"],
                record_type=atom["record_type"],
                atom_name=atom["atom_name"],
                altloc=atom["altloc"],
                residue_name=atom["residue_name"],
                chain=chain,
                residue_number=atom["residue_number"],
                insertion_code=atom["insertion_code"],
                element=atom["element"],
                x=atom["x"],
                y=atom["y"],
                z=atom["z"],
                radius_class=cls,
                explicit_assignment=explicit,
                base_radius=r,
                expanded_radius=R,
                weight=w,
                group=group,
            )
        )

    if not records:
        raise ValueError("No selected atoms remained after PDB filtering.")

    fallback_atoms = [a for a in records if not a.explicit_assignment]
    heavy_fallbacks = [a for a in fallback_atoms if a.element != "H"]

    if strict_heavy_fallbacks and heavy_fallbacks:
        preview = ", ".join(
            f"{a.chain}:{a.residue_name}{a.residue_number}:{a.atom_name}"
            for a in heavy_fallbacks[:20]
        )
        raise ValueError(
            f"{len(heavy_fallbacks)} selected heavy atoms lack explicit "
            f"Cazals residue/atom assignments. First cases: {preview}"
        )

    points = np.asarray([[a.x, a.y, a.z] for a in records], dtype=float)
    base = np.asarray([a.base_radius for a in records], dtype=float)
    expanded = np.asarray([a.expanded_radius for a in records], dtype=float)
    weights = np.asarray([a.weight for a in records], dtype=float)
    group_a = np.asarray([a.index for a in records if a.group == "A"], dtype=int)
    group_b = np.asarray([a.index for a in records if a.group == "B"], dtype=int)

    if not len(group_a) or not len(group_b):
        raise ValueError(
            f"Selected groups are empty after filtering: "
            f"A={len(group_a)}, B={len(group_b)}"
        )

    result = PreparedCazalsPDB(
        source=pdb_path,
        points=points,
        base_radii=base,
        expanded_radii=expanded,
        weights=weights,
        group_a=group_a,
        group_b=group_b,
        atoms=records,
        probe_radius=float(probe_radius),
        skipped_altlocs=skipped_altlocs,
        skipped_hydrogens=skipped_h,
        skipped_unselected_chains=skipped_chains,
        fallback_atoms=fallback_atoms,
    )
    if verbose:
        print_preparation_summary(result)
    return result


def print_preparation_summary(data: PreparedCazalsPDB) -> None:
    heavy_fallbacks = [a for a in data.fallback_atoms if a.element != "H"]
    print()
    print("=" * 68)
    print("CAZALS PDB PREPARATION")
    print("=" * 68)
    print(f"Input:                          {data.source}")
    print(f"Selected atoms:                 {len(data.atoms):,}")
    print(f"Group A atoms:                  {len(data.group_a):,}")
    print(f"Group B atoms:                  {len(data.group_b):,}")
    print(f"Explicit Cazals assignments:    {len(data.atoms)-len(data.fallback_atoms):,}")
    print(f"All fallback assignments:       {len(data.fallback_atoms):,}")
    print(f"Heavy-atom fallbacks:           {len(heavy_fallbacks):,}")
    print(f"Alternate conformers removed:   {data.skipped_altlocs:,}")
    print(f"Hydrogens removed:              {data.skipped_hydrogens:,}")
    print(f"Atoms in other chains removed:  {data.skipped_unselected_chains:,}")
    print(f"Probe radius:                   {data.probe_radius:.3f} A")
    print("Expanded radius:                R = r + probe")
    print("Power weight:                   w = R^2")
    if heavy_fallbacks:
        print("\nHeavy-atom fallbacks:")
        for atom in heavy_fallbacks[:30]:
            print(
                f"  idx={atom.index:<5d} serial={atom.pdb_serial:<5d} "
                f"{atom.chain}:{atom.residue_name}{atom.residue_number}:"
                f"{atom.atom_name:<4s} element={atom.element:<2s} "
                f"class={atom.radius_class}"
            )
    print("=" * 68)
    print()


def export_npz(data: PreparedCazalsPDB, path: str | Path) -> Path:
    """
    Export both the original standalone-script schema and additional audit arrays.

    Existing geometry scripts can simply read:
        points, radii, group_a, group_b
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        points=data.points,
        radii=data.expanded_radii,
        base_radii=data.base_radii,
        expanded_radii=data.expanded_radii,
        weights=data.weights,
        group_a=data.group_a,
        group_b=data.group_b,
        pdb_serial=np.asarray([a.pdb_serial for a in data.atoms], dtype=int),
        chain=np.asarray([a.chain for a in data.atoms]),
        residue_name=np.asarray([a.residue_name for a in data.atoms]),
        residue_number=np.asarray([a.residue_number for a in data.atoms], dtype=int),
        atom_name=np.asarray([a.atom_name for a in data.atoms]),
        element=np.asarray([a.element for a in data.atoms]),
        radius_class=np.asarray([a.radius_class for a in data.atoms]),
        probe_radius=np.asarray(data.probe_radius),
    )
    return path


def export_audit_csv(data: PreparedCazalsPDB, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "index", "pdb_serial", "group", "chain", "residue_name",
        "residue_number", "insertion_code", "atom_name", "element",
        "radius_class", "explicit_assignment", "x", "y", "z",
        "base_radius_A", "expanded_radius_A", "weight_A2",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for a in data.atoms:
            writer.writerow(
                {
                    "index": a.index,
                    "pdb_serial": a.pdb_serial,
                    "group": a.group,
                    "chain": a.chain,
                    "residue_name": a.residue_name,
                    "residue_number": a.residue_number,
                    "insertion_code": a.insertion_code,
                    "atom_name": a.atom_name,
                    "element": a.element,
                    "radius_class": a.radius_class,
                    "explicit_assignment": a.explicit_assignment,
                    "x": a.x, "y": a.y, "z": a.z,
                    "base_radius_A": a.base_radius,
                    "expanded_radius_A": a.expanded_radius,
                    "weight_A2": a.weight,
                }
            )
    return path
