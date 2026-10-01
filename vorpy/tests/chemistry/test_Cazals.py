"""
Tests for vorpy.src.chemistry.Cazals.

Expected repository layout
--------------------------
vorpy/
├── data/
├── src/
│   └── chemistry/
│       └── Cazals.py
└── tests/
    └── chemistry/
        └── test_Cazals.py

These tests deliberately include:
1. Exact tests of the published Cazals/Chothia group assignments encoded by
   Cazals.py.
2. Mathematical consistency tests for probe expansion and power weights.
3. Integration tests against real PDB files under vorpy/data.
4. An optional 2KAI-specific audit when a 2KAI PDB is present.

Run from the repository root with, for example:

    pytest -q vorpy/tests/chemistry/test_Cazals.py -s

The -s flag is useful because the data audit prints a classification summary.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import pytest

from vorpy.src.chemistry.Cazals import (
    ALIPHATIC_C,
    CHARGED_N,
    NEUTRAL_N,
    OTHER_RADIUS,
    OXYGEN,
    PROBE_RADIUS,
    SULFUR,
    TRIGONAL_C,
    base_radius,
    element_radii,
    expanded_radius,
    power_weight,
    special_radii,
)


# ---------------------------------------------------------------------------
# Repository/data discovery
# ---------------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()


def _find_data_dir() -> Path:
    """
    Locate the VorPy data directory without assuming pytest's working directory.

    Starting from this test file:
        .../vorpy/vorpy/tests/chemistry/test_Cazals.py

    the expected data directory is:
        .../vorpy/vorpy/data
    """
    for parent in THIS_FILE.parents:
        candidate = parent / "data"
        if candidate.is_dir():
            return candidate

    pytest.fail(
        "Could not locate VorPy's data directory by walking upward from "
        f"{THIS_FILE}"
    )


def _pdb_files(data_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in data_dir.rglob("*")
        if path.is_file() and path.suffix.lower() == ".pdb"
    )


def _find_2kai(data_dir: Path) -> Path | None:
    candidates = [
        path
        for path in _pdb_files(data_dir)
        if "2kai" in path.name.lower() or "2kai" in str(path.parent).lower()
    ]
    return candidates[0] if candidates else None


# ---------------------------------------------------------------------------
# Minimal PDB reader for the radii integration tests
# ---------------------------------------------------------------------------

def _infer_element(atom_name: str, element_field: str) -> str:
    """
    Infer an element only when the PDB element column is empty.

    For ordinary protein ATOM records this is sufficient for the present
    Cazals heavy-atom audit. The residue/atom-specific table remains the
    preferred source for C/N chemical classification.
    """
    element = element_field.strip().upper()
    if element:
        return element

    name = atom_name.strip().upper()
    if not name:
        return ""

    # Strip a leading digit used by some hydrogen naming conventions.
    while name and name[0].isdigit():
        name = name[1:]

    if not name:
        return ""

    # Protein atom names such as CA, CB, CG are carbons, not calcium.
    if name[0] in {"C", "N", "O", "S", "H", "P"}:
        return name[0]

    return name[:2]


def _read_pdb_atoms(path: Path) -> list[dict]:
    """
    Read ATOM/HETATM records using fixed PDB columns.

    This test intentionally does not depend on VorPy's full molecule parser:
    if a chemistry-radii test fails, parser behavior should not obscure the
    radius assignment being tested.
    """
    atoms = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[0:6].strip().upper()
            if record not in {"ATOM", "HETATM"}:
                continue

            atom_name = line[12:16].strip().upper()
            residue_name = line[17:20].strip().upper()
            chain = line[21:22].strip()
            residue_number = line[22:26].strip()
            element = _infer_element(atom_name, line[76:78] if len(line) >= 78 else "")

            atoms.append(
                {
                    "record": record,
                    "atom_name": atom_name,
                    "residue_name": residue_name,
                    "chain": chain,
                    "residue_number": residue_number,
                    "element": element,
                }
            )

    return atoms


# ---------------------------------------------------------------------------
# Constants and table structure
# ---------------------------------------------------------------------------

def test_cazals_constants():
    assert ALIPHATIC_C == pytest.approx(1.87)
    assert TRIGONAL_C == pytest.approx(1.76)
    assert NEUTRAL_N == pytest.approx(1.65)
    assert CHARGED_N == pytest.approx(1.50)
    assert OXYGEN == pytest.approx(1.40)
    assert SULFUR == pytest.approx(1.85)
    assert OTHER_RADIUS == pytest.approx(2.00)
    assert PROBE_RADIUS == pytest.approx(1.40)


def test_all_twenty_standard_amino_acids_are_present():
    expected = {
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    }
    assert expected == set(special_radii)


def test_every_standard_residue_has_backbone_assignments():
    for residue, atoms in special_radii.items():
        assert atoms["N"] == pytest.approx(NEUTRAL_N), residue
        assert atoms["CA"] == pytest.approx(ALIPHATIC_C), residue
        assert atoms["C"] == pytest.approx(TRIGONAL_C), residue
        assert atoms["O"] == pytest.approx(OXYGEN), residue


# ---------------------------------------------------------------------------
# Chemically important classification tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("residue", "atom", "expected"),
    [
        # Backbone
        ("ALA", "CA", ALIPHATIC_C),
        ("ALA", "C", TRIGONAL_C),
        ("ALA", "N", NEUTRAL_N),
        ("ALA", "O", OXYGEN),

        # Aliphatic side chains
        ("VAL", "CB", ALIPHATIC_C),
        ("VAL", "CG1", ALIPHATIC_C),
        ("LEU", "CD1", ALIPHATIC_C),
        ("ILE", "CG2", ALIPHATIC_C),

        # Aromatic/trigonal carbons
        ("PHE", "CG", TRIGONAL_C),
        ("PHE", "CZ", TRIGONAL_C),
        ("TYR", "CE1", TRIGONAL_C),
        ("TRP", "CD2", TRIGONAL_C),

        # Carbonyl/carboxyl carbons
        ("ASN", "CG", TRIGONAL_C),
        ("ASP", "CG", TRIGONAL_C),
        ("GLN", "CD", TRIGONAL_C),
        ("GLU", "CD", TRIGONAL_C),

        # Neutral nitrogens
        ("ASN", "ND2", NEUTRAL_N),
        ("GLN", "NE2", NEUTRAL_N),
        ("HIS", "ND1", NEUTRAL_N),
        ("TRP", "NE1", NEUTRAL_N),

        # Charged/basic nitrogens encoded by the model
        ("LYS", "NZ", CHARGED_N),
        ("ARG", "NE", CHARGED_N),
        ("ARG", "NH1", CHARGED_N),
        ("ARG", "NH2", CHARGED_N),

        # Heteroatoms
        ("SER", "OG", OXYGEN),
        ("THR", "OG1", OXYGEN),
        ("TYR", "OH", OXYGEN),
        ("CYS", "SG", SULFUR),
        ("MET", "SD", SULFUR),
    ],
)
def test_representative_atom_classifications(residue, atom, expected):
    assert base_radius(residue, atom) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# Radius/weight mathematics
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("residue", "atom", "element"),
    [
        ("ALA", "CA", "C"),
        ("ALA", "C", "C"),
        ("LYS", "NZ", "N"),
        ("ARG", "NH1", "N"),
        ("ASP", "OD1", "O"),
        ("CYS", "SG", "S"),
    ],
)
def test_probe_expansion_is_base_plus_1p4(residue, atom, element):
    r = base_radius(residue, atom, element)
    R = expanded_radius(residue, atom, element)

    assert R == pytest.approx(r + PROBE_RADIUS)


@pytest.mark.parametrize(
    ("residue", "atom", "element"),
    [
        ("ALA", "CA", "C"),
        ("PHE", "CZ", "C"),
        ("LYS", "NZ", "N"),
        ("SER", "OG", "O"),
        ("MET", "SD", "S"),
    ],
)
def test_power_weight_is_squared_expanded_radius(residue, atom, element):
    R = expanded_radius(residue, atom, element)
    w = power_weight(residue, atom, element)

    assert w == pytest.approx(R * R)
    assert w > 0.0
    assert math.isfinite(w)


def test_custom_probe_radius():
    probe = 2.25
    r = base_radius("ALA", "CA", "C")

    assert expanded_radius("ALA", "CA", "C", probe) == pytest.approx(r + probe)
    assert power_weight("ALA", "CA", "C", probe) == pytest.approx((r + probe) ** 2)


# ---------------------------------------------------------------------------
# Fallback behavior
# ---------------------------------------------------------------------------

def test_known_oxygen_element_fallback():
    assert base_radius("UNKNOWN", "OX", "O") == pytest.approx(OXYGEN)


def test_known_sulfur_element_fallback():
    assert base_radius("UNKNOWN", "SX", "S") == pytest.approx(SULFUR)


def test_ambiguous_carbon_fallback_does_not_guess_class():
    assert element_radii["C"] == pytest.approx(OTHER_RADIUS)
    assert base_radius("UNKNOWN", "CX", "C") == pytest.approx(OTHER_RADIUS)


def test_ambiguous_nitrogen_fallback_does_not_guess_charge():
    assert element_radii["N"] == pytest.approx(OTHER_RADIUS)
    assert base_radius("UNKNOWN", "NX", "N") == pytest.approx(OTHER_RADIUS)


def test_completely_unknown_atom_uses_other_radius():
    assert base_radius("UNKNOWN", "XX", "XX") == pytest.approx(OTHER_RADIUS)


# ---------------------------------------------------------------------------
# Integration against VorPy's real data files
# ---------------------------------------------------------------------------

def test_vorpy_data_contains_pdb_files():
    data_dir = _find_data_dir()
    pdbs = _pdb_files(data_dir)

    assert pdbs, f"No PDB files were found under {data_dir}"


def test_real_vorpy_protein_atoms_receive_valid_radii():
    """
    Audit real protein ATOM records in VorPy's data directory.

    We deliberately test only atoms whose residue+atom pair is explicitly
    represented by special_radii. This separates:
      * correct assignments,
      * unsupported/extra atom names,
      * and fallback behavior.

    The test also requires that at least one real protein atom is exercised.
    """
    data_dir = _find_data_dir()
    pdbs = _pdb_files(data_dir)

    tested = 0

    for pdb in pdbs:
        for atom in _read_pdb_atoms(pdb):
            residue = atom["residue_name"]
            name = atom["atom_name"]
            element = atom["element"]

            if residue not in special_radii:
                continue
            if name not in special_radii[residue]:
                continue

            r = base_radius(residue, name, element)
            R = expanded_radius(residue, name, element)
            w = power_weight(residue, name, element)

            assert math.isfinite(r) and r > 0.0, (pdb, atom)
            assert math.isfinite(R) and R > 0.0, (pdb, atom)
            assert math.isfinite(w) and w > 0.0, (pdb, atom)

            assert R == pytest.approx(r + PROBE_RADIUS), (pdb, atom)
            assert w == pytest.approx(R ** 2), (pdb, atom)

            tested += 1

    assert tested > 0, (
        "PDB files were found, but no standard protein atoms matched "
        "Cazals.special_radii."
    )


def test_real_vorpy_data_cazals_classification_audit():
    """
    Print a useful audit of standard-protein atoms in the repository.

    This is intentionally diagnostic rather than requiring 100% coverage:
    hydrogens, terminal aliases, alternate protonation states, and unusual
    residue naming may require explicit decisions before an exact Cazals
    reproduction.

    Run with `pytest -s` to see the report.
    """
    data_dir = _find_data_dir()
    pdbs = _pdb_files(data_dir)

    explicit = Counter()
    unsupported = Counter()
    total_standard = 0

    for pdb in pdbs:
        for atom in _read_pdb_atoms(pdb):
            residue = atom["residue_name"]
            name = atom["atom_name"]

            if residue not in special_radii:
                continue

            total_standard += 1

            if name in special_radii[residue]:
                explicit[(residue, name)] += 1
            else:
                unsupported[(residue, name, atom["element"])] += 1

    assert total_standard > 0

    explicit_count = sum(explicit.values())
    unsupported_count = sum(unsupported.values())
    coverage = 100.0 * explicit_count / total_standard

    print("\n")
    print("=" * 72)
    print("CAZALS RADII AUDIT: vorpy/data")
    print("=" * 72)
    print(f"Data directory:               {data_dir}")
    print(f"PDB files scanned:            {len(pdbs):,}")
    print(f"Standard-protein atoms:       {total_standard:,}")
    print(f"Explicitly classified atoms:  {explicit_count:,}")
    print(f"Unclassified atom records:    {unsupported_count:,}")
    print(f"Explicit coverage:             {coverage:.2f}%")

    if unsupported:
        print("\nMost common unclassified residue/atom combinations:")
        for (residue, name, element), count in unsupported.most_common(30):
            print(
                f"  {residue:>3s} {name:<5s} element={element:<2s} "
                f"count={count:,}"
            )

    print("=" * 72)


# ---------------------------------------------------------------------------
# Optional 2KAI-specific audit
# ---------------------------------------------------------------------------

def test_2kai_cazals_assignment_audit():
    """
    Audit 2KAI specifically when it exists in vorpy/data.

    This test is intended to become increasingly strict as the 17-degree
    Cazals reproduction is developed. For now it verifies:
      * 2KAI can be found/read,
      * standard protein heavy atoms are present,
      * explicitly mapped atoms have mathematically valid Cazals weights,
      * and it reports any heavy atoms lacking an explicit residue/atom map.

    It skips cleanly if 2KAI is not currently in the data directory.
    """
    data_dir = _find_data_dir()
    pdb = _find_2kai(data_dir)

    if pdb is None:
        pytest.skip("No 2KAI PDB was found under vorpy/data")

    atoms = _read_pdb_atoms(pdb)
    assert atoms, f"No ATOM/HETATM records found in {pdb}"

    standard_heavy = []
    unsupported_heavy = Counter()

    for atom in atoms:
        residue = atom["residue_name"]
        name = atom["atom_name"]
        element = atom["element"]

        if residue not in special_radii:
            continue
        if element == "H" or name.startswith("H"):
            continue

        standard_heavy.append(atom)

        if name not in special_radii[residue]:
            unsupported_heavy[(residue, name, element)] += 1
            continue

        r = base_radius(residue, name, element)
        R = expanded_radius(residue, name, element)
        w = power_weight(residue, name, element)

        assert r > 0.0 and math.isfinite(r), atom
        assert R == pytest.approx(r + PROBE_RADIUS), atom
        assert w == pytest.approx(R ** 2), atom

    assert standard_heavy, f"No standard protein heavy atoms found in {pdb}"

    print("\n")
    print("=" * 72)
    print("2KAI CAZALS RADII AUDIT")
    print("=" * 72)
    print(f"File:                         {pdb}")
    print(f"All PDB atom records:         {len(atoms):,}")
    print(f"Standard protein heavy atoms: {len(standard_heavy):,}")
    print(
        "Unclassified heavy atoms:    "
        f"{sum(unsupported_heavy.values()):,}"
    )

    if unsupported_heavy:
        print("\nUnclassified 2KAI heavy atoms:")
        for (residue, name, element), count in unsupported_heavy.most_common():
            print(
                f"  {residue:>3s} {name:<5s} element={element:<2s} "
                f"count={count:,}"
            )

    print("=" * 72)

    # This is deliberately not yet:
    #
    #     assert not unsupported_heavy
    #
    # First inspect the audit. Once all 2KAI heavy-atom naming/protonation
    # cases have been resolved, make that assertion strict. That will give
    # us a clean radii gate before attempting the published ~17 degree
    # interface-curvature reproduction.


@pytest.mark.parametrize(
    "residue",
    [
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    ],
)
def test_terminal_oxt_is_oxygen(residue):
    """Every standard amino acid may occur at the C terminus with OXT."""
    assert base_radius(residue, "OXT", "O") == pytest.approx(OXYGEN)

