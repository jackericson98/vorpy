"""Preflight tests: identity matching must never substitute row matching."""
from pathlib import Path

import pytest

from vorpy.src.analyze.comparative_geometry_study import (
    audit_input, biological_groups, canonical_identity, construct_views,
    discover_inputs, heavy_atoms, infer_chain_map, match_atoms, read_pdb, write_csv,
)


def pdb_atom(serial, name, resname, chain, resid, element, *, icode="", x=1.):
    return (f"ATOM  {serial:5d} {name:^4} {resname:>3} {chain or ' '}{resid:4d}{icode or ' '}   "
            f"{x:8.3f}{2.:8.3f}{3.:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2}  ")


def load(tmp_path, name, lines):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path, *read_pdb(path)


def test_atom_classification_does_not_use_atom_hetatm_for_solvent(tmp_path):
    _, atoms, _ = load(tmp_path, "solvent.pdb", [
        pdb_atom(1, "CA", "ALA", "A", 1, ""),
        pdb_atom(2, "HA", "ALA", "A", 1, ""),
        pdb_atom(3, "OW", "SOL", "", 9999, ""),
        pdb_atom(4, "HW1", "SOL", "", 9999, ""),
        pdb_atom(5, "OW", "SOL", "", 0, ""),
        pdb_atom(6, "OW", "SOL", "", 9999, ""),
        pdb_atom(7, "NA", "NA", "", 1, ""),
    ])
    assert atoms[0].element == "C"
    assert atoms[1].hydrogen
    assert len(heavy_atoms(atoms)) == 1
    assert len({a.residue_occurrence for a in atoms if a.category == "water"}) == 3
    assert atoms[-1].category == "ion"
    assert atoms[-1].element == "NA"


def test_identity_match_is_order_independent_and_preserves_insertion_codes(tmp_path):
    _, reference, _ = load(tmp_path, "crystal.pdb", [
        pdb_atom(1, "CD1", "ILE", "A", 16, "C"),
        pdb_atom(2, "CA", "SER", "B", 95, "C", icode="Y"),
        pdb_atom(3, "CA", "GLY", "I", 57, "C"),
    ])
    _, md, _ = load(tmp_path, "frame.pdb", [
        pdb_atom(50, "CA", "GLY", "C", 57, "C"),
        pdb_atom(70, "OXT", "GLY", "C", 57, "O"),
        pdb_atom(80, "CA", "SER", "B", 95, "C", icode="Y"),
        pdb_atom(20, "CD", "ILE", "A", 16, "C"),
    ])
    mapping, _ = infer_chain_map(md, reference)
    assert mapping == {"A": "A", "B": "B", "C": "I"}
    match = match_atoms(md, reference, mapping)
    assert match["matched"] == 3
    assert match["atom_aliases_applied"] == 1
    assert match["missing"] == []
    assert match["extra"] == [("I", 57, "", "GLY", "OXT")]
    assert canonical_identity(md[2], mapping)[2] == "Y"


def test_missing_duplicate_and_blank_chain_inputs_are_not_silently_harmonized(tmp_path):
    _, reference, _ = load(tmp_path, "ref.pdb", [pdb_atom(1, "CA", "ALA", "A", 1, "C")])
    path, md, headers = load(tmp_path, "dup.pdb", [pdb_atom(2, "CA", "ALA", "A", 1, "C")] * 2)
    row, match = audit_input(path, "2KAI", "MD_frame", md, headers, reference)
    assert match["duplicates"]
    assert not row["all_atoms_needed_for_group_matching_present"]
    assert construct_views(path, md, reference, row, tmp_path / "out") == ([], [])
    _, blank, _ = load(tmp_path, "blank.pdb", [pdb_atom(2, "CA", "ALA", "", 1, "C")])
    assert infer_chain_map(blank, reference)[0] == {}
    assert match_atoms(blank, reference, {})["missing"]


def test_paired_views_preserve_coordinates_and_solvent_columns(tmp_path):
    _, reference, _ = load(tmp_path, "ref.pdb", [pdb_atom(1, "CA", "ALA", "I", 1, "C")])
    path, md, headers = load(tmp_path, "frame.pdb", [
        pdb_atom(9, "CA", "ALA", "C", 1, "C", x=72.123),
        pdb_atom(10, "HA", "ALA", "C", 1, "H"),
        pdb_atom(11, "OXT", "ALA", "C", 1, "O"),
        pdb_atom(12, "OW", "SOL", "", 1, "O"),
    ])
    row, _ = audit_input(path, "2KAI", "MD_frame", md, headers, reference)
    manifest, matches = construct_views(path, md, reference, row, tmp_path / "out")
    by_tier = {r["tier"]: r for r in manifest}
    assert {k: r["atoms"] for k, r in by_tier.items()} == {"H": 1, "P": 3, "S": 4, "HS_context": 2}
    h, _ = read_pdb(by_tier["H"]["path"])
    hs, _ = read_pdb(by_tier["HS_context"]["path"])
    solvated, _ = read_pdb(by_tier["S"]["path"])
    assert h[0].identity == hs[0].identity == reference[0].identity
    assert h[0].xyz == hs[0].xyz == md[0].xyz
    assert all(a.line[30:54] == b.line[30:54] for a, b in zip(md, solvated))
    assert hs[-1].chain == ""
    assert hs[-1].xyz == md[-1].xyz
    assert matches[0]["source_id"][0] == "C"


def test_biological_groups_follow_author_assembly_not_all_molecule_copies():
    headers = [
        "COMPND    MOL_ID: 1; MOLECULE: CYCLOPHILIN; CHAIN: A, B;",
        "COMPND    MOL_ID: 2; MOLECULE: CAPSID; CHAIN: C, D;",
        "REMARK 350 BIOMOLECULE: 1",
        "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, D",
        "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000  0.00000",
        "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000  0.00000",
        "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000  0.00000",
    ]
    groups, evidence = biological_groups(headers)
    assert groups == {"A": ["A"], "B": ["D"]}
    assert "assembly 1" in evidence
    headers[-1] = headers[-1].replace("1.000000", "2.000000")
    assert biological_groups(headers)[0] == {"A": [], "B": []}


def test_discovery_does_not_drop_missing_systems_or_unusual_frame_names(tmp_path):
    frames = tmp_path / "frames"
    frames.mkdir()
    for path in (tmp_path / "2KAI.pdb", tmp_path / "2KAI_MD_centered.pdb", frames / "unexpected.pdb"):
        path.touch()
    found, missing = discover_inputs(tmp_path, frames, ["2KAI", "1AK4", "1BRC"])
    assert len(found) == 3
    assert missing == ["1AK4", "1BRC"]


def test_nonfinite_coordinates_and_csv_metrics_fail_explicitly(tmp_path):
    with pytest.raises(ValueError, match="Nonfinite coordinates"):
        load(tmp_path, "bad.pdb", [pdb_atom(1, "CA", "ALA", "A", 1, "C", x=float("nan"))])
    with pytest.raises(ValueError, match="Nonfinite CSV"):
        write_csv(tmp_path / "bad.csv", [{"area": float("nan")}])


def test_supplied_static_partner_definitions():
    data = Path(__file__).resolve().parents[2] / "data"
    for system, expected in [("1AK4", {"A": ["A"], "B": ["D"]}),
                             ("1BRC", {"A": ["E"], "B": ["I"]})]:
        path = data / f"{system}.pdb"
        if not path.exists():
            pytest.skip(f"Optional input unavailable: {path}")
        _, headers = read_pdb(path)
        assert biological_groups(headers)[0] == expected
