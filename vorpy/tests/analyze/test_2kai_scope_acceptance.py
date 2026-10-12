"""No-solve input gates for the 2KAI solvent-scope comparison."""
from collections import Counter
from pathlib import Path

from vorpy.src.analyze.comparative_geometry_study import audit_input, read_pdb


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
CORRECTED = DATA / "2KAI_MD_centered_no_sol_corrected.pdb"


def _raw_identity_and_coordinates(line):
    """Compare source records without depending on Python line numbers."""
    return (
        int(line[6:11]),
        line[12:16],
        line[17:20],
        line[21],
        line[22:26],
        line[26],
        line[30:54],
        line[54:66],
    )


def _source_atom_signature(atom):
    """Fields required to remain identical when chain labels are restored."""
    return (
        int(atom.line[6:11]), atom.atom_name, atom.altloc, atom.residue_name,
        atom.residue_number, atom.element, atom.xyz, atom.line[54:66],
    )


def test_2kai_inputs_record_presence_and_do_not_share_scope():
    solvated, _ = read_pdb(DATA / "2KAI_MD_centered.pdb")
    no_sol, _ = read_pdb(DATA / "2KAI_MD_centered_no_sol.pdb")

    assert len(solvated) == 44520
    assert sum(a.category == "protein" for a in solvated) == 4340
    assert sum(a.category == "water" for a in solvated) == 40086
    assert sum(a.category == "ion" for a in solvated) == 94
    assert Counter(a.residue_name for a in solvated if a.category == "ion") == {
        "NA": 53, "CL": 41,
    }
    assert not any(a.residue_name == "MG" for a in solvated)
    assert len(no_sol) == 3101
    assert all(a.category == "protein" for a in no_sol)
    assert all(a.chain == "" for a in solvated + no_sol)

    solvated_protein = [a for a in solvated if a.category == "protein"]
    assert [_raw_identity_and_coordinates(a.line) for a in no_sol] == [
        _raw_identity_and_coordinates(a.line) for a in solvated_protein[1239:]
    ]
    assert solvated_protein[0].residue_name == "ILE"
    assert no_sol[0].residue_name == "ALA"


def test_2kai_no_solvent_file_is_not_eligible_for_deposited_ab_vs_i_scope():
    reference, _ = read_pdb(DATA / "2KAI.pdb")
    reference_chains = {a.chain for a in reference if a.category == "protein"}
    assert reference_chains == {"A", "B", "I"}

    reference_first = {
        chain: next(a for a in reference if a.category == "protein" and a.chain == chain)
        for chain in ("A", "B", "I")
    }
    no_sol, _ = read_pdb(DATA / "2KAI_MD_centered_no_sol.pdb")
    # A blank-chain input cannot prove the deposited A+B versus I partition;
    # its exact suffix relationship also shows that the A prefix is absent.
    assert {a.chain for a in no_sol} == {""}
    solvated, _ = read_pdb(DATA / "2KAI_MD_centered.pdb")
    assert len(no_sol) < len([a for a in solvated if a.category == "protein"])
    assert (no_sol[0].residue_name, no_sol[0].residue_number, no_sol[0].atom_name) == (
        reference_first["B"].residue_name,
        reference_first["B"].residue_number,
        reference_first["B"].atom_name,
    )
    assert (solvated[0].residue_name, solvated[0].residue_number, solvated[0].atom_name) == (
        reference_first["A"].residue_name,
        reference_first["A"].residue_number,
        reference_first["A"].atom_name,
    )


def test_corrected_no_solvent_view_preserves_all_protein_atoms_and_coordinates():
    source, _ = read_pdb(DATA / "2KAI_MD_centered.pdb")
    corrected, _ = read_pdb(CORRECTED)
    source_protein = [a for a in source if a.category == "protein"]

    assert len(corrected) == len(source_protein) == 4340
    assert all(a.category == "protein" for a in corrected)
    assert Counter(a.chain for a in corrected) == Counter({"A": 1239, "B": 2243, "I": 858})
    assert len({a.identity for a in corrected}) == 4340
    assert [_source_atom_signature(a) for a in corrected] == [
        _source_atom_signature(a) for a in source_protein
    ]


def test_corrected_view_passes_authoritative_group_preflight():
    corrected_path = CORRECTED
    corrected, headers = read_pdb(corrected_path)
    reference, _ = read_pdb(DATA / "2KAI.pdb")
    row, match = audit_input(
        corrected_path, "2KAI", "corrected_no_solvent", corrected, headers, reference
    )

    assert row["group_A_chains"] == ["A", "B"]
    assert row["group_B_chains"] == ["I"]
    assert row["all_atoms_needed_for_group_matching_present"] is True
    assert match["missing"] == []
    assert match["duplicates"] == []
    assert match["extra"] == [
        ("A", 95, "B", "SER", "OXT"),
        ("I", 57, "", "GLY", "OXT"),
    ]
