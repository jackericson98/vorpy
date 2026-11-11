"""No-solve input gates for the 2KAI solvent-scope comparison."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest

from vorpy.src.boundary import SolventMode, prepare_network_geometry
from vorpy.src.analyze.comparative_geometry_study import audit_input, read_pdb
from vorpy.src.system.system import System


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
SOURCE = DATA / "2KAI_MD_centered.pdb"
REFERENCE = DATA / "2KAI.pdb"
CORRECTED = DATA / "2KAI_MD_centered_no_sol_corrected.pdb"
CHAIN_RESTORED = DATA / "2KAI_MD_centered_chainrestored_solvated.pdb"
FRAME_0100PS = DATA / "2KAI_Frames" / "frame_0100ps.pdb"
SIX_RUN_MANIFEST = ROOT.parent / "docs" / "performance" / "2KAI_six_run_profile.json"


EXPECTED_SHA256 = {
    SOURCE.name: "3eff4f2e5eaa3cc2d8e507ce50534d97efc7ba18bc24ee4e5faa1b9f91f228fc",
    REFERENCE.name: "d195dbaf0b446265cfb72261362414b96b9daa361d7fb81f92ad10dee7edf4cf",
    CORRECTED.name: "53cc21ad3d24ef6e89c8d8fb7a7bf67de9ef52eafd33bc92a30826de80596cb4",
    CHAIN_RESTORED.name: "2e869bee8330cffd4ddffe62f1e6063726fd1f8ebee59f5834bf04411d275e6a",
    FRAME_0100PS.name: "d15ce739150bc47d08d444211b8546af2f9da956840aec9c69d0f933f06185cd",
}


def _require_fixture_paths(*paths):
    """Skip only integration checks when their untracked research inputs are absent."""
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        pytest.skip("missing 2KAI integration fixture/path(s): " + ", ".join(missing))


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


def _frame_identity_signature(atom):
    """Compare frame identity without treating chain relabeling as geometry."""
    return (
        int(atom.line[6:11]), atom.atom_name, atom.altloc, atom.residue_name,
        atom.residue_number, atom.element,
    )


def _kabsch_rmsd(source, target):
    source = source - source.mean(axis=0)
    target = target - target.mean(axis=0)
    left, _, right = np.linalg.svd(source.T @ target)
    if np.linalg.det(left @ right) < 0:
        right[-1, :] *= -1
    rotation = left @ right
    return float(np.sqrt(np.mean(np.sum((source @ rotation - target) ** 2, axis=1))))


def test_synthetic_kabsch_is_rigid_transform_invariant_and_detects_deformation():
    source = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0],
    ])
    rotation = np.asarray([
        [0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0],
    ])
    rigid_target = source @ rotation + np.asarray([17.0, -3.0, 4.0])
    deformed_target = rigid_target.copy()
    deformed_target[3] += np.asarray([0.0, 0.0, 0.25])

    assert _kabsch_rmsd(source, rigid_target) < 1.0e-12
    assert _kabsch_rmsd(source, deformed_target) > 0.05


def test_2kai_canonical_fingerprints_and_authoritative_topology_are_stable():
    _require_fixture_paths(SOURCE, REFERENCE, CORRECTED, CHAIN_RESTORED, FRAME_0100PS)
    assert {
        path.name: sha256(path.read_bytes()).hexdigest()
        for path in (SOURCE, REFERENCE, CORRECTED, CHAIN_RESTORED, FRAME_0100PS)
    } == EXPECTED_SHA256

    reference, _ = read_pdb(REFERENCE)
    assert {a.chain for a in reference if a.category == "protein"} == {"A", "B", "I"}


def test_six_run_manifest_matches_chainrestored_scope_and_stays_scientifically_gated():
    _require_fixture_paths(SIX_RUN_MANIFEST)
    manifest = json.loads(SIX_RUN_MANIFEST.read_text(encoding="utf-8"))
    manifest_paths = [
        ROOT.parent / source["path"] for source in manifest["inputs"].values()
    ]
    _require_fixture_paths(*manifest_paths)
    assert manifest["enabled"] is False
    assert manifest["automatic_execution"] is False
    assert manifest["requires_manual_scientific_approval"] is True

    dry = manifest["inputs"]["no_solvent_corrected"]
    solvated = manifest["inputs"]["solvated_chainrestored"]
    assert dry["sha256"] == EXPECTED_SHA256[CORRECTED.name]
    assert solvated["path"].endswith(CHAIN_RESTORED.name)
    assert solvated["sha256"] == EXPECTED_SHA256[CHAIN_RESTORED.name]
    assert solvated["chain_selection"]["selection_status"] == "PASS"
    assert solvated["chain_selection"]["canonical_group_A"] == ["A", "B"]
    assert solvated["chain_selection"]["canonical_group_B"] == ["I"]
    assert solvated["solvent_ion_policy"] == {
        "mode": "SOLVENT_CONTEXT",
        "water_present": True,
        "water_atoms": 40086,
        "water_participates": True,
        "ions_present": True,
        "ion_species": ["NA", "CL"],
        "ion_atoms": 94,
        "ions_participate": True,
        "generator_counts": {
            "protein": 4340,
            "water": 40086,
            "ions": 94,
            "participating_total": 44520,
        },
        "solvent_cell_status": "NOT_CERTIFIED",
        "retention_rule": "retain explicit water and NA/CL only under explicit solvent policy",
    }

    runs = manifest["runs"]
    assert len(runs) == 6
    assert {(run["scheme"], run["solvent_mode"]) for run in runs} == {
        (scheme, mode)
        for scheme in ("prm", "pow", "aw")
        for mode in ("MOLECULE_ONLY", "SOLVENT_CONTEXT")
    }
    assert all(run["selection_status"] == "PASS" for run in runs)
    assert all(run["scientific_status"] == "PENDING_SCIENTIFIC_APPROVAL" for run in runs)
    assert all(run["group_selectors"] == {"A": "c A and c B", "B": "c I"} for run in runs)
    assert all(run["solvent_mode"] != "FULL_SOLVENT" for run in runs)
    assert "GLOBAL_COMPLETE" not in json.dumps(manifest)

    metric_columns = set(manifest["recording"]["canonical_log_columns"])
    assert {"quantity", "value", "units", "status", "provenance"} <= metric_columns
    assert manifest["recording"]["status_policy"]["zero_substitution"] is False


def test_protein_only_generator_sets_match_in_molecule_only_policy():
    _require_fixture_paths(SOURCE, CORRECTED, CHAIN_RESTORED)
    source_system = System(file=str(SOURCE), name="2kai-source-preflight", make_dir=False)
    corrected_system = System(file=str(CORRECTED), name="2kai-corrected-preflight", make_dir=False)
    restored_system = System(file=str(CHAIN_RESTORED), name="2kai-restored-preflight", make_dir=False)
    source_geometry = prepare_network_geometry(source_system, SolventMode.MOLECULE_ONLY)
    corrected_geometry = prepare_network_geometry(corrected_system, SolventMode.MOLECULE_ONLY)
    restored_geometry = prepare_network_geometry(restored_system, SolventMode.MOLECULE_ONLY)

    assert source_geometry.system_indices == corrected_geometry.system_indices == restored_geometry.system_indices == tuple(range(4340))
    assert source_geometry.stable_ids == corrected_geometry.stable_ids == restored_geometry.stable_ids
    assert np.array_equal(source_geometry.locations, corrected_geometry.locations)
    assert np.array_equal(source_geometry.locations, restored_geometry.locations)
    assert np.array_equal(source_geometry.radii, corrected_geometry.radii)
    assert np.array_equal(source_geometry.radii, restored_geometry.radii)
    assert np.array_equal(source_geometry.masses, corrected_geometry.masses)
    assert np.array_equal(source_geometry.masses, restored_geometry.masses)
    assert source_geometry.solvent_generator_ids == corrected_geometry.solvent_generator_ids == restored_geometry.solvent_generator_ids == ()

    context_geometry = prepare_network_geometry(source_system, SolventMode.SOLVENT_CONTEXT)
    full_geometry = prepare_network_geometry(source_system, SolventMode.FULL_SOLVENT)
    assert len(context_geometry.solvent_generator_ids) == 40086 + 94
    assert len(full_geometry.solvent_generator_ids) == 40086 + 94
    assert context_geometry.solvent_cells_complete is False
    assert full_geometry.solvent_cells_complete is True


def test_chainrestored_solvated_preserves_source_records_and_authoritative_mapping():
    _require_fixture_paths(SOURCE, CHAIN_RESTORED, CORRECTED, REFERENCE)
    source, _ = read_pdb(SOURCE)
    restored, restored_headers = read_pdb(CHAIN_RESTORED)
    corrected, _ = read_pdb(CORRECTED)
    reference, _ = read_pdb(REFERENCE)
    source_protein = [a for a in source if a.category == "protein"]
    restored_protein = [a for a in restored if a.category == "protein"]
    corrected_protein = [a for a in corrected if a.category == "protein"]

    assert len(restored) == len(source) == 44520
    assert [a.category for a in restored] == [a.category for a in source]
    assert [_source_atom_signature(a) for a in restored] == [
        _source_atom_signature(a) for a in source
    ]
    assert [a.line for a in restored if a.category != "protein"] == [
        a.line for a in source if a.category != "protein"
    ]
    assert [(a.chain, a.insertion_code) for a in restored_protein] == [
        (a.chain, a.insertion_code) for a in corrected_protein
    ]
    assert Counter(a.chain for a in restored_protein) == Counter({"A": 1239, "B": 2243, "I": 858})
    assert len(restored_protein) == len(source_protein) == len(corrected_protein) == 4340
    assert sum(a.category == "water" for a in restored) == 40086
    assert Counter(a.residue_name for a in restored if a.category == "ion") == {
        "NA": 53, "CL": 41,
    }

    row, match = audit_input(
        CHAIN_RESTORED, "2KAI", "chainrestored_solvated", restored,
        restored_headers, reference,
    )
    assert row["group_A_chains"] == ["A", "B"]
    assert row["group_B_chains"] == ["I"]
    assert row["all_atoms_needed_for_group_matching_present"] is True
    assert match["missing"] == []
    assert match["duplicates"] == []


def test_0100ps_frame_is_not_a_rigid_transform_of_canonical_source():
    _require_fixture_paths(SOURCE, FRAME_0100PS)
    source, source_headers = read_pdb(SOURCE)
    frame, frame_headers = read_pdb(FRAME_0100PS)
    source_protein = [a for a in source if a.category == "protein"]
    frame_protein = [a for a in frame if a.category == "protein"]
    source_heavy = [a for a in source_protein if not a.hydrogen]
    frame_heavy = [a for a in frame_protein if not a.hydrogen]

    assert len(frame) == len(source) == 44520
    assert len(frame_protein) == len(source_protein) == 4340
    assert [
        _frame_identity_signature(a) for a in source_protein
    ] == [
        _frame_identity_signature(a) for a in frame_protein
    ]
    assert Counter(a.chain for a in frame_protein) == Counter({"A": 1239, "B": 2243, "C": 858})
    assert any(h.startswith("REMARK    GENERATED BY TRJCONV") for h in frame_headers)
    assert any("t= 100.00000 step= 50000" in h for h in frame_headers)
    source_box = next(h for h in source_headers if h.startswith("CRYST1"))
    frame_box = next(h for h in frame_headers if h.startswith("CRYST1"))
    assert source_box != frame_box
    source_xyz = np.asarray([a.xyz for a in source_heavy])
    frame_xyz = np.asarray([a.xyz for a in frame_heavy])
    assert not np.array_equal(source_xyz, frame_xyz)
    assert _kabsch_rmsd(source_xyz, frame_xyz) > 1.0


def test_2kai_inputs_record_presence_and_do_not_share_scope():
    _require_fixture_paths(SOURCE, DATA / "2KAI_MD_centered_no_sol.pdb")
    solvated, _ = read_pdb(SOURCE)
    no_sol, _ = read_pdb(DATA / "2KAI_MD_centered_no_sol.pdb")

    assert len(solvated) == 44520
    assert sum(a.category == "protein" for a in solvated) == 4340
    assert sum(a.category == "water" for a in solvated) == 40086
    assert len({a.residue_occurrence for a in solvated if a.category == "water"}) == 13362
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
    _require_fixture_paths(REFERENCE, DATA / "2KAI_MD_centered_no_sol.pdb", SOURCE)
    reference, _ = read_pdb(REFERENCE)
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
    solvated, _ = read_pdb(SOURCE)
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
    _require_fixture_paths(SOURCE, CORRECTED)
    source, _ = read_pdb(SOURCE)
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
    _require_fixture_paths(CORRECTED, REFERENCE)
    corrected_path = CORRECTED
    corrected, headers = read_pdb(corrected_path)
    reference, _ = read_pdb(REFERENCE)
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
