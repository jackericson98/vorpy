"""Bounded, no-solve checks for the six-run 2KAI output manifest."""
import hashlib
import json
from pathlib import Path

import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PACKAGE_ROOT.parent
MANIFEST = REPOSITORY_ROOT / "docs" / "performance" / "2KAI_six_run_profile.json"


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _synthetic_manifest(tmp_path):
    reference = tmp_path / "reference.pdb"
    dry = tmp_path / "dry.pdb"
    solvated = tmp_path / "solvated_chainrestored.pdb"
    reference.write_text("REFERENCE\n", encoding="ascii")
    dry.write_text("DRY\n", encoding="ascii")
    solvated.write_text("SOLVATED\n", encoding="ascii")

    def source(path, scope, mode, water_present, ions_present):
        return {
            "path": str(path),
            "sha256": _sha256(path),
            "bytes": path.stat().st_size,
            "scope": scope,
            "preflight_status": "PASS",
            "chain_selection": {
                "raw_group_A": ["A", "B"],
                "raw_group_B": ["I"],
                "canonical_group_A": ["A", "B"],
                "canonical_group_B": ["I"],
                "chain_map": {"A": "A", "B": "B", "I": "I"},
                "selection_status": "PASS",
                "selection_source": "synthetic fixture",
            },
            "solvent_ion_policy": {
                "mode": mode,
                "water_present": water_present,
                "water_participates": water_present,
                "ions_present": ions_present,
                "ions_participate": ions_present,
            },
        }

    inputs = {
        "no_solvent_corrected": source(dry, "dry", "MOLECULE_ONLY", False, False),
        "solvated_chainrestored": source(
            solvated, "solvated", "SOLVENT_CONTEXT", True, True
        ),
    }
    runs = []
    for scheme in ("prm", "pow", "aw"):
        for suffix, input_key, mode, boundary in (
            ("no_solvent", "no_solvent_corrected", "MOLECULE_ONLY", "none"),
            ("solvent", "solvated_chainrestored", "SOLVENT_CONTEXT", "explicit"),
        ):
            run_id = f"synthetic_{scheme}_{suffix}"
            runs.append({
                "id": run_id,
                "scheme": scheme,
                "solvent_mode": mode,
                "boundary": boundary,
                "input_key": input_key,
                "producer_identity": f"producer:synthetic:{scheme}:{inputs[input_key]['sha256']}",
                "archive_identity": None,
                "dependency_validity": None,
                "lineage": None,
                "input_scope": inputs[input_key]["scope"],
                "run_status": "PENDING_SCIENTIFIC_APPROVAL",
                "selection_status": "PASS",
                "scientific_status": "PENDING_SCIENTIFIC_APPROVAL",
                "output_directory": f"output/synthetic/{run_id}",
                "archive": f"output/synthetic/{run_id}/{run_id}.vpy",
                "profile": f"output/synthetic/{run_id}/cache_export_profile.json",
                "canonical_results_log": f"output/synthetic/{run_id}/interface_*/logs.csv",
                "group_selectors": {"A": "c A and c B", "B": "c I"},
                "measurements": {
                    "wall_seconds": None,
                    "cpu_seconds": None,
                    "peak_rss_bytes": None,
                    "output_bytes": None,
                },
            })
    return {
        "format": "vorpy.2kai_six_run_profile",
        "version": 2,
        "enabled": False,
        "automatic_execution": False,
        "requires_manual_scientific_approval": True,
        "input_reference": {
            "path": str(reference),
            "sha256": _sha256(reference),
            "chain_groups": {"A": ["A", "B"], "B": ["I"]},
        },
        "inputs": inputs,
        "runs": runs,
        "recording": {
            "lineage_field": "provenance.details.geometry_lineage",
            "identity_fields": ["identity.source_network_id"],
            "status_policy": {"zero_substitution": False},
            "archive_timing_fields": ["total"],
            "canonical_log_columns": ["interface_id", "status", "provenance"],
            "windows_memory_limit": "Windows peak RSS is unavailable",
        },
    }


def _assert_six_run_schema(manifest):
    assert manifest["format"] == "vorpy.2kai_six_run_profile"
    assert manifest["version"] == 2
    assert manifest["enabled"] is False
    assert manifest["automatic_execution"] is False
    assert manifest["requires_manual_scientific_approval"] is True
    assert set(manifest["inputs"]) == {
        "no_solvent_corrected",
        "solvated_chainrestored",
    }
    assert manifest["inputs"]["no_solvent_corrected"]["scope"] == "dry"
    assert manifest["inputs"]["solvated_chainrestored"]["scope"] == "solvated"
    for source in manifest["inputs"].values():
        assert source["preflight_status"] == "PASS"
        assert source["chain_selection"]["raw_group_A"] == ["A", "B"]
        assert source["chain_selection"]["raw_group_B"] == ["I"]
        assert source["chain_selection"]["selection_status"] == "PASS"

    runs = manifest["runs"]
    assert len(runs) == 6
    assert {(run["scheme"], run["solvent_mode"]) for run in runs} == {
        (scheme, mode)
        for scheme in ("prm", "pow", "aw")
        for mode in ("MOLECULE_ONLY", "SOLVENT_CONTEXT")
    }
    assert len({run["output_directory"] for run in runs}) == 6
    assert len({run["archive"] for run in runs}) == 6
    assert len({run["producer_identity"] for run in runs}) == 6

    for run in runs:
        source = manifest["inputs"][run["input_key"]]
        assert source["sha256"] in run["producer_identity"]
        assert run["input_scope"] == source["scope"]
        assert run["solvent_mode"] == source["solvent_ion_policy"]["mode"]
        assert run["group_selectors"] == {"A": "c A and c B", "B": "c I"}
        assert run["selection_status"] == "PASS"
        assert run["run_status"] == "PENDING_SCIENTIFIC_APPROVAL"
        assert run["boundary"] == (
            "none" if run["solvent_mode"] == "MOLECULE_ONLY" else "explicit"
        )
        assert run["archive_identity"] is None
        assert run["dependency_validity"] is None
        assert run["lineage"] is None
        assert run["scientific_status"] == "PENDING_SCIENTIFIC_APPROVAL"
        assert run["measurements"] == {
            "wall_seconds": None,
            "cpu_seconds": None,
            "peak_rss_bytes": None,
            "output_bytes": None,
        }
        assert run["archive"].startswith(run["output_directory"] + "/")
        assert run["profile"].startswith(run["output_directory"] + "/")
        assert run["canonical_results_log"].startswith(run["output_directory"] + "/")

    recording = manifest["recording"]
    assert "provenance.details.geometry_lineage" == recording["lineage_field"]
    assert "identity.source_network_id" in recording["identity_fields"]
    assert recording["status_policy"]["zero_substitution"] is False
    assert recording["archive_timing_fields"]
    assert recording["canonical_log_columns"][-2:] == ["status", "provenance"]
    assert "Windows" in recording["windows_memory_limit"]


def _integration_manifest_or_skip(manifest_path=MANIFEST):
    if not manifest_path.is_file():
        pytest.skip(
            "fixture-dependent 2KAI integration skipped: disabled manifest is absent "
            "from this fresh checkout"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    required = [REPOSITORY_ROOT / manifest["input_reference"]["path"]]
    required.extend(
        REPOSITORY_ROOT / source["path"] for source in manifest["inputs"].values()
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        pytest.skip(
            "fixture-dependent 2KAI integration skipped; missing untracked fixture(s): "
            + ", ".join(missing)
        )
    return manifest


def test_2kai_manifest_schema_is_self_contained(tmp_path):
    manifest = _synthetic_manifest(tmp_path)
    _assert_six_run_schema(manifest)

    reference = Path(manifest["input_reference"]["path"])
    assert _sha256(reference) == manifest["input_reference"]["sha256"]
    for source in manifest["inputs"].values():
        path = Path(source["path"])
        assert path.is_file()
        assert path.stat().st_size == source["bytes"]
        assert _sha256(path) == source["sha256"]


def test_2kai_manifest_integration_skip_is_clear_without_manifest(tmp_path):
    with pytest.raises(pytest.skip.Exception, match="fixture-dependent 2KAI integration skipped"):
        _integration_manifest_or_skip(tmp_path / "absent-2kai-manifest.json")


def test_2kai_manifest_integration_skip_is_clear_without_pdb_fixtures(tmp_path):
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps({
        "input_reference": {"path": "missing-reference.pdb"},
        "inputs": {
            "dry": {"path": "missing-dry.pdb"},
            "solvated": {"path": "missing-solvated.pdb"},
        },
    }), encoding="utf-8")
    with pytest.raises(pytest.skip.Exception, match="missing untracked fixture"):
        _integration_manifest_or_skip(manifest_path)


def test_2kai_manifest_is_disabled_and_covers_unique_scheme_solvent_runs():
    manifest = _integration_manifest_or_skip()

    _assert_six_run_schema(manifest)
    reference = REPOSITORY_ROOT / manifest["input_reference"]["path"]
    assert _sha256(reference) == manifest["input_reference"]["sha256"]
    for key, source in manifest["inputs"].items():
        path = REPOSITORY_ROOT / source["path"]
        assert path.stat().st_size == source["bytes"], key
        assert _sha256(path) == source["sha256"]
        assert source["preflight_status"] == "PASS"

    assert "2KAI_Frames" not in json.dumps(manifest["inputs"])
    assert manifest["inputs"]["no_solvent_corrected"]["scope"] == "dry"
    solvated = manifest["inputs"]["solvated_chainrestored"]
    assert solvated["scope"] == "solvated"
    assert solvated["path"].endswith("2KAI_MD_centered_chainrestored_solvated.pdb")
    assert solvated["sha256"] == (
        "2e869bee8330cffd4ddffe62f1e6063726fd1f8ebee59f5834bf04411d275e6a"
    )
    assert solvated["chain_selection"] == {
        "raw_group_A": ["A", "B"],
        "raw_group_B": ["I"],
        "canonical_group_A": ["A", "B"],
        "canonical_group_B": ["I"],
        "chain_map": {"A": "A", "B": "B", "I": "I"},
        "selection_status": "PASS",
        "selection_source": (
            "GAUSS chain-restored mapping from corrected protein blocks to "
            "authoritative 2KAI topology"
        ),
    }
    assert solvated["provenance"]["source_sha256"] == (
        "3eff4f2e5eaa3cc2d8e507ce50534d97efc7ba18bc24ee4e5faa1b9f91f228fc"
    )
