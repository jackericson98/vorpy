"""Round 8 geometry-profile export dependency regressions."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from vorpy.src.command.calculation_plan import validate_export_dependencies


def _write_four_atom_source(path):
    path.write_text(
        """ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA B   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA C   1       0.000   2.000   0.000  1.00  0.00           C
ATOM      4  CA  ALA D   1       0.000   0.000   2.000  1.00  0.00           C
END
""",
        encoding="ascii",
    )


def _run_cli(source, output, profile, *export_tokens):
    root = Path(__file__).resolve().parents[2]
    command = [
        sys.executable,
        "-m",
        "vorpy",
        str(source),
        "-s",
        "nt",
        "aw",
        "and",
        "mv",
        "5",
        "--boundary",
        "shell",
        "--profile",
        profile,
        "-e",
        *export_tokens,
        "-e",
        "dir",
        str(output),
        "--fail-if-exists",
    ]
    return subprocess.run(
        command,
        cwd=source.parent,
        capture_output=True,
        text=True,
        timeout=120,
        env={**os.environ, "PYTHONPATH": str(root)},
    )


def _checkpoint_root(output, source):
    return output / source.stem


def test_geometry_profile_rejects_analysis_preset_after_checkpoint(tmp_path):
    source = tmp_path / "geometry_profile.pdb"
    _write_four_atom_source(source)
    output = tmp_path / "small"

    completed = _run_cli(source, output, "geometry", "small")

    assert completed.returncode != 0
    combined = completed.stdout + "\n" + completed.stderr
    assert "geometry-only" in combined
    assert "analysis-dependent" in combined
    assert "KeyError" not in combined
    assert list(_checkpoint_root(output, source).glob("*_checkpoint.vpy"))


@pytest.mark.parametrize("preset", ("medium", "all"))
def test_geometry_profile_rejects_other_analysis_presets_cleanly(tmp_path, preset):
    source = tmp_path / f"geometry_{preset}.pdb"
    _write_four_atom_source(source)
    output = tmp_path / preset

    completed = _run_cli(source, output, "geometry", preset)

    assert completed.returncode != 0
    combined = completed.stdout + "\n" + completed.stderr
    assert "Use --profile full" in combined
    assert "Traceback" not in combined
    assert list(_checkpoint_root(output, source).glob("*_checkpoint.vpy"))


def test_geometry_profile_none_keeps_checkpoint_and_no_optional_meshes(tmp_path):
    source = tmp_path / "geometry_none.pdb"
    _write_four_atom_source(source)
    output = tmp_path / "none"

    completed = _run_cli(source, output, "geometry", "none")

    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    root = _checkpoint_root(output, source)
    assert list(root.glob("*_checkpoint.vpy"))
    assert not list(root.rglob("*.off"))


def test_geometry_profile_explicit_mesh_component_remains_available(tmp_path):
    source = tmp_path / "geometry_edges.pdb"
    _write_four_atom_source(source)
    output = tmp_path / "edges"

    completed = _run_cli(source, output, "geometry", "only", "edges")

    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    root = _checkpoint_root(output, source)
    assert list(root.glob("*_checkpoint.vpy"))
    assert (root / "system" / "edges.off").is_file()


def test_full_profile_keeps_small_preset_available(tmp_path):
    source = tmp_path / "full_small.pdb"
    _write_four_atom_source(source)
    output = tmp_path / "full"

    completed = _run_cli(source, output, "full", "small")

    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    root = _checkpoint_root(output, source)
    assert list(root.glob("*_checkpoint.vpy"))
    assert (root / "system" / "system_info.txt").is_file()


def test_profile_dependency_contract_is_explicit():
    with pytest.raises(ValueError, match="analysis-dependent"):
        validate_export_dependencies("geometry", ("edges", "info"))

    validate_export_dependencies("geometry", ("edges", "surfaces", "vertices"))
    validate_export_dependencies("geometry", ("info", "surfaces"), scope="interface")
    validate_export_dependencies("full", ("info", "logs"))
