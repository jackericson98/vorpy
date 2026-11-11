"""Round 7 bounded public export and checkpoint regressions."""

import os
import subprocess
import sys
from pathlib import Path

import pytest


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


def _run_cli(source, output, method, *export_tokens):
    root = Path(__file__).resolve().parents[2]
    command = [
        sys.executable,
        "-m",
        "vorpy",
        str(source),
        "-s",
        "nt",
        method,
        "and",
        "mv",
        "5",
        "--boundary",
        "shell",
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


@pytest.mark.parametrize("method", ("aw", "pow", "prm"))
def test_supported_small_cli_exports_are_method_specific_and_checkpointed(tmp_path, method):
    source = tmp_path / "four_groups.pdb"
    _write_four_atom_source(source)
    output = tmp_path / method / "small"

    completed = _run_cli(source, output, method, "small")

    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    assert "Export stage failed" not in completed.stdout
    result_root = output / source.stem
    system_root = result_root / "system"
    assert list(result_root.glob("*_checkpoint.vpy"))
    assert (system_root / "edges.off").is_file()
    assert (system_root / "surfs.off").is_file()
    assert (system_root / "system_info.txt").is_file()
    assert (system_root / "group_info.txt").is_file()
    assert (system_root / f"{method}_logs.csv").is_file()


def test_medium_and_none_cli_exports_keep_the_checkpoint_contract(tmp_path):
    source = tmp_path / "four_groups.pdb"
    _write_four_atom_source(source)

    medium = tmp_path / "medium"
    completed = _run_cli(source, medium, "aw", "medium")
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    medium_root = medium / source.stem
    assert list(medium_root.glob("*_checkpoint.vpy"))
    assert (medium_root / "system" / "verts.off").is_file()

    none = tmp_path / "none"
    completed = _run_cli(source, none, "aw", "none")
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    none_root = none / source.stem
    assert list(none_root.glob("*_checkpoint.vpy"))
    assert not list(none_root.glob("*.off"))
