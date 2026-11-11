import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from vorpy.src.command.commands import print_help


@pytest.mark.parametrize("flag", ["-h", "--help"])
def test_module_help_from_dataset_free_directory(flag, tmp_path):
    repo_root = Path(__file__).parents[3]
    environment = dict(os.environ)
    environment["PYTHONPATH"] = (
        str(repo_root) + os.pathsep + environment.get("PYTHONPATH", "")
    )
    result = subprocess.run(
        [sys.executable, "-m", "vorpy", flag],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "Input file not found" not in result.stdout
    assert "Solving" not in result.stdout
    assert "Usage:" in result.stdout
    assert result.stderr == ""


def test_installed_console_entry_point_help_from_dataset_free_directory(tmp_path):
    console = shutil.which("vorpy")
    if console is None:
        pytest.skip("vorpy console entry point is not installed")
    repo_root = Path(__file__).parents[3]
    environment = dict(os.environ)
    environment["PYTHONPATH"] = (
        str(repo_root) + os.pathsep + environment.get("PYTHONPATH", "")
    )
    result = subprocess.run(
        [console, "--help"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "Input file not found" not in result.stdout
    assert "Solving" not in result.stdout
    assert "Usage:" in result.stdout
    assert result.stderr == ""


def test_help_lists_current_cli_options(capsys):
    print_help()
    output = capsys.readouterr().out

    for option in (
        "--profile",
        "--save-network",
        "--fail-if-exists",
        "--load-network",
        "--profile-cache-export",
        "--boundary",
        "--solvent-mode",
        "--diagnose-edges",
        "--power-interface",
        "--aw-alpha-experimental",
        "--export-dual-visualization",
    ):
        assert option in output

    assert "current CLI execution default for mv is 5 A" in output
    assert "0.5-5000 A" in output
    assert "Common examples:" in output
    assert "vorpy --help" in output
