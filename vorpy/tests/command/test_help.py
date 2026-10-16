import subprocess
import sys
from pathlib import Path

import pytest

from vorpy.src.command.commands import print_help


@pytest.mark.parametrize("flag", ["-h", "--help"])
def test_module_help_exits_cleanly_without_resolving_input(flag):
    result = subprocess.run(
        [sys.executable, "-m", "vorpy", flag],
        cwd=Path(__file__).parents[2],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "Input file not found" not in result.stdout
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
