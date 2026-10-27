"""Round-2 acceptance checks for the canonical geometry checkpoint writer."""

from pathlib import Path

import pytest

from vorpy.src.output import checkpoint as checkpoint_writer
from vorpy.src.output.checkpoint import write_checkpoint
from vorpy.workbench.services.vorpy_backend import (
    VorPyBackend,
    VorPySolveSettings,
)


TINY_PDB = """\
ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C
ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C
END
"""


def test_atomic_checkpoint_replace_preserves_previous_valid_archive(tmp_path, monkeypatch):
    source = tmp_path / "tiny.pdb"
    source.write_text(TINY_PDB, encoding="ascii")
    solved = VorPyBackend(
        VorPySolveSettings(network_type="aw", boundary_mode="shell", max_vertices=5)
    ).solve(source, lambda *args: None, lambda: False)
    system = solved.export_group.net.sys
    destination = tmp_path / "checkpoint" / "tiny_checkpoint.vpy"
    destination.parent.mkdir()

    write_checkpoint(system, destination)
    previous = destination.read_bytes()

    def interrupted_replace(*args, **kwargs):
        raise OSError("synthetic interruption before checkpoint replacement")

    monkeypatch.setattr(checkpoint_writer.os, "replace", interrupted_replace)
    with pytest.raises(OSError, match="synthetic interruption"):
        write_checkpoint(system, destination)

    assert destination.read_bytes() == previous
    assert not list(destination.parent.glob(f".{destination.name}.*.tmp"))
