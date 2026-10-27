"""Bounded AW/POW/PRM solve and source-independent replay acceptance tests."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from vorpy.src.output.checkpoint import write_checkpoint
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


TINY_PDB = """\
ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C
ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C
ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C
END
"""


@pytest.mark.parametrize("scheme", ("aw", "pow", "prm"))
@pytest.mark.timeout(90)
def test_round4_solver_checkpoint_replay_and_export_without_build(tmp_path, scheme):
    """Each supported scheme survives a fresh-process, source-independent replay."""
    source = tmp_path / "tiny.pdb"
    source.write_text(TINY_PDB, encoding="ascii")

    result = VorPyBackend(
        VorPySolveSettings(
            network_type=scheme,
            boundary_mode="shell",
            max_vertices=5,
            typed_group_selections={"tiny": (0,)},
            source_network_id=f"producer:round4:{scheme}",
        )
    ).solve(source, lambda *args: None, lambda: False)
    network = result.export_group.net
    expected = {
        "scheme": (network.settings or {}).get("net_type", scheme),
        "source_network_id": network.source_network_id,
        "counts": {
            name: len(getattr(network, name))
            for name in ("balls", "verts", "edges", "surfs")
        },
    }

    archive = tmp_path / f"{scheme}.vpy"
    write_checkpoint(network.sys, archive)
    archive_hash = archive.read_bytes()
    source.unlink()
    replay_cwd = tmp_path / "replay-cwd"
    replay_cwd.mkdir()
    output = tmp_path / "replay-output"
    root = Path(__file__).resolve().parents[2]
    child = r'''
import json
import sys
import warnings
from pathlib import Path

warnings.simplefilter("ignore", ResourceWarning)
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.network import Network
from vorpy.src.output.checkpoint import _load_checkpoint

archive = Path(sys.argv[1])
output = Path(sys.argv[2])
original_build = Network.build
def forbidden_build(*args, **kwargs):
    raise AssertionError("checkpoint replay called Network.build()")
Network.build = forbidden_build

container = _load_checkpoint(archive)
document = container.document
loaded = Command()._load_canonical_checkpoint(archive)
sys.argv = ["vorpy", str(archive), "-e", "dir", str(output), "-e", "medium"]
Command().run()
print("NEXUS_ROUND4=" + json.dumps({
    "scheme": (loaded.settings or {}).get("net_type"),
    "source_network_id": loaded.source_network_id,
    "counts": {
        name: len(getattr(loaded, name))
        for name in ("balls", "verts", "edges", "surfs")
    },
    "source_bytes": len(document.source.content),
    "group_ids": sorted({group_id for atom in document.atoms for group_id in atom.group_ids}),
    "export_files": sum(1 for path in output.rglob("*") if path.is_file()),
}))
'''
    env = {**os.environ, "PYTHONPATH": str(root)}
    replay = subprocess.run(
        [sys.executable, "-c", child, str(archive), str(output)],
        cwd=replay_cwd,
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )
    assert replay.returncode == 0, replay.stdout + "\n" + replay.stderr
    marker = next(
        line for line in replay.stdout.splitlines() if line.startswith("NEXUS_ROUND4=")
    )
    observed = json.loads(marker.split("=", 1)[1])
    assert observed["scheme"] == expected["scheme"]
    assert observed["source_network_id"] == expected["source_network_id"]
    assert observed["counts"] == expected["counts"]
    assert observed["source_bytes"] > 0
    assert observed["group_ids"] == ["tiny"]
    assert observed["export_files"] > 0
    assert archive.read_bytes() == archive_hash
