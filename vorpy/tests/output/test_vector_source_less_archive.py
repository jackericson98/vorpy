import json
import sys

import pytest

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.io import load_network, save_network
from vorpy.src.network import Network
from vorpy.src.output.checkpoint import (
    CheckpointSourceUnavailable,
    write_checkpoint,
)
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


def _solve(tmp_path):
    source = tmp_path / "legacy_source.pdb"
    source.write_text(
        "ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C\n"
        "ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C\n"
        "END\n",
        encoding="ascii",
    )
    result = VorPyBackend(
        VorPySolveSettings(network_type="aw", boundary_mode="shell", max_vertices=5)
    ).solve(source, lambda *args: None, lambda: False)
    return source, result.export_group.net


def test_legacy_source_less_cache_export_is_explicit_and_build_free(tmp_path, monkeypatch):
    source, network = _solve(tmp_path)
    archive = save_network(network, tmp_path / "legacy.vpy")
    archive_before = archive.read_bytes()
    source.unlink()
    output = tmp_path / "cache-output"

    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("cache-only export rebuilt the network"),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["vorpy", str(archive), "-e", "dir", str(output), "-e", "tiny"],
    )

    command = Command()
    command.run()

    output_root = output / network.sys.name
    status_path = output_root / "system" / "checkpoint_status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    assert status["status"] == "NOT_WRITTEN"
    assert status["reason"] == "SOURCE_BYTES_UNAVAILABLE"
    assert status["source"]["exact_bytes_available"] is False
    assert status["source"]["sha256"] is None
    assert len(status["source"]["original_recorded_sha256"]) == 64
    assert status["scientific_archive_reloadable"] is True
    assert status["checkpoint"]["written"] is False
    assert not (output_root / "legacy_source_checkpoint.vpy").exists()
    assert not (output_root / "legacy_source.pdb").exists()
    assert list(output_root.rglob("group_info.txt"))
    assert archive.read_bytes() == archive_before

    restored = load_network(archive)
    assert len(restored.verts) == len(network.verts)
    assert len(restored.edges) == len(network.edges)
    assert len(restored.surfs) == len(network.surfs)


def test_source_loss_does_not_replace_existing_checkpoint(tmp_path):
    source, network = _solve(tmp_path)
    output = tmp_path / "checkpoint-output"
    output.mkdir()
    network.sys.files["dir"] = str(output)
    checkpoint = write_checkpoint(network.sys)
    before = checkpoint.read_bytes()

    source.unlink()
    with pytest.raises(CheckpointSourceUnavailable, match="Exact original source bytes"):
        write_checkpoint(network.sys)

    assert checkpoint.read_bytes() == before
