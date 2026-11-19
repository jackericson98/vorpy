from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.atlas import cache_export


def _fake_network():
    balls = pd.DataFrame([
        {"loc": [0, 0, 0], "rad": 1.0, "system_num": 0, "stable_id": "system:0", "complete": True},
        {"loc": [2, 0, 0], "rad": 1.0, "system_num": 1, "stable_id": "system:1", "complete": True},
    ])
    vertices = pd.DataFrame([
        {"loc": [1, 1, 0], "balls": [0, 1], "edges": [0], "surfs": [0]},
        {"loc": [1, -1, 0], "balls": [0, 1], "edges": [0], "surfs": [0]},
    ])
    edges = pd.DataFrame([
        {"balls": [0, 1], "verts": [0, 1], "surfs": [0], "points": [[1, 1, 0], [1, -1, 0]], "length": 2.0},
    ])
    surfaces = pd.DataFrame([
        {"balls": [0, 1], "verts": [0, 1], "edges": [0],
         "points": [[1, 1, 0], [1, -1, 0], [1, 0, 1]], "tris": [[0, 1, 2]],
         "sa": 1.0, "func": None, "int_gauss_curv": 0.0},
    ])
    source_balls = pd.DataFrame([
        {"loc": [0, 0, 0], "rad": 1.0, "name": "CA", "res_name": "ALA", "res_seq": 1,
         "chain_name": "A", "element": "C"},
        {"loc": [2, 0, 0], "rad": 1.0, "name": "CA", "res_name": "ALA", "res_seq": 2,
         "chain_name": "B", "element": "C"},
    ])
    return SimpleNamespace(
        archive_id="fake-archive", source_network_id=None, settings={"net_type": "pow"},
        balls=balls, verts=vertices, edges=edges, surfs=surfaces,
        generator_stable_ids=("system:0", "system:1"), solvent_local_indices=(),
        solvent_cells_complete=False,
        sys=SimpleNamespace(balls=source_balls, groups=[SimpleNamespace(sa=1.0, vol=2.0, int_gauss_curv=0.0)]),
    )


def test_cache_export_does_not_rebuild_and_refuses_existing_output(tmp_path, monkeypatch):
    archive = tmp_path / "network.vpy"
    archive.write_bytes(b"synthetic archive identity")
    network = _fake_network()
    calls = []
    monkeypatch.setattr(cache_export, "load_network", lambda path: calls.append(path) or network)

    output = cache_export.export_atlas_archive(archive, tmp_path / "atlas")
    assert output.is_dir()
    assert len(calls) == 1
    for name in ("surfaces.csv", "edges.csv", "vertices.csv", "contacts.csv",
                 "interface_summary.csv", "validation_report.json", "atlas_manifest.json"):
        assert (output / name).is_file()
    assert (output / "surfaces.csv").read_text(encoding="utf-8").count("\n") == 2

    with pytest.raises(FileExistsError):
        cache_export.export_atlas_archive(archive, output)
