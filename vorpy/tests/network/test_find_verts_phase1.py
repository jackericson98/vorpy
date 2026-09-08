import time
import numpy as np

from vorpy.src.network import find_verts as module


def test_find_verts_phase1_normalizes_inputs_and_records_timing(monkeypatch):
    seed = {"balls": [0, 1, 2, 3], "loc": np.array([0.0, 0.0, 0.0]), "rad": 1.0}
    monkeypatch.setattr(module, "find_v0", lambda **kwargs: seed)
    monkeypatch.setattr(module, "find_site_container", lambda **kwargs: None)
    timing = {}
    result = module.find_verts(
        locs=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
        rads=[1, 1, 1, 1],
        max_vert=40,
        net_type="aw",
        check_ndxs=[0, 1, 2, 3],
        my_group=[0],
        timing=timing, start_time=time.perf_counter(),
    )
    assert result is not None
    vert_ndxs, vlocs, vrads, _, _, remaining, _ = result
    assert vert_ndxs == [[0, 1, 2, 3]]
    assert len(vlocs) == len(vrads) == 1
    assert remaining == [0, 1, 2, 3]
    assert timing["edge_search_calls"] == 4
    assert timing["accepted_vertices"] == 0
    assert timing["site_container"] >= 0.0


def test_find_verts_phase1_preserves_vertex_state_shape(monkeypatch):
    seed = {"balls": [0, 1, 2, 3], "loc": np.zeros(3), "rad": 1.0, "loc2": np.ones(3), "rad2": 0.5}
    monkeypatch.setattr(module, "find_v0", lambda **kwargs: seed)
    monkeypatch.setattr(module, "find_site_container", lambda **kwargs: None)
    result = module.find_verts(
        locs=np.zeros((4, 3)), rads=np.ones(4), max_vert=40, net_type="aw",
        check_ndxs=list(range(4)), my_group=[0], timing={}, start_time=time.perf_counter()
    )
    assert result[2] == [1.0]
    np.testing.assert_allclose(result[3][0], np.ones(3))
    assert result[4] == [0.5]
