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
        max_vert=5,
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
        locs=np.zeros((4, 3)), rads=np.ones(4), max_vert=5, net_type="aw",
        check_ndxs=list(range(4)), my_group=[0], timing={}, start_time=time.perf_counter()
    )
    assert result[2] == [1.0]
    np.testing.assert_allclose(result[3][0], np.ones(3))
    assert result[4] == [0.5]


def test_edge_spatial_query_cache_reuses_box_and_candidate_queries(monkeypatch):
    from vorpy.src.network import fast

    calls = {"box": 0, "balls": 0}

    def fake_box_search(loc):
        calls["box"] += 1
        return tuple(loc)

    def fake_get_balls(cells, dist):
        calls["balls"] += 1
        return [1, 2]

    monkeypatch.setattr(fast, "box_search", fake_box_search)
    monkeypatch.setattr(fast, "get_balls", fake_get_balls)
    locs = np.zeros((3, 3))
    cache = {}
    first = fast._edge_spatial_query([2, 1, 0], locs, 0.45, cache)
    second = fast._edge_spatial_query([0, 2, 1], locs, 0.45, cache)
    assert first == second
    assert calls == {"box": 3, "balls": 1}


def test_squared_distance_preserves_distance_ordering():
    from vorpy.src.network.fast import _squared_distance

    reference = np.array([1.0, -2.0, 0.5])
    points = [np.array([4.0, -2.0, 0.5]), np.array([1.0, 0.0, 0.5]), np.array([1.0, -2.0, 2.5])]
    squared_order = sorted(range(len(points)), key=lambda i: _squared_distance(points[i], reference))
    euclidean_order = sorted(range(len(points)), key=lambda i: np.linalg.norm(points[i] - reference))
    assert squared_order == euclidean_order


def test_edge_surrounding_query_cache_reuses_arrays(monkeypatch):
    from vorpy.src.network import fast

    calls = {"balls": 0}

    monkeypatch.setattr(fast, "_edge_spatial_query",
                        lambda edge, locs, dist, cache: (
                            [], [0, 1, 2]
                        ))
    original = fast._edge_spatial_query
    # Count through a wrapper while retaining deterministic candidate output.
    def counted(edge, locs, dist, cache):
        calls["balls"] += 1
        return original(edge, locs, dist, cache)
    monkeypatch.setattr(fast, "_edge_spatial_query", counted)
    locs = np.zeros((3, 3))
    rads = np.ones(3)
    cache = {}
    first = fast._edge_surrounding_query([2, 1, 0], locs, rads, 5.0, cache)
    second = fast._edge_surrounding_query([0, 2, 1], locs, rads, 5.0, cache)
    assert first[0] == second[0] == [0, 1, 2]
    assert first[1] is second[1]
    assert first[2] is second[2]
    assert first[3] is second[3]
    assert calls["balls"] == 1


def test_cached_geometry_evaluates_each_key_once():
    from vorpy.src.network.fast import _cached_geometry

    calls = []
    cache = {}

    def calculate():
        calls.append(True)
        return (1, 2, 3)

    first = _cached_geometry(cache, ("aw", (0, 1, 2, 3)), calculate)
    second = _cached_geometry(cache, ("aw", (0, 1, 2, 3)), calculate)
    other = _cached_geometry(cache, ("flat", True, (0, 1, 2, 3)), calculate)
    assert first == second == other
    assert len(calls) == 2
