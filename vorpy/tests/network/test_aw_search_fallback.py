import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd

from vorpy.src.network.find_net_verts import (
    _load_seed_vertex_state,
    _vertex_seed_dependency,
)


def test_bulk_seed_adjacency_is_sorted_for_binary_vertex_lookup():
    net = SimpleNamespace(
        balls=pd.DataFrame({
            'loc': [[float(index), 0., 0.] for index in range(5)],
            'rad': [1.] * 5,
        }),
        generator_stable_ids=tuple(f'atom-{index}' for index in range(5)),
        settings={'net_type': 'aw', 'max_vert': 5., 'box_size': 1.},
        group=[0], iface_grps=None,
    )
    net.seed_vertex_state = [
        {'balls': [4, 3, 2, 0], 'loc': [0., 0., 0.], 'rad': 1.},
        {'balls': [4, 3, 1, 0], 'loc': [1., 0., 0.], 'rad': 1.},
    ]
    net.seed_vertex_state_dependency = _vertex_seed_dependency(net)

    vert_ndxs, _, _, _, _, b_verts = _load_seed_vertex_state(net)

    assert vert_ndxs == [[0, 2, 3, 4], [0, 1, 3, 4]]
    assert b_verts[0] == [1, 0]
    assert [vert_ndxs[index] for index in b_verts[0]] == sorted(vert_ndxs)


def test_vertex_search_traverses_supplied_seed_before_finding_v0(monkeypatch):
    module = importlib.import_module('vorpy.src.network.find_verts')
    monkeypatch.setattr(module, 'find_v0', lambda *args, **kwargs:
                        (_ for _ in ()).throw(AssertionError('unexpected seed search')))
    monkeypatch.setattr(module, 'find_site_container', lambda *args, **kwargs: None)
    result = module.find_verts(
        locs=np.zeros((4, 3)), rads=np.ones(4), max_vert=5., net_type='pow',
        check_ndxs=[3], my_group=[3], vert_ndxs=[[0, 1, 2, 3]],
        vlocs=[[0., 0., 0.]], vrads=[0.5], vloc2s=[[None, None, None]],
        vrad2s=[None], b_verts=[[0], [0], [0], [0]], traverse_existing=True,
    )
    assert result[0] == [[0, 1, 2, 3]]
    assert result[5] == []


def test_aw_search_checks_remaining_candidates_after_hull_neighbors_fail(monkeypatch):
    fast = importlib.import_module("vorpy.src.network.fast")
    locs = np.array([
        [0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
        [1., 1., 1.], [1., 1., 2.], [1., 1., 3.], [1., 1., -1.],
    ])
    monkeypatch.setattr(fast, "_edge_spatial_query",
                        lambda *args: ([], [3, 4, 5, 6]))
    monkeypatch.setattr(fast, "calc_vert_numba_indexed",
                        lambda centers, radii, balls: (centers[balls[-1]], 0.5, None, None))
    monkeypatch.setattr(fast, "calc_circ", lambda *args: (np.zeros(3), 0.0))
    checked = []

    def choose(candidate, *args, **kwargs):
        fourth = candidate["balls"][-1]
        checked.append(fourth)
        return (candidate, None) if fourth == 4 else (None, fourth)

    monkeypatch.setattr(fast, "choose_vert", choose)
    result, _ = fast.find_site_aw(
        edge_balls=[0, 1, 2], locs=locs, rads=np.ones(7),
        b_verts=[[] for _ in range(7)], vert_ndxs=[],
        max_vert=40, mv_inc=4.5, check_ndxs=False, surr_balls=None,
        my_boxes=[], invalid_ndxs=[], vn_1=[0, 1, 2, 7],
        vn_1_loc=np.zeros(3),
    )
    assert checked == [6, 3, 4]
    assert result["balls"] == [0, 1, 2, 4]
