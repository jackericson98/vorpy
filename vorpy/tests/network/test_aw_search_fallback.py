import importlib

import numpy as np


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
