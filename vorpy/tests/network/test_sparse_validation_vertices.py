from types import SimpleNamespace

import numpy as np
import pandas as pd

from vorpy.src.network.find_net_verts import _enumerate_sparse_aw_vertices


def test_sparse_seven_generator_box_fallback_finds_eight_vertices():
    locs = np.array([
        [0, 0, 0], [20, 0, 0], [-20, 0, 0],
        [0, 20, 0], [0, -20, 0], [0, 0, 20], [0, 0, -20],
    ], dtype=float)
    net = SimpleNamespace(
        settings={"net_type": "aw", "max_vert": 40.0},
        balls=pd.DataFrame({"loc": list(locs), "rad": [1.0] * 7}),
        group=list(range(7)),
        iface_grps=None,
    )
    state = _enumerate_sparse_aw_vertices(net)
    assert state is not None
    vertices = np.asarray(state[1], dtype=float)
    assert len(vertices) == 8
    assert np.unique(np.round(vertices, 8), axis=0).shape[0] == 8
    assert np.allclose(np.unique(np.abs(vertices), axis=0), [[10, 10, 10]])

