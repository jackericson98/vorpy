from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.calculations import edge_curvature
from vorpy.src.calculations.edge_geometry import LineEdgeGeometry


@pytest.mark.parametrize('verbose', [False, True])
def test_network_curvature_preserves_short_regular_edges(monkeypatch, verbose):
    length = 5e-7
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 4., 0.]])
    geometry = LineEdgeGeometry([2., 2., 0.], [2., 2., length])
    resolved = SimpleNamespace(
        geometry=geometry, ball_indices=(0, 1, 2), locations=locations,
    )
    net = SimpleNamespace(
        settings={'net_type': 'aw', 'verbose': verbose},
        group=[0, 1, 2],
        edges=pd.DataFrame([{'surfs': [0, 1, 2]}]),
        surfs=pd.DataFrame({'balls': [(0, 1), (0, 2), (1, 2)]}),
        update_progress=lambda *args: None,
    )
    monkeypatch.setattr(edge_curvature, 'get_aw_edge_geometry_cache',
                        lambda *args, **kwargs: ({0: resolved}, False, 0.))

    mean, gaussian = edge_curvature.calculate_aw_network_edge_curvatures(net)

    # Equal-radius generators give planar faces and a straight trisector.
    assert mean[0] == pytest.approx(
        {0: np.pi * length / 4, 1: np.pi * length / 8, 2: np.pi * length / 8},
        rel=1e-12, abs=1e-20,
    )
    assert gaussian[0] == {0: 0., 1: 0., 2: 0.}
