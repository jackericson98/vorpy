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


def _collapsed_network(length=0.0):
    return SimpleNamespace(
        settings={'net_type': 'aw'}, group=[0, 1, 2],
        balls=pd.DataFrame({'loc': [[0., 0., 0.], [4., 0., 0.], [0., 4., 0.]],
                            'rad': [1., 1., 1.]}),
        verts=pd.DataFrame({'loc': [[2., 2., 0.], [2., 2., length]]}),
        edges=pd.DataFrame([{'balls': [0, 1, 2], 'verts': [0, 1], 'surfs': [0, 1, 2],
                             'points': [[2., 2., 0.], [2., 2., length]]}]),
        surfs=pd.DataFrame({'balls': [(0, 1), (0, 2), (1, 2)]}),
        update_progress=lambda *args: None,
    )


@pytest.mark.parametrize('length', [0.0, 3.7e-16])
def test_collapsed_straight_edge_has_explicit_zero_line_integrals(length):
    from vorpy.src.calculations.edge_resolution_cache import get_aw_edge_geometry_cache
    net = _collapsed_network(length)
    mean, gaussian = edge_curvature.calculate_aw_network_edge_curvatures(net)
    assert mean == [{0: 0., 1: 0., 2: 0.}]
    assert gaussian == [{0: 0., 1: 0., 2: 0.}]
    cache = net._aw_fused_edge_curvature_cache
    assert cache['gaussian_by_face'] == [{(i, j): 0. for i in range(3) for j in range(3) if i != j}]
    # A zero line integral does not manufacture a tangent for vertex angles.
    assert get_aw_edge_geometry_cache(net)[0][0] is None


def test_noncollapsed_resolution_failure_is_not_silently_zeroed():
    net = _collapsed_network(1.0)
    net.edges.at[0, 'points'] = []
    with pytest.raises(ValueError, match='Unable to resolve AW edge 0.*missing_samples'):
        edge_curvature.calculate_aw_network_edge_curvatures(net)


def test_coincident_unequal_radius_edge_is_not_assumed_empty(monkeypatch):
    from vorpy.src.calculations import edge_resolution_cache
    net = _collapsed_network()
    net.balls.loc[0, 'rad'] = 2.0
    def failed(*args, **kwargs):
        raise ValueError('unresolved curve')
    monkeypatch.setattr(edge_resolution_cache, 'resolve_aw_network_edge', failed)
    with pytest.raises(ValueError, match='Unable to resolve AW edge 0.*unresolved curve'):
        edge_curvature.calculate_aw_network_edge_curvatures(net)
