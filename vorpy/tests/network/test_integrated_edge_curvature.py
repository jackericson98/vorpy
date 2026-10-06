from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.network.network import Network
from vorpy.src.calculations import edge_curvature, edge_resolution_cache
from vorpy.src.calculations.edge_geometry import AdditivelyWeightedTrisectorBranch, LineEdgeGeometry


@pytest.mark.parametrize('curved', [False, True])
@pytest.mark.parametrize('verbose', [False, True])
def test_build_integrates_curvature_and_reuses_resolved_geometry(monkeypatch, curved, verbose):
    import vorpy.src.network.network as network_module
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    radii = np.array([1., 1.5, 2.]) if curved else np.ones(3)
    geometry = (AdditivelyWeightedTrisectorBranch(locations, radii, 2., 4., branch=1)
                if curved else LineEdgeGeometry([2., 2.5, -2.], [2., 2.5, 2.]))
    samples = [geometry.point(t) for t in np.linspace(geometry.t_min, geometry.t_max, 17)]
    net = SimpleNamespace(
        settings={'net_type': 'aw', 'surf_res': .2, 'verbose': verbose}, group=[0, 1],
        balls=pd.DataFrame({'loc': list(locations), 'rad': radii}),
        verts=pd.DataFrame({'loc': [samples[-1], samples[0]], 'dub': [0, 0]}, index=[8, 3]),
        edges=pd.DataFrame({'balls': [[0, 1, 2]], 'verts': [[8, 3]], 'surfs': [[0, 1, 2]]}, index=[5]),
        surfs=pd.DataFrame({'balls': [(0, 1), (0, 2), (1, 2)]}),
        update_progress=lambda *args: None,
    )
    monkeypatch.setattr(network_module, 'build_edge', lambda **kwargs: (samples, {}))
    separate = deepcopy(net)
    Network.build_edges(separate, compute_curvature=False)
    expected = edge_curvature.calculate_aw_network_edge_curvatures(separate)
    calls = []
    original = edge_resolution_cache.resolve_aw_network_edge
    def resolve(*args, **kwargs):
        calls.append(kwargs.get('edge_data'))
        return original(*args, **kwargs)
    monkeypatch.setattr(edge_resolution_cache, 'resolve_aw_network_edge', resolve)

    Network.build_edges(net)
    assert len(calls) == 1 and calls[0] is not None
    assert net.edges['int_mean_curv_by_ball'].iloc[0] == pytest.approx(expected[0][0], rel=1e-12)
    assert net.edges['int_gauss_curv_by_ball'].iloc[0] == pytest.approx(expected[1][0], rel=1e-12)
    assert net.edges['int_gauss_curv_by_face'].iloc[0] == pytest.approx(
        separate._aw_fused_edge_curvature_cache['gaussian_by_face'][0], rel=1e-12)
    assert edge_resolution_cache.get_aw_edge_geometry_cache(net)[1] is False
    assert len(calls) == 1
    Network.build_edge_mean_curvature(net)
    Network.build_edge_gaussian_curvature(net)
    assert len(calls) == 1
    assert net.edges['length'].iloc[0] == separate.edges['length'].iloc[0]
    # A rejected surface is removed and references renumbered by meshing.
    # Refresh cell G from cached face integrals, without repeating quadrature.
    net.surfs.drop(index=1, inplace=True)
    net.surfs.reset_index(drop=True, inplace=True)
    net.edges.at[5, 'surfs'] = [0, 1]
    def no_reintegration(*args, **kwargs):
        raise AssertionError('Surface rejection must not repeat edge quadrature')
    monkeypatch.setattr(edge_curvature, 'aw_edge_curvature_measures', no_reintegration)
    Network.build_edge_mean_curvature(net)
    Network.build_edge_gaussian_curvature(net)
    assert set(net.edges['int_gauss_curv_by_ball'].iloc[0]) == {1}
    assert net.edges['int_gauss_curv_by_ball'].iloc[0][1] == pytest.approx(expected[1][0][1])
    assert len(calls) == 1


def test_verbose_integration_keeps_compiled_kernel(monkeypatch):
    from vorpy.src.calculations import edge_geometry, edge_curvature_kernel
    original = edge_curvature_kernel.integrate_edge_samples
    calls = []
    def integrate(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(edge_curvature_kernel, 'integrate_edge_samples', integrate)
    timing = {}
    edge_geometry.aw_edge_curvature_measures(
        LineEdgeGeometry([2., 2.5, -2.], [2., 2.5, 2.]),
        [0, 1, 2], [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
        [(0, 1)], timing=timing,
    )
    assert calls == [1]
    assert timing['compiled_accumulation'] > 0


@pytest.mark.parametrize('reverse', [False, True])
def test_batched_branch_geometry_matches_scalar_reference(reverse):
    from vorpy.src.calculations.edge_geometry import AffineEdgeGeometry, aw_edge_curvature_measures
    locations = np.array([[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]])
    branch = AdditivelyWeightedTrisectorBranch(locations, [1., 1.5, 2.], 2., 4., branch=1)
    edge = AffineEdgeGeometry(branch, 4. if reverse else 2., 2. if reverse else 4.)
    parameters = np.linspace(.01, .99, 32)
    points, firsts, seconds = edge.samples(parameters)
    for actual, method in [(points, edge.point), (firsts, edge.tangent), (seconds, edge.second_derivative)]:
        np.testing.assert_allclose(actual, [method(t) for t in parameters], rtol=1e-12, atol=1e-12)
    pairs = [(i, j) for i in range(3) for j in range(3) if i != j]
    actual = aw_edge_curvature_measures(edge, [0, 1, 2], locations, pairs, timing={})
    expected = aw_edge_curvature_measures(edge, [0, 1, 2], locations, pairs, reference=True)
    for measure in ('mean', 'gaussian'):
        assert actual[measure] == pytest.approx(expected[measure], rel=1e-12, abs=1e-12)
