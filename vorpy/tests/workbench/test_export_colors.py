from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from matplotlib import colormaps
from matplotlib.colors import to_rgb

from vorpy.workbench.domain import AnalysisResult
from vorpy.workbench.services.vorpy_backend import _layers_from_network
from vorpy.workbench.services.export_colors import make_color_provider
from vorpy.workbench.services.display_colors import DisplayColors
from vorpy.src.output.surfs import prepare_surfs
from vorpy.src.output.edges import prepare_edges
from vorpy.src.output.verts import prepare_verts


def result_with_geometry():
    net = SimpleNamespace(
        settings={'net_type': 'pow', 'surf_scheme': 'solid'},
        edges=pd.DataFrame({'points': [np.array([[0., 0., 0.], [1., 0., 0.]])]}),
        verts=pd.DataFrame({'loc': [np.array([0., 0., 0.])]}),
        surfs=pd.DataFrame({'points': [np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])],
                            'tris': [np.array([[0, 1, 2]])],
                            'mean_tri_curvs': [np.array([.25])], 'gauss_tri_curvs': [np.array([.125])]}),
    )
    group = SimpleNamespace(net=net, layer_surfs=[[0]], layer_edges=[[0]], layer_verts=[[0]])
    return AnalysisResult(None, 'sample', layers=_layers_from_network(net, [0], [0], [0]), export_group=group)


def test_export_solid_colors_follow_each_viewer_layer():
    result = result_with_geometry()
    for index, layer in enumerate(result.layers):
        layer.color = ['#112233', '#445566', '#778899', '#aabbcc', '#123456', '#abcdef'][index]
    provider = make_color_provider(result, {})
    net = result.export_group.net
    net._export_color_provider = provider
    for kind, prepare in [('surfaces', prepare_surfs), ('edges', prepare_edges), ('vertices', prepare_verts)]:
        for mode in ('magnitude', 'boundary', 'cell'):
            mesh = prepare(net, [0], color_scheme='solid', color_mode=mode)
            layer = next(layer for layer in result.layers if layer.kind == kind and
                         layer.interpretation == ('boundary' if mode == 'boundary' else 'magnitude'))
            assert mesh.face_colors == pytest.approx(np.tile(to_rgb(layer.color), (len(mesh.triangles), 1)))


def test_scalar_colors_and_overrides_do_not_modify_viewer():
    result = result_with_geometry()
    for layer in result.layers:
        layer.color_scheme = 'integrated_mean_curvature'
        layer.color_map = 'viridis'
        layer.scale_mode = 'log100'
        layer.cell_scalars[layer.color_scheme][:] = -.5
    provider = make_color_provider(result, {})
    for layer in result.layers:
        values, (low, high) = DisplayColors(result.layers)._display_scalars(layer, layer.cell_scalars[layer.color_scheme])
        expected = colormaps[layer.color_map]((values[0] - low) / (high - low))[:3]
        assert provider(layer.kind, 0, layer.interpretation, 1)[0] == pytest.approx(expected)
    override = make_color_provider(result, {'edges': '#ff0000'})
    assert override('edges', 0, 'boundary', 1)[0] == pytest.approx([1., 0., 0.])
    assert all(layer.color_scheme == 'integrated_mean_curvature' for layer in result.layers)
