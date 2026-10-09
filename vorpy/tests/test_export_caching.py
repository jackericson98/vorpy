import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.output.curvature_colors import curvature_color_limit, export_color_cache
from vorpy.src.output.draw import draw_joint
from vorpy.src.output.output import _group_export_cache


@pytest.mark.parametrize('preset', ['micro', 'tiny', 'medium', 'large', 'all'])
def test_every_export_preset_saves_network_archive(tmp_path, monkeypatch, preset):
    output = importlib.import_module('vorpy.src.output.output')
    io = importlib.import_module('vorpy.src.io')
    saved = []
    monkeypatch.setattr(io, 'save_network', lambda network, path: saved.append((network, path)))
    monkeypatch.setattr(output, '_move_vert_file', lambda *args: None)
    monkeypatch.setattr(output, '_export_nonpolar_geometry', lambda **kwargs: None)
    group = SimpleNamespace(name='group', net=SimpleNamespace(), settings={}, exports=lambda **kwargs: None)
    system = SimpleNamespace(name='example', groups=[group], ifaces=[], files={'dir': str(tmp_path)},
                             exports=lambda **kwargs: None, update_progress=lambda **kwargs: None)
    output.export_preset(system, preset)
    assert saved == [(system, str(tmp_path / 'example.vpy'))]
    assert 'network archive' in system.export_timing
    saved.clear()
    system._export_skip_archive = True
    output.export_preset(system, preset)
    assert saved == []


@pytest.mark.parametrize('preset', ['small', 'medium', 'all'])
def test_interface_log_presets_also_write_canonical_results_log(tmp_path, monkeypatch, preset):
    output = importlib.import_module('vorpy.src.output.output')
    written = []
    monkeypatch.setattr(output, '_export_canonical_interface_log', lambda iface: written.append(iface))
    monkeypatch.setattr(output, '_export_dual_geometry_without_launcher', lambda **kwargs: None)
    monkeypatch.setattr(output, '_export_nonpolar_geometry', lambda **kwargs: None)

    interface = SimpleNamespace(
        name='A_B', dir=str(tmp_path / 'interface_A_B'),
        net=SimpleNamespace(settings={'net_type': 'aw'}),
        export=lambda **kwargs: None,
    )
    system = SimpleNamespace(
        name='example', groups=[], ifaces=[interface], files={'dir': str(tmp_path)},
        exports=lambda **kwargs: None, update_progress=lambda **kwargs: None,
        _export_skip_archive=True,
    )

    output.export_preset(system, preset)

    assert written == [interface]


@pytest.mark.parametrize('commands, expected', [
    ([['large'], ['medium']], 1),
    ([['large'], ['no_archive']], 0),
    ([['none']], 0),
    ([['none'], ['logs']], 0),
    ([['only', 'logs']], 0),
])
def test_archive_cli_policy_and_deduplication(tmp_path, monkeypatch, commands, expected):
    output = importlib.import_module('vorpy.src.output.output')
    command = importlib.import_module('vorpy.src.command.command_export')
    io = importlib.import_module('vorpy.src.io')
    saved = []
    monkeypatch.setattr(io, 'save_network', lambda network, path: saved.append(path))
    monkeypatch.setattr(output, '_move_vert_file', lambda *args: None)
    monkeypatch.setattr(output, '_export_nonpolar_geometry', lambda **kwargs: None)
    group = SimpleNamespace(name='group', net=SimpleNamespace(), settings={}, exports=lambda **kwargs: None)
    system = SimpleNamespace(name='example', groups=[group], ifaces=[], files={'dir': str(tmp_path)},
                             exports=lambda **kwargs: None, update_progress=lambda **kwargs: None)
    command.argv_export(system, commands)
    assert len(saved) == expected
    assert not hasattr(system, '_export_skip_archive')
    assert not hasattr(system, '_export_archive_written')


def test_surface_parts_combine_without_retriangulation_and_refresh_cache(tmp_path):
    from vorpy.src.output.surfs import prepare_surfs, write_surfs
    net = _colored_network()
    net.surfs = pd.concat([net.surfs, net.surfs], ignore_index=True)
    with export_color_cache(net):
        mesh = prepare_surfs(net, [0, 1], target_cells=[0], color_limit=3.)
        assert mesh.triangles.tolist() == [[0, 1, 2], [3, 4, 5]]
        assert len(mesh.points) == 6
        assert set(net._export_surface_rows) == {0, 1}
        path = write_surfs(net, [0, 1], 'cached', directory=tmp_path, target_cells=[0], color_limit=3.)
    assert not hasattr(net, '_export_surface_rows')
    reference = write_surfs(net, [0, 1], 'reference', directory=tmp_path, target_cells=[0], color_limit=3.)
    assert path.read_bytes() == reference.read_bytes()
    net.surfs.at[0, 'points'] = [[2., 0., 0.], [3., 0., 0.], [2., 1., 0.]]
    with export_color_cache(net):
        mesh = prepare_surfs(net, [0], target_cells=[0], color_limit=3., include_face_data=False)
        assert mesh.points[0].tolist() == [2., 0., 0.]
        assert mesh.face_data == {}


def test_export_color_scales_are_refreshed_between_plans():
    net = SimpleNamespace(
        surfs=pd.DataFrame([{"balls": [0, 1], "int_mean_curv_by_ball": {0: 2., 1: -2.}}]),
        edges=pd.DataFrame(), verts=pd.DataFrame(),
    )
    with export_color_cache(net):
        assert curvature_color_limit(net, "int_mean_curv", [0]) == 2.
        net.surfs.at[0, "int_mean_curv_by_ball"] = {0: 7., 1: -7.}
        assert curvature_color_limit(net, "int_mean_curv", [0, 0]) == 2.
        assert curvature_color_limit(net, "int_mean_curv", [1]) == 7.
    assert not hasattr(net, "_export_color_limits")
    assert curvature_color_limit(net, "int_mean_curv", [0]) == 7.


def test_export_summary_is_reused_and_cache_cleared_on_failure(monkeypatch):
    module = importlib.import_module("vorpy.src.group.group")
    group = SimpleNamespace(net=SimpleNamespace())
    calls = []
    monkeypatch.setattr(module, "get_info", lambda grp: calls.append(grp))
    with pytest.raises(RuntimeError):
        with _group_export_cache(group):
            module.Group.get_info(group)
            module.Group.get_info(group)
            assert len(calls) == 1
            raise RuntimeError("Export failed")
    assert not hasattr(group, "_export_info_ready")
    assert not hasattr(group.net, "_export_color_limits")
    module.Group.get_info(group)
    assert len(calls) == 2


@pytest.mark.parametrize("subdivisions", [0, 1, 2])
def test_joint_template_preserves_radius_translation_and_independent_results(subdivisions):
    points, triangles = draw_joint([0., 0., 0.], radius=1., subdivisions=subdivisions)
    moved, moved_triangles = draw_joint([2., 3., 4.], radius=0.5, subdivisions=subdivisions)
    assert np.asarray(moved) == pytest.approx(np.asarray(points) * 0.5 + [2., 3., 4.])
    assert np.linalg.norm(points, axis=1) == pytest.approx(np.ones(len(points)))
    assert len(triangles) == 20 * 4 ** subdivisions
    assert triangles == moved_triangles
    triangles[0][0] = -1
    assert moved_triangles[0][0] >= 0


def _colored_network():
    return SimpleNamespace(
        settings={'net_type': 'aw', 'surf_scheme': 'int_mean_curv',
                  'surf_col': 'coolwarm', 'scheme_factor': 'log'},
        surfs=pd.DataFrame({'balls': [[0, 2]], 'points': [[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]],
                            'tris': [[[0, 1, 2]]], 'int_mean_curv_by_ball': [{0: 2., 2: -2.}]}),
        edges=pd.DataFrame({'balls': [[0, 2, 3]], 'points': [[[0., 0., 0.], [1., 0., 0.]]],
                            'int_mean_curv_by_ball': [{0: 1., 2: 2., 3: 3.}]}),
        verts=pd.DataFrame({'loc': [[0., 0., 0.], [1., 0., 0.]], 'edges': [[0], [0]]}),
    )


def test_target_selection_is_consumed_once_in_vertex_loop():
    from vorpy.src.output.verts import prepare_verts
    class SingleUse:
        used = False
        def __iter__(self):
            assert not self.used, 'Repeated construction of target-cell selection'
            self.used = True
            return iter(['0', 0, 1])
    mesh = prepare_verts(_colored_network(), [0, 1], target_cells=SingleUse(), color_limit=3.)
    assert len(mesh.triangles) == 40


def test_edge_vertex_row_cache_refresh_and_off_fields():
    from vorpy.src.output.edges import prepare_edges
    from vorpy.src.output.verts import prepare_verts
    net = _colored_network()
    with export_color_cache(net):
        # Vertex colors read edges before tube drawing columns are populated.
        vertex_mesh = prepare_verts(net, [0], target_cells=[0], color_limit=3., include_face_data=False)
        edge_mesh = prepare_edges(net, [0], target_cells=[0], color_limit=3., include_face_data=False)
        assert vertex_mesh.face_data == edge_mesh.face_data == {}
        row = net._export_component_rows[('edges', 0)]
        assert len(row['draw_points']) > 0
        repeated = prepare_edges(net, [0], target_cells=[0], color_limit=3.)
        np.testing.assert_array_equal(repeated.points, edge_mesh.points)
        assert set(repeated.face_data) == {'edge_index', 'geometry_kind'}
    assert not hasattr(net, '_export_component_rows')
    net.verts.at[0, 'loc'] = [4., 0., 0.]
    with export_color_cache(net):
        shifted = prepare_verts(net, [0], target_cells=[0], color_limit=3.)
        np.testing.assert_allclose(shifted.points, vertex_mesh.points + [4., 0., 0.])


def test_incident_edge_color_cache_is_scoped_to_plan_and_selection():
    from vorpy.src.output.curvature_colors import mean_vertex_display_value
    net = _colored_network()
    vertex = net.verts.iloc[0]
    with export_color_cache(net):
        assert mean_vertex_display_value(net, vertex, [0]) == 1.
        net.edges.at[0, 'int_mean_curv_by_ball'] = {0: 7., 2: 2., 3: 3.}
        assert mean_vertex_display_value(net, vertex, [0, 0]) == 1.
        assert mean_vertex_display_value(net, vertex, [2]) == 2.
    assert not hasattr(net, '_export_edge_color_values')
    assert mean_vertex_display_value(net, vertex, [0]) == 7.


@pytest.mark.parametrize('mode', ['boundary', 'magnitude', 'cell'])
@pytest.mark.parametrize('file_type', ['off', 'ply', 'vtp'])
def test_optimized_coloring_preserves_mesh_file_contents(tmp_path, monkeypatch, mode, file_type):
    from vorpy.src.output import curvature_colors as colors
    from vorpy.src.output import verts as vertices
    from vorpy.src.output import edges
    from vorpy.src.output import surfs
    net = _colored_network()
    writers = [('vertices', vertices.write_verts, [0, 1]),
               ('edges', edges.write_edges, [0]), ('surfaces', surfs.write_surfs, [0])]
    actual = {}
    with export_color_cache(net):
        for name, writer, indices in writers:
            path = writer(net, indices, name, directory=tmp_path, file_type=file_type,
                          target_cells=[0], color_mode=mode, color_limit=3.)
            actual[name] = path.read_bytes()
    # Reference behavior rebuilt a plain integer set on every scalar lookup.
    def legacy_targets(targets):
        return None if targets is None else {int(value) for value in targets}
    monkeypatch.setattr(colors, '_target_set', legacy_targets)
    monkeypatch.setattr(vertices, '_target_set', legacy_targets)
    monkeypatch.setattr(edges, '_target_set', legacy_targets)
    for name, writer, indices in writers:
        path = writer(net, indices, name + '_reference', directory=tmp_path, file_type=file_type,
                      target_cells=[0], color_mode=mode, color_limit=3.)
        assert path.read_bytes() == actual[name]


def test_export_timings_include_mesh_preparation_and_writing(tmp_path):
    from vorpy.src.output.output import ExportProgress, _run_export
    from vorpy.src.output.verts import write_verts
    net = _colored_network()
    class Owner:
        def __init__(self):
            self.net = net
        def export(self):
            return write_verts(net, [0], 'vertex', directory=tmp_path, color_limit=3., target_cells=[0])
    sys = SimpleNamespace(verbose=True, update_progress=lambda **kwargs: None)
    progress = ExportProgress(1, sys)
    with export_color_cache(net):
        _run_export(progress, 'vertex', Owner().export)
    progress.finish()
    assert sys.export_mesh_timing['vertex']['prepare'] > 0
    assert sys.export_mesh_timing['vertex']['write'] > 0


def test_group_indices_refresh_between_plans():
    from vorpy.src.group.export import _group_topology_indices
    net = SimpleNamespace(balls=pd.DataFrame({'num': [0, 1], 'system_num': [10, 11]}))
    group = SimpleNamespace(net=net, ball_ndxs=[11])
    with _group_export_cache(group):
        indices = _group_topology_indices(group)
        assert indices == [1]
        indices.clear()
        assert _group_topology_indices(group) == [1]
    assert not hasattr(group, '_export_topology_indices')
    group.ball_ndxs = [10]
    with _group_export_cache(group):
        assert _group_topology_indices(group) == [0]
