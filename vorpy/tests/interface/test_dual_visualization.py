"""Normal AW interface exports consume the same solved/cached geometry."""
import csv
import importlib
import json
from dataclasses import replace
from io import StringIO
from types import ModuleType, SimpleNamespace
import sys

import numpy as np
import pytest

from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.filtrations.alpha import AlphaSimplex
from vorpy.src.geometry.interfaces import build_alpha_interface
from vorpy.src.geometry.visualization.interface import export_interface_dual_visualization, write_interface_dual_summary
from vorpy.src.interface.interface import Interface
from vorpy.tests.geometry.test_alpha_interfaces import _network, _filtration


def cached_interface(path, *, unresolved=False, reverse=False):
    net = _network('aw')
    net.surfs['points'] = [np.asarray([[1., 0., 0.], [1., 3., 0.], [1., 3., 1.], [1., 0., 1.]])] * 3
    net.surfs['tris'] = [np.asarray([[0, 1, 2], [0, 2, 3]])] * 3
    net.surfs['int_gauss_curv'] = [.125, .25, .5]
    net.edges['points'] = [np.asarray([[1., 1., -1.], [1., 1., 1.]])]
    dual = build_dual(net)
    filtration = _filtration(net, 'aw')
    if unresolved:
        pair = (10, 16)
        filtration.records[1][pair] = AlphaSimplex(1, pair, None, None, 'unresolved_aw', False, notes='restricted birth unresolved')
        filtration.blocked.append(filtration.records[1][pair])
    a, b = ({13, 16}, {10}) if reverse else ({10}, {13, 16})
    rep = build_alpha_interface(net, dual, filtration, 0., a, b, calculate_curvature=False)
    rep.surfaces = [replace(row, smooth_H_A=-.25) if row.selected else row for row in rep.surfaces]
    iface = Interface.__new__(Interface)
    iface.net, iface.dir, iface.name = net, str(path), 'A_B'
    iface.group1, iface.group2 = SimpleNamespace(name='A', group_id='A'), SimpleNamespace(name='B', group_id='B')
    iface._geometry_dual, iface._geometry_filtration, iface._geometry_representation = dual, filtration, rep
    selected = [list(record.generator_tuple) for record in filtration.interface_at(0., a, b)]
    iface.geometry_analysis = SimpleNamespace(metadata={'alpha_value': 0.}, alpha_selection={
        'selected_generator_pairs': selected, 'eligible_bicolor_pairs': 2})
    return iface


def rows(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def forbid_recalculation(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Visualization attempted geometry/filtration/curvature recalculation')
    monkeypatch.setattr('vorpy.src.network.network.Network.build', forbidden)
    monkeypatch.setattr('vorpy.src.geometry.duals.build_dual', forbidden)
    monkeypatch.setattr('vorpy.src.geometry.filtrations.alpha.build_alpha_filtration', forbidden)
    monkeypatch.setattr('vorpy.src.geometry.interfaces.build_alpha_interface', forbidden)
    monkeypatch.setattr(Interface, 'analyze_geometry', forbidden)


def test_normal_export_uses_cached_indices_selection_and_curvature(tmp_path, monkeypatch):
    iface = cached_interface(tmp_path)
    before_points = [points.copy() for points in iface.net.surfs['points']]
    forbid_recalculation(monkeypatch)
    iface.export(dual=True)
    root = tmp_path / 'dual'
    mapping = rows(root / 'dual_interface_mapping.csv')
    assert len(mapping) == 1
    record = mapping[0]
    assert (record['generator_1'], record['generator_2'], record['physical_surface_id']) == ('10', '13', '20')
    assert record['dual_feature_id'] == '1:10,13'
    assert record['stable_atom_1'] == 'A|1|GLY|CA'
    assert record['stable_atom_2'] == 'B|2|ALA|CB'
    atom = (root / 'molecule_A.pdb').read_text().splitlines()[0]
    assert atom[12:16].strip() == 'CA' and atom[17:20] == 'GLY' and atom[21] == 'A'
    identity = rows(root / 'generator_identity.csv')
    assert identity[0]['stable_atom_id'] == 'A|1|GLY|CA'
    assert float(record['smooth_integrated_H']) == -.25
    assert float(record['smooth_integrated_K']) == .125
    selected = rows(root / 'apollonius_selected_edges.csv')
    assert {row['simplex_id'] for row in selected} == {row['dual_feature_id'] for row in mapping}
    vertices = rows(root / 'aw' / 'dual' / 'dual_full_vertices.csv')
    assert {int(row['generator_0']) for row in vertices} == set(iface.net.balls.index)
    contacts = json.loads((root / 'contacts.json').read_text())
    np.testing.assert_array_equal(contacts['1:10,13']['points'], iface.net.surfs.loc[20, 'points'])
    np.testing.assert_array_equal(contacts['1:10,13']['triangles'], iface.net.surfs.loc[20, 'tris'])
    summary = json.loads((root / 'interface_dual_summary.json').read_text())
    assert summary['selected_pairs'] == summary['mapped_selected_surfaces'] == 1
    assert summary['certified_rejected_pairs'] == 1 and summary['unresolved_pairs'] == 0
    assert summary['additional_network_solves'] == 0
    for before, after in zip(before_points, iface.net.surfs['points']):
        np.testing.assert_array_equal(before, after)


def test_unresolved_aw_pair_is_never_exported_as_rejected(tmp_path, monkeypatch):
    iface = cached_interface(tmp_path, unresolved=True)
    forbid_recalculation(monkeypatch)
    result = export_interface_dual_visualization(iface)
    assert result['unresolved_pairs'] == 1 and result['certified_rejected_pairs'] == 0
    statuses = rows(tmp_path / 'dual' / 'dual_pair_status.csv')
    unresolved = next(row for row in statuses if row['generator_2'] == '16')
    assert unresolved['restriction_status'] == 'unresolved' and unresolved['birth_value'] == ''
    assert (tmp_path / 'dual' / 'apollonius_unresolved_edges.off').exists()


def test_cached_lower_bound_exclusion_does_not_invent_a_birth(tmp_path):
    iface = cached_interface(tmp_path, unresolved=True)
    export_interface_dual_visualization(iface)
    assert (tmp_path / 'dual' / 'apollonius_unresolved_edges.off').exists()
    pair = (10, 16)
    record = replace(iface._geometry_filtration.records[1][pair], diagnostics={
        'restriction_status': 'rejected', 'restriction_certificate': 'global_pair_minimum_lower_bound',
        'birth_lower_bound_A': 1.})
    iface._geometry_filtration.records[1][pair] = record
    result = export_interface_dual_visualization(iface)
    assert result['certified_rejected_pairs'] == 1 and result['unresolved_pairs'] == 0
    assert not (tmp_path / 'dual' / 'apollonius_unresolved_edges.off').exists()
    row = next(row for row in rows(tmp_path / 'dual' / 'dual_pair_status.csv') if row['generator_2'] == '16')
    assert row['birth_value'] == '' and row['birth_lower_bound'] == '1.0'
    assert row['certification_status'] == 'certified_lower_bound_exclusion'


@pytest.mark.parametrize('preset, expected', [('micro', False), ('tiny', False), ('medium', False), ('large', True), ('all', True)])
def test_presets_include_dual_only_in_large_and_all(tmp_path, monkeypatch, preset, expected):
    output = importlib.import_module('vorpy.src.output.output')
    iface = cached_interface(tmp_path)
    iface.dir = None
    calls = []
    iface.export = lambda **kwargs: calls.append(kwargs)
    monkeypatch.setattr(output, '_export_nonpolar_geometry', lambda **kwargs: None)
    system = SimpleNamespace(groups=[], ifaces=[iface], exports=lambda **kwargs: None,
        update_progress=lambda **kwargs: None, _export_skip_archive=True)
    output.export_preset(system, preset)
    assert any(call.get('dual') for call in calls) is expected


def test_export_refuses_missing_or_disagreeing_cache_without_rebuilding(tmp_path, monkeypatch):
    iface = cached_interface(tmp_path)
    forbid_recalculation(monkeypatch)
    iface.geometry_analysis.alpha_selection['selected_generator_pairs'] = []
    with pytest.raises(ValueError, match='disagree'):
        iface.export(dual=True)
    iface._geometry_representation = None
    with pytest.raises(ValueError, match='cached geometry analysis'):
        iface.export(dual=True)


def test_wrong_physical_surface_mapping_is_rejected(tmp_path):
    iface = cached_interface(tmp_path)
    iface._geometry_representation.selection_mappings[0]['surface_id'] = 21
    with pytest.raises(ValueError, match='wrong physical AW surface'):
        iface.export(dual=True)


def test_group_reversal_preserves_unoriented_topology_and_mapping(tmp_path):
    records = []
    for reverse in (False, True):
        iface = cached_interface(tmp_path / str(reverse), reverse=reverse)
        export_interface_dual_visualization(iface)
        records.append(rows(tmp_path / str(reverse) / 'dual' / 'dual_interface_mapping.csv'))
    for column in ('dual_feature_id', 'generator_1', 'generator_2', 'physical_surface_id', 'stable_atom_1', 'stable_atom_2'):
        assert [row[column] for row in records[0]] == [row[column] for row in records[1]]
    assert records[0][0]['group_1'] != records[1][0]['group_1']


def test_generated_scene_loads_off_and_isolates_generic_contact(tmp_path, monkeypatch):
    iface = cached_interface(tmp_path)
    export_interface_dual_visualization(iface)
    loaded, calls = {}, []
    class Commands:
        def load_cgo(self, geometry, name):
            loaded[name] = geometry
        def __getattr__(self, name):
            return lambda *args, **kwargs: calls.append((name, args))
    pymol = ModuleType('pymol'); pymol.cmd = Commands()
    cgo = ModuleType('pymol.cgo')
    for index, name in enumerate(('BEGIN', 'END', 'TRIANGLES', 'COLOR', 'VERTEX', 'CYLINDER')):
        setattr(cgo, name, index)
    monkeypatch.setitem(sys.modules, 'pymol', pymol)
    monkeypatch.setitem(sys.modules, 'pymol.cgo', cgo)
    namespace = {}
    exec((tmp_path / 'dual' / 'apollonius_interface.py').read_text(), namespace)
    assert {'apollonius_full', 'apollonius_selected', 'aw_interface'} <= set(loaded)
    assert loaded['apollonius_full'] != loaded['apollonius_selected']
    namespace['show_contact']('1:10,13')
    assert {'contact_dual_edge', 'contact_aw_surface'} <= set(loaded)
    assert ('extend', ('show_contact', namespace['show_contact'])) in calls


def test_info_summary_and_existing_exports_do_not_create_dual_tree(tmp_path, monkeypatch):
    iface = cached_interface(tmp_path)
    stream = StringIO()
    write_interface_dual_summary(stream, iface)
    assert 'Certified rejected: 1' in stream.getvalue()
    iface.export()
    assert not (tmp_path / 'dual').exists()
    iface.net.settings['net_type'] = 'pow'
    with pytest.raises(ValueError, match='AW interface'):
        iface.export(dual=True)


def test_cli_dual_and_apollonius_export_from_interface_only(tmp_path):
    output = importlib.import_module('vorpy.src.output.output')
    command = importlib.import_module('vorpy.src.command.command_export')
    iface = cached_interface(tmp_path)
    system = SimpleNamespace(groups=[], ifaces=[iface], update_progress=lambda **kwargs: None)
    command.argv_export(system, [['only', 'dual']])
    assert (tmp_path / 'dual' / 'apollonius_interface.pml').exists()
    output.other_exports(system, 'apollonius')
    assert system.export_timing['total'] > 0
    system.ifaces = []
    with pytest.raises(ValueError, match='two-group -i'):
        output.other_exports(system, 'dual')


def test_visualize_export_is_cache_only_and_reports_missing_layers(tmp_path, monkeypatch):
    output = importlib.import_module('vorpy.src.output.output')
    iface = cached_interface(tmp_path)
    iface.dir = None
    iface.group1.ball_ndxs = []
    iface.group2.ball_ndxs = []
    system = SimpleNamespace(
        groups=[], ifaces=[iface], files={'dir': str(tmp_path)},
        update_progress=lambda **kwargs: None,
        verbose=False,
    )
    iface.sys = system
    forbid_recalculation(monkeypatch)

    output.other_exports(system, 'visualize')

    launcher = tmp_path / 'interface_A_B' / 'visualization' / 'interface_A_B_visualization.pml'
    manifest = tmp_path / 'interface_A_B' / 'visualization' / 'visualization_manifest.json'
    assert launcher.exists()
    data = json.loads(manifest.read_text(encoding='utf-8'))
    assert data['layers']['partition_interface']['state'] == 'AVAILABLE'
    assert data['layers']['dual_contact_complex']['state'] == 'AVAILABLE'
    assert data['layers']['molecular_contact_surface:A']['state'] == 'NOT_CALCULATED'


def test_visualize_reuses_cached_molecular_contact_exporter(tmp_path, monkeypatch):
    output = importlib.import_module('vorpy.src.output.output')
    from vorpy.tests.geometry.test_molecular_contact_visualization import _surface

    iface = cached_interface(tmp_path)
    iface.dir = None
    iface.group1.ball_ndxs = []
    iface.group2.ball_ndxs = []
    iface.representation_caches = {'molecular_contact_surface:A': _surface('A')}
    system = SimpleNamespace(
        groups=[], ifaces=[iface], files={'dir': str(tmp_path)},
        update_progress=lambda **kwargs: None,
        verbose=False,
    )
    iface.sys = system
    forbid_recalculation(monkeypatch)

    output.other_exports(system, 'visualize')

    root = tmp_path / 'interface_A_B' / 'molecular_contact_surface_A'
    assert (root / 'molecular_contact_surface_A.obj').exists()
    manifest = json.loads(
        (tmp_path / 'interface_A_B' / 'visualization' / 'visualization_manifest.json').read_text()
    )
    assert manifest['layers']['molecular_contact_surface:A']['state'] == 'PARTIAL'


def test_compact_interface_bundle_has_one_log_and_physical_surface(tmp_path):
    from vorpy.src.output.visualization import export_compact_visualization_bundle

    iface = cached_interface(tmp_path)
    iface.dir = None
    iface.group1.ball_ndxs = []
    iface.group2.ball_ndxs = []
    system = SimpleNamespace(
        groups=[], ifaces=[iface], files={'dir': str(tmp_path)},
        update_progress=lambda **kwargs: None, verbose=False, name='fixture',
        balls=None,
    )
    iface.sys = system

    export_compact_visualization_bundle(system)

    root = tmp_path / 'interface_A_B'
    assert (root / 'aw_surfs.off').stat().st_size > 0
    assert (root / 'aw_edges.off').stat().st_size > 0
    assert (root / 'aw_verts.pdb').stat().st_size > 0
    assert (root / 'logs.csv').read_text(encoding='utf-8').startswith(
        'interface_id,representation,side,quantity,value,units,status,provenance\n'
    )
    log_rows = rows(root / 'logs.csv')
    assert log_rows
    assert all(row['status'] != 'AVAILABLE' for row in log_rows)
    assert {row['representation'] for row in log_rows} <= {
        'physical_partition_interface', 'molecular_contact_surface',
        'alpha_selection', 'water_analysis', 'timing',
    }
    assert not (root / 'interface.pml').exists()
    assert not (root / 'dual_mapping_edges.off').exists()


def test_compact_interface_make_net_does_not_create_legacy_directory(tmp_path, monkeypatch):
    module = importlib.import_module('vorpy.src.interface.interface')

    class FakeNetwork:
        def __init__(self, **kwargs):
            self.balls = kwargs['locs']
            self.settings = kwargs['settings']

    system = SimpleNamespace(
        files={'dir': str(tmp_path)}, _compact_interface_workflow=True,
        boundary_mode='shell',
    )
    iface = Interface.__new__(Interface)
    iface.sys = system
    iface.dir = None
    iface.name = 'A_B_interface'
    iface.settings = {'net_type': 'aw'}
    iface.group1 = SimpleNamespace(ball_ndxs=[1], name='A', group_id='A')
    iface.group2 = SimpleNamespace(ball_ndxs=[2], name='B', group_id='B')
    monkeypatch.setattr(module, 'Network', FakeNetwork)
    monkeypatch.setattr(
        module, 'network_geometry',
        lambda _sys: ([], [], [], set(), None),
    )
    monkeypatch.setattr(Interface, '_update_group_metadata', lambda *_args, **_kwargs: None)

    iface.make_net()

    assert iface.dir is None
    assert not (tmp_path / iface.name).exists()
