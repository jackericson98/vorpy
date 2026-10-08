"""Cached scientific object contracts, with calculation guards after reload."""
from dataclasses import asdict, replace
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np
import pandas as pd
import pytest

from vorpy.src.io import Session, load_session, save_session, save_network, load_network, ArchiveError
from vorpy.src.io.capabilities import CapabilityState
from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis, reverse_perspective


ROOT = Path(__file__).resolve().parents[3]


def guard_calculations(monkeypatch):
    import importlib
    from vorpy.src.network import Network
    def forbidden(*args, **kwargs):
        raise AssertionError('Scientific calculation attempted after reload')
    for name in ('__init__', 'build', 'build_edges', 'build_surfaces', 'find_verts', 'solve',
                 'build_vertex_gaussian_curvature', 'build_edge_mean_curvature',
                 'build_edge_gaussian_curvature', 'build_surface_mean_curvature'):
        if hasattr(Network, name):
            monkeypatch.setattr(Network, name, forbidden)
    targets = {
        'vorpy.src.geometry.duals': ('build_dual',),
        'vorpy.src.geometry.duals.builders': ('build_dual',),
        'vorpy.src.geometry.filtrations': ('build_alpha_filtration',),
        'vorpy.src.geometry.filtrations.alpha': ('build_alpha_filtration',),
        'vorpy.src.geometry.interfaces': ('build_alpha_interface', 'build_full_interface'),
        'vorpy.src.geometry.interfaces.alpha': ('build_alpha_interface', '_curvature_for_selected_interface'),
        'vorpy.src.geometry.interfaces.builders': ('build_full_interface',),
        'vorpy.src.interface.geometry_analysis': ('analyze_interface_geometry', 'summarize_physical_interface'),
        'vorpy.src.analyze.aw_interface_curvature': ('analyze_aw_interface_curvature',),
        'vorpy.src.analyze.open_interface_curvature': ('analyze_surface_topology', 'integrated_mean_curvature',
            'polygon_gaussian_terms', '_aw_supported_gauss_bonnet'),
        'vorpy.src.network.triangulate': ('triangulate_2D_Surface',),
        'vorpy.src.network.build_surf': ('triangulate_2D_Surface',),
        'vorpy.src.calculations.curvature': ('calc_surf_tri_curvs', 'calc_surf_tri_curvs_both', 'calc_avg_surface_curvature'),
        'vorpy.src.calculations.edge_mean_curvature': ('calculate_aw_network_edge_mean_curvatures',),
        'vorpy.src.calculations.edge_gaussian_curvature': ('calculate_aw_network_edge_gaussian_curvatures',),
        'vorpy.src.calculations.surface_mean_curvature': ('calculate_aw_network_surface_mean_curvatures',),
        'vorpy.src.geometry.molecular_contact': (
            'build_molecular_contact_surface', 'build_molecular_contact_surface_pair',
            '_build_carrier', '_pair_cap_area', '_face_area', '_arc_integral',
            '_circle_intersections',
        ),
    }
    for module, names in targets.items():
        target = importlib.import_module(module)
        for name in names:
            monkeypatch.setattr(target, name, forbidden)
    from vorpy.src.geometry.aw_alpha.complex import AWAlphaFiltration
    from vorpy.src.analyze.apollonius import ApolloniusComplex
    from vorpy.src.geometry.duals.model import DualComplex
    from vorpy.src.geometry.filtrations.alpha import AlphaFiltration
    from vorpy.src.geometry.interfaces.model import InterfaceRepresentation
    from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
    from vorpy.src.geometry.molecular_contact import (
        MolecularContactSurface, SphericalPatch, CircularArc, JunctionVertex,
        RefinedSeamPiece, TopologicalCut,
    )
    for cls in (AWAlphaFiltration, ApolloniusComplex, DualComplex, AlphaFiltration,
                InterfaceRepresentation, InterfaceGeometryAnalysis,
                MolecularContactSurface, SphericalPatch, CircularArc, JunctionVertex,
                RefinedSeamPiece, TopologicalCut):
        monkeypatch.setattr(cls, '__init__', forbidden)
    monkeypatch.setattr(AWAlphaFiltration, 'calculate_births', forbidden)
    import scipy.spatial
    monkeypatch.setattr(scipy.spatial, 'Delaunay', forbidden)


def attach_interface(net, a, b, name='cached interface'):
    from vorpy.src.group import Group
    from vorpy.src.interface import Interface
    groups = []
    for label, indices in (('A', a), ('B', b)):
        group = Group.__new__(Group)
        group.sys, group.net, group.name, group.group_id = net.sys, None, label, label
        group.ball_ndxs, group.settings, group.interfaces = sorted(indices), dict(net.settings), {}
        group.parent_interface = None
        groups.append(group)
    iface = Interface.__new__(Interface)
    iface.sys, iface.net, iface.name = net.sys, net, name
    iface.interface_id = name
    iface.group1, iface.group2 = groups
    iface.group1_name, iface.group2_name = 'A', 'B'
    iface.group1_indices, iface.group2_indices = set(a), set(b)
    iface.settings = dict(net.settings)
    iface.water_topology = iface.water_geometries = None
    iface.dir = None
    net.sys.groups.extend(groups)
    net.sys.ifaces.append(iface)
    for group in groups:
        group.interfaces[name] = iface
    return iface


@pytest.fixture
def cached_session():
    from vorpy.src.system import System
    from vorpy.src.network import Network
    from vorpy.src.geometry.duals import build_dual
    from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex
    from vorpy.src.geometry.interfaces import build_alpha_interface
    from vorpy.src.boundary import BoundaryConfig, BoundaryGenerator, BoundaryMode
    system = System.__new__(System)
    system.name, system.files = 'fixture', {'base_file': None, 'dir': None}
    system.groups, system.ifaces, system.residues, system.chains, system.sol = [], [], [], [], None
    system.balls = pd.DataFrame({'loc': [(0.,0.,0.), (2.,0.,0.), (0.,2.,0.), (0.,0.,2.)],
        'rad': [1.] * 4, 'num': range(4), 'res': [None] * 4, 'chn': [None] * 4,
        'mass': [12.] * 4, 'name': ['C'] * 4, 'element': ['C'] * 4,
        'res_name': ['ALA'] * 4, 'res_seq': ['1','2','3','4'], 'chain_name': ['A','B','B','C']})
    system.boundary_config = BoundaryConfig(BoundaryMode.SHELL)
    system.boundary_generators = (BoundaryGenerator((9.,9.,9.)),)
    system.boundary_mode = 'shell'
    net = Network.__new__(Network)
    net.sys, net.group, net.group_name, net.iface_grps = system, [0], 'fixture', ({0}, {1,2})
    net.settings, net.metrics, net.box = {'net_type': 'aw'}, {'tot': 1.}, {}
    net.balls = system.balls.copy()
    net.balls['complete'] = True
    net.verts = pd.DataFrame({'loc': [(1.,1.,1.), (1.,1.,-1.)], 'rad': [1.,1.],
        'balls': [[0,1,2,3], [0,1,2,3]], 'edges': [[0],[0]], 'surfs': [[0,1],[0,1]], 'complete': [True,True]})
    net.edges = pd.DataFrame({'balls': [[0,1,2]], 'verts': [[0,1]], 'surfs': [[0,1]],
        'points': [np.array([[1.,1.,1.],[1.,1.,-1.]])], 'complete': [True]})
    points = np.array([[1.,1.,1.],[1.,1.,-1.],[0.,0.,0.]])
    net.surfs = pd.DataFrame({'balls': [[0,1],[0,2]], 'verts': [[0,1],[0,1]], 'edges': [[0],[0]],
        'points': [points, points], 'tris': [np.array([[0,1,2]]),np.array([[0,1,2]])],
        'sa': [3.,4.], 'complete': [True,True], 'int_gauss_curv': [None,None]})
    net.boundary_indices, net.boundary_config = (3,), system.boundary_config
    iface = attach_interface(net, {0}, {1,2})
    dual = build_dual(net)
    records = {d: {} for d in range(4)}
    for pair, birth in (((0,1), -1.), ((0,2), 2.)):
        records[1][pair] = AlphaSimplex(1, pair, birth, birth, 'fixture', True)
    blocked = AlphaSimplex(2, (0,1,2), None, None, 'unresolved', False, notes='birth not calculated')
    filtration = AlphaFiltration(net, 'aw', records, {'native_alpha_units': 'A', 'filtration_kind': 'pair-only'}, blocked=[blocked])
    from vorpy.src.geometry.aw_alpha.complex import AWAlphaFiltration, AWAlphaRecord
    source = AWAlphaFiltration(net)
    source.records[1] = {pair: AWAlphaRecord(1, pair, row.geometric_birth, row.filtration_birth,
        row.birth_source, supported=row.supported) for pair, row in records[1].items()}
    source._max_dimension_built = 1
    source.issues = ['higher-dimensional births unresolved']
    filtration._aw_source = source
    rep = build_alpha_interface(net, dual, filtration, 0., {0}, {1,2}, calculate_curvature=False)
    from vorpy.src.geometry.interfaces.curvature import CurvatureSummary
    from dataclasses import fields
    values = {item.name: None for item in fields(CurvatureSummary)}
    values.update(scheme='aw', selection_mode='alpha', status='partial', source='cached fixture',
                  area_A2=rep.area_A2, smooth_H_A=-54.150161)
    rep.curvature = CurvatureSummary(**values)
    iface._geometry_dual, iface._geometry_filtration, iface._geometry_representation = dual, filtration, rep
    side = {'partial_H': -54.150161, 'total_H': None, 'total_K': None,
        'mean_certified': False, 'gaussian_total_certified': False, 'area': rep.area_A2,
        'topology_vertices': 2, 'topology_edges': 1, 'topology_faces': 1, 'euler_characteristic': 2,
        'boundary_loops': None, 'unresolved_reasons': ['Gaussian curvature unresolved']}
    iface.geometry_analysis = InterfaceGeometryAnalysis(metadata={'alpha_value': 0., 'scheme': 'aw',
        'solvent_context': False, 'additional_network_solves': 0},
        alpha_selection={'selected_generator_pairs': [[0,1]], 'selected_physical_pairs': 1,
            'eligible_bicolor_pairs': 2, 'status': 'partial'}, voronoi_side_1=side,
        voronoi_side_2=reverse_perspective(side), coverage={'mean_coverage': .5}, unresolved=['missing births'])
    system._interface_phase1_dry_cache = {'fixture': iface.geometry_analysis}
    system.interface_geometry_cache = {'surfs': {'fixture': {'points': points}}, 'verts': {}, 'edges': {}}
    return Session.from_network(net)


def scientific_snapshot(session):
    net = session.active_network
    iface = next(iface for iface in net.sys.ifaces if iface.net is net)
    dual, filtration, rep = iface._geometry_dual, iface._geometry_filtration, iface._geometry_representation
    return json.loads(json.dumps({
        'counts': {key: len(getattr(net, key)) for key in ('balls','verts','edges','surfs')},
        'network_id': net.archive_id, 'system_id': net.sys.archive_id, 'scope': net.network_scope,
        'generator_provenance': net.generator_provenance,
        'analysis': iface.geometry_analysis.to_dict(), 'analysis_id': iface.geometry_analysis.archive_id,
        'cache_provenance': {name: getattr(getattr(iface, name), 'archive_provenance', {}) for name in (
            'geometry_analysis', '_geometry_dual', '_geometry_filtration', '_geometry_representation')},
        'dual_id': dual.archive_id, 'filtration_id': filtration.archive_id,
        'representation_id': rep.archive_id, 'representation': rep.summary(),
        'surfaces': [asdict(row) for row in rep.surfaces], 'edges': [asdict(row) for row in rep.edges],
        'vertices': [asdict(row) for row in rep.vertices], 'mappings': rep.selection_mappings,
        'incidence_audit': rep.alpha_incidence_audit,
        'dual': [[asdict(row) for row in records.values()] for records in dual.simplices.values()],
        'dual_mapping': [(list(key), value) for key, value in dual.feature_to_simplices.items()],
        'dual_generators': dual.generators, 'dual_metadata': dual.metadata,
        'filtration': [[asdict(row) for row in records.values()] for records in filtration.records.values()],
        'blocked': [asdict(row) for row in filtration.blocked], 'filtration_metadata': filtration.metadata,
    }, sort_keys=True, allow_nan=False))


def export_cached(session, directory):
    from vorpy.src.interface.selected_export import export_selected_interface
    from vorpy.src.geometry.visualization.interface import export_interface_dual_visualization
    iface = next(iface for iface in session.active_network.sys.ifaces if iface.net is session.active_network)
    iface.dir = str(directory)
    export_selected_interface(iface)
    export_interface_dual_visualization(iface)


def compare_exports(before, after):
    paths = {path.relative_to(before) for path in before.rglob('*') if path.is_file()}
    assert paths == {path.relative_to(after) for path in after.rglob('*') if path.is_file()}
    for path in paths:
        left, right = (before / path).read_text(encoding='utf-8'), (after / path).read_text(encoding='utf-8')
        # Absolute launcher paths are presentation only; all numerical tables,
        # mappings, meshes, selection/status metadata must match exactly.
        for root, text in ((before, left), (after, right)):
            text = text.replace(str(root).replace('\\', '\\\\'), '<EXPORT>')
            text = text.replace(str(root), '<EXPORT>').replace(root.as_posix(), '<EXPORT>')
            if root == before:
                left = text
            else:
                right = text
        assert left == right, str(path)
    return len(paths)


def test_scientific_session_and_network_round_trip(cached_session, tmp_path, monkeypatch):
    archive = save_session(cached_session, tmp_path / 'session.vpy')
    assert cached_session.active_network.sys.archive_timing['total'] > 0
    expected = scientific_snapshot(cached_session)
    save_network(cached_session.active_network, tmp_path / 'network.vpy')
    guard_calculations(monkeypatch)
    restored = load_session(archive)
    assert scientific_snapshot(restored) == expected
    assert scientific_snapshot(load_session(tmp_path / 'network.vpy')) == expected
    iface = restored.active_network.sys.ifaces[0]
    assert iface._geometry_dual.network is restored.active_network
    assert iface._geometry_filtration.network is restored.active_network
    assert iface.sys._interface_phase1_dry_cache['fixture'] is iface.geometry_analysis
    assert iface.net.boundary_config is iface.sys.boundary_config
    assert iface.net.boundary_config.mode.value == 'shell'
    assert iface.sys.boundary_generators[0].position == (9.,9.,9.)
    assert iface.geometry_analysis.voronoi_side_1['partial_H'] == -54.150161
    assert iface.geometry_analysis.voronoi_side_1['total_H'] is None
    assert iface._geometry_filtration.simplex_birth((0,1,2)) is None
    assert iface._geometry_filtration._aw_source.network is iface.net
    assert iface._geometry_filtration._aw_source.incidence.network is iface.net
    assert iface._geometry_filtration._aw_source._max_dimension_built == 1
    assert iface._geometry_filtration._aw_source.records[2] == {}
    cache_points = restored.active_network.sys.interface_geometry_cache['surfs']['fixture']['points']
    np.testing.assert_array_equal(cache_points, np.array([[1., 1., 1.], [1., 1., -1.], [0., 0., 0.]]))
    assert restored.capabilities(result='mean_curvature')[0].state is CapabilityState.PARTIAL
    assert restored.capabilities(result='gaussian_curvature')[0].state is CapabilityState.UNRESOLVED
    assert restored.capabilities(environment='wet', partition='Power')[0].state is CapabilityState.NOT_CALCULATED
    assert restored.capabilities(representation='Dual A')[0].state is CapabilityState.NOT_SUPPORTED
    before, after = tmp_path / 'before', tmp_path / 'after'
    export_cached(cached_session, before)
    export_cached(restored, after)
    assert compare_exports(before, after) > 10


def test_cached_interface_fresh_process_export_fidelity(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'cached-interface.vpy')
    expected = scientific_snapshot(cached_session)
    export_cached(cached_session, tmp_path / 'before')
    code = '''from pathlib import Path
import json, sys, pytest
from vorpy.tests.io.test_scientific_session import guard_calculations, scientific_snapshot, export_cached
from vorpy.src.io import load_session
with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    session = load_session(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(scientific_snapshot(session), sort_keys=True), encoding='utf-8')
    export_cached(session, Path(sys.argv[3]))
'''
    completed = subprocess.run([sys.executable, '-c', code, str(archive), str(tmp_path / 'after.json'),
                                str(tmp_path / 'after')], cwd=ROOT,
                               env=dict(os.environ, PYTHONPATH=str(ROOT)),
                               capture_output=True, text=True, timeout=120)
    assert completed.returncode == 0, completed.stderr
    assert json.loads((tmp_path / 'after.json').read_text(encoding='utf-8')) == expected
    assert compare_exports(tmp_path / 'before', tmp_path / 'after') > 10


def test_session_vocabulary_normalizes_old_labels(cached_session, tmp_path):
    from vorpy.src.io.scientific_adapters import canonical_environment, canonical_network_scope
    net = cached_session.active_network
    iface = net.sys.ifaces[0]
    net.network_scope = 'interface/dedicated'
    iface.geometry_analysis.metadata['environment'] = 'wet'
    archive = save_session(cached_session, tmp_path / 'vocabulary.vpy')

    # Make the saved payload look like a historical archive. The loader must
    # accept the old labels and expose only the canonical session vocabulary.
    with zipfile.ZipFile(archive) as source:
        members = {name: source.read(name) for name in source.namelist()}
    state = json.loads(members['state.json'])
    for record in state['objects']:
        if record['kind'] == 'network':
            record['fields']['network_scope'] = 'interface/dedicated'
        if record['kind'] == 'interface_analysis':
            for pair in record['fields']['metadata']['map']:
                if pair[0] == 'environment':
                    pair[1] = 'wet'
            provenance = record['fields']['archive_provenance']
            for pair in provenance['map']:
                if pair[0] == 'environment':
                    pair[1] = 'wet'
    members['state.json'] = json.dumps(state).encode()
    with zipfile.ZipFile(archive, 'w') as destination:
        for name, content in members.items():
            destination.writestr(name, content)

    restored = load_session(archive)
    analysis = restored.active_network.sys.ifaces[0].geometry_analysis
    assert canonical_network_scope('system network') == 'system'
    assert canonical_environment('wet') == 'solvent_competing'
    assert restored.active_network.network_scope == 'interface'
    assert analysis.metadata['environment'] == 'solvent_competing'
    assert analysis.archive_provenance['environment'] == 'solvent_competing'
    assert restored.capabilities(result='analysis')[0].environment == 'solvent_competing'
    assert restored.capabilities(result='analysis')[0].network_scope == 'interface'


def test_missing_molecular_contact_surfaces_are_not_calculated(cached_session):
    for side in ('A', 'B'):
        rows = cached_session.capabilities(
            representation=f'molecular_contact_surface:{side}', result='geometry')
        assert len(rows) == 1
        assert rows[0].state is CapabilityState.NOT_CALCULATED


def _adapter_surface(contact_selection_source=None):
    from vorpy.src.geometry.molecular_contact import MolecularContactSurface

    return MolecularContactSurface(
        side='A', convention='expanded_union_of_balls', definition_version='fixture-v1',
        patches=(), arcs=(), junctions=(), carrier_atom_ids=(), selected_contact_ids=(),
        total_area_A2=1.0, naive_pairwise_area_A2=1.0, components=1, boundary_loops=0,
        euler_characteristic=1, nonmanifold_vertices=0, unresolved_contact_ids=(),
        area_method='adapter fixture', area_error_A2=0.0, status='PARTIAL',
        contact_selection_id='fixture:selected_contact_cache',
        contact_selection_source=contact_selection_source,
    )


def test_cached_surface_adapter_uses_fresh_surface_selection_source(cached_session):
    from vorpy.src.output.visualization import _cached_surface_result

    iface = cached_session.active_network.sys.ifaces[0]
    surface = _adapter_surface('fresh_surface_field')

    result = _cached_surface_result(iface, surface)

    assert result.identity.contact_selection_source == 'fresh_surface_field'
    assert result.provenance.details['contact_selection_source'] == 'fresh_surface_field'


def test_cached_surface_adapter_falls_back_to_explicit_archive_source(cached_session):
    from vorpy.src.output.visualization import _cached_surface_result

    iface = cached_session.active_network.sys.ifaces[0]
    surface = _adapter_surface()
    object.__setattr__(surface, 'archive_provenance', {
        'contact_selection_source': 'explicit_archive_field',
    })

    result = _cached_surface_result(iface, surface)

    assert result.identity.contact_selection_source == 'explicit_archive_field'
    assert result.provenance.details['contact_selection_source'] == 'explicit_archive_field'


def test_cached_surface_adapter_preserves_selection_source_after_reload(cached_session, tmp_path):
    from vorpy.src.output.visualization import _cached_surface_result

    iface = cached_session.active_network.sys.ifaces[0]
    surface = _adapter_surface('archived_surface_field')
    iface.representation_caches = {'molecular_contact_surface:A': surface}
    iface.representation_links = {'molecular_contact_surface:A': {'state': 'PARTIAL'}}

    archive = save_session(cached_session, tmp_path / 'surface-source.vpy')
    restored_iface = load_session(archive).active_network.sys.ifaces[0]
    restored_surface = restored_iface.representation_caches['molecular_contact_surface:A']
    result = _cached_surface_result(restored_iface, restored_surface)

    assert restored_surface.contact_selection_source == 'archived_surface_field'
    assert result.identity.contact_selection_source == 'archived_surface_field'


def test_cached_surface_adapter_rejects_missing_selection_provenance(cached_session):
    from vorpy.src.output.visualization import _cached_surface_result

    iface = cached_session.active_network.sys.ifaces[0]
    with pytest.raises(ValueError, match='identity'):
        _cached_surface_result(iface, _adapter_surface())


def test_cached_surface_adapter_prefers_surface_over_conflicting_archive_metadata(cached_session):
    from vorpy.src.output.visualization import _cached_surface_result

    iface = cached_session.active_network.sys.ifaces[0]
    surface = _adapter_surface('authoritative_surface_field')
    object.__setattr__(surface, 'archive_provenance', {
        'contact_selection_source': 'stale_archive_field',
        'source': 'conflict fixture',
    })

    result = _cached_surface_result(iface, surface)

    assert result.identity.contact_selection_source == 'authoritative_surface_field'
    assert result.provenance.details['contact_selection_source'] == 'authoritative_surface_field'
    assert result.provenance.details['source'] == 'conflict fixture'


def test_archive_profile_reports_stages_and_archive_size(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'profiled.vpy')
    profile = cached_session.active_network.sys.archive_timing

    for name in (
        'archive_preparation', 'object_collection', 'object_table_traversal',
        'array_conversion', 'array_serialization', 'metadata_serialization',
        'compression', 'zip_archive_write', 'archive_member_write',
        'filesystem_replace', 'total',
    ):
        assert profile[name] >= 0.0
    assert profile['archive_bytes'] == archive.stat().st_size
    assert profile['object_count'] > 0
    assert profile['array_count'] >= 0
    assert profile['table_count'] >= 0


def test_archive_profile_reports_object_type_diagnostics(cached_session, tmp_path):
    save_session(cached_session, tmp_path / 'object-types.vpy')
    profile = cached_session.active_network.sys.archive_timing
    diagnostics = profile['object_type_timings']

    assert diagnostics
    assert sum(item['count'] for item in diagnostics.values()) == profile['object_count']
    assert all(item['seconds'] >= 0.0 for item in diagnostics.values())
    assert all(all(value >= 0.0 for value in item['fields'].values())
               for item in diagnostics.values())
    tables = profile['table_diagnostics']
    assert tables
    assert all(item['rows'] >= 0 and item['columns'] for item in tables.values())
    assert all(all(value['seconds'] >= 0.0 for value in item['columns'].values())
               for item in tables.values())


def test_archive_timing_manifest_persists_bounded_diagnostics(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'timing-manifest.vpy')
    manifest_path = archive.with_name(archive.name + '.timing.json')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    profile = manifest['timing']
    with zipfile.ZipFile(archive) as stored:
        assert not any(name.endswith('.timing.json') for name in stored.namelist())

    assert manifest['format'] == 'vorpy.archive_timing'
    assert manifest['manifest_version'] == 1
    assert profile['archive_bytes'] == archive.stat().st_size
    assert profile['member_count'] == sum(item['count'] for item in profile['size_histograms']['members'].values())
    assert profile['array_count'] == sum(item['count'] for item in profile['size_histograms']['arrays'].values())
    assert profile['array_compressed_bytes'] <= profile['compressed_member_bytes']
    assert profile['array_uncompressed_bytes'] <= profile['uncompressed_member_bytes']
    for name in ('json_staging', 'member_compression_write', 'zip_finalization', 'total'):
        assert profile[name] >= 0.0
    assert profile['timing_semantics']['json_staging'].startswith('exclusive')


def test_archive_progress_checkpoints_are_opt_in(cached_session, tmp_path):
    events = []
    save_session(cached_session, tmp_path / 'progress.vpy',
                 progress=lambda stage, completed, total: events.append((stage, completed, total)))
    assert events
    assert events[-1][0] == 'zip_finalization'


def test_molecular_contact_curvature_round_trip_and_legacy_absence(cached_session, tmp_path):
    from vorpy.src.geometry.molecular_contact import MolecularContactSurface, MolecularContactSurfaceCurvature

    curvature = MolecularContactSurfaceCurvature(
        1.25, -2.5, -1.25, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
        1.6, 1.7, -0.1, 'PARTIAL', 'CERTIFIED', 'PARTIAL',
        'curvature fixture', ('missing corner audit',), 2, 3, 4,
    )
    object.__setattr__(curvature, 'archive_provenance', {'source': 'curvature-fixture'})
    surface = MolecularContactSurface(
        'A', 'expanded_union_of_balls', 'fixture-v1', (), (), (), (), (),
        1.0, 1.0, 1, 0, 1, 0, (), 'fixture', 0.0, 'PARTIAL',
        curvature=curvature, contact_selection_source='selected_contact_cache',
    )
    iface = cached_session.active_network.sys.ifaces[0]
    object.__setattr__(surface, 'archive_provenance', {
        'source': 'surface-fixture', 'network_id': iface.net.archive_id,
    })
    iface.representation_caches = {'molecular_contact_surface:A': surface}
    iface.representation_links = {'molecular_contact_surface:A': {'state': 'PARTIAL'}}

    archive = save_session(cached_session, tmp_path / 'curvature.vpy')
    loaded = load_session(archive).active_network.sys.ifaces[0]
    restored = loaded.representation_caches['molecular_contact_surface:A']
    assert asdict(restored.curvature) == asdict(curvature)
    assert restored.curvature.archive_provenance == curvature.archive_provenance

    with zipfile.ZipFile(archive) as source:
        members = {name: source.read(name) for name in source.namelist()}
    state = json.loads(members['state.json'])
    record = next(item for item in state['objects'] if item['kind'] == 'molecular_contact_surface')
    record['fields'].pop('curvature')
    record['fields'].pop('contact_selection_source')
    members['state.json'] = json.dumps(state, separators=(',', ':')).encode()
    legacy = tmp_path / 'legacy-no-curvature.vpy'
    with zipfile.ZipFile(legacy, 'w') as destination:
        for name, content in members.items():
            destination.writestr(name, content)
    legacy_surface = load_session(legacy).active_network.sys.ifaces[0].representation_caches[
        'molecular_contact_surface:A']
    assert legacy_surface.curvature is None
    assert legacy_surface.contact_selection_source is None


def _typed_result_snapshot(session):
    from vorpy.src.output.visualization import _cached_interface_result, _cached_surface_result

    net = session.active_network
    iface = next(iface for iface in net.sys.ifaces if iface.net is net)
    physical = _cached_interface_result(iface)
    molecular = _cached_surface_result(
        iface, iface.representation_caches['molecular_contact_surface:A']
    )
    return json.loads(json.dumps({
        'physical': physical.to_dict(),
        'molecular_A': molecular.to_dict(),
    }, sort_keys=True, allow_nan=False))


def test_typed_results_survive_fresh_process_vpy_round_trip(cached_session, tmp_path):
    from vorpy.src.geometry.molecular_contact import (
        MolecularContactSurface, MolecularContactSurfaceCurvature,
    )

    curvature = MolecularContactSurfaceCurvature(
        1.25, -2.5, -1.25, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
        1.6, 1.7, -0.1, 'PARTIAL', 'CERTIFIED', 'PARTIAL',
        'curvature fixture', ('missing corner audit',), 2, 3, 4,
    )
    surface = MolecularContactSurface(
        'A', 'expanded_union_of_balls', 'fixture-v1', (), (), (), (), (),
        1.0, 1.0, 1, 0, 1, 0, (), 'fixture', 0.0, 'PARTIAL',
        contact_selection_id='fixture:selected_contact_cache',
        contact_selection_source='selected_contact_cache', curvature=curvature,
    )
    iface = cached_session.active_network.sys.ifaces[0]
    object.__setattr__(surface, 'archive_provenance', {
        'source': 'typed-results-round-trip',
        'network_id': iface.net.archive_id,
    })
    iface.representation_caches = {'molecular_contact_surface:A': surface}
    iface.representation_links = {'molecular_contact_surface:A': {'state': 'PARTIAL'}}

    # Seed the same complete cache provenance that archive preparation preserves.
    from hashlib import sha256
    from vorpy.src.io.scientific_adapters import canonical_environment, canonical_network_scope, scope, stable_id
    analysis = iface.geometry_analysis
    alpha = analysis.metadata['alpha_value']
    object.__setattr__(analysis, 'archive_provenance', {
        'network_id': stable_id(iface.net), 'system_id': stable_id(iface.sys),
        'network_scope': canonical_network_scope(getattr(iface.net, 'network_scope', None) or scope(iface.net)),
        'partition': iface.net.settings.get('net_type'),
        'interface_id': stable_id(iface),
        'alpha': alpha,
        'environment': canonical_environment(analysis.metadata.get('environment', 'dry')),
        'representation': 'Voronoi',
        'alpha_convention': analysis.metadata.get('alpha_convention'),
        'alpha_units': analysis.metadata.get('alpha_units'),
        'query_id': sha256(f'{stable_id(iface._geometry_filtration) if iface._geometry_filtration else iface.net.archive_id}:{alpha!r}'.encode()).hexdigest(),
    })

    archive = save_session(cached_session, tmp_path / 'typed-results.vpy')
    assert surface.archive_provenance['contact_selection_source'] == surface.contact_selection_source
    expected = _typed_result_snapshot(cached_session)
    physical = expected['physical']
    molecular = expected['molecular_A']
    assert physical['identity']['representation'] == 'voronoi'
    assert physical['identity']['side'] == 'shared'
    assert physical['metrics']['area']['units'] == 'Å²'
    assert physical['mean_curvature']['complete_total']['value'] is None
    assert physical['mean_curvature']['complete_total']['status'] == 'UNRESOLVED'
    assert molecular['identity']['representation'] == 'molecular_contact_surface'
    assert molecular['identity']['side'] == 'A'
    assert molecular['metrics']['area']['value'] == 1.0
    assert molecular['metrics']['area']['units'] == 'Å²'
    assert molecular['mean_curvature']['boundary_contribution']['value'] is None
    assert molecular['mean_curvature']['boundary_contribution']['status'] == 'NOT_SUPPORTED'
    assert molecular['provenance']['details']['source'] == 'typed-results-round-trip'

    assert archive.is_file() and archive.suffix == '.vpy'
    output = tmp_path / 'typed-results-after.json'
    code = '''
import json
import sys
from pathlib import Path
import pytest
from vorpy.tests.io.test_scientific_session import guard_calculations, _typed_result_snapshot
from vorpy.src.io import load_session

with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    restored = load_session(sys.argv[1])
    Path(sys.argv[2]).write_text(
        json.dumps(_typed_result_snapshot(restored), sort_keys=True),
        encoding='utf-8',
    )
'''
    completed = subprocess.run(
        [sys.executable, '-c', code, str(archive), str(output)],
        cwd=ROOT, env=dict(os.environ, PYTHONPATH=str(ROOT)),
        capture_output=True, text=True, timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(output.read_text(encoding='utf-8')) == expected


def _load_2kai_molecular_surface(side):
    from vorpy.src.geometry.molecular_contact import (
        MolecularContactSurface, SphericalPatch, CircularArc, JunctionVertex,
        RefinedSeamPiece, TopologicalCut,
    )
    source = ROOT / 'output/dual_side_incidence_audit/2kai_power/molecular_contact_surface' / f'surface_{side}.json'
    payload = json.loads(source.read_text(encoding='utf-8'))
    summary = payload['summary']
    with (ROOT / '2kai_cazals_atoms.csv').open(newline='', encoding='utf-8') as stream:
        atom_rows = {int(row['index']): row for row in csv.DictReader(stream)}
    pairs = set()
    with (ROOT / '2kai_alpha0_edges.csv').open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            pairs.add(tuple(sorted((int(row['facet_1_i']), int(row['facet_1_j'])))))
            pairs.add(tuple(sorted((int(row['facet_2_i']), int(row['facet_2_j'])))))

    def atom_stable_id(row):
        return (f"{row['index']}|pdb={row['pdb_serial']}|group={row['group']}|"
                f"{row['chain']}:{row['residue_name']}{row['residue_number']}:{row['atom_name']}")

    selected_ids = tuple(f'power-alpha0:{a}:{b}' for a, b in sorted(pairs))
    carrier_indices = {
        (a if atom_rows[a]['group'] == 'A' else b) if side == 'A'
        else (b if atom_rows[a]['group'] == 'A' else a)
        for a, b in pairs
    }
    carrier_atom_ids = tuple(atom_stable_id(atom_rows[index]) for index in sorted(carrier_indices))

    def record(cls, value, tuple_fields):
        value = dict(value)
        for name in tuple_fields:
            value[name] = tuple(value[name])
        return cls(**value)

    patches = tuple(record(SphericalPatch, value, (
        'partner_atom_ids', 'self_occluding_atom_ids', 'source_contact_ids',
        'boundary_arc_ids', 'junction_ids')) for value in payload['patches'])
    arcs = tuple(record(CircularArc, value, (
        'start_xyz', 'end_xyz', 'circle_normal', 'source_atom_ids', 'source_contact_ids'))
                 for value in payload['arcs'])
    junctions = tuple(record(JunctionVertex, value, (
        'xyz', 'unit_direction', 'source_circle_ids', 'incident_arc_ids',
        'source_carrier_atom_ids', 'incident_edge_ids')) for value in payload['junctions'])
    refined_seams = tuple(record(RefinedSeamPiece, value, (
        'carrier_atom_ids', 'circle_center', 'circle_normal', 'parameter_interval',
        'endpoint_junction_ids', 'start_xyz', 'end_xyz', 'incident_patch_ids',
        'source_arc_ids', 'source_contact_ids')) for value in payload.get('refined_seams', ()))
    for seam in refined_seams:
        object.__setattr__(seam, 'orientation_on_patches',
                           tuple(tuple(item) for item in seam.orientation_on_patches))
    topological_cuts = tuple(record(TopologicalCut, value, (
        'cycle_indices', 'endpoint_junction_ids', 'start_xyz', 'end_xyz'))
                             for value in payload.get('topological_cuts', ()))
    return MolecularContactSurface(
        side=summary['side'], convention=summary['convention'],
        definition_version=summary['definition_version'], patches=patches,
        arcs=arcs, junctions=junctions,
        carrier_atom_ids=carrier_atom_ids,
        selected_contact_ids=selected_ids,
        total_area_A2=summary['area_A2'], naive_pairwise_area_A2=summary['naive_pairwise_area_A2'],
        components=summary['components'], boundary_loops=summary['boundary_loops'],
        euler_characteristic=summary['euler_characteristic'],
        nonmanifold_vertices=summary['nonmanifold_vertices'],
        unresolved_contact_ids=tuple(summary['unresolved_contacts']),
        area_method=summary['area_method'], area_error_A2=summary['area_error_A2'],
        status=summary['status'], diagnostics=tuple(summary['diagnostics']),
        contact_selection_source=summary.get('contact_selection_source'),
        contact_selection_id=summary['contact_selection_id'],
        refined_seams=refined_seams,
        global_vertex_count=summary.get('global_vertices', len(junctions)),
        global_edge_count=summary.get('global_edges', 0),
        global_face_count=summary.get('global_faces', len(patches)),
        manifold_vertices=summary.get('manifold_vertices', 0),
        singular_vertices=summary.get('singular_vertices', 0),
        junction_link_status=tuple(sorted(summary.get('junction_link_status', {}).items())),
        orientability=summary.get('orientability'),
        genus_by_component=(tuple(summary['genus_by_component'])
                            if summary.get('genus_by_component') is not None else None),
        area_status=summary.get('area_status', 'CERTIFIED'),
        topological_cuts=topological_cuts,
        topological_face_count=summary.get('topological_faces', len(patches)),
        spherical_patch_record_count=summary.get('spherical_patch_records', len(patches)),
        face_boundary_cycle_counts=tuple(sorted(summary.get('face_boundary_cycle_counts', {}).items())),
        geometric_edge_count=summary.get('geometric_edges', summary.get('global_edges', 0)),
        betti_numbers=(tuple(summary['betti_numbers'])
                       if summary.get('betti_numbers') is not None else None),
        incidence_status=summary.get('incidence_status', 'UNRESOLVED'),
    )


def _molecular_snapshot(surface):
    return json.loads(json.dumps({
        'scientific': asdict(surface),
        'archive_id': getattr(surface, 'archive_id', None),
        'archive_provenance': getattr(surface, 'archive_provenance', {}),
        'patch_archive_ids': [getattr(value, 'archive_id', None) for value in surface.patches],
        'arc_archive_ids': [getattr(value, 'archive_id', None) for value in surface.arcs],
        'junction_archive_ids': [getattr(value, 'archive_id', None) for value in surface.junctions],
    }, sort_keys=True))


def _molecular_result_identity(surface):
    from vorpy.src.results import ResultIdentity
    source_network_id = surface.patches[0].source_network_id
    source_interface_id = surface.patches[0].source_interface_id
    return ResultIdentity(
        result_kind='molecular_contact_surface', system='2KAI', frame=None,
        network_id=source_network_id, network_scope='interface', environment='dry',
        partition='power', representation='molecular_contact_surface',
        radius_configuration_id='cazals-expanded-v1', system_id='2KAI',
        interface_id=source_interface_id, group_A='A', group_B='B',
        orientation='A -> B', side=surface.side,
        partner_id='B' if surface.side == 'A' else 'A',
        contact_selection_source=surface.contact_selection_source,
        contact_selection_id=surface.contact_selection_id,
        molecular_surface_convention=surface.convention,
        definition_version=surface.definition_version,
        source_network_id=source_network_id,
        source_interface_id=source_interface_id,
        interaction_id='2KAI:AB', orientation_convention='outward_molecular_normal',
    )


def test_2kai_molecular_surface_fresh_process_round_trip(cached_session, tmp_path, monkeypatch):
    source = ROOT / 'output/dual_side_incidence_audit/2kai_power/molecular_contact_surface'
    if not all((source / f'surface_{side}.json').exists() for side in ('A', 'B')):
        pytest.skip('Real 2KAI molecular-contact cache is not present')
    surfaces = {side: _load_2kai_molecular_surface(side) for side in ('A', 'B')}
    assert surfaces['A'].summary()['carrier_atoms'] == 131
    assert surfaces['A'].summary()['patches'] == 111
    assert surfaces['A'].summary()['arcs'] == 1169
    assert surfaces['A'].summary()['junction_vertices'] == 670
    assert surfaces['A'].summary()['components'] == 1
    assert surfaces['A'].summary()['boundary_loops'] == 5
    assert surfaces['A'].summary()['global_vertices'] == 670
    assert surfaces['A'].summary()['global_edges'] == 786
    assert surfaces['A'].summary()['geometric_edges'] == 784
    assert surfaces['A'].summary()['topological_cuts'] == 2
    assert surfaces['A'].summary()['euler_characteristic'] == -5
    assert surfaces['A'].summary()['betti_numbers'] == (1, 6, 0)
    assert surfaces['A'].summary()['orientability'] is True
    assert surfaces['A'].summary()['genus_by_component'] == (1,)
    assert surfaces['A'].status == 'CERTIFIED'
    assert surfaces['A'].total_area_A2 == pytest.approx(638.9827835, abs=1e-6)
    assert surfaces['B'].summary()['carrier_atoms'] == 65
    assert surfaces['B'].summary()['patches'] == 70
    assert surfaces['B'].summary()['arcs'] == 1067
    assert surfaces['B'].summary()['junction_vertices'] == 684
    assert surfaces['B'].summary()['components'] == 1
    assert surfaces['B'].summary()['boundary_loops'] == 5
    assert surfaces['B'].summary()['global_vertices'] == 684
    assert surfaces['B'].summary()['global_edges'] == 757
    assert surfaces['B'].summary()['geometric_edges'] == 754
    assert surfaces['B'].summary()['topological_cuts'] == 3
    assert surfaces['B'].summary()['euler_characteristic'] == -3
    assert surfaces['B'].summary()['betti_numbers'] == (1, 4, 0)
    assert surfaces['B'].summary()['orientability'] is True
    assert surfaces['B'].summary()['genus_by_component'] == (0,)
    assert surfaces['B'].status == 'CERTIFIED'
    assert surfaces['B'].total_area_A2 == pytest.approx(774.8019194, abs=1e-6)
    # The checked-in JSON cache predates explicit selection-source provenance.
    assert surfaces['A'].contact_selection_source is None
    assert surfaces['B'].contact_selection_source is None

    iface = cached_session.active_network.sys.ifaces[0]
    iface.representation_caches = {
        f'molecular_contact_surface:{side}': surfaces[side] for side in ('A', 'B')
    }
    iface.representation_links = {
        f'molecular_contact_surface:{side}': {
            'state': ('AVAILABLE' if surfaces[side].status == 'CERTIFIED' else 'PARTIAL'),
            'reason': ('' if surfaces[side].status == 'CERTIFIED' else 'Topology remains uncertified'),
            'dependencies': {'selected_ab_contacts': 'AVAILABLE'},
        }
        for side in ('A', 'B')
    }
    archive = save_session(cached_session, tmp_path / '2kai_molecular.vpy')
    archive_bytes = archive.read_bytes()
    expected = {side: _molecular_snapshot(surfaces[side]) for side in ('A', 'B')}
    (tmp_path / 'expected.json').write_text(json.dumps(expected, sort_keys=True), encoding='utf-8')
    with zipfile.ZipFile(archive) as payload:
        assert not any(name.lower().endswith(('.obj', '.ply', '.stl')) for name in payload.namelist())

    code = '''import json, sys
from pathlib import Path
import pytest
from vorpy.tests.io.test_scientific_session import (
    guard_calculations, _molecular_result_identity, _molecular_snapshot,
)
from vorpy.src.io import load_session
from vorpy.src.io.capabilities import CapabilityState
from vorpy.src.geometry.visualization import export_molecular_contact_bundle
from vorpy.src.results import Provenance, MolecularContactSurfaceResult

with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    session = load_session(sys.argv[1])
    iface = session.active_network.sys.ifaces[0]
    surfaces = {side: iface.representation_caches[f"molecular_contact_surface:{side}"]
                for side in ("A", "B")}
    expected = json.loads(Path(sys.argv[2]).read_text())
    assert {side: _molecular_snapshot(surfaces[side]) for side in surfaces} == expected
    for side, surface in surfaces.items():
        row = session.capabilities(
            representation=f"molecular_contact_surface:{side}", result="geometry")[0]
        expected_state = CapabilityState.AVAILABLE if surface.status == "CERTIFIED" else CapabilityState.PARTIAL
        assert row.state is expected_state
        if expected_state is CapabilityState.PARTIAL:
            assert row.reason == "Topology remains uncertified"
        if surface.contact_selection_source is None:
            with pytest.raises(ValueError):
                _molecular_result_identity(surface)
            continue
        result = surface.to_contract(
            _molecular_result_identity(surface),
            Provenance("fresh-process cached molecular contact surface"),
        )
        assert isinstance(result, MolecularContactSurfaceResult)
        assert result.geometric_counts is not None
        assert result.geometric_counts.quantities["geometric_edge_count"].value == surface.geometric_edge_count
        assert result.geometric_counts.quantities["spherical_patch_count"].value == len(surface.patches)
        assert result.topology is not None
        topology = result.topology.quantities
        assert topology["vertex_count"].value == surface.global_vertex_count
        assert topology["edge_count"].value == surface.global_edge_count
        assert topology["topological_edge_count"].value == surface.global_edge_count
        assert topology["topological_cut_count"].value == len(surface.topological_cuts)
        assert topology["patch_attachment_record_count"].value == surface.spherical_patch_record_count
        assert topology["face_count"].value == surface.global_face_count
        assert topology["boundary_component_count"].value == surface.boundary_loops
        assert topology["euler_characteristic"].value == surface.euler_characteristic
        assert topology["orientable"].value is surface.orientability
        assert topology["genus"].value == surface.genus_by_component[0]
        assert topology["manifold_vertex_count"].value == surface.manifold_vertices
        assert topology["singular_vertex_count"].value == surface.singular_vertices
        assert topology["chain_complex_valid"].value is True
        assert topology["boundary_of_boundary_zero"].value is True
        assert topology["beta0"].value == surface.betti_numbers[0]
        assert topology["beta1"].value == surface.betti_numbers[1]
        assert topology["beta2"].value == surface.betti_numbers[2]
    export_root = Path(sys.argv[4])
    manifest = export_molecular_contact_bundle(
        {f"molecular_contact_surface:{side}": value for side, value in surfaces.items()},
        export_root, interaction="2KAI_AB",
    )
    for side, surface in surfaces.items():
        expected = "AVAILABLE" if surface.status == "CERTIFIED" else "PARTIAL"
        assert manifest["capabilities"][f"molecular_contact_surface:{side}"]["state"] == expected
    Path(sys.argv[3]).write_text(json.dumps({side: _molecular_snapshot(value)
        for side, value in surfaces.items()}, sort_keys=True), encoding="utf-8")
'''
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    completed = subprocess.run(
        [sys.executable, '-c', code, str(archive), str(tmp_path / 'expected.json'),
         str(tmp_path / 'after.json'), str(tmp_path / 'export')],
        cwd=ROOT, env=env, capture_output=True, text=True, timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    assert archive.read_bytes() == archive_bytes
    assert (tmp_path / 'export/interfaces/2KAI_AB/visualization/molecular_contact_visualization_manifest.json').exists()
    print(f'2KAI molecular archive: {archive.stat().st_size} bytes; canonical A/B round-trip and cache-only export passed')


def test_refined_seam_and_topological_cut_round_trip(cached_session, tmp_path):
    from vorpy.src.geometry.molecular_contact import RefinedSeamPiece, TopologicalCut

    surface = _load_2kai_molecular_surface('A')
    patch = surface.patches[0]
    arc = surface.arcs[0]
    seam = RefinedSeamPiece(
        'refined-test-seam', (patch.carrier_atom_id, 'same-group-carrier'), arc.circle_id,
        patch.carrier_center, patch.carrier_radius, arc.circle_normal, (0.0, 0.25),
        (arc.start_vertex_id or 'start', arc.end_vertex_id or 'end'),
        arc.start_xyz, arc.end_xyz, (patch.patch_id,), (arc.arc_id,),
        arc.source_contact_ids, ((patch.patch_id, 1),), 'internal', 'PARTIAL',
    )
    cut = TopologicalCut(
        'topological-test-cut', patch.patch_id, patch.carrier_atom_id, (0, 1),
        (arc.start_vertex_id or 'start', arc.end_vertex_id or 'end'),
        arc.start_xyz, arc.end_xyz,
    )
    cached = replace(surface, refined_seams=(seam,), topological_cuts=(cut,))
    iface = cached_session.active_network.sys.ifaces[0]
    iface.representation_caches = {'molecular_contact_surface:A': cached}
    iface.representation_links = {'molecular_contact_surface:A': {
        'state': 'PARTIAL', 'reason': 'Topology remains uncertified',
    }}
    archive = save_session(cached_session, tmp_path / 'refined.vpy')
    restored = load_session(archive).active_network.sys.ifaces[0]
    loaded = restored.representation_caches['molecular_contact_surface:A']
    assert asdict(loaded.refined_seams[0]) == asdict(seam)
    assert asdict(loaded.topological_cuts[0]) == asdict(cut)


def test_old_version_two_archive_loads_with_new_reader(cached_session, tmp_path, monkeypatch):
    archive = save_session(cached_session, tmp_path / 'v2.vpy')
    with zipfile.ZipFile(archive) as source:
        members = {name: source.read(name) for name in source.namelist()}
    metadata = json.loads(members['metadata.json'])
    metadata['format_version'] = 2
    members['metadata.json'] = json.dumps(metadata).encode()
    with zipfile.ZipFile(archive, 'w') as destination:
        for name, content in members.items():
            destination.writestr(name, content)
    guard_calculations(monkeypatch)
    restored = load_session(archive)
    assert restored.active_network is not None


def test_version_two_rejects_reserved_representation_fields(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'v2-with-future-field.vpy')
    with zipfile.ZipFile(archive) as source:
        members = {name: source.read(name) for name in source.namelist()}
    metadata = json.loads(members['metadata.json'])
    metadata['format_version'] = 2
    state = json.loads(members['state.json'])
    next(record for record in state['objects'] if record['kind'] == 'interface')['fields'][
        'representation_links'] = {'map': []}
    members['metadata.json'] = json.dumps(metadata).encode()
    members['state.json'] = json.dumps(state).encode()
    with zipfile.ZipFile(archive, 'w') as destination:
        for name, content in members.items():
            destination.writestr(name, content)
    with pytest.raises(ArchiveError, match='Unknown fields'):
        load_session(archive)


def test_unknown_representation_links_are_metadata_only(cached_session, tmp_path, monkeypatch):
    from vorpy.src.io.scientific_adapters import stable_id
    net = cached_session.active_network
    iface = net.sys.ifaces[0]
    iface.representation_links = {
        'future_representation:A': {
            'state': 'PARTIAL', 'network_id': stable_id(net), 'interface_id': stable_id(iface)},
        'future_representation:B': {
            'state': 'UNRESOLVED', 'network_id': net.archive_id, 'interface_id': iface.archive_id},
        'future_representation:C': {
            'dependencies': {
                'selected_ab_contacts': 'AVAILABLE',
                'spherical_clipping': 'NOT_CALCULATED',
            },
            'network_id': net.archive_id, 'interface_id': iface.archive_id},
    }
    archive = save_session(cached_session, tmp_path / 'links.vpy')
    guard_calculations(monkeypatch)
    restored = load_session(archive)
    for key, state in (
        ('future_representation:A', CapabilityState.PARTIAL),
        ('future_representation:B', CapabilityState.UNRESOLVED),
        ('future_representation:C', CapabilityState.NOT_CALCULATED),
    ):
        row = restored.capabilities(representation=key, result='geometry')[0]
        assert row.state is state
    assert (restored.capabilities(
        representation='future_representation:C', result='geometry')[0].reason
        == 'clipping geometry unavailable')


def test_unregistered_representation_adapter_is_rejected(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'unknown-adapter.vpy')
    with zipfile.ZipFile(archive) as source:
        members = {name: source.read(name) for name in source.namelist()}
    state = json.loads(members['state.json'])
    next(record for record in state['objects'] if record['kind'] == 'dual_complex')['kind'] = \
        'unregistered_representation'
    members['state.json'] = json.dumps(state).encode()
    with zipfile.ZipFile(archive, 'w') as destination:
        for name, content in members.items():
            destination.writestr(name, content)
    with pytest.raises(ArchiveError, match='Invalid object identity/type'):
        load_session(archive)


def test_cache_query_reuses_restored_analysis(cached_session, tmp_path, monkeypatch):
    from vorpy.src.interface.geometry_analysis import analyze_interface_geometry
    archive = save_session(cached_session, tmp_path / 'session.vpy')
    guard_calculations(monkeypatch)
    restored = load_session(archive)
    iface = restored.active_network.sys.ifaces[0]
    # A previously bound public function may inspect its key and return the
    # original cached object; all calculation helpers remain forbidden.
    assert analyze_interface_geometry(iface) is iface.geometry_analysis


def test_wrong_network_cache_rejected(cached_session, tmp_path):
    from copy import copy
    iface = cached_session.active_network.sys.ifaces[0]
    iface._geometry_dual.network = copy(iface.net)
    with pytest.raises(ArchiveError, match='another network'):
        save_session(cached_session, tmp_path / 'bad.vpy')


def test_empty_and_unsolved_workbench_session(tmp_path):
    from vorpy.workbench.domain import AnalysisResult, Atom
    from vorpy.workbench.project import Project
    from vorpy.workbench.session import session_from_workbench, workbench_from_session
    source = tmp_path / 'missing.pdb'
    result = AnalysisResult(source, 'unsolved', atoms=[Atom(0,1,'C','C',(0.,0.,0.))])
    session = session_from_workbench(Project(view_state={'camera_position': [[1,2,3],[0,0,0],[0,1,0]]}),
        {source: result}, source, {source: ({0}, {'A': (0,)}, {})})
    saved = save_session(session, tmp_path / 'unsolved.vpy')
    loaded = load_session(saved)
    project, results, active, states = workbench_from_session(loaded, saved)
    assert active == source and results[source].export_group is None
    assert states[source][0] == {0} and project.view_state['camera_position'][0] == [1,2,3]
    with pytest.raises(ArchiveError, match='no network'):
        load_network(saved)


def test_workbench_scientific_reopen(cached_session, tmp_path, monkeypatch):
    from vorpy.workbench.domain import AnalysisResult, GeometryLayer
    from vorpy.workbench.project import Project
    from vorpy.workbench.session import session_from_workbench, workbench_from_session
    net = cached_session.active_network
    group = net.sys.groups[0]
    group.net = net
    source = tmp_path / 'deleted.pdb'
    result = AnalysisResult(source, 'cached', export_group=group,
        layers=[GeometryLayer('selected', 'surfaces', points=np.array([[1.,2.,3.]]), color='#123456')])
    session = session_from_workbench(Project(view_state={'export_config': {'atom_format': 'xyz'}}),
        {source: result}, source, trajectory_results={(source, 1): result})
    archive = save_session(session, tmp_path / 'workbench.vpy')
    guard_calculations(monkeypatch)
    loaded = load_session(archive)
    project, results, active, _ = workbench_from_session(loaded, archive)
    assert active == source
    assert results[source].export_group.net is loaded.active_network
    assert results[source].layers[0].color == '#123456'
    assert loaded.workbench['trajectory_results'][(source, 1)] is results[source]
    assert project.view_state['export_config']['atom_format'] == 'xyz'


def test_old_network_has_honest_missing_capabilities(tmp_path, monkeypatch):
    source = ROOT / 'data/dry_interface_regression/water_callback_roundtrip.vpy'
    if not source.exists():
        pytest.skip('Local legacy archive is not present')
    guard_calculations(monkeypatch)
    session = load_session(source)
    assert session.active_network is not None
    assert all(row.state is CapabilityState.NOT_CALCULATED for row in session.capabilities(result='analysis'))


def test_version_one_archive_loads_without_reconstructing_caches(cached_session, tmp_path, monkeypatch):
    from vorpy.src.io.network_archive import V1_FIELDS, LINKS
    path = save_network(cached_session.active_network, tmp_path / 'legacy.vpy')
    with zipfile.ZipFile(path) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    state = json.loads(members['state.json'])
    records = [record for record in state['objects'] if record['kind'] in V1_FIELDS]
    identities = {record['id']: index for index, record in enumerate(records)}
    def remap(value):
        if isinstance(value, dict):
            if 'ref' in value:
                return {'ref': identities[value['ref']]}
            return {key: remap(item) for key, item in value.items()}
        if isinstance(value, list):
            return [remap(item) for item in value]
        return value
    for record in records:
        record['id'] = identities[record['id']]
        record.pop('schema_version')
        record['fields'] = remap({key: value for key, value in record['fields'].items()
                                 if key in V1_FIELDS[record['kind']] + LINKS[record['kind']]})
    metadata = json.loads(members['metadata.json'])
    metadata.update(format_version=1, root_network=identities[metadata['root_network']])
    metadata.pop('document_kind')
    metadata.pop('document_root')
    members['metadata.json'] = json.dumps(metadata).encode()
    members['state.json'] = json.dumps({'objects': records}).encode()
    for name in list(members):
        if name.startswith('tables/'):
            members[name] = json.dumps(remap(json.loads(members[name]))).encode()
    with zipfile.ZipFile(path, 'w') as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    guard_calculations(monkeypatch)
    session = load_session(path)
    iface = session.active_network.sys.ifaces[0]
    assert not hasattr(iface, 'geometry_analysis')
    assert session.capabilities(result='analysis')[0].state is CapabilityState.NOT_CALCULATED
    assert session.capabilities(representation='dual', result='geometry')[0].state is CapabilityState.NOT_CALCULATED


def test_multiple_systems_and_network_scopes(cached_session, tmp_path, monkeypatch):
    from copy import deepcopy
    second = deepcopy(cached_session)
    # Before the first save only the root network has received a stable ID.
    del second.active_network.archive_id
    second.active_network_id = None
    second.networks[0].iface_grps = None
    second.networks[0].network_mode = 'network'
    merged = Session(cached_session.systems + second.systems, cached_session.networks + second.networks,
        active_network_id=cached_session.active_network_id)
    archive = save_session(merged, tmp_path / 'multiple.vpy')
    guard_calculations(monkeypatch)
    restored = load_session(archive)
    assert len(restored.systems) == len(restored.networks) == 2
    assert [net.network_scope for net in restored.networks] == ['interface', 'system']
    for net in restored.networks:
        iface = net.sys.ifaces[0]
        assert iface._geometry_dual.network is net
        assert iface._geometry_filtration.network is net
        assert iface._geometry_representation.archive_provenance['network_id'] == net.archive_id
    assert restored.networks[0].archive_id != restored.networks[1].archive_id


def test_adapter_version_rejected(cached_session, tmp_path):
    path = save_session(cached_session, tmp_path / 'session.vpy')
    with zipfile.ZipFile(path) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
    state = json.loads(members['state.json'])
    next(record for record in state['objects'] if record['kind'] == 'dual_complex')['schema_version'] = 99
    members['state.json'] = json.dumps(state).encode()
    with zipfile.ZipFile(path, 'w') as archive:
        for name, data in members.items():
            archive.writestr(name, data)
    with pytest.raises(ArchiveError, match='schema version'):
        load_session(path)


def test_stale_and_inconsistent_state_is_not_repaired(cached_session, tmp_path, monkeypatch):
    from vorpy.src.interface.geometry_analysis import geometry_cache_key
    iface = cached_session.active_network.sys.ifaces[0]
    iface._geometry_analysis_key = geometry_cache_key(iface, 0.)
    iface.net.surfs = iface.net.surfs.copy()
    iface.geometry_analysis.alpha_selection['selected_physical_pairs'] = 99
    saved = save_session(cached_session, tmp_path / 'stale.vpy')
    guard_calculations(monkeypatch)
    session = load_session(saved)
    restored = session.active_network.sys.ifaces[0]
    assert restored.geometry_analysis.alpha_selection['selected_physical_pairs'] == 99
    assert len(restored._geometry_representation.selected_surfaces) == 1
    capabilities = session.capabilities(result='analysis')
    assert capabilities[0].state is CapabilityState.STALE
    assert capabilities[0].provenance_issues


def test_gui_session_save_open_keeps_scientific_results(cached_session, tmp_path, monkeypatch):
    from vorpy.tests.workbench.test_main_window import make_window
    from vorpy.workbench.ui import main_window
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.project import InterfaceDefinition
    net = cached_session.active_network
    group = net.sys.groups[0]
    group.net = net
    source = tmp_path / 'missing.pdb'
    window = make_window(monkeypatch)
    result = AnalysisResult(source, 'cached result', export_group=group)
    window.source = source
    window._display_result(result)
    window._groups = {'A': (0,), 'B': (1,2)}
    window._interfaces = {'Contact': ('A','B')}
    window.viewer.plotter.camera_position = [(3.,4.,5.),(0.,0.,0.),(0.,1.,0.)]
    window.project_file = tmp_path / 'session.vpy'
    window.export_panel.config.color_overrides = {'surfaces': '#123456'}
    assert window.save_project()
    restored = make_window(monkeypatch)
    monkeypatch.setattr(main_window.QFileDialog, 'getOpenFileName', lambda *args: (str(window.project_file), ''))
    guard_calculations(monkeypatch)
    restored.open_project()
    assert restored.current_result.export_group.net is restored._scientific_session.active_network
    assert restored.current_result.export_group.net.sys.ifaces[0].geometry_analysis.voronoi_side_1['total_H'] is None
    assert restored._groups == {'A': (0,), 'B': (1,2)}
    assert restored.export_panel.config.color_overrides == {'surfaces': '#123456'}
    assert restored.viewer.plotter.camera_position[0] == [3.,4.,5.]
    assert restored.export_panel.export_button.isEnabled()
    window.close()
    restored.close()


def test_gui_legacy_project_open_and_save_migrates(tmp_path, monkeypatch):
    from vorpy.tests.workbench.test_main_window import make_window
    from vorpy.workbench.ui import main_window
    from vorpy.workbench.project import Project, save_project
    legacy = tmp_path / 'legacy.vpyworkbench.json'
    save_project(Project(name='Legacy', view_state={'export_config': {'atom_format': 'xyz'}}), legacy)
    original = legacy.read_bytes()
    window = make_window(monkeypatch)
    monkeypatch.setattr(main_window.QFileDialog, 'getOpenFileName', lambda *args: (str(legacy), ''))
    window.open_project()
    assert window.export_panel.config.atom_format == 'xyz'
    assert window.save_project()
    assert window.project_file == tmp_path / 'legacy.vpy'
    assert legacy.read_bytes() == original
    assert load_session(window.project_file).workbench['project'].name == 'Legacy'
    window.close()


def prepare_validated_2kai():
    """Materialize the existing validated checkpoints BEFORE the guarded reload.

    No birth or curvature calculation: load saved births/components as-is,
    build incidence/selection once to form the pre-save reference objects.
    """
    import csv
    from dataclasses import replace
    from vorpy.src.geometry.duals import build_dual
    from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex
    from vorpy.src.geometry.interfaces import build_alpha_interface
    from vorpy.src.geometry.interfaces.curvature import CurvatureSummary
    def read(path):
        with path.open(newline='', encoding='utf-8-sig') as stream:
            return list(csv.DictReader(stream))
    cache = ROOT / 'output/comparative_study/matched_probe_control/2KAI'
    net = load_network(ROOT / 'output/comparative_study/cases/2KAI/H/full_aw.vpy')
    atoms = read(cache / 'input_identity_audit.csv')
    for row in atoms:
        generator = int(row['generator_id'])
        np.testing.assert_allclose(net.balls.loc[generator, 'loc'],
            [float(row[key]) for key in ('x_A', 'y_A', 'z_A')], rtol=0., atol=1e-12)
        assert net.balls.loc[generator, 'rad'] == float(row['base_radius_A'])
    a = frozenset(int(row['generator_id']) for row in atoms if row['chain'] in {'A','B'})
    b = frozenset(int(row['generator_id']) for row in atoms if row['chain'] == 'I')
    checkpoint = read(cache / 'aw_pair_birth_checkpoint.csv')
    thresholds = {float(row['threshold']) for row in checkpoint}
    assert len(thresholds) == 1
    alpha = thresholds.pop()
    records = {d: {} for d in range(4)}
    for row in checkpoint:
        pair = tuple(sorted((int(row['generator_i']), int(row['generator_j']))))
        birth = float(row['surface_birth_A']) if row['surface_birth_A'] else None
        records[1][pair] = AlphaSimplex(1, pair, birth, birth, row['birth_source'], birth is not None,
            notes=row['birth_reason'], diagnostics={
                'restriction_status': row['restriction_status'],
                'restriction_certificate': row['restriction_certificate'],
                'birth_lower_bound_A': float(row['birth_lower_bound_A']) if row['birth_lower_bound_A'] else None})
    filtration = AlphaFiltration(net, 'aw', records, {
        'native_alpha_units': 'A', 'filtration_kind': 'persisted experimental AW restricted-cell pair births',
        'birth_dimensions_calculated': 1, 'source': str(cache / 'aw_pair_birth_checkpoint.csv'),
        'stable_atom_metadata': {int(row['generator_id']): {
            key: row[key] for key in ('chain','residue_number','insertion_code','residue_name','atom_name','stable_atom_id','element')}
            for row in atoms}}, tolerance=1e-6)
    dual = build_dual(net)
    rep = build_alpha_interface(net, dual, filtration, alpha, a, b, calculate_curvature=False)
    curvature_root = ROOT / 'output/2KAI_curvature_showcase'
    surfaces = {int(row['surface_id']): row for row in read(curvature_root / 'aw_surface_curvature_components.csv')}
    rep.surfaces = [replace(surface, smooth_H_A=(float(surfaces[surface.feature_id]['integral_H_dA_A'])
        if surface.feature_id in surfaces and surfaces[surface.feature_id]['support_status'] == 'included'
        and surfaces[surface.feature_id]['integral_H_dA_A'] else None)) for surface in rep.surfaces]
    measured_edges = {int(row['edge_id']): row for row in read(curvature_root / 'aw_edge_curvature_components.csv')}
    rep.edges = [replace(edge,
        length_A=float(measured_edges[edge.feature_id]['length_A']) if measured_edges.get(edge.feature_id, {}).get('length_A') else None,
        signed_integral_A_rad=float(measured_edges[edge.feature_id]['raw_signed_contribution_A_rad'])
            if measured_edges.get(edge.feature_id, {}).get('raw_signed_contribution_A_rad') else None,
        curvature_status=measured_edges.get(edge.feature_id, {}).get('support_status', 'unresolved')) for edge in rep.edges]
    component_rows = read(curvature_root / 'open_interface_curvature/aw_curvature_components.csv')
    components = {row['term']: float(row['value']) if row['value'] else None for row in component_rows}
    topology = read(curvature_root / 'open_interface_curvature/aw_topology_audit.csv')
    selected_topology = next(row for row in topology
        if row['scope'] == 'selected_physical_interface' and row['component'] == 'TOTAL')
    selected_components = [row for row in topology
        if row['scope'] == 'selected_physical_interface' and row['component'] != 'TOTAL']
    side = {
        'area': components['area_A2'],
        'topology_vertices': int(selected_topology['V_physical']),
        'topology_edges': int(selected_topology['E_physical']),
        'topology_faces': int(selected_topology['F_physical']),
        'euler_characteristic': int(selected_topology['chi']),
        'connected_components': len(selected_components), 'boundary_loops': int(selected_topology['boundary_loops']),
        'smooth_H': components['smooth_face_H'], 'edge_H': components['internal_edge_H_half_raw'],
        'raw_signed_edge': components['internal_crease_raw_signed_beta_ds'],
        'boundary_integrated_mean_curvature': components['free_boundary_H'],
        'partial_H': components['partial_H_sum'], 'total_H': components['total_H'],
        'smooth_gaussian_curvature': components['smooth_K'],
        'intrinsic_seam_gaussian_curvature': components['internal_seam_face_geodesic_K'],
        'boundary_geodesic_gaussian': components['free_boundary_face_geodesic_K'],
        'interior_vertex_gaussian_curvature': components['interior_vertex_defect_K'],
        'boundary_corner_gaussian': components['boundary_vertex_corner_K'],
        'gauss_bonnet_residual': components['supported_subcomplex_GB_residual_order_256'],
        'gauss_bonnet_lhs': components['supported_subcomplex_GB_LHS_order_256'],
        'gauss_bonnet_expected': components['supported_subcomplex_two_pi_chi'],
        'total_K': components['total_Gaussian_GB'], 'mean_certified': False,
        'gaussian_total_certified': False,
    }
    from dataclasses import fields
    values = {item.name: None for item in fields(CurvatureSummary)}
    values.update(scheme='aw', selection_mode='alpha', area_A2=side['area'], smooth_H_A=side['smooth_H'],
        raw_signed_edge_A_rad=side['raw_signed_edge'], conventional_edge_H_A=side['edge_H'],
        combined_H_A=None, status='partial', source=str(curvature_root / 'open_interface_curvature/aw_curvature_components.csv'))
    rep.curvature = CurvatureSummary(**values)
    iface = attach_interface(net, a, b, '2KAI_cached_AB_I')
    selected = sorted(row.generator_tuple for row in filtration.interface_at(alpha, a, b))
    iface._geometry_dual, iface._geometry_filtration, iface._geometry_representation = dual, filtration, rep
    iface.geometry_analysis = InterfaceGeometryAnalysis(metadata={
        'alpha_value': alpha, 'alpha_units': 'A', 'scheme': 'aw', 'solvent_context': False,
        'additional_network_solves': 0, 'source': 'existing validated 2KAI checkpoints',
        'component_table_records': component_rows, 'topology_audit_records': topology},
        alpha_selection={'selected_generator_pairs': [list(pair) for pair in selected],
            'selected_physical_pairs': len(rep.selected_surfaces),
            'eligible_bicolor_pairs': sum(surface.supported and surface.complete and surface.bounded
                and surface.area_A2 is not None for surface in rep.candidate_surfaces), 'status': 'partial'},
        voronoi_side_1=side, voronoi_side_2=reverse_perspective(side),
        coverage={'mean_coverage': components['mean_surface_coverage_fraction'],
            'internal_edge_coverage': components['mean_internal_edge_coverage_fraction']},
        unresolved=[row['support'] for row in component_rows if 'UNRESOLVED' in row['support']])
    return Session.from_network(net)


@pytest.mark.slow
@pytest.mark.timeout(600)
def test_validated_2kai_fresh_process_export_fidelity(tmp_path):
    source = ROOT / 'output/comparative_study/cases/2KAI/H/full_aw.vpy'
    if not source.exists():
        pytest.skip('Local validated 2KAI scientific checkpoint is not present')
    session = prepare_validated_2kai()
    archive = save_session(session, tmp_path / '2kai.vpy')
    expected = scientific_snapshot(session)
    side = expected['analysis']['voronoi_side_1']
    assert [side[key] for key in ('topology_vertices','topology_edges','topology_faces',
        'euler_characteristic','connected_components','boundary_loops')] == [584,924,339,-1,2,5]
    assert len(session.active_network.sys.ifaces[-1]._geometry_representation.selected_surfaces) == 339
    assert side['partial_H'] == pytest.approx(-54.150160911300546, abs=1e-12)
    assert side['area'] == pytest.approx(787.0925777403216, abs=1e-12)
    (tmp_path / 'before.json').write_text(json.dumps(expected, sort_keys=True), encoding='utf-8')
    export_cached(session, tmp_path / 'before')
    del session
    code = '''from pathlib import Path
import json, sys, pytest
from vorpy.tests.io.test_scientific_session import guard_calculations, scientific_snapshot, export_cached
from vorpy.src.io import load_session
with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    session = load_session(sys.argv[1])
    snapshot = scientific_snapshot(session)
    assert session.active_network.sys.ifaces[-1].geometry_analysis.voronoi_side_1['total_H'] is None
    Path(sys.argv[2]).write_text(json.dumps(snapshot, sort_keys=True), encoding='utf-8')
    export_cached(session, Path(sys.argv[3]))
'''
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    completed = subprocess.run([sys.executable, '-c', code, str(archive), str(tmp_path / 'after.json'),
        str(tmp_path / 'after')], cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
    assert completed.returncode == 0, completed.stderr
    assert json.loads((tmp_path / 'after.json').read_text()) == expected
    count = compare_exports(tmp_path / 'before', tmp_path / 'after')
    print(f'2KAI: {count} export files equivalent; network {expected["counts"]}; '
          f'partial H={expected["analysis"]["voronoi_side_1"]["partial_H"]}; complete H/K unresolved')


def test_visualize_export_fresh_process_is_cache_only(cached_session, tmp_path):
    from vorpy.src.output.visualization import export_visualization_bundle

    before_root = tmp_path / 'before_visualize'
    after_root = tmp_path / 'after_visualize'
    cached_session.active_network.sys.type = 'pdb'
    for iface in cached_session.active_network.sys.ifaces:
        iface.group1.ball_ndxs = []
        iface.group2.ball_ndxs = []
    cached_session.active_network.sys.files['dir'] = str(before_root)
    before = export_visualization_bundle(cached_session.active_network.sys)
    archive = save_session(cached_session, tmp_path / 'visualize.vpy')

    code = '''
import json, sys
from pathlib import Path
import pytest
from vorpy.tests.io.test_scientific_session import guard_calculations
from vorpy.src.io import load_session
from vorpy.src.output.visualization import export_visualization_bundle

with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    session = load_session(sys.argv[1])
    session.active_network.sys.type = 'pdb'
    session.active_network.sys.files['dir'] = sys.argv[2]
    items = export_visualization_bundle(session.active_network.sys)
    item = items[0]
    assert Path(item['launcher']).is_file()
    assert Path(item['manifest']).is_file()
    Path(sys.argv[3]).write_text(json.dumps({
        name: layer['state'] for name, layer in item['layers'].items()
    }, sort_keys=True), encoding='utf-8')
'''
    state_path = tmp_path / 'after_visualize.json'
    completed = subprocess.run(
        [sys.executable, '-c', code, str(archive), str(after_root), str(state_path)],
        cwd=ROOT,
        env=dict(os.environ, PYTHONPATH=str(ROOT)),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    after = json.loads(state_path.read_text(encoding='utf-8'))
    assert before[0]['layers']['partition_interface']['state'] == after['partition_interface']
    assert before[0]['layers']['dual_contact_complex']['state'] == after['dual_contact_complex']


def test_compact_visualize_export_fresh_process_is_cache_only(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'compact.vpy')
    code = '''
import csv
import sys
from pathlib import Path
import pytest
from vorpy.tests.io.test_scientific_session import guard_calculations
from vorpy.src.io import load_session
from vorpy.src.output.visualization import export_compact_visualization_bundle

with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    session = load_session(sys.argv[1])
    system = session.active_network.sys
    for iface in system.ifaces:
        iface.group1.ball_ndxs = []
        iface.group2.ball_ndxs = []
    system.files['dir'] = sys.argv[2]
    system._loaded_from_archive = True
    item = export_compact_visualization_bundle(system)[0]
    assert item['launcher'] is None
    root = Path(sys.argv[2])
    assert not (root / 'interface_A_B' / 'interface.pml').exists()
    log_path = root / 'interface_A_B' / 'logs.csv'
    assert log_path.is_file()
    with log_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    assert reader.fieldnames == ['interface_id', 'representation', 'side', 'quantity', 'value', 'units', 'status', 'provenance']
    assert rows and all(row['status'] != 'AVAILABLE' for row in rows)
    assert {'physical_partition_interface', 'molecular_contact_surface'} <= {
        row['representation'] for row in rows
    }
    assert (root / 'interface_A_B' / 'aw_surfs.off').stat().st_size > 0
    assert not (root / 'groups').exists()
    assert not (root / 'interfaces').exists()
'''
    completed = subprocess.run(
        [sys.executable, '-c', code, str(archive), str(tmp_path / 'compact_after')],
        cwd=ROOT,
        env=dict(os.environ, PYTHONPATH=str(ROOT)),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr


def test_small_preset_fresh_process_is_cache_only(cached_session, tmp_path):
    archive = save_session(cached_session, tmp_path / 'preset.vpy')
    code = '''
import csv
import json
import sys
from pathlib import Path
import pytest
from vorpy.tests.io.test_scientific_session import guard_calculations
from vorpy.src.io import load_session
from vorpy.src.output.output import export_preset

with pytest.MonkeyPatch.context() as guard:
    guard_calculations(guard)
    session = load_session(sys.argv[1])
    iface = next(iface for iface in session.active_network.sys.ifaces
                 if iface.net is session.active_network)
    archive_provenance = dict(getattr(iface.geometry_analysis, 'archive_provenance', {}) or {})
    assert archive_provenance
    from vorpy.src.output.visualization import _cached_interface_result
    cached_result = _cached_interface_result(iface)
    assert dict(cached_result.provenance.details) == archive_provenance
    system = session.active_network.sys
    for network in session.networks:
        network.settings.setdefault('surf_scheme', 'int_mean_curv')
        network.settings.setdefault('surf_col', 'coolwarm')
        network.settings.setdefault('scheme_factor', 'log')
        network.settings.setdefault('edge_col', 'yellow')
        network.settings.setdefault('vert_col', 'red')
    root = Path(sys.argv[2])
    root.mkdir(parents=True, exist_ok=True)
    system.files['dir'] = str(root)
    system._loaded_from_archive = True
    export_preset(system, 'small')
    assert (root / 'fixture.vpy').is_file()
    assert (root / 'info.txt').is_file()
    assert list(root.rglob('*.off'))
    assert list(root.rglob('logs.csv'))
    assert not list(root.rglob('*.pdb'))
    assert not list(root.rglob('*.pml'))
    assert not list(root.rglob('waters'))
    log_path = next(root.rglob('logs.csv'))
    with log_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ['interface_id', 'representation', 'side', 'quantity', 'value', 'units', 'status', 'provenance']
        rows = list(reader)
        assert any(json.loads(row['provenance']).get('details') == archive_provenance
                   for row in rows if row['provenance'])
        assert all(row['status'] != 'AVAILABLE' for row in rows)
'''
    completed = subprocess.run(
        [sys.executable, '-c', code, str(archive), str(tmp_path / 'preset_after')],
        cwd=ROOT,
        env=dict(os.environ, PYTHONPATH=str(ROOT)),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
