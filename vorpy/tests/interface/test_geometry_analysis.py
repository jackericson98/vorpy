import csv
import json
from dataclasses import replace
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.interface.interface import Interface
from vorpy.src.interface.geometry_analysis import (
    analyze_interface_geometry, reverse_perspective, solvent_environment,
    summarize_physical_interface, write_geometry_info, write_geometry_log,
)
from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.interfaces import build_alpha_interface
from vorpy.tests.geometry.test_alpha_interfaces import _network, _filtration
from vorpy.tests.inputs.test_log_column_order import write_log
from vorpy.src.log_columns import LogWriter
from vorpy.src.inputs.logs import read_logs
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.analyze.open_interface_curvature import gauss_bonnet_accounting, integrated_mean_curvature
from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis


def iface_for(network):
    system_balls = network.balls.copy()
    system_balls['num'] = system_balls.index
    system_balls['res_name'] = 'ALA'
    iface = Interface.__new__(Interface)
    iface.net = network
    # The fixture's num is an intentionally different parent-system identifier.
    iface.group1_indices = {101}
    iface.group2_indices = {102, 103}
    iface.interface_id = 'A_B'
    iface.geometry_analysis = None
    iface.sys = SimpleNamespace(balls=system_balls, sol=None)
    return iface


def install_filtration(monkeypatch, network, scheme='pow', **kwargs):
    filtration = _filtration(network, scheme, **kwargs)
    calls = []
    def build(*args, **options):
        calls.append(options)
        return filtration
    monkeypatch.setattr('vorpy.src.geometry.filtrations.alpha.build_alpha_filtration', build)
    return calls


def test_cache_selection_and_planar_partial_gaussian(monkeypatch):
    network = _network()
    calls = install_filtration(monkeypatch, network)
    iface = iface_for(network)
    result = iface.analyze_geometry()
    assert iface.geometry_analysis is result
    assert iface.analyze_geometry() is result
    assert len(calls) == 1
    assert result.alpha_selection['selected_pairs'] == 1
    assert result.alpha_selection['selected_physical_pairs'] == 1
    assert result.voronoi_side_1['area'] == 3.
    assert result.voronoi_side_1['smooth_H'] == 0.
    assert result.voronoi_side_1['smooth_K_partial'] == 0.
    assert result.voronoi_side_1['total_K'] is None
    assert result.voronoi_side_1['gaussian_total_certified'] is False
    assert result.metadata['solvent_context'] is False
    assert not hasattr(result, 'alpha_side_1')
    assert result.metadata['molecular_contact_patch_status'] == 'not_implemented'


def test_static_matched_probe_2kai_curvature_regression_fixture():
    fixture_path = Path(__file__).parent / 'fixtures' / 'phase1_2kai_matched_probe_curvature.json'
    fixture = json.loads(fixture_path.read_text(encoding='utf-8'))
    assert fixture['selected_faces'] == 350
    assert fixture['topology_vertices'] - fixture['topology_edges'] + fixture['selected_faces'] == fixture['euler_characteristic'] == 0
    assert fixture['components'] == 2 and fixture['boundary_loops'] == 4
    assert fixture['area_A2'] == pytest.approx(898.1186, abs=5e-5)
    mean = integrated_mean_curvature(
        [fixture['smooth_mean_curvature_A']], [fixture['internal_crease_mean_curvature_A']],
    )
    assert mean['total_H'] == pytest.approx(-64.450496)
    assert mean['certified'] is fixture['certified_mean_total'] is True
    gb = gauss_bonnet_accounting(
        smooth_K=0., face_boundary_geodesic=0., vertex_corner_terms=0.,
        euler_characteristic=fixture['euler_characteristic'],
    )
    assert gb['certified'] is fixture['certified_gaussian_total'] is True

    aw = fixture['aw']
    assert aw['selected_faces'] == 339
    assert aw['topology_vertices'] - aw['topology_edges'] + aw['selected_faces'] == aw['euler_characteristic'] == -1
    assert aw['components'] == 2 and aw['boundary_loops'] == 5
    assert aw['smooth_mean_curvature_A'] == pytest.approx(1.621963)
    assert aw['internal_crease_mean_curvature_A'] == pytest.approx(-55.772124)
    assert aw['smooth_mean_curvature_A'] + aw['internal_crease_mean_curvature_A'] == pytest.approx(aw['partial_mean_curvature_A'])
    assert aw['supported_gauss_bonnet_residual'] == pytest.approx(-0.0006654952459843599)
    aw_mean = integrated_mean_curvature(
        [aw['smooth_mean_curvature_A']], [aw['internal_crease_mean_curvature_A']],
        unresolved_faces=aw['unsupported_surfaces'],
        unresolved_internal_edges=aw['unresolved_internal_edges'],
    )
    assert aw_mean['total_H'] is None and not aw_mean['certified']

    phase1 = InterfaceGeometryAnalysis(voronoi_side_1={
        'topology_vertices': fixture['topology_vertices'], 'topology_edges': fixture['topology_edges'],
        'topology_faces': fixture['selected_faces'], 'euler_characteristic': fixture['euler_characteristic'],
        'components': fixture['components'], 'boundary_loops': fixture['boundary_loops'],
        'area': fixture['area_A2'], 'smooth_H': fixture['smooth_mean_curvature_A'],
        'edge_H': fixture['internal_crease_mean_curvature_A'], 'total_H': fixture['total_mean_curvature_A'],
        'smooth_mean_curvature': fixture['smooth_mean_curvature_A'],
        'internal_crease_mean_curvature': fixture['internal_crease_mean_curvature_A'],
        'boundary_mean_curvature': fixture['boundary_mean_curvature_A'],
        'total_mean_curvature': fixture['total_mean_curvature_A'],
        'certified_mean_total': True, 'mean_certified': True,
        'smooth_gaussian_curvature': 0., 'intrinsic_vertex_gaussian': 0.,
        'boundary_geodesic_gaussian': 0., 'boundary_corner_gaussian': 0.,
        'gauss_bonnet_lhs': fixture['gauss_bonnet_lhs'],
        'gauss_bonnet_expected': fixture['gauss_bonnet_expected'],
        'gauss_bonnet_residual': fixture['gauss_bonnet_residual'],
        'total_gaussian_curvature': 0., 'total_K': 0.,
        'certified_gaussian_total': True, 'gaussian_total_certified': True,
    })
    reverse = reverse_perspective(phase1.voronoi_side_1)
    assert reverse['total_mean_curvature'] == -fixture['total_mean_curvature_A']
    assert reverse['euler_characteristic'] == fixture['euler_characteristic']
    assert reverse['total_gaussian_curvature'] == 0.
    phase1.voronoi_side_2 = reverse
    record = phase1.log_record()
    assert record['vor_topology_vertices'] == 605
    assert record['vor_g1_to_g2_certified_mean_total'] is True
    assert record['vor_g1_to_g2_total_gaussian_curvature'] == 0.
    assert record['vor_g2_to_g1_total_mean_curvature'] == 64.450496


def test_aw_reuses_existing_values_and_reverse_geometry(monkeypatch):
    network = _network('aw')
    network.surfs['int_gauss_curv'] = [2.25, 3., 4.]
    install_filtration(monkeypatch, network, 'aw')
    calls = []
    def analyzer(*args, **kwargs):
        calls.append(args)
        return SimpleNamespace(selected_surfaces=[SimpleNamespace(surface_id=20, integrated_mean_curvature=.75)], selected_edges=[])
    monkeypatch.setattr('vorpy.src.analyze.aw_interface_curvature.analyze_aw_interface_curvature', analyzer)
    result = analyze_interface_geometry(iface_for(network))
    forward, reverse = result.voronoi_side_1, result.voronoi_side_2
    assert len(calls) == 1
    assert forward['smooth_H'] == .75
    assert reverse['smooth_H'] == -.75
    assert reverse['smooth_H_partial'] == -.75
    assert forward['total_H'] is reverse['total_H'] is None
    assert forward['smooth_mean_curvature'] == .75
    assert reverse['smooth_mean_curvature'] == -.75
    assert not forward['certified_mean_total']
    assert forward['area'] == reverse['area'] == 3.
    assert forward['unsigned_edge'] == reverse['unsigned_edge'] == 0.
    assert forward['smooth_K_partial'] == reverse['smooth_K_partial'] == 2.25
    assert not forward['gaussian_total_certified']


def test_unresolved_seam_never_becomes_certified_zero():
    network = _network()
    representation = build_alpha_interface(network, build_dual(network), _filtration(network, second_pair_birth=-1.), 0., {10}, {13, 16})
    result = summarize_physical_interface(representation, network)
    assert result['unresolved_edges'] == 1
    assert result['supported_edges'] == 0
    assert result['raw_signed_edge'] is None
    assert result['edge_H'] is None
    assert result['total_H'] is None
    assert not result['mean_certified']
    assert result['mean_coverage'] == pytest.approx(2 / 3)


def test_signed_edge_reversal_and_pointwise_positive_negative():
    network = _network()
    rep = build_alpha_interface(network, build_dual(network), _filtration(network, second_pair_birth=-1.), 0., {10}, {13, 16})
    rep.edges[0] = replace(rep.edges[0], signed_integral_A_rad=-2., unsigned_integral_A_rad=6., curvature_status='integrated')
    rep.curvature = SimpleNamespace(positive_edge_A_rad=2., negative_edge_A_rad=-4.)
    side = summarize_physical_interface(rep, network)
    reverse = reverse_perspective(side)
    assert side['total_H'] is None
    assert side['internal_crease_mean_curvature'] == -1.
    assert reverse['internal_crease_mean_curvature'] == 1.
    assert reverse['total_H'] is None
    assert reverse['partial_H'] == 1.
    assert side['unsigned_edge'] == reverse['unsigned_edge'] == 6.
    assert reverse['positive_edge'] == 4.
    assert reverse['negative_edge'] == -2.


def test_solvent_detection_excludes_virtual_shell(monkeypatch):
    network = _network()
    iface = iface_for(network)
    iface.sys.balls = pd.DataFrame([
        {'num': 101, 'res_name': 'ALA', 'res_seq': 1},
        {'num': 102, 'res_name': 'ALA', 'res_seq': 2},
        {'num': 103, 'res_name': 'ALA', 'res_seq': 3},
        {'num': 104, 'res_name': 'HOH', 'res_seq': 4, 'is_boundary_generator': True},
    ], index=[101, 102, 103, 104]).fillna(False)
    network.balls.loc[19] = network.balls.loc[16].copy()
    network.balls.loc[19, 'num'] = 104
    network.balls['is_boundary_generator'] = [False, False, False, True]
    environment, ids = solvent_environment(iface)
    assert not environment['solvent_detected']
    assert not ids
    iface.sys.balls.loc[104, 'is_boundary_generator'] = False
    network.balls.loc[19, 'is_boundary_generator'] = False
    environment, ids = solvent_environment(iface)
    assert environment['solvent_detected']
    assert environment['solvent_atom_count'] == 1
    assert environment['solvent_residue_count'] == 1
    assert ids == {19}
    install_filtration(monkeypatch, network)
    result = analyze_interface_geometry(iface)
    assert result.metadata['solvent_context']
    assert not result.dry_solvent_comparison['available']
    assert result.metadata['additional_network_solves'] == 0


@pytest.mark.parametrize('reader', [read_logs, read_logs2])
def test_appended_log_summary_preserves_standard_sections(tmp_path, monkeypatch, reader):
    network = _network()
    install_filtration(monkeypatch, network)
    result = analyze_interface_geometry(iface_for(network))
    path = tmp_path / 'interface.csv'
    write_log(path, LogWriter)
    original = reader(str(path), return_dict=True)
    with path.open('a', newline='') as stream:
        write_geometry_log(LogWriter(stream), result)
    parsed = reader(str(path), return_dict=True)
    extension = parsed.pop('interface geometry')
    assert parsed == original
    assert extension['vor_g1_to_g2_smooth_H'] == 0.
    assert extension['vor_total_K'] is None
    assert extension['vor_certified_gaussian'] is False
    assert extension['alpha_selected_pairs'] == 1
    assert extension == result.log_record()
    text = StringIO()
    write_geometry_info(text, result)
    assert 'ALPHA-SELECTED INTERFACE GEOMETRY' in text.getvalue()
    assert 'PARTIAL / UNRESOLVED' in text.getvalue()
    assert 'Gauss-Bonnet' in text.getvalue()


def test_interface_build_populates_analysis(monkeypatch):
    network = _network()
    network.build = lambda: None
    install_filtration(monkeypatch, network)
    iface = iface_for(network)
    iface.sys.cache_interface_geometry = lambda _iface: None
    iface.water_topology = None
    iface.name = 'example'
    iface._update_group_metadata = lambda **kwargs: None
    monkeypatch.setattr('vorpy.src.interface.interface.analyze_interface_waters', lambda **kwargs: [])
    monkeypatch.setattr('vorpy.src.interface.interface.build_buried_water_groups', lambda _iface: [])
    iface.build()
    assert iface.geometry_analysis is not None
    assert iface.geometry_analysis.voronoi_side_1['area'] == 3.


def test_interface_build_can_skip_optional_water_pipeline(monkeypatch):
    network = _network()
    network.build = lambda: None
    install_filtration(monkeypatch, network)
    iface = iface_for(network)
    iface.sys.cache_interface_geometry = lambda _iface: None
    iface._update_group_metadata = lambda **kwargs: None
    monkeypatch.setattr(
        'vorpy.src.interface.interface.analyze_interface_waters',
        lambda **kwargs: (_ for _ in ()).throw(AssertionError('water analysis ran')),
    )
    monkeypatch.setattr(
        'vorpy.src.interface.interface.build_buried_water_groups',
        lambda _iface: (_ for _ in ()).throw(AssertionError('buried solve ran')),
    )
    iface.build(analyze_waters=False)
    assert iface.geometry_analysis is not None
    assert iface.water_topology is None
    assert iface.buried_water_groups == []


def test_missing_dependency_is_reported_not_zero_filled(monkeypatch):
    def unavailable(*args, **kwargs):
        raise RuntimeError('GUDHI unavailable')
    monkeypatch.setattr('vorpy.src.geometry.filtrations.alpha.build_alpha_filtration', unavailable)
    result = analyze_interface_geometry(iface_for(_network()))
    assert result.voronoi_side_1['total_H'] is None
    assert not result.voronoi_side_1['mean_certified']
    assert 'GUDHI unavailable' in result.unresolved[0]


def test_interface_selection_can_skip_optional_curvature_and_topology(monkeypatch):
    network = _network()
    install_filtration(monkeypatch, network)
    result = analyze_interface_geometry(iface_for(network), calculate_curvature=False)
    assert result.metadata['curvature_requested'] is False
    assert result.voronoi_side_1['area'] == 3.
    assert result.voronoi_side_1['total_H'] is None
    assert result.voronoi_side_1['euler_characteristic'] is None


def test_primitive_planar_behavior(monkeypatch):
    network = _network('prm')
    install_filtration(monkeypatch, network, 'prm')
    result = analyze_interface_geometry(iface_for(network))
    assert result.voronoi_side_1['smooth_H'] == 0.
    assert result.voronoi_side_1['smooth_K_partial'] == 0.
    assert result.metadata['alpha_units'] == 'A^2'


def test_comparison_requires_matching_genuine_dry_cache(monkeypatch):
    network = _network()
    install_filtration(monkeypatch, network)
    dry_iface = iface_for(network)
    dry = analyze_interface_geometry(dry_iface)
    assert not dry.voronoi_side_1['mean_certified']
    wet_network = _network()
    wet_network.balls.loc[19] = wet_network.balls.loc[16].copy()
    wet_network.balls.at[19, 'num'] = 104
    wet_network.balls.at[19, 'loc'] = (5., 5., 5.)
    wet_iface = iface_for(wet_network)
    wet_iface.sys = dry_iface.sys
    wet_iface.sys.balls.loc[104, 'res_name'] = 'WAT'
    wet_iface.sys.balls.loc[104, 'res_seq'] = 4
    install_filtration(monkeypatch, wet_network)
    wet = analyze_interface_geometry(wet_iface)
    comparison = wet.dry_solvent_comparison
    assert comparison['available']
    assert comparison['pair_jaccard'] == 1.
    assert comparison['delta_H'] is None
    assert comparison['delta_K'] is None
    assert wet.metadata['additional_network_solves'] == 0
    # Changed molecular radii cannot reuse the previous dry analysis.
    wet_network.balls.at[10, 'rad'] = 1.1
    mismatch = analyze_interface_geometry(wet_iface, refresh=True)
    assert not mismatch.dry_solvent_comparison['available']


def test_cached_aw_surface_partial_does_not_certify_missing_support():
    network = _network('aw')
    rep = build_alpha_interface(network, build_dual(network), _filtration(network, 'aw'), 0., {10}, {13, 16}, calculate_curvature=False)
    network.surfs['int_mean_curv_by_ball'] = [{10: -.75, 13: .75}, {}, {}]
    result = summarize_physical_interface(rep, network)
    assert result['smooth_H_partial'] == -.75
    assert result['smooth_H'] is None
    assert result['total_H'] is None
    assert not result['mean_certified']
    assert reverse_perspective(result)['smooth_H_partial'] == .75


def test_cache_is_available_before_later_water_solver_failure(monkeypatch):
    network = _network()
    network.build = lambda: None
    install_filtration(monkeypatch, network)
    iface = iface_for(network)
    iface.sys.cache_interface_geometry = lambda _iface: None
    iface.water_topology = None
    iface.name = 'example'
    monkeypatch.setattr('vorpy.src.interface.interface.analyze_interface_waters', lambda **kwargs: [])
    def fail(_iface):
        raise RuntimeError('water solver failure')
    monkeypatch.setattr('vorpy.src.interface.interface.build_buried_water_groups', fail)
    with pytest.raises(RuntimeError, match='water solver failure'):
        iface.build()
    assert iface.geometry_analysis.voronoi_side_1['area'] == 3.
