"""Cached adapters for alpha selection and validated open-interface accounting.

The mathematical kernel remains independent; exporters only serialize this cache.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
import hashlib
import json
import math

from vorpy.src.boundary import SOLVENT_RESIDUES


PHYSICAL_FIELDS = (
    'interface_atoms', 'components', 'area', 'area_scope', 'area_selection_complete',
    'smooth_H', 'smooth_H_partial', 'raw_signed_edge', 'edge_H', 'partial_H', 'total_H',
    'unsigned_edge', 'positive_edge', 'negative_edge', 'edge_candidates', 'boundary_edges',
    'supported_edges', 'unresolved_edges', 'supported_surfaces', 'unresolved_surfaces',
    'mean_coverage', 'mean_certified', 'coverage_basis', 'smooth_K_partial', 'smooth_K_coverage',
    'total_K', 'gaussian_total_certified', 'euler_incidence', 'boundary_policy', 'unresolved_reasons',
    'selected_pairs',
    'topology_vertices', 'topology_edges', 'topology_faces', 'euler_characteristic',
    'boundary_loops', 'smooth_mean_curvature', 'internal_crease_mean_curvature',
    'boundary_mean_curvature', 'total_mean_curvature', 'certified_mean_total',
    'smooth_gaussian_curvature', 'intrinsic_vertex_gaussian',
    'boundary_geodesic_gaussian', 'boundary_corner_gaussian', 'gauss_bonnet_lhs',
    'gauss_bonnet_expected', 'gauss_bonnet_residual', 'total_gaussian_curvature',
    'certified_gaussian_total',
)
OPEN_INTERFACE_FIELDS = (
    'connected_components', 'genus_by_component', 'manifold_incidence_certified',
    'smooth_integrated_mean_curvature', 'internal_seam_integrated_mean_curvature',
    'boundary_integrated_mean_curvature', 'partial_integrated_mean_curvature',
    'total_integrated_mean_curvature', 'mean_curvature_certified',
    'mean_supported_faces', 'mean_unsupported_faces', 'mean_supported_internal_seams',
    'mean_unresolved_internal_seams', 'smooth_integrated_gaussian_curvature',
    'intrinsic_seam_gaussian_curvature', 'interior_vertex_gaussian_curvature',
    'intrinsic_integrated_gaussian_curvature', 'boundary_geodesic_curvature',
    'boundary_corner_turning', 'gaussian_curvature_certified',
    'gaussian_accounting_scope', 'gaussian_supported_faces', 'gaussian_unsupported_faces',
    'gaussian_coverage', 'supported_subcomplex_gauss_bonnet',
)
PHYSICAL_FIELDS += OPEN_INTERFACE_FIELDS


@dataclass
class InterfaceGeometryAnalysis:
    metadata: dict = field(default_factory=dict)
    alpha_selection: dict = field(default_factory=dict)
    voronoi_side_1: dict = field(default_factory=dict)
    voronoi_side_2: dict = field(default_factory=dict)
    dry_solvent_comparison: dict = field(default_factory=dict)
    coverage: dict = field(default_factory=dict)
    unresolved: list[str] = field(default_factory=list)

    def to_dict(self):
        return asdict(self)

    def log_record(self):
        m, a, v, reverse = self.metadata, self.alpha_selection, self.voronoi_side_1, self.voronoi_side_2
        record = {key: m.get(key) for key in (
            'schema_version', 'certification_definition', 'interface_id', 'scheme', 'alpha_convention', 'alpha_value', 'alpha_units',
            'solvent_context', 'solvent_detected', 'solvent_atom_count', 'solvent_residue_count',
            'solvent_definition', 'solvent_competing_generator_count', 'additional_network_solves', 'molecular_contact_patch_status')}
        record.update({f'alpha_{key}': value for key, value in a.items()})
        record.update({f'vor_{key}': value for key, value in v.items() if key != 'orientation'})
        for prefix, side in (('vor_g1_to_g2', v), ('vor_g2_to_g1', reverse)):
            for name in OPEN_INTERFACE_FIELDS:
                record[f'{prefix}_{name}'] = side.get(name)
            for term, name in (('smooth_H', 'smooth_H'), ('edge_H', 'edge_H'), ('total_H', 'total_H'),
                               ('smooth_mean_curvature', 'smooth_mean_curvature'),
                               ('internal_crease_mean_curvature', 'internal_crease_mean_curvature'),
                               ('boundary_mean_curvature', 'boundary_mean_curvature'),
                               ('total_mean_curvature', 'total_mean_curvature'),
                               ('certified_mean_total', 'certified_mean_total'),
                               ('smooth_gaussian_curvature', 'smooth_gaussian_curvature'),
                               ('intrinsic_vertex_gaussian', 'intrinsic_vertex_gaussian'),
                               ('boundary_geodesic_gaussian', 'boundary_geodesic_gaussian'),
                               ('boundary_corner_gaussian', 'boundary_corner_gaussian'),
                               ('gauss_bonnet_lhs', 'gauss_bonnet_lhs'),
                               ('gauss_bonnet_expected', 'gauss_bonnet_expected'),
                               ('gauss_bonnet_residual', 'gauss_bonnet_residual'),
                               ('total_gaussian_curvature', 'total_gaussian_curvature'),
                               ('certified_gaussian_total', 'certified_gaussian_total')):
                record[f'{prefix}_{name}'] = side.get(term)
        record['vor_total_K'] = v.get('total_K')
        record['vor_certified_gaussian'] = v.get('gaussian_total_certified', False)
        record['vor_certified_mean'] = v.get('mean_certified', False)
        record.update({f'solv_{key}': value for key, value in self.dry_solvent_comparison.items()})
        record['unresolved_reasons'] = self.unresolved
        return record


def finite(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _parent_id(index, row):
    for name in ('system_num', 'num'):
        value = finite(row.get(name))
        if value is not None:
            return int(value)
    return int(index)


def _boundary_flag(row):
    return str(row.get('is_boundary_generator', False)).strip().lower() in {'true', '1'}


def _system_pairs(network, pairs):
    return [sorted(_parent_id(g, network.balls.loc[g]) for g in pair) for pair in pairs]


def solvent_environment(iface):
    """Use physical System solvent metadata, never the generic surrounding group."""
    system = iface.sys
    balls = getattr(system, 'balls', None)
    boundary = set(getattr(iface.net, 'boundary_indices', ()) or ())
    classified, residue_sets = set(), []
    boundary_parent = set()
    for index, row in iface.net.balls.iterrows():
        if int(index) in boundary or _boundary_flag(row):
            boundary_parent.add(_parent_id(index, row))
    sol = getattr(system, 'sol', None)
    for position, residue in enumerate(getattr(sol, 'residues', ()) or ()):
        ids = set(map(int, getattr(residue, 'atoms', ()) or ()))
        classified.update(ids)
        if ids:
            residue_sets.append(ids)
    recognized_keys = set()
    if balls is not None:
        for index, row in balls.iterrows():
            index = int(index)
            if index in boundary_parent or _boundary_flag(row):
                classified.discard(index)
                continue
            if str(row.get('res_name', '')).strip().upper() in SOLVENT_RESIDUES:
                classified.add(index)
                recognized_keys.add((str(row.get('chain_name', '')), str(row.get('res_seq', '')),
                                     str(row.get('res_name', ''))))
    # sys.sol residue grouping is authoritative when available.
    classified -= boundary_parent
    residue_count = sum(bool(ids & classified) for ids in residue_sets) if residue_sets else len(recognized_keys)
    present = set()
    for index, row in iface.net.balls.iterrows():
        if int(index) in boundary or _boundary_flag(row):
            continue
        parent = _parent_id(index, row)
        if parent in classified and parent not in (set(iface.group1_indices) | set(iface.group2_indices)):
            present.add(int(index))
    return {
        'solvent_detected': bool(classified), 'solvent_atom_count': len(classified),
        'solvent_residue_count': residue_count,
        'solvent_definition': 'System.sol residues and shared SOLVENT_RESIDUES; physical generators only',
        'solvent_competing_generator_count': len(present),
        'solvent_context': bool(present),
    }, present


def _network_groups(iface):
    groups = (set(map(int, iface.group1_indices)), set(map(int, iface.group2_indices)))
    mapped = (set(), set())
    for index, row in iface.net.balls.iterrows():
        parent = _parent_id(index, row)
        for source, target in zip(groups, mapped):
            if parent in source:
                target.add(int(index))
    if any(len(source) != len(target) for source, target in zip(groups, mapped)):
        raise ValueError('Interface side membership cannot be mapped uniquely to network generators')
    return mapped


def reverse_perspective(side):
    """Same patches, reversed normal: H changes sign, K does not."""
    result = dict(side)
    result['orientation'] = 'group2 -> group1'
    for key in ('smooth_H', 'smooth_H_partial', 'raw_signed_edge', 'edge_H', 'partial_H', 'total_H',
                'smooth_mean_curvature', 'internal_crease_mean_curvature',
                'boundary_mean_curvature', 'total_mean_curvature',
                'smooth_integrated_mean_curvature', 'internal_seam_integrated_mean_curvature',
                'boundary_integrated_mean_curvature', 'partial_integrated_mean_curvature',
                'total_integrated_mean_curvature'):
        result[key] = None if side.get(key) is None else -side[key]
    result['positive_edge'] = None if side.get('negative_edge') is None else -side['negative_edge']
    result['negative_edge'] = None if side.get('positive_edge') is None else -side['positive_edge']
    return result


def _open_interface_accounting(representation, network, smooth_partial_values,
                               supported_edges, unresolved_surfaces, unresolved_edges,
                               selection_complete):
    """Use the validated open-interface kernel on the solved physical incidence."""
    from vorpy.src.analyze.open_interface_curvature import (
        _aw_supported_gauss_bonnet, _network_face_cycle, analyze_surface_topology,
        gauss_bonnet_accounting, integrated_mean_curvature, polygon_gaussian_terms,
    )

    cycles = []
    face_records = []
    cycle_error = None
    for surface in representation.selected_surfaces:
        try:
            cycle = _network_face_cycle(network, network.surfs.loc[surface.feature_id])
            cycles.append(cycle)
            face_records.append((surface.feature_id, surface.generator_ids, cycle))
        except (ValueError, KeyError, TypeError, AttributeError) as exc:
            cycle_error = f'{type(exc).__name__}: {exc}'
    try:
        topology = analyze_surface_topology(cycles) if len(cycles) == len(representation.selected_surfaces) else None
    except ValueError as exc:
        topology, cycle_error = None, str(exc)

    mean = integrated_mean_curvature(
        [v for v in smooth_partial_values if v is not None],
        [0.5 * edge.signed_integral_A_rad for edge in supported_edges],
        unresolved_faces=unresolved_surfaces,
        unresolved_internal_edges=unresolved_edges,
    )
    mean_certified = bool(mean['certified'] and selection_complete and topology is not None
                          and topology['boundary_is_manifold'] and topology['orientable'])
    if not mean_certified:
        mean['total_H'] = None
        mean['certified'] = False

    gauss = None
    if topology is not None and topology['boundary_is_manifold'] and topology['orientable']:
        if representation.scheme in {'prm', 'pow'}:
            coordinates = {}
            for i, row in network.verts.iterrows():
                xyz = row.get('loc', row.get('position', row.get('xyz')))
                coordinates[int(i)] = tuple(map(float, xyz))
            gauss = polygon_gaussian_terms(cycles, coordinates)
        elif representation.scheme == 'aw':
            supported_ids = {
                int(s.feature_id) for s in representation.selected_surfaces
                if finite(network.surfs.loc[s.feature_id].get('int_gauss_curv')) is not None
                and finite(s.smooth_H_A) is not None
            }
            gauss = _aw_supported_gauss_bonnet(
                network, face_records, supported_ids, set(map(int, representation.group_a)), order=256,
            )
            full_support = len(supported_ids) == len(representation.selected_surfaces)
            gauss['complete'] = bool(
                full_support and not gauss['unresolved_edges'] and not gauss['unresolved_vertices']
                and not gauss['unresolved_faces'] and len(supported_ids) == topology['faces']
            )
            gauss['certified'] = bool(gauss['certified'] and gauss['complete'])
            gauss['accounting_scope'] = 'full_selected_interface' if full_support else 'supported_subcomplex'
            gauss['supported_face_count'] = len(supported_ids)
        else:
            gauss = gauss_bonnet_accounting(
                smooth_K=0.0, face_boundary_geodesic=0.0, vertex_corner_terms=0.0,
                euler_characteristic=topology['euler_characteristic'], complete=False,
            )
    if gauss is not None:
        gauss['certified'] = bool(gauss['certified'] and selection_complete)
    return topology, mean, gauss, cycle_error


def summarize_physical_interface(representation, network, selection_complete=True):
    """Read feature-level availability instead of trusting a combined partial."""
    surfaces = representation.selected_surfaces
    edges = representation.interface_edges
    reasons = []
    smooth_values = [finite(surface.smooth_H_A) for surface in surfaces]
    smooth_supported = sum(value is not None for value in smooth_values)
    partial_values = list(smooth_values)
    # Cell-completeness policy may exclude a measured patch from the certified
    # analyzer. Retain its cached analytic integral as an explicitly uncertified
    # partial, using the same cell-outward Group-1 orientation as that analyzer.
    if representation.scheme == 'aw' and smooth_supported == 0:
        for i, surface in enumerate(surfaces):
            if partial_values[i] is not None:
                continue
            row = network.surfs.loc[surface.feature_id]
            by_ball = row.get('int_mean_curv_by_ball')
            owner = next((g for g in surface.generator_ids if g in representation.group_a), None)
            if isinstance(by_ball, dict) and owner is not None:
                partial_values[i] = finite(by_ball.get(owner, by_ball.get(str(owner))))
    partial_supported = sum(value is not None for value in partial_values)
    smooth_partial = sum(v for v in partial_values if v is not None) if partial_supported or not surfaces else None
    smooth = sum(smooth_values) if smooth_supported == len(surfaces) else None
    if smooth_supported != len(surfaces):
        reasons.append('Selected smooth mean-curvature contributions are missing')
    seam_edges = [edge for edge in edges if edge.classification != 'boundary_selected']
    supported = [edge for edge in seam_edges if edge.classification == 'interior_selected'
                 and edge.curvature_status == 'integrated' and edge.supported
                 and edge.complete and edge.bounded
                 and finite(edge.signed_integral_A_rad) is not None]
    invalid_boundary = [edge for edge in edges if edge.classification == 'boundary_selected'
                        and not (edge.supported and edge.bounded and edge.complete)]
    missing_edges = {edge_id for surface in surfaces for edge_id in surface.boundary_edge_ids} - {edge.feature_id for edge in edges}
    # Free edges have no scalar H term. Their incidence still participates in
    # the independent manifold audit, but missing boundary curvature is not H.
    unresolved_edges = len(seam_edges) - len(supported)
    signed = sum(edge.signed_integral_A_rad for edge in supported) if supported or not seam_edges else None
    unsigned = (sum(edge.unsigned_integral_A_rad for edge in supported)
                if all(finite(edge.unsigned_integral_A_rad) is not None for edge in supported)
                and (supported or not seam_edges) else None)
    edge_half = None if signed is None else .5 * signed
    if unresolved_edges:
        reasons.append('Required interface edge turning contributions are unresolved')
    if missing_edges:
        reasons.append(f'Selected surface incidence references {len(missing_edges)} missing physical edges')
    for edge in seam_edges:
        if edge not in supported:
            reasons.append(f'Edge {edge.feature_id}: {edge.reason or edge.curvature_status}')
    if not selection_complete:
        reasons.append('Alpha selection or selected-pair physical mapping is incomplete')
    candidates = len(surfaces) + len(seam_edges)
    coverage = (smooth_supported + len(supported)) / candidates if candidates else (1. if selection_complete else 0.)
    certified = selection_complete and smooth_supported == len(surfaces) and unresolved_edges == 0
    partial = None if smooth_partial is None or edge_half is None else smooth_partial + edge_half
    gaussian_values = [finite(network.surfs.loc[s.feature_id].get('int_gauss_curv')) for s in surfaces]
    if representation.scheme in {'prm', 'pow'}:
        gaussian_values = [0. for _ in surfaces]  # Exact planar smooth term only.
    gaussian_supported = sum(value is not None for value in gaussian_values)
    gaussian_partial = sum(v for v in gaussian_values if v is not None) if gaussian_supported or not surfaces else None
    try:
        topology, validated_mean, validated_gauss, cycle_error = _open_interface_accounting(
            representation, network, partial_values, supported, len(surfaces) - smooth_supported,
            unresolved_edges, selection_complete,
        )
    except (ValueError, KeyError, TypeError, AttributeError, RuntimeError, IndexError) as exc:
        # A failed optional topology/Gaussian audit must not erase a valid
        # alpha selection, measured area or previously available partials.
        topology, validated_gauss, cycle_error = None, None, f'{type(exc).__name__}: {exc}'
        validated_mean = {'smooth_H': smooth_partial, 'internal_edge_H': edge_half,
                          'boundary_H': 0., 'total_H': None, 'certified': False}
    if cycle_error:
        reasons.append(f'Physical face topology unresolved: {cycle_error}')
    if topology is None:
        reasons.append('Validated physical face-cycle topology is unavailable')
    elif not topology['boundary_is_manifold'] or not topology['orientable']:
        reasons.append('Physical selected-surface incidence is nonmanifold or nonorientable')
    gaussian_certified = bool(validated_gauss and validated_gauss.get('certified'))
    if validated_gauss and not gaussian_certified:
        reasons.append(validated_gauss.get('status', 'Gauss-Bonnet validation did not certify'))
    smooth_k = finite(validated_gauss.get('smooth_K')) if validated_gauss else gaussian_partial
    vertex_gaussian = (finite(validated_gauss.get('interior_vertex_defect_K'))
                       if validated_gauss else None)
    seam_gaussian = (finite(validated_gauss.get('internal_seam_face_geodesic_K', 0.))
                     if validated_gauss else None)
    boundary_geodesic = (finite(validated_gauss.get('free_boundary_face_geodesic_K',
                                                  validated_gauss.get('boundary_geodesic_K')))
                         if validated_gauss else None)
    boundary_corner = (finite(validated_gauss.get('boundary_vertex_corner_K'))
                       if validated_gauss else finite(validated_gauss.get('boundary_corner_K'))
                       if validated_gauss else None)
    gauss_lhs = finite(validated_gauss.get('total_GB_left')) if validated_gauss else None
    gauss_expected = finite(validated_gauss.get('two_pi_chi')) if validated_gauss else None
    gauss_residual = finite(validated_gauss.get('residual')) if validated_gauss else None
    total_gaussian = (smooth_k + seam_gaussian + vertex_gaussian
                      if gaussian_certified and all(x is not None for x in (smooth_k, seam_gaussian, vertex_gaussian)) else None)
    scope = validated_gauss.get('accounting_scope', 'full_selected_interface') if validated_gauss else 'unavailable'
    subset = None
    if scope == 'supported_subcomplex':
        reasons.append('Gaussian partial terms describe a supported subcomplex, not the full selected interface')
        subset = {key: validated_gauss.get(key) for key in (
            'smooth_K', 'internal_seam_face_geodesic_K', 'interior_vertex_defect_K',
            'free_boundary_face_geodesic_K', 'boundary_vertex_corner_K',
            'total_GB_left', 'two_pi_chi', 'residual', 'certified', 'status')}
        subset['certified'] = False  # Never certify the full selection from its subset.
        subtop = validated_gauss['topology']
        subset['topology'] = {key: subtop[key] for key in (
            'vertices', 'edges', 'faces', 'euler_characteristic', 'components', 'boundary_components')}
        subset['quadrature_order'] = 256
        gauss_lhs = gauss_residual = None
        gauss_expected = 2 * math.pi * topology['euler_characteristic']
    return {
        'orientation': 'group1 -> group2', 'interface_atoms': len(representation.interface_atoms),
        'components': topology['components'] if topology else representation.component_count, 'area': representation.area_A2,
        'area_scope': 'retained complete alpha-selected physical surfaces; excluded mappings are reported separately',
        'area_selection_complete': selection_complete,
        'smooth_H': smooth, 'smooth_H_partial': smooth_partial,
        'raw_signed_edge': signed, 'edge_H': edge_half, 'partial_H': partial,
        'total_H': validated_mean['total_H'], 'unsigned_edge': unsigned,
        'positive_edge': finite(getattr(representation.curvature, 'positive_edge_A_rad', None)) if signed is not None else None,
        'negative_edge': finite(getattr(representation.curvature, 'negative_edge_A_rad', None)) if signed is not None else None,
        'edge_candidates': len(seam_edges) + len(invalid_boundary) + len(missing_edges),
        'boundary_edges': sum(e.classification == 'boundary_selected' for e in edges),
        'supported_edges': len(supported), 'unresolved_edges': unresolved_edges,
        'supported_surfaces': smooth_supported, 'unresolved_surfaces': len(surfaces) - smooth_supported,
        'mean_coverage': coverage, 'mean_certified': validated_mean['certified'],
        'coverage_basis': 'supported selected surfaces and required seam features / candidate features',
        'smooth_K_partial': gaussian_partial,
        'smooth_K_coverage': gaussian_supported / len(surfaces) if surfaces else 1.,
        'total_K': total_gaussian, 'gaussian_total_certified': gaussian_certified,
        'euler_incidence': topology['euler_characteristic'] if topology else representation.euler_characteristic,
        'boundary_policy': 'Mean curvature uses smooth patches plus interior seams; no free-boundary turning term',
        'unresolved_reasons': list(dict.fromkeys(reasons)),
        'topology_vertices': topology['vertices'] if topology else None,
        'topology_edges': topology['edges'] if topology else None,
        'topology_faces': topology['faces'] if topology else None,
        'euler_characteristic': topology['euler_characteristic'] if topology else None,
        'boundary_loops': topology['boundary_components'] if topology else None,
        'smooth_mean_curvature': validated_mean['smooth_H'],
        'internal_crease_mean_curvature': validated_mean['internal_edge_H'],
        'boundary_mean_curvature': validated_mean['boundary_H'],
        'total_mean_curvature': validated_mean['total_H'],
        'certified_mean_total': validated_mean['certified'],
        'smooth_gaussian_curvature': smooth_k,
        'intrinsic_vertex_gaussian': vertex_gaussian,
        'boundary_geodesic_gaussian': boundary_geodesic,
        'boundary_corner_gaussian': boundary_corner,
        'gauss_bonnet_lhs': gauss_lhs,
        'gauss_bonnet_expected': gauss_expected,
        'gauss_bonnet_residual': gauss_residual,
        'total_gaussian_curvature': total_gaussian,
        'certified_gaussian_total': gaussian_certified,
        'connected_components': topology['components'] if topology else None,
        'genus_by_component': [c.genus for c in topology['component_records']] if topology else None,
        'manifold_incidence_certified': bool(topology and topology['boundary_is_manifold'] and topology['orientable']),
        'smooth_integrated_mean_curvature': smooth_partial,
        'internal_seam_integrated_mean_curvature': edge_half,
        'boundary_integrated_mean_curvature': 0.,
        'partial_integrated_mean_curvature': partial,
        'total_integrated_mean_curvature': validated_mean['total_H'],
        'mean_curvature_certified': validated_mean['certified'],
        'mean_supported_faces': smooth_supported,
        'mean_unsupported_faces': len(surfaces)-smooth_supported,
        'mean_supported_internal_seams': len(supported),
        'mean_unresolved_internal_seams': unresolved_edges,
        'smooth_integrated_gaussian_curvature': smooth_k,
        'intrinsic_seam_gaussian_curvature': seam_gaussian,
        'interior_vertex_gaussian_curvature': vertex_gaussian,
        'intrinsic_integrated_gaussian_curvature': total_gaussian,
        'boundary_geodesic_curvature': boundary_geodesic,
        'boundary_corner_turning': boundary_corner,
        'gaussian_curvature_certified': gaussian_certified,
        'gaussian_accounting_scope': scope,
        'gaussian_supported_faces': validated_gauss.get('supported_face_count', gaussian_supported) if validated_gauss else gaussian_supported,
        'gaussian_unsupported_faces': len(surfaces)-(validated_gauss.get('supported_face_count', gaussian_supported) if validated_gauss else gaussian_supported),
        'gaussian_coverage': (validated_gauss.get('supported_face_count', gaussian_supported) if validated_gauss else gaussian_supported)/len(surfaces) if surfaces else 1.,
        'supported_subcomplex_gauss_bonnet': subset,
    }


def _comparison_key(iface, solvent_ids, alpha):
    """Exact non-solvent generator geometry, including numerical boundary sites."""
    rows = []
    for index, row in iface.net.balls.iterrows():
        if int(index) not in solvent_ids:
            rows.append((_parent_id(index, row),
                         tuple(map(float, row['loc'])), float(row['rad'])))
    payload = (sorted(rows), sorted(iface.group1_indices), sorted(iface.group2_indices),
               iface.net.settings.get('net_type'), alpha,
               {k: iface.net.settings.get(k) for k in ('surf_res', 'max_vert', 'box_size')})
    return hashlib.sha256(repr(payload).encode()).hexdigest()


def analyze_interface_geometry(iface, alpha=None, *, refresh=False):
    """Analyze once per solve/query. No additional network solve is performed."""
    from vorpy.src.geometry.duals import build_dual
    from vorpy.src.geometry.filtrations.alpha import build_alpha_filtration
    from vorpy.src.geometry.interfaces import build_alpha_interface

    settings = iface.net.settings
    explicit_alpha = alpha is not None
    alpha = float(settings.get('interface_alpha', 0.) if alpha is None else alpha)
    if not math.isfinite(alpha):
        raise ValueError('Interface alpha must be finite')
    scheme = settings.get('net_type', 'aw')
    key = (id(iface.net), id(iface.net.surfs), id(iface.net.edges), id(iface.net.verts), alpha, scheme, 2)
    if not refresh and getattr(iface, '_geometry_analysis_key', None) == key:
        return iface.geometry_analysis
    environment, solvent_ids = solvent_environment(iface)
    result = InterfaceGeometryAnalysis(metadata={
        'schema_version': 2, 'interface_id': iface.interface_id, 'scheme': scheme,
        'certification_definition': 'complete numerical accounting under the validated implementation; not interval-arithmetic proof',
        'alpha_value': alpha, 'alpha_convention': 'additive' if scheme == 'aw' else 'squared_radius' if scheme == 'prm' else 'power_distance',
        'alpha_units': 'A' if scheme == 'aw' else 'A^2',
        'alpha_value_source': 'explicit query' if explicit_alpha else 'interface_alpha setting or existing alpha-zero convention',
        'molecular_contact_patch_status': 'not_implemented',
        'additional_network_solves': 0, **environment,
    })
    try:
        groups = _network_groups(iface)
        source_key = (key[0:4], scheme)
        if refresh or getattr(iface, '_geometry_source_key', None) != source_key:
            iface._geometry_dual = build_dual(iface.net)
            # Phase 1 AW selection requires only validated pair births.
            iface._geometry_filtration = build_alpha_filtration(iface.net, max_dimension=1 if scheme == 'aw' else 3)
            iface._geometry_source_key = source_key
        filtration = iface._geometry_filtration
        representation = build_alpha_interface(iface.net, iface._geometry_dual, filtration, alpha, *groups)
        iface._geometry_representation = representation
        subcomplex = filtration.simplices_at(alpha)
        blocked_pairs = [r for r in subcomplex.blocked if r.dimension == 1
                         and ((r.generator_tuple[0] in groups[0] and r.generator_tuple[1] in groups[1])
                              or (r.generator_tuple[1] in groups[0] and r.generator_tuple[0] in groups[1]))]
        candidates = {tuple(s.generator_ids) for s in representation.candidate_surfaces}
        eligible = {tuple(s.generator_ids) for s in representation.candidate_surfaces
                    if s.supported and s.complete and s.bounded and s.area_A2 is not None
                    }
        selected = {tuple(sorted(r.generator_tuple)) for r in filtration.interface_at(alpha, *groups)}
        failures = sum(not r['selected'] for r in representation.selection_mappings)
        unresolved_pairs = {tuple(r.generator_tuple) for r in blocked_pairs}
        for pair in candidates:
            birth = filtration.simplex_birth(pair)
            if birth is None or not birth.supported or finite(birth.filtration_birth) is None:
                unresolved_pairs.add(pair)
        resolved = len(candidates - unresolved_pairs)
        complete = not unresolved_pairs and not failures
        result.alpha_selection = {
            'eligible_bicolor_pairs': len(eligible), 'candidate_bicolor_pairs': len(candidates),
            'selected_pairs': len(selected), 'selected_physical_pairs': len(representation.selected_surfaces),
            'excluded_pairs': len(eligible - selected), 'unresolved_pairs': len(unresolved_pairs),
            'mapping_failures': failures, 'coverage': resolved / len(candidates) if candidates else (0. if blocked_pairs else 1.),
            'certified': complete, 'status': 'available' if complete else 'partial',
            'selected_generator_pairs': [list(pair) for pair in sorted(selected)],
            'selected_system_pairs': _system_pairs(iface.net, sorted(selected)),
            'selected_physical_system_pairs': _system_pairs(iface.net, (s.generator_ids for s in representation.selected_surfaces)),
        }
        result.metadata['filtration_kind'] = filtration.metadata.get('filtration_kind')
        result.metadata['higher_dimension_policy'] = 'AW higher-dimensional births are not invented; pair-only Phase 1 selection'
        result.voronoi_side_1 = summarize_physical_interface(representation, iface.net, complete)
        result.voronoi_side_1['selected_pairs'] = len(selected)
        result.voronoi_side_2 = reverse_perspective(result.voronoi_side_1)
        result.coverage = {k: result.voronoi_side_1[k] for k in (
            'mean_coverage', 'smooth_K_coverage', 'mean_certified', 'gaussian_total_certified',
            'gaussian_coverage', 'manifold_incidence_certified',
            'mean_supported_faces', 'mean_supported_internal_seams')}
        result.unresolved = list(result.voronoi_side_1['unresolved_reasons'])
        result.unresolved.extend(
            f'Alpha-selected pair ({row["generator_i"]}, {row["generator_j"]}): {row["exclusion_reason"]}'
            for row in representation.selection_mappings if not row['selected']
        )
        if blocked_pairs:
            result.unresolved.extend(f'Alpha pair {r.generator_tuple}: {r.notes or r.birth_source}' for r in blocked_pairs)
        cache = getattr(iface.sys, '_interface_phase1_dry_cache', None)
        if cache is None:
            cache = iface.sys._interface_phase1_dry_cache = {}
        comparison_key = _comparison_key(iface, solvent_ids, alpha)
        dry = cache.get(comparison_key)
        if environment['solvent_context'] and dry is not None:
            dry_pairs = {tuple(p) for p in dry.alpha_selection['selected_physical_system_pairs']}
            wet_pairs = {tuple(p) for p in result.alpha_selection['selected_physical_system_pairs']}
            d, w = dry.voronoi_side_1, result.voronoi_side_1
            result.dry_solvent_comparison = {
                'available': True, 'status': 'matching genuine dry geometry cached',
                'pairs_retained': len(dry_pairs & wet_pairs), 'pairs_removed': len(dry_pairs - wet_pairs),
                'pairs_added': len(wet_pairs - dry_pairs),
                'pair_jaccard': len(dry_pairs & wet_pairs) / len(dry_pairs | wet_pairs) if dry_pairs | wet_pairs else 1.,
                'delta_area': w['area'] - d['area'] if complete and dry.alpha_selection['certified'] else None,
                'delta_H': w['total_H'] - d['total_H'] if w['mean_certified'] and d['mean_certified'] else None,
                'delta_K': None,
            }
        else:
            result.dry_solvent_comparison = {'available': False, 'status': 'no matching genuine dry solve cached' if environment['solvent_context'] else 'not applicable: no solvent competition',
                                            'delta_area': None, 'delta_H': None, 'delta_K': None}
        if not environment['solvent_context']:
            cache[comparison_key] = result
    except (RuntimeError, ValueError, KeyError, AttributeError, IndexError, TypeError) as error:
        # Optional analysis cannot make a solved interface unexportable.
        reason = f'{type(error).__name__}: {error}'
        result.alpha_selection = {k: None for k in (
            'eligible_bicolor_pairs', 'candidate_bicolor_pairs', 'selected_pairs', 'selected_physical_pairs',
            'excluded_pairs', 'unresolved_pairs', 'mapping_failures', 'selected_generator_pairs',
            'selected_system_pairs', 'selected_physical_system_pairs')}
        result.alpha_selection.update({'status': 'unresolved', 'coverage': 0., 'certified': False})
        result.voronoi_side_1 = {k: None for k in PHYSICAL_FIELDS}
        result.voronoi_side_1.update({'orientation': 'group1 -> group2', 'mean_certified': False, 'mean_coverage': 0.,
                                'total_H': None, 'total_K': None, 'gaussian_total_certified': False,
                                'mean_curvature_certified': False, 'gaussian_curvature_certified': False,
                                'manifold_incidence_certified': False, 'gaussian_coverage': 0.,
                                'gaussian_accounting_scope': 'unavailable',
                                'boundary_integrated_mean_curvature': 0.,
                                'unresolved_reasons': [reason]})
        result.voronoi_side_2 = reverse_perspective(result.voronoi_side_1)
        result.coverage = {'mean_coverage': 0., 'smooth_K_coverage': None,
                           'mean_certified': False, 'gaussian_total_certified': False}
        result.unresolved = [reason, 'Validated open-interface accounting could not be completed']
        result.dry_solvent_comparison = {'available': False, 'status': 'analysis unresolved'}
    iface.geometry_analysis = result
    iface._geometry_analysis_key = key
    return result


def parse_geometry_record(headers, row):
    """Typed extension records, including explicit nulls and diagnostic lists."""
    result = {}
    for key, raw in zip(headers, row):
        try:
            result[key] = json.loads(raw) if raw else None
        except (ValueError, TypeError):
            result[key] = raw
    return result


def write_geometry_log(writer, analysis):
    if analysis is None:
        return
    record = analysis.log_record()
    writer.writerow(['Interface Geometry'])
    writer.writerow(list(record))
    writer.writerow([json.dumps(value, allow_nan=False, separators=(',', ':')) for value in record.values()])


def _write_open_geometry_info(stream, side):
    def show(value):
        return 'UNRESOLVED' if value is None else f'{value:.6f}' if isinstance(value, float) else str(value)
    stream.write('\nTopology\n--------\n')
    stream.write('  V / E / F: ' + ' / '.join(show(side.get(k)) for k in
                 ('topology_vertices', 'topology_edges', 'topology_faces')) + '\n')
    for title, key in (('Euler characteristic', 'euler_characteristic'),
                       ('Components', 'connected_components'), ('Boundary loops', 'boundary_loops'),
                       ('Genus by component', 'genus_by_component'),
                       ('Manifold incidence certified', 'manifold_incidence_certified')):
        stream.write(f'  {title}: {show(side.get(key))}\n')
    mean_state = 'TOTAL / CERTIFIED' if side.get('mean_curvature_certified') else 'PARTIAL / UNRESOLVED'
    stream.write(f'\nIntegrated mean curvature [{mean_state}]\n' + '-' * 42 + '\n')
    for title, key in (('Smooth contribution (partial if uncertified)', 'smooth_integrated_mean_curvature'),
                       ('Internal-seam contribution', 'internal_seam_integrated_mean_curvature'),
                       ('Free-boundary contribution', 'boundary_integrated_mean_curvature'),
                       ('Partial sum', 'partial_integrated_mean_curvature'),
                       ('Total', 'total_integrated_mean_curvature')):
        stream.write(f'  {title} (A): {show(side.get(key))}\n')
    stream.write(f'  TOTAL / CERTIFIED: {"yes" if side.get("mean_curvature_certified") else "no"}\n')
    gaussian_state = 'TOTAL / CERTIFIED' if side.get('gaussian_curvature_certified') else 'PARTIAL / UNRESOLVED'
    stream.write(f'\nIntrinsic Gaussian curvature [{gaussian_state}]\n' + '-' * 42 + '\n')
    stream.write(f'  Accounting scope: {side.get("gaussian_accounting_scope", "unavailable")}\n')
    for title, key in (('Smooth contribution', 'smooth_integrated_gaussian_curvature'),
                       ('Intrinsic seam contribution', 'intrinsic_seam_gaussian_curvature'),
                       ('Interior-vertex contribution', 'interior_vertex_gaussian_curvature'),
                       ('Intrinsic total', 'intrinsic_integrated_gaussian_curvature')):
        stream.write(f'  {title}: {show(side.get(key))}\n')
    stream.write('\nGauss-Bonnet\n------------\n')
    stream.write(f'  Boundary-term scope: {side.get("gaussian_accounting_scope", "unavailable")}\n')
    for title, key in (('Boundary geodesic curvature', 'boundary_geodesic_curvature'),
                       ('Boundary corner turning', 'boundary_corner_turning'),
                       ('LHS', 'gauss_bonnet_lhs'), ('2pi*chi', 'gauss_bonnet_expected'),
                       ('Residual', 'gauss_bonnet_residual')):
        stream.write(f'  {title}: {show(side.get(key))}\n')
    stream.write(f'  TOTAL / CERTIFIED: {"yes" if side.get("gaussian_curvature_certified") else "no"}\n')
    subset = side.get('supported_subcomplex_gauss_bonnet')
    if subset:
        stream.write('  Supported-subcomplex diagnostic (not full-interface completion):\n')
        for key, value in subset.items():
            stream.write(f'    {key}: {show(value)}\n')


def write_geometry_info(stream, analysis):
    if analysis is None:
        return
    m, a, v = analysis.metadata, analysis.alpha_selection, analysis.voronoi_side_1
    def show(value):
        return 'UNRESOLVED' if value is None else f'{value:.6f}' if isinstance(value, float) else str(value)
    stream.write('\nALPHA-SELECTED INTERFACE GEOMETRY\n' + '=' * 40 + '\n')
    stream.write('Physical selected geometry: alpha_interface/selected_interface.pml (default Large export or -e surfs).\n')
    stream.write('Full physical interface: surfs.*; selected geometry uses this cached analysis, with zero additional solves.\n')
    for label, key in (('Scheme', 'scheme'), ('Alpha convention', 'alpha_convention'), ('Alpha value', 'alpha_value'), ('Alpha units', 'alpha_units'), ('Solvent context', 'solvent_context')):
        stream.write(f'{label}: {show(m.get(key))}\n')
    stream.write('\nAlpha selection\n---------------\n')
    for label, key in (('Eligible bicolor pairs', 'eligible_bicolor_pairs'), ('Selected pairs', 'selected_pairs'), ('Excluded pairs', 'excluded_pairs'), ('Selected physical pairs', 'selected_physical_pairs'), ('Selection status', 'status')):
        stream.write(f'{label}: {show(a.get(key))}\n')
    stream.write(f'Selection coverage: {100 * a.get("coverage", 0.):.1f}%\n')
    stream.write('\nPhysical Voronoi interface\n--------------------------\n')
    for label, key in (('Interface atoms', 'interface_atoms'), ('Components', 'components'), ('Area (A^2)', 'area')):
        stream.write(f'{label}: {show(v.get(key))}\n')
    for label, side in (('Group 1 -> Group 2', v), ('Group 2 -> Group 1', analysis.voronoi_side_2)):
        coverage = 100 * side.get('mean_coverage', 0.)
        stream.write(f'\n{label}\n')
        for title, key in (('Smooth integrated H', 'smooth_H'), ('Edge integrated H', 'edge_H'), ('Partial integrated H', 'partial_H'), ('Total integrated H', 'total_H')):
            value = show(side.get(key))
            if key == 'total_H' and not side.get('mean_certified'):
                value = f'UNRESOLVED ({coverage:.1f}% coverage)'
            stream.write(f'  {title} (A): {value}\n')
        stream.write(f'  Mean-curvature coverage: {coverage:.1f}%\n  Certified: {"yes" if side.get("mean_certified") else "no"}\n')
    for label, key in (('Raw signed edge integral', 'raw_signed_edge'), ('Unsigned edge curvature', 'unsigned_edge'), ('Positive edge contribution', 'positive_edge'), ('Negative edge contribution', 'negative_edge'), ('Edge candidates', 'edge_candidates'), ('Supported edges', 'supported_edges'), ('Unresolved edges', 'unresolved_edges')):
        stream.write(f'{label}: {show(v.get(key))}\n')
    _write_open_geometry_info(stream, v)
    stream.write(f'Molecular contact patches: {m["molecular_contact_patch_status"]}\n')
    stream.write(f'Dry/solvent comparison: {analysis.dry_solvent_comparison.get("status")}\nAdditional network solves: {m["additional_network_solves"]}\n')
    if analysis.dry_solvent_comparison.get('available'):
        for key, value in analysis.dry_solvent_comparison.items():
            stream.write(f'  {key}: {show(value)}\n')
    for reason in analysis.unresolved:
        stream.write(f'Unresolved: {reason}\n')
    stream.write('\n')
