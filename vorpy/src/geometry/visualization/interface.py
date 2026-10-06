"""Normal interface visualization from cached dual/filtration/primal records.

This module performs no network construction, birth calculation, or curvature
integration. Center-connected dual meshes are the existing abstract dual
realization, not a replacement for the curved physical AW surfaces.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from vorpy.src.geometry.interfaces.alpha import _stable_atom
from vorpy.src.geometry.visualization.export import (
    export_dual_visualization, _write_simplex_layers, _write_simplex_tables,
    _write_csv, _COLORS, pymol_off_loader_lines,
)
from vorpy.src.output.draw import DEFAULT_EDGE_RADIUS
from vorpy.src.output.pdb import make_pdb_line


def _cached_sources(iface):
    net = iface.net
    if net.settings.get('net_type', 'aw') != 'aw':
        raise ValueError('Apollonius interface visualization requires an AW interface.')
    dual = getattr(iface, '_geometry_dual', None)
    filtration = getattr(iface, '_geometry_filtration', None)
    representation = getattr(iface, '_geometry_representation', None)
    analysis = getattr(iface, 'geometry_analysis', None)
    if any(value is None for value in (dual, filtration, representation, analysis)):
        raise ValueError('Interface dual export requires the cached geometry analysis from Interface.build().')
    if dual.network is not net or filtration.network is not net or dual.scheme != 'aw' or filtration.scheme != 'aw':
        raise ValueError('Cached AW dual and filtration must belong to the interface network.')
    alpha = float(analysis.metadata['alpha_value'])
    if representation.scheme != 'aw' or representation.selection_mode != 'alpha' or representation.alpha != alpha:
        raise ValueError('Cached physical representation and geometry analysis use different alpha queries.')
    selected_raw = analysis.alpha_selection.get('selected_generator_pairs')
    if selected_raw is None:
        raise ValueError('Interface alpha selection is unresolved; no selected-pair result is available to visualize.')
    selected = {tuple(sorted(map(int, pair))) for pair in selected_raw}
    actual = {tuple(sorted(row.generator_tuple)) for row in filtration.interface_at(
        alpha, representation.group_a, representation.group_b)}
    if selected != actual:
        raise ValueError('Cached geometry_analysis selected pairs disagree with the AW filtration query.')
    mapped_pairs = {(int(row['generator_i']), int(row['generator_j'])) for row in representation.selection_mappings}
    if mapped_pairs != selected:
        raise ValueError('Cached dual-to-physical mappings disagree with geometry_analysis selected pairs.')
    for mapping in representation.selection_mappings:
        surface_id = mapping.get('surface_id')
        if surface_id is not None:
            if surface_id not in net.surfs.index or tuple(sorted(map(int, net.surfs.loc[surface_id, 'balls']))) != (int(mapping['generator_i']), int(mapping['generator_j'])):
                raise ValueError('Cached dual mapping points to the wrong physical AW surface.')
    return net, dual, filtration, representation, analysis, alpha, selected


def interface_dual_records(iface):
    """Read selection/certification and smooth curvature without recalculation."""
    net, dual, filtration, rep, analysis, alpha, selected = _cached_sources(iface)
    pairs = set(dual.simplices[1]) | set(filtration.records.get(1, {}))
    pairs.update(row.generator_tuple for row in filtration.blocked if row.dimension == 1)
    pairs.update(tuple(sorted(surface.generator_ids)) for surface in rep.candidate_surfaces)
    a, b = rep.group_a, rep.group_b
    pairs = sorted(pair for pair in pairs if (pair[0] in a and pair[1] in b) or (pair[1] in a and pair[0] in b))
    mappings = {(int(row['generator_i']), int(row['generator_j'])): row for row in rep.selection_mappings}
    surfaces = {surface.feature_id: surface for surface in rep.surfaces}
    records = []
    for pair in pairs:
        birth = filtration.simplex_birth(pair)
        resolved = bool(birth and birth.supported and birth.filtration_birth is not None
                        and math.isfinite(birth.filtration_birth))
        diagnostic = birth.diagnostics if birth else {}
        lower_bound = diagnostic.get('birth_lower_bound_A')
        excluded_by_bound = bool(
            diagnostic.get('restriction_certificate') == 'global_pair_minimum_lower_bound'
            and diagnostic.get('restriction_status') == 'rejected'
            and lower_bound is not None and math.isfinite(float(lower_bound))
            and float(lower_bound) > alpha + filtration.tolerance)
        status = 'selected' if pair in selected else 'certified_rejected' if excluded_by_bound or (resolved and birth.filtration_birth > alpha + filtration.tolerance) else 'unresolved'
        mapping = mappings.get(pair, {})
        surface = surfaces.get(mapping.get('surface_id'))
        atom_metadata = [_stable_atom(net, generator, a, b, filtration) for generator in pair]
        simplex = dual.simplex(1, pair)
        row = {
            'dual_feature_id': birth.simplex_id if birth else simplex.simplex_id if simplex else '1:' + ','.join(map(str, pair)),
            'generator_1': pair[0], 'generator_2': pair[1],
            'group_1': atom_metadata[0]['group'], 'group_2': atom_metadata[1]['group'],
            'group_A_identity': str(getattr(iface.group1, 'group_id', getattr(iface.group1, 'name', 'A'))),
            'group_B_identity': str(getattr(iface.group2, 'group_id', getattr(iface.group2, 'name', 'B'))),
            'stable_atom_1': mapping.get('stable_atom_i', atom_metadata[0]['stable_atom_id']),
            'stable_atom_2': mapping.get('stable_atom_j', atom_metadata[1]['stable_atom_id']),
            'restriction_status': status, 'alpha_value': alpha,
            'birth_value': float(birth.filtration_birth) if resolved else None,
            'birth_lower_bound': lower_bound,
            'birth_source': birth.birth_source if birth else 'missing AW birth',
            'birth_supported': bool(birth and birth.supported),
            'dual_supported': bool(simplex and simplex.supported),
            'physical_surface_id': surface.feature_id if surface else None,
            'physical_surface_area': surface.area_A2 if surface else None,
            'smooth_integrated_H': surface.smooth_H_A if surface else None,
            'smooth_integrated_K': None,
            'physical_selected': bool(surface and surface.selected),
            'mapping_status': 'mapped_selected' if surface and surface.selected else 'mapping_unresolved' if status == 'selected' else 'not_selected',
            'physical_supported': bool(surface and surface.supported),
            'physical_complete': bool(surface and surface.complete),
            'generator_cells_complete': bool(surface and surface.generator_cells_complete),
            'certification_status': 'certified_lower_bound_exclusion' if excluded_by_bound else 'resolved_numerical_birth' if resolved else 'unresolved',
            'reason': mapping.get('exclusion_reason', '') if status == 'selected' else birth.notes if birth else 'missing AW birth',
        }
        if surface and surface.feature_id in net.surfs.index:
            value = net.surfs.loc[surface.feature_id].get('int_gauss_curv')
            if value is not None and np.isscalar(value) and math.isfinite(float(value)):
                row['smooth_integrated_K'] = float(value)
        records.append(row)
    summary = {
        'full_dual_counts': {str(d): len(dual.simplices[d]) for d in range(4)},
        'eligible_bicolor_pairs': analysis.alpha_selection.get('eligible_bicolor_pairs'),
        'candidate_bicolor_pairs': len(records),
        'selected_pairs': sum(row['restriction_status'] == 'selected' for row in records),
        'certified_rejected_pairs': sum(row['restriction_status'] == 'certified_rejected' for row in records),
        'unresolved_pairs': sum(row['restriction_status'] == 'unresolved' for row in records),
        'mapped_selected_surfaces': len(rep.selected_surfaces), 'selected_area_A2': rep.area_A2,
        'alpha_value': alpha, 'alpha_units': 'A', 'additional_network_solves': 0,
    }
    return records, summary


def export_interface_dual_visualization(iface):
    net, dual, filtration, rep, analysis, alpha, selected = _cached_sources(iface)
    records, summary = interface_dual_records(iface)
    root = Path(iface.dir).resolve() / 'dual'
    root.mkdir(parents=True, exist_ok=True)
    # Existing exporter owns full/alpha complex meshes and primal mappings.
    metadata = export_dual_visualization(net, dual, root, filtration=filtration,
        alpha=alpha, max_dimension=1, interface_alpha=rep, edge_radius=DEFAULT_EDGE_RADIUS * .5)
    coordinates = {int(row['generator_id']): np.asarray(row['coordinates']) for row in dual.generators}
    layers = {}
    for status, prefix in [('selected', 'apollonius_selected'),
                           ('certified_rejected', 'apollonius_rejected'), ('unresolved', 'apollonius_unresolved')]:
        features = []
        for row in records:
            if row['restriction_status'] != status:
                continue
            pair = (row['generator_1'], row['generator_2'])
            feature = dual.simplex(1, pair) or filtration.simplex_birth(pair)
            if feature is not None and all(generator in coordinates for generator in pair):
                features.append(feature)
        layers[status] = _write_simplex_layers(root, prefix, features, coordinates, DEFAULT_EDGE_RADIUS * 2)
        # Re-exporting a query with an empty category must not leave its old
        # mesh masquerading as a current selected/rejected/unresolved layer.
        if layers[status]['edges_file'] is None:
            (root / (prefix + '_edges.off')).unlink(missing_ok=True)
        _write_simplex_tables(root, prefix, features, coordinates)
        vertex_ids = sorted({generator for feature in features for generator in feature.generator_tuple})
        _write_atoms(net, vertex_ids, root / (prefix + '_vertices.pdb'), filtration, rep.group_a, rep.group_b)
    _write_atoms(net, sorted(rep.group_a), root / 'molecule_A.pdb', filtration, rep.group_a, rep.group_b)
    _write_atoms(net, sorted(rep.group_b), root / 'molecule_B.pdb', filtration, rep.group_a, rep.group_b)
    _write_generator_identity(net, filtration, rep.group_a, rep.group_b, root)
    _write_csv(root / 'dual_interface_mapping.csv', [row for row in records if row['restriction_status'] == 'selected'])
    _write_csv(root / 'dual_pair_status.csv', records)
    contacts = {}
    for row in records:
        if row['restriction_status'] != 'selected':
            continue
        surface_id = row['physical_surface_id'] if row['physical_selected'] else None
        surface = net.surfs.loc[surface_id] if surface_id is not None else None
        contacts[row['dual_feature_id']] = {
            'generators': [coordinates[row['generator_1']].tolist(), coordinates[row['generator_2']].tolist()],
            'surface_id': surface_id,
            'points': [] if surface is None else np.asarray(surface.get('points', ())).tolist(),
            'triangles': [] if surface is None else np.asarray(surface.get('tris', ()), dtype=int).tolist(),
        }
    (root / 'contacts.json').write_text(json.dumps(contacts, allow_nan=False), encoding='utf-8')
    _write_scene(root, metadata, layers)
    summary.update({'pymol_script': str(root / 'apollonius_interface.pml'),
                    'dual_geometry': 'existing straight generator-center realization of the cached AW abstract dual',
                    'physical_geometry': 'cached alpha-selected curved AW surface meshes'})
    (root / 'interface_dual_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    iface.dual_visualization_summary = summary
    return summary


def _write_atoms(net, ids, path, filtration, group_a, group_b):
    """Use existing PDB formatting on the solved generator coordinates."""
    with path.open('w', encoding='utf-8', newline='\n') as stream:
        for serial, generator in enumerate(ids, 1):
            row = net.balls.loc[generator]
            atom = _stable_atom(net, generator, group_a, group_b, filtration)
            override = filtration.metadata.get('stable_atom_metadata', {}).get(generator, {})
            element = override.get('element') if isinstance(override, dict) else None
            stream.write(make_pdb_line(ser_num=serial % 100000, name=(atom['atom_name'] or 'C')[:4],
                res_name=(atom['residue_name'] or 'GEN')[:3],
                chain=atom['chain'][:1] or 'G', cfir=atom['insertion_code'][:1],
                res_seq=atom['residue_number'] or generator % 9999, elem=str(element or row.get('element', row.get('elem', 'C')))[:2],
                x=float(row['loc'][0]), y=float(row['loc'][1]), z=float(row['loc'][2])))
        stream.write('END\n')


def _write_generator_identity(net, filtration, group_a, group_b, root):
    serials = {generator: (name, serial) for name, group in [('molecule_A', group_a), ('molecule_B', group_b)]
               for serial, generator in enumerate(sorted(group), 1)}
    records = []
    for generator in net.balls.index:
        atom = _stable_atom(net, int(generator), group_a, group_b, filtration)
        name, serial = serials.get(int(generator), ('', None))
        serial = serial % 100000 if serial is not None else None
        records.append({'generator_id': int(generator), **atom, 'pymol_object': name, 'pdb_serial': serial})
    _write_csv(root / 'generator_identity.csv', records)


def write_interface_dual_summary(stream, iface):
    if iface.net.settings.get('net_type', 'aw') != 'aw' or getattr(iface, '_geometry_representation', None) is None:
        return
    try:
        _, summary = interface_dual_records(iface)
    except ValueError as error:
        stream.write(f'\nApollonius dual visualization: unavailable ({error})\n')
        return
    counts = summary['full_dual_counts']
    stream.write('\nApollonius full dual\n--------------------\n')
    stream.write('  Vertices / edges / faces / tetrahedra: ' + ' / '.join(str(counts[str(d)]) for d in range(4)) + '\n')
    stream.write('Restricted interface dual\n')
    for label, key in [('Eligible bicolor pairs', 'eligible_bicolor_pairs'), ('Selected pairs', 'selected_pairs'),
                       ('Certified rejected', 'certified_rejected_pairs'), ('Unresolved', 'unresolved_pairs'),
                       ('Mapped selected surfaces', 'mapped_selected_surfaces'), ('Selected area (A^2)', 'selected_area_A2')]:
        stream.write(f'  {label}: {summary[key]}\n')


def _write_scene(root, metadata, layers):
    full_dir = root / 'aw' / 'dual'
    physical = metadata['interface_layers']['interface_alpha']
    lines = ['from pymol import cmd', 'from pymol.cgo import BEGIN, END, TRIANGLES, COLOR, VERTEX, CYLINDER',
             'import os, json', *pymol_off_loader_lines()]
    for group, color in [('A', 'marine'), ('B', 'salmon')]:
        lines += [f"cmd.load({str(root / ('molecule_' + group + '.pdb'))!r}, 'molecule_{group}')",
                  f"cmd.hide('everything', 'molecule_{group}')", f"cmd.show('spheres', 'molecule_{group}')",
                  f"cmd.set('sphere_scale', 0.22, 'molecule_{group}')", f"cmd.color('{color}', 'molecule_{group}')"]
    lines += [f"_load_off_layer({str(full_dir / 'dual_full_edges.off')!r}, 'apollonius_full', {_COLORS['dual_full_edges']!r})",
              f"cmd.load({str(full_dir / 'generators.pdb')!r}, 'dual_vertices')", "cmd.disable('dual_vertices')"]
    for status, name, color in [('selected', 'apollonius_selected', _COLORS['alpha_edges']),
                                ('certified_rejected', 'rejected_dual', _COLORS['dual_full_edges']),
                                ('unresolved', 'unresolved_dual', (0.9, 0.25, 0.8))]:
        if layers[status]['edges_file']:
            lines.append(f"_load_off_layer({str(root / layers[status]['edges_file'])!r}, {name!r}, {color!r})")
            if status != 'selected':
                lines.append(f"cmd.disable({name!r})")
    if physical['surfaces']:
        lines += [f"_load_off_layer({str(Path(physical['directory']) / physical['surfaces'])!r}, 'aw_interface', {_COLORS['interface_alpha_surfaces']!r})",
                  "cmd.set('cgo_transparency', 0.5, 'aw_interface')"]
    for filename, name in [('dual_full_triangles.off', 'apollonius_faces'), ('dual_full_tetrahedra_faces.off', 'apollonius_tetrahedra')]:
        if (full_dir / filename).exists():
            lines += [f"_load_off_layer({str(full_dir / filename)!r}, {name!r}, {_COLORS['dual_full_triangles']!r})", f"cmd.disable({name!r})"]
    lines += ["cmd.group('apollonius_dual', 'apollonius_full apollonius_selected dual_vertices rejected_dual unresolved_dual apollonius_faces apollonius_tetrahedra')",
              "cmd.group('interface_geometry', 'aw_interface')", "cmd.zoom('molecule_A or molecule_B')",
              f"_contacts_path = {str(root / 'contacts.json')!r}",
              'def show_contact(dual_feature_id, second_generator=None):',
              "    if second_generator is not None: dual_feature_id = str(dual_feature_id).strip() + ',' + str(second_generator).strip()",
              "    with open(_contacts_path, encoding='utf-8') as stream: contacts = json.load(stream)",
              "    contact = contacts[str(dual_feature_id)]",
              "    a, b = contact['generators']",
              "    cmd.delete('contact_*')",
              "    cmd.load_cgo([CYLINDER, *a, *b, 0.15, 1., 0.55, 0.05, 1., 0.55, 0.05], 'contact_dual_edge')",
              "    cmd.pseudoatom('contact_generators', pos=a); cmd.pseudoatom('contact_generators', pos=b)",
              "    cmd.show('spheres', 'contact_generators')",
              "    obj=[BEGIN, TRIANGLES, COLOR, 1., 0.15, 0.25]",
              "    for triangle in contact['triangles']:",
              "        for index in triangle: obj.extend([VERTEX, *contact['points'][index]])",
              "    obj.append(END)",
              "    if contact['triangles']: cmd.load_cgo(obj, 'contact_aw_surface')",
              "    cmd.group('representative_contact', 'contact_dual_edge contact_generators contact_aw_surface')",
              "    cmd.disable('apollonius_dual'); cmd.disable('interface_geometry')",
              "    cmd.zoom('contact_generators')",
              "cmd.extend('show_contact', show_contact)",
              "# Use show_contact with a dual_feature_id from dual_interface_mapping.csv."]
    script = root / 'apollonius_interface.py'
    script.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    (root / 'apollonius_interface.pml').write_text(f'run "{script.as_posix()}"\n', encoding='utf-8')
