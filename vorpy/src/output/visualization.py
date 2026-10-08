"""Cache-only interaction visualization bundle.

Every writer here consumes an already-cached representation.  It never calls
the scientific builders; OBJ/OFF/PDB/PML files are disposable derivatives.
"""

from __future__ import annotations

import json
import re
import csv
import os
from collections.abc import Mapping
from dataclasses import fields
from time import perf_counter
from pathlib import Path

from vorpy.src.geometry.visualization import export_molecular_contact_bundle
from vorpy.src.geometry.visualization.export import (
    DEFAULT_EDGE_RADIUS,
    _physical_edge_mesh,
    _surface_mesh,
    _write_simplex_layers,
    export_dual_visualization,
    pymol_off_loader_lines,
)
from vorpy.src.interface.selected_export import export_selected_interface
from vorpy.src.output.mesh import write_mesh
from vorpy.src.output.pdb import write_pdb
from vorpy.src.io.scientific_adapters import (
    canonical_environment, canonical_network_scope, scope, stable_id,
)
from vorpy.src.output.layout import interface_folder_name
from vorpy.src.results import (
    INTERFACE_LOG_COLUMNS,
    INTERFACE_LOG_QUANTITY_UNITS,
    RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG,
    RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG,
    RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG,
    InterfaceLogRow,
    Provenance,
    ResultIdentity,
    Status,
    adapt_interface_analysis,
    unavailable_interface,
)


def _slug(value):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value)).strip('._') or 'interaction'


def _status(state, reason='', files=()):
    return {'state': state, 'reason': reason, 'files': [str(path) for path in files]}


def _write_molecule(sys, group, name, directory):
    indices = sorted(set(int(index) for index in (getattr(group, 'ball_ndxs', ()) or ())))
    if not indices:
        return _status('NOT_CALCULATED', f'{name} has no cached atom selection')
    directory.mkdir(parents=True, exist_ok=True)
    if not hasattr(sys, 'type'):
        sys.type = 'pdb'
    write_pdb(indices, name, sys, directory=str(directory))
    path = directory / f'{name}.pdb'
    return _status('AVAILABLE', files=(path,))


def _physical_interface(iface, directory):
    analysis = getattr(iface, 'geometry_analysis', None)
    representation = getattr(iface, '_geometry_representation', None)
    if analysis is None or representation is None:
        return _status('NOT_CALCULATED', 'Cached physical interface is absent')
    old_directory = getattr(iface, 'dir', None)
    try:
        directory.mkdir(parents=True, exist_ok=True)
        iface.dir = str(directory)
        export_selected_interface(iface)
    except (AttributeError, KeyError, TypeError, ValueError, RuntimeError) as error:
        return _status('UNRESOLVED', f'Cached physical interface is not exportable: {error}')
    finally:
        iface.dir = old_directory
    root = directory / 'alpha_interface'
    return _status('AVAILABLE', files=tuple(root.glob('*')))


def _dual(iface, directory):
    dual = getattr(iface, '_geometry_dual', None)
    if dual is None:
        return _status('NOT_CALCULATED', 'Cached dual is absent')
    if getattr(dual, 'network', None) is not getattr(iface, 'net', None):
        return _status('STALE', 'Cached dual references another network')
    analysis = getattr(iface, 'geometry_analysis', None)
    filtration = getattr(iface, '_geometry_filtration', None)
    alpha = None
    if analysis is not None:
        alpha = (getattr(analysis, 'metadata', None) or {}).get('alpha_value')
    if filtration is not None and alpha is None:
        return _status('UNRESOLVED', 'Cached filtration has no alpha query')
    representation = getattr(iface, '_geometry_representation', None)
    try:
        metadata = export_dual_visualization(
            iface.net,
            dual,
            directory,
            filtration=filtration if alpha is not None else None,
            alpha=alpha,
            max_dimension=1 if getattr(dual, 'scheme', None) == 'aw' else None,
            interface_alpha=representation,
        )
    except (AttributeError, KeyError, TypeError, ValueError, RuntimeError) as error:
        return _status('UNRESOLVED', f'Cached dual is not exportable: {error}')
    files = tuple(directory.rglob('*'))
    return _status('AVAILABLE', files=files) | {'metadata': metadata}


def _write_launcher(path, layers):
    def pymol_path(target):
        return os.path.relpath(Path(target), Path(path).parent).replace(os.sep, '/')

    lines = [
        '# Cache-only VorPy interaction visualization launcher.',
        'from pymol import cmd',
    ]
    for name in ('molecule_A', 'molecule_B'):
        layer = layers.get(name, {})
        files = layer.get('files', ())
        pdb = next((str(file) for file in files if str(file).lower().endswith('.pdb')), None)
        if pdb:
            lines.append(f"cmd.load({pymol_path(pdb)!r}, {name!r})")
    physical = layers.get('partition_interface', {})
    if physical.get('state') == 'AVAILABLE':
        script = Path(path).parent.parent / 'partition_interface' / 'alpha_interface' / 'selected_interface.py'
        if script.exists():
            lines.append(f'run {pymol_path(script)!r}')
            lines.append("cmd.group('partition_interface', 'alpha_interface alpha_edges alpha_vertices')")
    dual = layers.get('dual_contact_complex', {})
    dual_script = (dual.get('metadata') or {}).get('pymol_script')
    if dual_script:
        script = Path(dual_script).with_suffix('.py')
        if script.exists():
            lines.append(f'run {pymol_path(script)!r}')
            lines.append("cmd.group('dual_contact_complex', 'dual* *_dual* aw_* pow_* prm_*')")
    for side in ('A', 'B'):
        key = f'molecular_contact_surface:{side}'
        layer = layers.get(key, {})
        loaded = []
        for filename, object_name in (
            (f'molecular_contact_surface_{side}.obj', f'molecular_contact_surface_{side}'),
            (f'molecular_contact_surface_{side}_boundary.obj', f'molecular_contact_surface_{side}_boundary'),
            (f'molecular_contact_surface_{side}_junctions.obj', f'molecular_contact_surface_{side}_junctions'),
        ):
            files = [str(file) for file in layer.get('files', ())]
            candidate = next((file for file in files if file.endswith(filename)), None)
            if candidate:
                lines.append(f"cmd.load({pymol_path(candidate)!r}, {object_name!r})")
                loaded.append(object_name)
        if loaded:
            lines.append(f"cmd.group('molecular_contact_surface_{side}', {' '.join(loaded)!r})")
    lines.append("cmd.zoom('molecule_A or molecule_B')")
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def _progress(sys, process, progress):
    updater = getattr(sys, 'update_progress', None)
    if updater is not None and hasattr(sys, 'run_start'):
        updater(process=process, progress=float(progress))


def export_visualization_bundle(sys):
    """Export all representations already cached on each solved interface."""
    root = Path((getattr(sys, 'files', {}) or {}).get('dir') or Path.cwd()).expanduser().resolve()
    interfaces = [iface for iface in (getattr(sys, 'ifaces', None) or []) if getattr(iface, 'net', None) is not None]
    results = []
    for iface in interfaces:
        if getattr(iface, 'dir', None) is None and hasattr(iface, 'set_dir'):
            iface.set_dir()
        interaction = interface_folder_name(getattr(iface, 'group1', None), getattr(iface, 'group2', None))
        interaction_root = Path(getattr(iface, 'dir', None) or root / interaction)
        iface.dir = str(interaction_root)
        visualization_root = interaction_root / 'visualization'
        visualization_root.mkdir(parents=True, exist_ok=True)
        timings = {}
        _progress(sys, f'Export | Preparing visualization: {interaction}', 0.0)
        started = perf_counter()
        layers = {
            'molecule_A': _write_molecule(sys, iface.group1, 'molecule_A', visualization_root),
            'molecule_B': _write_molecule(sys, iface.group2, 'molecule_B', visualization_root),
            'partition_interface': _physical_interface(iface, interaction_root / 'partition_interface'),
            'dual_contact_complex': _dual(iface, interaction_root / 'dual'),
            'solvent': _status('NOT_CALCULATED', 'No cached solvent visualization is available'),
        }
        timings['cached interface and dual preparation'] = perf_counter() - started

        caches = getattr(iface, 'representation_caches', None) or {}
        surfaces = {}
        from vorpy.src.geometry.molecular_contact import MolecularContactSurface
        for key in ('molecular_contact_surface:A', 'molecular_contact_surface:B'):
            value = caches.get(key)
            if isinstance(value, MolecularContactSurface):
                surfaces[key] = value
            else:
                layers[key] = _status('NOT_CALCULATED', 'Molecular contact surface cache is absent')
        _progress(
            sys,
            f'Export | {"Writing" if surfaces else "Recording unavailable"} '
            f'molecular contact surfaces: {interaction}',
            50.0,
        )
        started = perf_counter()
        molecular_manifest = export_molecular_contact_bundle(
            surfaces,
            root,
            interaction=interaction,
            layers=layers,
            directory=interaction_root,
        )
        timings['molecular contact rendering'] = perf_counter() - started
        for key, record in molecular_manifest.get('layers', {}).items():
            if key.startswith('molecular_contact_surface:'):
                directory = interaction_root / record['directory']
                layers[key] = {
                    **_status(
                    molecular_manifest['capabilities'][key]['state'],
                    molecular_manifest['capabilities'][key].get('reason', ''),
                    files=tuple(directory.glob('*')),
                    ),
                    'directory': record['directory'],
                    'metadata': record.get('metadata', {}),
                }
        launcher = visualization_root / f'{interaction}_visualization.pml'
        _progress(sys, f'Export | Writing visualization launcher: {interaction}', 90.0)
        started = perf_counter()
        _write_launcher(launcher, layers)
        timings['launcher and manifest writing'] = perf_counter() - started
        molecular_manifest['combined_pymol_launcher'] = str(launcher)
        molecular_manifest['layers'] = layers
        molecular_manifest.setdefault('capabilities', {})['solvent'] = layers['solvent']
        manifest_path = visualization_root / 'visualization_manifest.json'
        manifest_path.write_text(json.dumps(molecular_manifest, indent=2, sort_keys=True, default=str) + '\n', encoding='utf-8')
        results.append({
            'interaction': interaction,
            'launcher': str(launcher),
            'manifest': str(manifest_path),
            'layers': layers,
            'timings': timings,
        })
    if not interfaces:
        return [{'interaction': None, 'launcher': None, 'manifest': None, 'layers': {}}]
    return results


_TOPOLOGY_LOG_NAMES = {
    'vertices': 'vertex_count',
    'edges': 'edge_count',
    'faces': 'face_count',
    'connected_components': 'connected_component_count',
    'boundary_loops': 'boundary_component_count',
}
_ALPHA_LOG_NAMES = {
    'eligible_pairs': 'alpha_candidate_count',
    'selected_features': 'alpha_selected_count',
    'mapped_surfaces': 'alpha_selected_physical_system_pairs',
}


def _cached_result_identity(iface, *, result_kind='interface', representation='voronoi',
                            side='shared', surface=None):
    net = iface.net
    sys = iface.sys
    analysis = getattr(iface, 'geometry_analysis', None)
    metadata = getattr(analysis, 'metadata', {}) or {}
    filtration_metadata = getattr(getattr(iface, '_geometry_filtration', None), 'metadata', {}) or {}
    alpha_value = metadata.get('alpha_value')
    alpha_units = metadata.get('alpha_units') or filtration_metadata.get('native_alpha_units')
    alpha_method = metadata.get('alpha_method') or metadata.get('alpha_convention')
    if alpha_value is not None:
        alpha_units = alpha_units or 'A'
        alpha_method = alpha_method or 'cached_alpha'
    scheme = str(getattr(getattr(iface, '_geometry_representation', None), 'scheme', None)
                 or getattr(net, 'settings', {}).get('net_type', 'aw')).lower()
    partition = 'power' if scheme in {'power', 'pow', 'prm'} else 'aw'
    archive_provenance = getattr(surface, 'archive_provenance', None) if surface is not None else None
    source = archive_provenance if isinstance(archive_provenance, Mapping) else {}
    surface_selection_source = (
        getattr(surface, 'contact_selection_source', None) if surface is not None else None
    )
    contact_selection_source = (
        surface_selection_source
        if surface_selection_source is not None
        else source.get('contact_selection_source')
    )
    environment = canonical_environment(source.get('environment', metadata.get(
        'environment', 'solvent_competing' if metadata.get('solvent_context') else 'dry')))
    return ResultIdentity(
        result_kind=result_kind,
        system=str(getattr(sys, 'name', 'system')),
        frame=None,
        network_id=stable_id(net),
        network_scope=canonical_network_scope(getattr(net, 'network_scope', None) or scope(net)),
        environment=environment,
        partition=partition,
        representation=representation,
        radius_configuration_id=str(metadata.get('radius_configuration_id', 'cached')),
        interface_id=stable_id(iface),
        group_A=str(getattr(iface.group1, 'name', 'A')),
        group_B=str(getattr(iface.group2, 'name', 'B')),
        side=side,
        orientation=('A -> B' if result_kind == 'molecular_contact_surface'
                     else 'group1 -> group2'),
        alpha_method=alpha_method,
        alpha_value=alpha_value,
        alpha_units=alpha_units,
        system_id=stable_id(sys),
        contact_selection_source=contact_selection_source,
        contact_selection_id=getattr(surface, 'contact_selection_id', None),
        molecular_surface_convention=getattr(surface, 'convention', None),
        definition_version=getattr(surface, 'definition_version', None),
        partner_id='B' if side == 'A' else 'A' if side == 'B' else None,
        source_network_id=stable_id(net),
        source_interface_id=stable_id(iface),
        interaction_id=f'{getattr(sys, "name", "system")}:{getattr(iface, "name", "interface")}',
        orientation_convention=(
            'outward_molecular_normal' if result_kind == 'molecular_contact_surface'
            else 'physical_interface_A_to_B'
        ),
    )


def _cached_interface_result(iface):
    analysis = getattr(iface, 'geometry_analysis', None)
    identity = _cached_result_identity(iface)
    archive_provenance = getattr(analysis, 'archive_provenance', None) if analysis is not None else None
    details = dict(archive_provenance) if isinstance(archive_provenance, Mapping) else {}
    provenance = Provenance('cached interface scientific state', scientific_state_id=stable_id(analysis)
                            if analysis is not None else None, details=details)
    return (adapt_interface_analysis(analysis, identity, provenance=provenance)
            if analysis is not None else unavailable_interface(identity, provenance=provenance))


def _cached_surface_result(iface, surface):
    identity = _cached_result_identity(
        iface, result_kind='molecular_contact_surface',
        representation='molecular_contact_surface', side=surface.side, surface=surface,
    )
    archive_provenance = getattr(surface, 'archive_provenance', None)
    details = dict(archive_provenance) if isinstance(archive_provenance, Mapping) else {}
    surface_selection_source = getattr(surface, 'contact_selection_source', None)
    if surface_selection_source is not None:
        details['contact_selection_source'] = surface_selection_source
    provenance = Provenance('cached molecular contact scientific state',
                            scientific_state_id=stable_id(surface),
                            details=details)
    return surface.to_contract(identity, provenance)


def _iter_result_quantities(result):
    for name, quantity in result.metrics.items():
        yield f'metrics.{name}', quantity
    for name, quantity in result.selection.items():
        yield f'selection.{name}', quantity
    if result.topology is not None:
        for name, quantity in result.topology.quantities.items():
            yield f'topology.{name}', quantity
    if result.mean_curvature is not None:
        for field in fields(result.mean_curvature):
            if field.name != 'convention':
                yield f'mean_curvature.{field.name}', getattr(result.mean_curvature, field.name)
    if result.gaussian_curvature is not None:
        for field in fields(result.gaussian_curvature):
            if field.name != 'convention':
                yield f'gaussian_curvature.{field.name}', getattr(result.gaussian_curvature, field.name)
    for block_name in ('molecular_selection', 'geometry', 'geometric_counts'):
        block = getattr(result, block_name, None)
        for name, quantity in (getattr(block, 'quantities', {}) or {}).items():
            yield f'{block_name}.{name}', quantity


def _compact_result_rows(iface):
    rows = []

    def add(result, representation, side, path, quantity, name=None):
        log_name = name or RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG.get(
            path,
            RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG.get(
                path,
                RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG.get(path, quantity.name),
            ),
        )
        if log_name not in INTERFACE_LOG_QUANTITY_UNITS and log_name not in _TOPOLOGY_LOG_NAMES.values():
            if not (log_name.startswith(('alpha_', 'timing_', 'water_'))
                    or log_name in {
                        'smooth_integrated_mean_curvature',
                        'physical_seam_integrated_mean_curvature',
                        'total_integrated_mean_curvature',
                        'boundary_scalar_mean_curvature',
                        'contact_count', 'selected_ab_contact_count',
                        'direct_contact_atom_count', 'contributing_surface_atom_count',
                        'occluding_same_group_atom_count', 'partner_clipping_atom_count',
                        'unresolved_contact_count', 'carrier_atom_count',
                        'geometric_edge_count', 'spherical_patch_count', 'circular_arc_count',
                        'junction_record_count', 'refined_seam_count', 'boundary_arc_count',
                        'topological_edge_count', 'topological_cut_count',
                        'patch_attachment_record_count', 'connected_component_count',
                        'boundary_component_count', 'euler_characteristic', 'orientable',
                        'genus', 'nonorientable_genus', 'manifold_vertex_count',
                        'singular_vertex_count', 'chain_complex_valid',
                        'boundary_of_boundary_zero', 'beta0', 'beta1', 'beta2',
                    }):
                return
        rows.append(InterfaceLogRow(
            interface_id=result.identity.interface_id or result.identity.result_id,
            representation=representation,
            side=side,
            quantity=log_name,
            value=quantity.value,
            units=quantity.units,
            status=quantity.status,
            provenance=quantity.provenance or result.provenance,
        ).as_csv_row())

    result = _cached_interface_result(iface)
    for path, quantity in _iter_result_quantities(result):
        name = _TOPOLOGY_LOG_NAMES.get(quantity.name)
        if path.startswith('selection.'):
            name = _ALPHA_LOG_NAMES.get(quantity.name)
        add(result, 'physical_partition_interface', 'shared', path, quantity, name=name)

    identity = result.identity
    if identity.alpha_value is not None:
        rows.append(InterfaceLogRow(
            identity.interface_id or identity.result_id, 'alpha_selection', 'shared',
            'alpha_value', identity.alpha_value, identity.alpha_units or '1',
            result.status, provenance=result.provenance,
        ).as_csv_row())

    metadata = getattr(getattr(iface, 'geometry_analysis', None), 'metadata', {}) or {}
    for quantity, value in (metadata.get('timings_seconds', {}) or {}).items():
        if value is not None:
            rows.append(InterfaceLogRow(
                identity.interface_id or identity.result_id, 'timing', 'shared',
                f'timing_{quantity}', value, 's', Status.CERTIFIED,
                provenance=result.provenance,
            ).as_csv_row())
    if metadata.get('timing_total_seconds') is not None:
        rows.append(InterfaceLogRow(
            identity.interface_id or identity.result_id, 'timing', 'shared',
            'timing_total', metadata['timing_total_seconds'], 's', Status.CERTIFIED,
            provenance=result.provenance,
        ).as_csv_row())

    topology = getattr(iface, 'water_topology', None) or {}
    if isinstance(topology, dict):
        for quantity, value in (topology.get('summary', {}) or {}).items():
            if value is not None:
                rows.append(InterfaceLogRow(
                    identity.interface_id or identity.result_id, 'water_analysis', 'shared',
                    f'water_{quantity}', value, '1', Status.PARTIAL,
                    provenance=result.provenance,
                ).as_csv_row())

    caches = getattr(iface, 'representation_caches', None) or {}
    for side in ('A', 'B'):
        surface = caches.get(f'molecular_contact_surface:{side}')
        if surface is None:
            surface_result = unavailable_interface(
                _cached_result_identity(iface), provenance=result.provenance)
        else:
            try:
                surface_result = _cached_surface_result(iface, surface)
            except (AttributeError, KeyError, TypeError, ValueError):
                surface_result = unavailable_interface(
                    _cached_result_identity(iface), provenance=result.provenance)
        for path, quantity in _iter_result_quantities(surface_result):
            add(surface_result, 'molecular_contact_surface', side, path, quantity,
                name=_TOPOLOGY_LOG_NAMES.get(quantity.name))

    return rows


def _compact_log(iface, path):
    """Write validated canonical Results with provenance preserved."""
    rows = sorted(
        _compact_result_rows(iface),
        key=lambda row: tuple(str(row.get(column, '')) for column in (
            'interface_id', 'representation', 'side', 'quantity',
        )),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=INTERFACE_LOG_COLUMNS)
        writer.writeheader()
        writer.writerows(row for row in rows)


def _write_vertex_pdb(path, vertices):
    rows = []
    for serial, vertex in enumerate(vertices, start=1):
        xyz = getattr(vertex, 'xyz', None)
        if xyz is None and hasattr(vertex, 'get'):
            xyz = vertex.get('loc', vertex.get('point'))
        if xyz is None:
            continue
        x, y, z = (float(value) for value in xyz)
        rows.append(
            f'HETATM{serial:5d}  V   VTX A{serial:4d}    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          V\n'
        )
    if rows:
        path.write_text(''.join(rows) + 'END\n', encoding='ascii')


def _compact_physical_interface(iface, directory):
    network = getattr(iface, 'net', None)
    representation = getattr(iface, '_geometry_representation', None)
    surfaces = getattr(network, 'surfs', None)
    if network is None or surfaces is None or len(surfaces) == 0:
        return _status('NOT_CALCULATED', 'Cached physical interface is absent')
    directory.mkdir(parents=True, exist_ok=True)
    scheme = str((getattr(network, 'settings', {}) or {}).get('net_type',
                  getattr(representation, 'scheme', 'aw'))).lower()
    prefix = 'power' if scheme in {'pow', 'power', 'prm'} else 'aw'
    names = {}

    # The compact physical layer is the complete cached network partition.
    # Alpha-selected surfaces remain scientific analysis, not the default mesh.
    surface_mesh = _surface_mesh(network, list(range(len(surfaces))))
    surface_target = directory / f'{prefix}_surfs.off'
    if surface_mesh is not None:
        write_mesh(surface_mesh, surface_target.name, directory=directory)
        names['surfaces'] = surface_target

    edges = getattr(network, 'edges', None)
    edge_mesh = (_physical_edge_mesh(network, list(range(len(edges))), DEFAULT_EDGE_RADIUS)
                 if edges is not None and len(edges) else None)
    edge_target = directory / f'{prefix}_edges.off'
    if edge_mesh is not None:
        write_mesh(edge_mesh, edge_target.name, directory=directory)
        names['physical_edges'] = edge_target

    vertex_target = directory / f'{prefix}_verts.pdb'
    vertices = getattr(network, 'verts', None)
    _write_vertex_pdb(
        vertex_target,
        vertices.to_dict('records') if hasattr(vertices, 'to_dict') else (),
    )
    if vertex_target.exists():
        names['physical_vertices'] = vertex_target
    files = tuple(path for path in names.values() if path.exists() and path.stat().st_size)
    network_area = (float(surfaces['sa'].sum())
                    if hasattr(surfaces, 'columns') and 'sa' in surfaces else None)
    return _status('AVAILABLE', files=files) | {
        'scheme': scheme,
        'area_A2': network_area,
        'alpha_area_A2': getattr(representation, 'area_A2', None),
        'selected_surface_count': len(getattr(representation, 'selected_surfaces', ()) or ()),
        'surface_count': len(surfaces),
        'selection_mode': 'full_network',
    }


def _compact_dual(iface, directory):
    dual = getattr(iface, '_geometry_dual', None)
    if dual is None:
        return _status('NOT_CALCULATED', 'Cached dual is absent')
    directory.mkdir(parents=True, exist_ok=True)
    coordinates = {
        int(row['generator_id']): row['coordinates']
        for row in dual.generators
    }
    prefix = f'{"power" if dual.scheme in {"pow", "prm"} else "aw"}_dual'
    result = _write_simplex_layers(
        directory, prefix, list(dual.simplices[1].values()),
        coordinates, DEFAULT_EDGE_RADIUS,
    )
    path = directory / result['edges_file'] if result.get('edges_file') else None
    if path is None or not path.exists() or not path.stat().st_size:
        return _status('UNRESOLVED', 'Cached dual has no renderable edges')
    return _status('AVAILABLE', files=(path,))


def _write_compact_launcher(path, layers):
    def pymol_path(target):
        return os.path.relpath(Path(target), Path(path).parent).replace(os.sep, '/')

    lines = [
        '# Compact cache-only VorPy interface launcher.',
        'import os',
        'from pymol import cmd',
        '_VORPY_INTERFACE_DIR = os.path.dirname(os.path.abspath(__file__))',
        'def _vorpy_file(name): return os.path.join(_VORPY_INTERFACE_DIR, name)',
        *pymol_off_loader_lines(),
    ]
    for name in ('molecule_A', 'molecule_B'):
        pdb = next((Path(file) for file in layers.get(name, {}).get('files', ())
                    if str(file).lower().endswith('.pdb')), None)
        if pdb and pdb.exists():
            lines += [f'cmd.load(_vorpy_file({pymol_path(pdb)!r}), {name!r})',
                      f"cmd.show('cartoon', {name!r})",
                      f"cmd.show('sticks', {name!r})",
                      f"cmd.set('cartoon_transparency', 35, {name!r})"]
    physical = layers.get('partition_interface', {})
    for suffix, name, color in (
        ('_surfs.off', 'partition_interface', (1.0, 0.15, 0.25)),
        ('_edges.off', 'partition_interface_edges', (1.0, 0.8, 0.1)),
    ):
        candidate = next((Path(item) for item in physical.get('files', ())
                          if str(item).endswith(suffix)), None)
        if candidate and candidate.exists() and candidate.stat().st_size:
            lines.append(f'_load_off_layer(_vorpy_file({pymol_path(candidate)!r}), {name!r}, {color!r})')
    if physical.get('state') == 'AVAILABLE':
        lines.append("cmd.group('partition_interface', 'partition_interface partition_interface_edges')")
    dual = layers.get('dual_contact_complex', {})
    dual_file = next((Path(item) for item in dual.get('files', ())
                      if Path(item).exists() and Path(item).stat().st_size), None)
    if dual_file:
        lines.append(
            f'_load_off_layer(_vorpy_file({pymol_path(dual_file)!r}), "dual_contact_complex", '
            f'{(0.45, 0.55, 0.75)!r})'
        )
    for side in ('A', 'B'):
        key = f'molecular_contact_surface:{side}'
        loaded = []
        for file in layers.get(key, {}).get('files', ()):
            if str(file).lower().endswith(('.obj', '.off')) and Path(file).exists() and Path(file).stat().st_size:
                object_name = f'molecular_contact_surface_{side}'
                lines.append(f'cmd.load(_vorpy_file({pymol_path(file)!r}), {object_name!r})')
                loaded.append(object_name)
                break
        if loaded:
            lines.append(f"cmd.group('molecular_contact_surface_{side}', '{' '.join(loaded)}')")
    lines.append("cmd.zoom('molecule_A or molecule_B')")
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def export_compact_visualization_bundle(sys):
    """Write the small interface bundle from cached state only."""
    timings = {}
    total_started = perf_counter()
    root = Path((getattr(sys, 'files', {}) or {}).get('dir') or Path.cwd()).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    try:
        sys.files['dir'] = str(root)
        sys.exports(pdb=True, set_atoms=True, info=True)
    except (AttributeError, OSError, TypeError, ValueError):
        (root / 'info.txt').write_text(f'VorPy system: {getattr(sys, "name", "system")}\n', encoding='utf-8')
    timings['system atom/info export'] = perf_counter() - started

    results = []
    interfaces = [iface for iface in (getattr(sys, 'ifaces', None) or [])
                  if getattr(iface, 'net', None) is not None]
    for iface in interfaces:
        if getattr(iface, 'dir', None) is None and hasattr(iface, 'set_dir'):
            iface.set_dir()
        interface_dir = Path(getattr(iface, 'dir', None) or root / interface_folder_name(
            getattr(iface, 'group1', None), getattr(iface, 'group2', None)))
        interface_dir.mkdir(parents=True, exist_ok=True)
        iface.dir = str(interface_dir)
        layers = {}
        started = perf_counter()
        for group, name in ((getattr(iface, 'group1', None), 'molecule_A'),
                            (getattr(iface, 'group2', None), 'molecule_B')):
            if group is not None and name not in layers:
                atom_name = 'group_1_atoms' if name == 'molecule_A' else 'group_2_atoms'
                layers[name] = _write_molecule(sys, group, atom_name, interface_dir)
        try:
            from vorpy.src.interface.export import export_info
            export_info(iface, directory=str(interface_dir))
        except (AttributeError, KeyError, TypeError, ValueError):
            network = getattr(iface, 'net', None)
            vertices = getattr(network, 'verts', None)
            edges = getattr(network, 'edges', None)
            surfaces = getattr(network, 'surfs', None)
            representation = getattr(iface, '_geometry_representation', None)
            (interface_dir / 'info.txt').write_text(
                f"{getattr(iface, 'name', interface_dir.name)}\n\n"
                "Interface network topology:\n"
                f"  Network vertices: {0 if vertices is None else len(vertices)}\n"
                f"  Network edges: {0 if edges is None else len(edges)}\n"
                f"  Network surfaces: {0 if surfaces is None else len(surfaces)}\n\n"
                "Surface classification:\n"
                "  Physical export selection: full_network (all cached network surfaces).\n"
                f"  Alpha metadata selection: {getattr(representation, 'selection_mode', 'not recorded')}\n"
                "  Detailed statistics: unavailable for this cached interface.\n",
                encoding='utf-8',
            )
        timings[f'{interface_dir.name}: atom/info export'] = perf_counter() - started
        started = perf_counter()
        layers['partition_interface'] = _compact_physical_interface(iface, interface_dir)
        layers['dual_contact_complex'] = _compact_dual(iface, interface_dir)
        timings[f'{interface_dir.name}: physical interface writing'] = perf_counter() - started
        started = perf_counter()
        _compact_log(iface, interface_dir / 'logs.csv')
        timings[f'{interface_dir.name}: canonical log writing'] = perf_counter() - started

        caches = getattr(iface, 'representation_caches', None) or {}
        surfaces = {}
        from vorpy.src.geometry.molecular_contact import MolecularContactSurface
        for key in ('molecular_contact_surface:A', 'molecular_contact_surface:B'):
            value = caches.get(key)
            if isinstance(value, MolecularContactSurface):
                surfaces[key] = value
            else:
                layers[key] = _status('NOT_CALCULATED', 'Molecular contact surface cache is absent')
        if surfaces:
            started = perf_counter()
            molecular_manifest = export_molecular_contact_bundle(
                surfaces, root, interaction=interface_dir.name, layers=layers,
                directory=interface_dir,
            )
            for key, record in molecular_manifest.get('layers', {}).items():
                if not key.startswith('molecular_contact_surface:'):
                    continue
                directory = interface_dir / record['directory']
                layers[key] = {
                    **_status(
                        molecular_manifest['capabilities'][key]['state'],
                        molecular_manifest['capabilities'][key].get('reason', ''),
                        files=tuple(directory.glob('*')),
                    ),
                    'directory': record['directory'],
                    'metadata': record.get('metadata', {}),
                }
            timings[f'{interface_dir.name}: molecular contact rendering'] = perf_counter() - started

        # The compact derivative keeps only the system launcher.  Remove a
        # stale interface-local launcher from older exports so the output
        # cannot be mistaken for a second system visualization entry point.
        launcher = interface_dir / 'interface.pml'
        if launcher.exists():
            launcher.unlink()
        results.append({
            'interaction': interface_dir.name,
            'launcher': None,
            'manifest': str(interface_dir / 'visualization_manifest.json'),
            'layers': layers,
        })

    archive = root / f'{getattr(sys, "name", "system")}.vpy'
    started = perf_counter()
    if getattr(sys, '_loaded_from_archive', False):
        archive = archive if archive.exists() else None
    else:
        try:
            from vorpy.src.io import save_network
            save_network(sys, archive)
            archive_profile = getattr(sys, 'archive_timing', {}) or {}
            timings.update({f'archive: {key}': value for key, value in archive_profile.items()
                            if isinstance(value, (int, float))})
        except (AttributeError, OSError, TypeError, ValueError, RuntimeError) as error:
            archive = None
            print(f'Archive save unavailable: {error}')
    timings['archive writing'] = perf_counter() - started
    generated_log = root / 'aw_verts.txt'
    if generated_log.exists():
        generated_log.unlink()
    for item in results:
        payload = {
            'archive': None if archive is None else str(archive),
            'layers': item['layers'],
            'launcher': item['launcher'],
            'canonical_state': 'cached scientific objects; rendering files are derivatives',
            'timings': {**timings, 'total': perf_counter() - total_started},
        }
        Path(item['manifest']).write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str) + '\n',
            encoding='utf-8',
        )
        item['archive'] = None if archive is None else str(archive)
        item['timings'] = payload['timings']
    return results or [{'interaction': None, 'launcher': None, 'manifest': None,
                        'layers': {}, 'archive': None,
                        'timings': {**timings, 'total': perf_counter() - total_started}}]
