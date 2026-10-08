"""Cache-only visualization for canonical molecular-contact geometry.

Spherical patches, arcs, junctions, and refined seams are scientific input.
OBJ/CSV/JSON files written here are disposable rendering derivatives.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Mapping


_EPS = 1e-10


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _scale(a, value):
    return tuple(value * x for x in a)


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a):
    length = _norm(a)
    return _scale(a, 1.0 / length) if length > _EPS else None


def _distance(a, b):
    return _norm(_sub(a, b))


def _basis(normal):
    axis = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.9 else (0.0, 1.0, 0.0)
    first = _unit(_cross(normal, axis)) or (0.0, 0.0, 1.0)
    return first, _unit(_cross(normal, first)) or (0.0, 1.0, 0.0)


def _arc_points(patch, arc, samples):
    center = tuple(patch.carrier_center)
    normal = _unit(tuple(arc.circle_normal)) or (0.0, 0.0, 1.0)
    radius = float(patch.carrier_radius)
    offset = float(arc.circle_offset)
    circle_center = _add(center, _scale(normal, radius * offset))
    circle_radius = radius * math.sqrt(max(0.0, 1.0 - offset * offset))
    if circle_radius <= _EPS:
        return [tuple(arc.start_xyz), tuple(arc.end_xyz)]

    start = _unit(_sub(tuple(arc.start_xyz), circle_center))
    end = _unit(_sub(tuple(arc.end_xyz), circle_center))
    if start is None:
        start, _ = _basis(normal)
    tangent = _unit(_cross(normal, start)) or _basis(normal)[1]
    if arc.start_vertex_id is None and arc.end_vertex_id is None:
        delta = math.tau
    elif end is None:
        delta = math.copysign(min(abs(float(arc.angular_span)), math.tau), float(arc.angular_span) or 1.0)
    else:
        observed = math.atan2(_dot(end, tangent), _dot(end, start))
        magnitude = min(abs(float(arc.angular_span)), math.tau)
        if magnitude <= _EPS:
            delta = observed
        else:
            candidates = (magnitude, -magnitude)
            delta = min(candidates, key=lambda value: abs(math.atan2(math.sin(value - observed), math.cos(value - observed))))

    points = []
    for index in range(max(1, int(samples)) + 1):
        fraction = index / max(1, int(samples))
        angle = delta * fraction
        direction = _add(_scale(start, math.cos(angle)), _scale(tangent, math.sin(angle)))
        points.append(_add(circle_center, _scale(direction, circle_radius)))
    if arc.start_vertex_id is not None:
        points[0] = tuple(arc.start_xyz)
    if arc.end_vertex_id is not None:
        points[-1] = tuple(arc.end_xyz)
    return points


def _seam_points(seam, samples):
    normal = _unit(tuple(seam.circle_normal)) or (0.0, 0.0, 1.0)
    first, second = _basis(normal)
    begin, end = seam.parameter_interval
    count = max(1, int(samples))
    points = []
    for index in range(count + 1):
        angle = begin + (end - begin) * index / count
        direction = _add(_scale(first, math.cos(angle)), _scale(second, math.sin(angle)))
        points.append(_add(tuple(seam.circle_center), _scale(direction, seam.circle_radius)))
    points[0] = tuple(seam.start_xyz)
    points[-1] = tuple(seam.end_xyz)
    return points


def _dedupe(points):
    result = []
    for point in points:
        if not result or _distance(result[-1], point) > _EPS:
            result.append(tuple(point))
    if len(result) > 1 and _distance(result[0], result[-1]) <= _EPS:
        result.pop()
    return result


def _write_obj(path, vertices, groups, lines=()):
    text = ['# Derived visualization; canonical source is the .vpy cache\n']
    text.extend('v %.12g %.12g %.12g\n' % tuple(point) for point in vertices)
    for group in groups:
        text.append(f"g {group['name']}\n")
        text.extend('f %d %d %d\n' % face for face in group.get('faces', ()))
    for line in lines:
        text.append('l ' + ' '.join(str(index) for index in line) + '\n')
    path.write_text(''.join(text), encoding='utf-8')


def _write_mapping(path, rows):
    fields = ('kind', 'scientific_id', 'source_patch_id', 'source_arc_id',
              'source_junction_id', 'source_seam_id', 'vertex_start',
              'vertex_count', 'face_count', 'status')
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({field: row.get(field, '') for field in fields} for row in rows)


def _surface_manifest(surface, representation_key, files, counts):
    patches = tuple(getattr(surface, 'patches', ()))
    source_network_ids = sorted({patch.source_network_id for patch in patches if patch.source_network_id})
    source_interface_ids = sorted({patch.source_interface_id for patch in patches if patch.source_interface_id})
    return {
        'archive_id': getattr(surface, 'archive_id', None),
        'representation_key': representation_key,
        'side': surface.side,
        'source_network_ids': source_network_ids,
        'source_interface_ids': source_interface_ids,
        'contact_selection_id': surface.contact_selection_id,
        'molecular_surface_convention': surface.convention,
        'definition_version': surface.definition_version,
        'geometry_status': getattr(surface, 'area_status', surface.status),
        'topology_status': getattr(surface, 'incidence_status', 'UNRESOLVED'),
        'scientific_status': surface.status,
        'topology': {
            'vertex_count': getattr(surface, 'global_vertex_count', None),
            'topological_edge_count': getattr(surface, 'global_edge_count', None),
            'geometric_edge_count': getattr(surface, 'geometric_edge_count', None),
            'face_count': getattr(surface, 'global_face_count', None),
            'topological_cut_count': len(getattr(surface, 'topological_cuts', ())),
            'connected_component_count': getattr(surface, 'components', None),
            'boundary_component_count': getattr(surface, 'boundary_loops', None),
            'euler_characteristic': surface.euler_characteristic,
            'betti_numbers': getattr(surface, 'betti_numbers', None),
            'orientability': getattr(surface, 'orientability', None),
            'genus_by_component': getattr(surface, 'genus_by_component', None),
            'manifold_vertex_count': getattr(surface, 'manifold_vertices', None),
            'singular_vertex_count': getattr(surface, 'singular_vertices', None),
            'chain_complex_status': getattr(surface, 'incidence_status', 'UNRESOLVED'),
            'diagnostics': list(surface.diagnostics),
        },
        'canonical': ['SphericalPatch', 'CircularArc', 'JunctionVertex'] + (
            ['RefinedSeamPiece'] if getattr(surface, 'refined_seams', ()) else []),
        'derived': ['rendering_triangles', 'OBJ', 'boundary_lines', 'junction_points'],
        'generated_rendering_files': files,
        'counts': counts,
    }


def export_molecular_contact_surface(surface, output_dir, *, representation_key=None,
                                     samples_per_arc=8):
    """Render one loaded canonical surface without scientific reconstruction."""
    if surface.side not in {'A', 'B'}:
        raise ValueError('Molecular contact surface side must be A or B')
    samples_per_arc = max(1, int(samples_per_arc))
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    key = representation_key or f'molecular_contact_surface:{surface.side}'
    arcs = {arc.arc_id: arc for arc in surface.arcs}
    seams = tuple(getattr(surface, 'refined_seams', ()))
    mesh_vertices = []
    mesh_groups = []
    boundary_vertices = []
    boundary_lines = []
    junction_vertices = [tuple(junction.xyz) for junction in surface.junctions]
    mapping = []
    rendered_patches = 0

    def add_mesh_group(name, points, faces, row):
        start = len(mesh_vertices) + 1
        mesh_vertices.extend(points)
        mesh_groups.append({'name': name, 'faces': [tuple(start + index for index in face) for face in faces]})
        row.update(vertex_start=start, vertex_count=len(points), face_count=len(faces))

    for patch in surface.patches:
        row = {'kind': 'patch', 'scientific_id': patch.patch_id,
               'source_patch_id': patch.patch_id, 'status': patch.status}
        ring = []
        for arc_id in patch.boundary_arc_ids:
            arc = arcs.get(arc_id)
            if arc is None:
                continue
            points = _arc_points(patch, arc, samples_per_arc)
            if ring and points and _distance(ring[-1], points[0]) <= _EPS:
                points = points[1:]
            ring.extend(points)
            boundary_start = len(boundary_vertices) + 1
            boundary_vertices.extend(points)
            boundary_lines.append(tuple(range(boundary_start, boundary_start + len(points))))
            mapping.append({'kind': 'boundary_arc', 'scientific_id': arc.arc_id,
                            'source_patch_id': patch.patch_id, 'source_arc_id': arc.arc_id,
                            'vertex_start': boundary_start, 'vertex_count': len(points),
                            'status': arc.status})
        ring = _dedupe(ring)
        if len(ring) >= 3:
            direction = _unit(tuple(sum(point[index] - patch.carrier_center[index] for point in ring) for index in range(3)))
            interior = _add(tuple(patch.carrier_center), _scale(direction or (0.0, 0.0, 1.0), patch.carrier_radius))
            points = [interior] + ring
            faces = [(0, index + 1, (index + 1) % len(ring) + 1) for index in range(len(ring))]
            add_mesh_group(f'patch_{patch.patch_id}', points, faces, row)
            rendered_patches += 1
        elif patch.region_kind == 'full_sphere':
            center, radius = tuple(patch.carrier_center), patch.carrier_radius
            points = [center, _add(center, (radius, 0.0, 0.0)), _add(center, (-radius, 0.0, 0.0)),
                      _add(center, (0.0, radius, 0.0)), _add(center, (0.0, -radius, 0.0)),
                      _add(center, (0.0, 0.0, radius)), _add(center, (0.0, 0.0, -radius))]
            faces = [(0, 1, 5), (0, 5, 2), (0, 2, 6), (0, 6, 1),
                     (0, 3, 5), (0, 4, 3), (0, 6, 4), (0, 5, 3)]
            add_mesh_group(f'patch_{patch.patch_id}', points, faces, row)
            rendered_patches += 1
        mapping.append(row)

    for junction in surface.junctions:
        mapping.append({'kind': 'junction', 'scientific_id': junction.vertex_id,
                        'source_junction_id': junction.vertex_id, 'status': junction.status})
    for seam in seams:
        points = _seam_points(seam, samples_per_arc)
        start = len(boundary_vertices) + 1
        boundary_vertices.extend(points)
        boundary_lines.append(tuple(range(start, start + len(points))))
        mapping.append({'kind': 'refined_seam', 'scientific_id': seam.seam_id,
                        'source_seam_id': seam.seam_id, 'vertex_start': start,
                        'vertex_count': len(points), 'status': seam.status})

    prefix = f'molecular_contact_surface_{surface.side}'
    mesh_path = output_dir / f'{prefix}.obj'
    boundary_path = output_dir / f'{prefix}_boundary.obj'
    junction_path = output_dir / f'{prefix}_junctions.obj'
    mapping_path = output_dir / f'{prefix}_mapping.csv'
    metadata_path = output_dir / f'{prefix}_metadata.json'
    _write_obj(mesh_path, mesh_vertices, mesh_groups)
    _write_obj(boundary_path, boundary_vertices, (), boundary_lines)
    _write_obj(junction_path, junction_vertices, (), ((index,) for index in range(1, len(junction_vertices) + 1)))
    _write_mapping(mapping_path, mapping)
    files = [path.name for path in (mesh_path, boundary_path, junction_path, mapping_path, metadata_path)]
    metadata = _surface_manifest(surface, key, files, {
        'patches': len(surface.patches), 'rendered_patches': rendered_patches,
        'arcs': len(surface.arcs), 'refined_seams': len(seams),
        'topological_cuts': len(getattr(surface, 'topological_cuts', ())),
        'junctions': len(surface.junctions), 'mesh_vertices': len(mesh_vertices),
        'mesh_triangles': sum(len(group['faces']) for group in mesh_groups),
    })
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    return metadata


def _layer_metadata(value, root):
    if value is None:
        return {'state': 'NOT_CALCULATED', 'reason': 'Representation is absent', 'files': []}
    if isinstance(value, Mapping):
        result = dict(value)
        result.setdefault('state', 'AVAILABLE')
        result.setdefault('files', [])
        return result
    path = Path(value)
    try:
        path = path.resolve().relative_to(root.resolve())
    except ValueError:
        path = Path(str(value))
    return {'state': 'AVAILABLE', 'files': [str(path)]}


def export_molecular_contact_bundle(surfaces, output_root, *, interaction='interaction',
                                    layers=None, samples_per_arc=8, directory=None):
    """Write independently toggleable A/B layers and a capability manifest."""
    from vorpy.src.geometry.molecular_contact import MolecularContactSurface

    root = (Path(directory) if directory is not None else
            Path(output_root).expanduser().resolve() / 'interfaces' / str(interaction))
    root = root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    surfaces = surfaces or {}
    by_side = {}
    for key, surface in surfaces.items():
        if not isinstance(surface, MolecularContactSurface):
            raise TypeError(f'{key} is not a MolecularContactSurface cache')
        if key not in {surface.side, f'molecular_contact_surface:{surface.side}'}:
            raise ValueError(f'Representation key does not match side: {key}')
        by_side[surface.side] = surface

    manifest = {
        'interaction': str(interaction),
        'canonical': ['SphericalPatch', 'CircularArc', 'JunctionVertex', 'RefinedSeamPiece'],
        'derived': ['rendering_triangles', 'OBJ', 'boundary_lines', 'junction_points'],
        'layers': {}, 'capabilities': {},
    }
    for side in ('A', 'B'):
        key = f'molecular_contact_surface:{side}'
        surface = by_side.get(side)
        if surface is None:
            manifest['capabilities'][key] = _layer_metadata(None, root)
            continue
        directory = root / f'molecular_contact_surface_{side}'
        metadata = export_molecular_contact_surface(
            surface, directory, representation_key=key, samples_per_arc=samples_per_arc)
        manifest['layers'][key] = {'directory': str(directory.relative_to(root)), 'metadata': metadata}
        state = ('AVAILABLE' if metadata['scientific_status'] == 'CERTIFIED'
                 and metadata['topology_status'] == 'CERTIFIED'
                 else metadata['scientific_status'])
        manifest['capabilities'][key] = {'state': state, 'reason': metadata['topology']['diagnostics']}

    supplied_layers = layers or {}
    for key in ('molecule_A', 'physical_partition_interface', 'molecule_B', 'dual_contact_complex'):
        manifest['layers'][key] = _layer_metadata(supplied_layers.get(key), root)
        manifest['capabilities'][key] = manifest['layers'][key]
    visualization_dir = root / 'visualization'
    visualization_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = visualization_dir / 'molecular_contact_visualization_manifest.json'
    manifest['manifest'] = str(manifest_path.relative_to(root))
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    return manifest
