"""Measure an assembled VorPy selected-cell union; never substitute targets."""
from __future__ import annotations

import argparse
import contextlib
import json
from pathlib import Path

import numpy as np


def components(adjacency):
    unseen = set(adjacency)
    result = []
    while unseen:
        pending = [min(unseen)]
        reached = set()
        while pending:
            node = pending.pop()
            if node in unseen:
                unseen.remove(node)
                reached.add(node)
                pending.extend(adjacency[node] & unseen)
        result.append(sorted(reached))
    return result


def boundary_checks(net, boundary):
    """Independently inspect edge incidence, vertex links and face orientation."""
    edge_faces, vertex_faces, vertex_edges = {}, {}, {}
    face_directions = {}
    invalid_cycles = []
    for face_id in boundary:
        face = net.surfs.iloc[face_id]
        local = {}
        for raw_edge in face['edges']:
            edge = int(raw_edge)
            edge_faces.setdefault(edge, set()).add(face_id)
            ends = list(map(int, net.edges.iloc[edge]['verts']))
            if len(ends) != 2:
                invalid_cycles.append(face_id)
                continue
            a, b = ends
            local.setdefault(a, []).append((b, edge))
            local.setdefault(b, []).append((a, edge))
            for vertex in ends:
                vertex_faces.setdefault(vertex, set()).add(face_id)
                vertex_edges.setdefault(vertex, set()).add(edge)
        if not local or any(len(value) != 2 for value in local.values()):
            invalid_cycles.append(face_id)
            continue
        start = min(local)
        current, previous, directions = start, None, {}
        for _ in range(len(local)):
            choices = [(v, e) for v, e in local[current] if e != previous]
            nxt, edge = choices[0]
            directions[edge] = 1 if current < nxt else -1
            current, previous = nxt, edge
            if current == start:
                break
        if current != start or len(directions) != len(local):
            invalid_cycles.append(face_id)
        face_directions[face_id] = directions
    bad_edges = sorted(e for e, faces in edge_faces.items() if len(faces) != 2)
    bad_vertices = []
    for vertex, faces in vertex_faces.items():
        link = {f: set() for f in faces}
        for edge in vertex_edges[vertex]:
            touching = edge_faces[edge] & faces
            for f in touching:
                link[f].update(touching - {f})
        if any(len(n) != 2 for n in link.values()) or len(components(link)) != 1:
            bad_vertices.append(vertex)
    # Assign a winding to each polygon such that the two uses of every edge
    # oppose each other. Unlike Euler parity, this tests orientability itself.
    adjacency = {f: [] for f in boundary}
    for edge, faces in edge_faces.items():
        if len(faces) != 2:
            continue
        a, b = sorted(faces)
        da, db = face_directions.get(a, {}).get(edge), face_directions.get(b, {}).get(edge)
        if da is not None and db is not None:
            relation = -da * db
            adjacency[a].append((b, relation))
            adjacency[b].append((a, relation))
    signs, conflicts = {}, []
    for root in boundary:
        if root in signs:
            continue
        signs[root] = 1
        pending = [root]
        while pending:
            a = pending.pop()
            for b, relation in adjacency[a]:
                wanted = signs[a] * relation
                if b not in signs:
                    signs[b] = wanted
                    pending.append(b)
                elif signs[b] != wanted:
                    conflicts.append([a, b])
    return {
        'bad_edges': bad_edges, 'bad_vertices': sorted(bad_vertices),
        'invalid_face_cycles': sorted(set(invalid_cycles)),
        'orientation_conflicts': conflicts,
        'orientable_by_face_winding': bool(boundary) and not (
            bad_edges or bad_vertices or invalid_cycles or conflicts
        ),
    }


def _oriented_polygon(group, face_id):
    """Order a planar face toward the excluded generator (solid outward)."""
    net = group.net
    face = net.surfs.iloc[face_id]
    if not face['flat']:
        raise ValueError('Polygon orientation/export requires planar faces.')
    body = set(map(int, group.layer_net_atoms[0]))
    inside, = body.intersection(map(int, face['balls']))
    outside, = set(map(int, face['balls'])) - body
    normal = np.asarray(net.balls.iloc[outside]['loc']) - net.balls.iloc[inside]['loc']
    normal /= np.linalg.norm(normal)
    neighbors = {}
    for e in face['edges']:
        a, b = map(int, net.edges.iloc[int(e)]['verts'])
        neighbors.setdefault(a, []).append(b)
        neighbors.setdefault(b, []).append(a)
    if any(len(v) != 2 for v in neighbors.values()):
        raise ValueError('Face does not have a closed polygon cycle.')
    order = [min(neighbors)]
    previous = None
    while True:
        nxt = next(v for v in neighbors[order[-1]] if v != previous)
        if nxt == order[0]:
            break
        previous = order[-1]
        order.append(nxt)
        if len(order) > len(neighbors):
            raise ValueError('Invalid polygon cycle.')
    if len(order) != len(neighbors):
        raise ValueError('Disconnected polygon cycle.')
    xyz = np.array([net.verts.iloc[v]['loc'] for v in order])
    area_vector = np.cross(xyz - xyz.mean(axis=0),
                           np.roll(xyz, -1, axis=0) - xyz.mean(axis=0)).sum(axis=0)
    if np.dot(area_vector, normal) < 0:
        order.reverse()
        xyz = xyz[::-1]
    return order, xyz, normal


def boundary_component_reports(group):
    """Partition actual faces and sum VorPy's measured curvature terms."""
    net = group.net
    boundary = set(map(int, group.layer_surfs[0]))
    edge_faces = {}
    for f in boundary:
        for e in net.surfs.iloc[f]['edges']:
            edge_faces.setdefault(int(e), set()).add(f)
    adjacency = {f: set() for f in boundary}
    for faces in edge_faces.values():
        for f in faces:
            adjacency[f].update(faces - {f})
    reports = []
    for faces in components(adjacency):
        edges = {int(e) for f in faces for e in net.surfs.iloc[f]['edges']}
        vertices = {int(v) for e in edges for v in net.edges.iloc[e]['verts']}
        checks = boundary_checks(net, faces)
        gf = sum(float(net.surfs.iloc[f]['int_gauss_curv']) for f in faces)
        ge = sum(group.boundary_edge_gaussian_contributions[e] for e in edges)
        gv = sum(group.boundary_vertex_gaussian_contributions[v] for v in vertices)
        report = {
            'faces': faces, 'V': len(vertices), 'E': len(edges), 'F': len(faces),
            'chi': len(vertices) - len(edges) + len(faces),
            'closed': not checks['bad_edges'] and not checks['invalid_face_cycles'],
            'manifold': not (checks['bad_edges'] or checks['bad_vertices']
                             or checks['invalid_face_cycles']),
            'orientable': checks['orientable_by_face_winding'],
            'gaussian_face': gf, 'gaussian_edge': ge, 'gaussian_vertex': gv,
            'integrated_gaussian_curvature': gf + ge + gv,
        }
        if all(net.surfs.iloc[f]['flat'] for f in faces):
            volume, radial_dots, directed_edges = 0., [], {}
            origin = np.mean([net.verts.iloc[v]['loc'] for v in vertices], axis=0)
            report['orientation_reference_point'] = origin.tolist()
            for f in faces:
                order, xyz, normal = _oriented_polygon(group, f)
                radial_dots.append(float(np.dot(normal, xyz.mean(axis=0) - origin)))
                for a, b in zip(order, order[1:] + order[:1]):
                    directed_edges.setdefault(tuple(sorted((a, b))), []).append(1 if a < b else -1)
                for j in range(1, len(xyz) - 1):
                    volume += np.dot(xyz[0], np.cross(xyz[j], xyz[j + 1])) / 6
            report.update({
                'signed_volume': float(volume),
                'solid_outward_winding_consistent': all(
                    len(signs) == 2 and sum(signs) == 0 for signs in directed_edges.values()),
                'normal_dot_position_min': min(radial_dots),
                'normal_dot_position_max': max(radial_dots),
            })
        reports.append(report)
    # Exterior first, identified geometrically by its positive signed volume.
    reports.sort(key=lambda r: -r.get('signed_volume', 0.))
    for index, report in enumerate(reports):
        report['component'] = index
    return reports


def export_planar_boundary(group, reports, directory):
    """Export the assembled faces, with cavity winding preserved, to OBJ/PyMOL."""
    obj, script, offset = [], [
        'from pymol import cmd', 'from pymol.cgo import BEGIN, END, TRIANGLES, NORMAL, VERTEX, COLOR',
    ], 1
    for report in reports:
        name = f"boundary_{report['component']}"
        obj.append(f'o {name}')
        cgo = []
        color = (0.3, 0.6, 0.9) if report['signed_volume'] > 0 else (1., 0.45, 0.2)
        script.append(f'mesh = [BEGIN, TRIANGLES, COLOR, {color[0]}, {color[1]}, {color[2]}]')
        for f in report['faces']:
            _, xyz, normal = _oriented_polygon(group, f)
            obj.extend('v ' + ' '.join(map(str, point)) for point in xyz)
            obj.append('f ' + ' '.join(str(offset + j) for j in range(len(xyz))))
            offset += len(xyz)
            for j in range(1, len(xyz) - 1):
                cgo.append('NORMAL, ' + ', '.join(map(str, normal)))
                cgo.extend('VERTEX, ' + ', '.join(map(str, xyz[k])) for k in (0, j, j + 1))
        script.extend(['mesh.extend([' + ',\n'.join(cgo) + ', END])',
                       f'cmd.load_cgo(mesh, "{name}")'])
        if report['signed_volume'] > 0:
            script.append(f'cmd.set("cgo_transparency", 0.7, "{name}")')
    script.extend(['cmd.bg_color("white")', 'cmd.orient("boundary_*")'])
    (directory / 'boundary.obj').write_text('\n'.join(obj) + '\n')
    (directory / 'boundary.py').write_text('\n'.join(script) + '\n')
    (directory / 'boundary.pml').write_text(f'run "{directory / "boundary.py"}"\n')


def inspect_group(group, loop_cells=(), junction_positions=()):
    net = group.net
    selected = set(map(int, group.ball_ndxs))
    topo_to_system, selected_balls = {}, {}
    for row, ball in net.balls.iterrows():
        topo = int(ball['num'])
        system = int(ball.get('system_num', topo))
        topo_to_system[topo] = system
        if system in selected:
            selected_balls[system] = ball
    adjacency = {s: set() for s in selected}
    interfaces, external = [], set()
    for face_id, face in enumerate(net.surfs.to_dict('records')):
        pair = [topo_to_system[int(b)] for b in face['balls']]
        if len(pair) != 2:
            continue
        count = sum(b in selected for b in pair)
        if count == 1:
            external.add(face_id)
        if count == 2:
            area = float(face['sa'] or 0.0)
            record = {'face': face_id, 'cells': pair, 'area': area}
            interfaces.append(record)
            if np.isfinite(area) and area > 1e-10:
                a, b = pair
                adjacency[a].add(b)
                adjacency[b].add(a)
    boundary = set(map(int, group.layer_surfs[0]))
    checks = boundary_checks(net, boundary)
    selected_components = components(adjacency)
    junctions = []
    for index, position in enumerate(junction_positions):
        a, b = set(loop_cells[index]), set(loop_cells[index + 1])
        shared = a & b
        cross = [f for f in interfaces if (
            (f['cells'][0] in a and f['cells'][1] in b)
            or (f['cells'][1] in a and f['cells'][0] in b)
        )]
        local = shared | {c for f in cross for c in f['cells']}
        local_faces = {f for f in boundary if any(
            topo_to_system[int(c)] in local for c in net.surfs.iloc[f]['balls']
        )}
        local_edges = {int(e) for f in local_faces for e in net.surfs.iloc[f]['edges']}
        local_vertices = {int(v) for e in local_edges for v in net.edges.iloc[e]['verts']}
        local_bad_edges = sorted(local_edges & set(checks['bad_edges']))
        local_bad_vertices = sorted(local_vertices & set(checks['bad_vertices']))
        junctions.append({
            'junction': index, 'position': list(position),
            'handle_a_cells': sorted(a), 'handle_b_cells': sorted(b),
            'shared_junction_cells': sorted(shared),
            'shared_selected_selected_faces': cross,
            'internal_faces_retained_on_boundary': [f['face'] for f in cross if f['face'] in boundary],
            'union_connected': any(a | b <= set(c) for c in selected_components),
            'external_boundary_manifold': bool(local_faces) and not (local_bad_edges or local_bad_vertices),
            'local_bad_edges': local_bad_edges, 'local_bad_vertices': local_bad_vertices,
        })
    chi = int(group.boundary_vertex_count - group.boundary_edge_count + group.boundary_face_count)
    applicable = (group.boundary_is_closed and group.boundary_is_manifold
                  and checks['orientable_by_face_winding'] and group.boundary_component_count == 1)
    return {
        'selected_cells': len(selected),
        'complete_selected_cells': sum(bool(b['complete']) for b in selected_balls.values()),
        'incomplete_selected_cells': sorted(s for s in selected if s not in selected_balls or not selected_balls[s]['complete']),
        'selected_cell_components': len(selected_components),
        'selected_cell_component_members': selected_components,
        'selected_selected_faces': interfaces,
        'boundary_vertices': int(group.boundary_vertex_count),
        'boundary_edges': int(group.boundary_edge_count),
        'boundary_faces': int(group.boundary_face_count),
        'euler_characteristic': chi,
        'boundary_components': int(group.boundary_component_count),
        'complete': bool(group.boundary_is_complete and len(selected_balls) == len(selected)
                         and all(bool(b['complete']) for b in selected_balls.values())),
        'boundary_complete': bool(group.boundary_is_complete),
        'closed': bool(group.boundary_is_closed),
        'manifold': bool(group.boundary_is_manifold),
        'orientable': bool(group.boundary_is_orientable),
        'inferred_genus': 1 - chi / 2 if applicable else None,
        'integrated_gaussian_curvature': float(group.int_gauss_curv),
        'gaussian_face': float(group.int_gauss_curv_face),
        'gaussian_edge': float(group.int_gauss_curv_edge),
        'gaussian_vertex': float(group.int_gauss_curv_vertex),
        'boundary_missing_faces': int(group.boundary_missing_faces),
        'boundary_missing_edges': int(group.boundary_missing_edges),
        'boundary_missing_vertices': int(group.boundary_missing_vertices),
        'boundary_face_symmetric_difference': sorted(boundary ^ external),
        'independent_boundary_checks': checks,
        'junctions': junctions,
    }


def acceptance_errors(report, genus, atol=1e-8):
    """Compare measured results with requested acceptance criteria only."""
    expected = {
        'selected_cell_components': 1, 'boundary_components': 1,
        'complete': True, 'closed': True, 'manifold': True, 'orientable': True,
        'euler_characteristic': 2 - 2 * genus, 'inferred_genus': genus,
    }
    errors = [f'{key}: measured {report[key]!r}, expected {value!r}'
              for key, value in expected.items() if report[key] != value]
    if report['complete_selected_cells'] != report['selected_cells']:
        errors.append('Not all selected cells are complete.')
    if not report['independent_boundary_checks']['orientable_by_face_winding']:
        errors.append('Independent face-winding check failed.')
    if report['boundary_face_symmetric_difference']:
        errors.append('Assembled boundary differs from selected/exterior interfaces.')
    expected_g = 4 * np.pi * (1 - genus)
    if not np.isclose(report['integrated_gaussian_curvature'], expected_g, rtol=0, atol=atol):
        errors.append(f'Gaussian curvature differs from {expected_g:.12g}.')
    for junction in report['junctions']:
        if not junction['union_connected'] or not junction['external_boundary_manifold']:
            errors.append(f'Junction {junction["junction"]} is disconnected or nonmanifold.')
        if junction['internal_faces_retained_on_boundary']:
            errors.append(f'Junction {junction["junction"]} retains internal boundary faces.')
        faces = junction['shared_selected_selected_faces']
        if not faces or not all(np.isfinite(f['area']) and f['area'] > 1e-10 for f in faces):
            errors.append(f'Junction {junction["junction"]} lacks positive-area interfaces.')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shape', choices=['multitorus', 'sphere_with_cavities'], default='multitorus')
    parser.add_argument('--n-cavities', type=int, default=1)
    parser.add_argument('--outer-radius', type=float, default=None,
                        help='Omit to grow the exterior safely with cavity count.')
    parser.add_argument('--cavity-radius', type=float, default=3.)
    parser.add_argument('--resolution', type=int, default=32)
    parser.add_argument('--torus-count', type=int)
    parser.add_argument('--interior-count', type=int, default=16)
    parser.add_argument('--cross-section-resolution', type=int, default=12)
    parser.add_argument('--major-radius', type=float, default=10)
    parser.add_argument('--minor-radius', type=float, default=3)
    parser.add_argument('--junction-offset', type=float, default=None,
                        help='Use zero when auditing a previously exported unperturbed PDB.')
    parser.add_argument('--input', type=Path, help='Inspect an existing PDB instead of generating it.')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expect-genus', type=int,
                        help='Assert acceptance criteria after measurement; never changes geometry/results.')
    args = parser.parse_args()
    if args.shape == 'sphere_with_cavities':
        if args.input or args.expect_genus is not None:
            parser.error('Cavity diagnostic generates its own input and checks cavity targets.')
        run_cavity_diagnostic(args)
        return
    if args.torus_count is None:
        parser.error('--torus-count is required for multitorus.')
    # Network construction may change the working directory during export.
    args.output = args.output.resolve()
    if args.input:
        args.input = args.input.resolve()
    from vorpy.src.geometry.validation_shapes import multitorus
    from vorpy.src.command.set import sett
    from vorpy.src.group.group import Group
    from vorpy.src.system.system import System

    args.output.mkdir(parents=True, exist_ok=True)
    shape = multitorus(args.major_radius, args.minor_radius, args.torus_count,
                      args.interior_count, args.cross_section_resolution,
                      junction_offset=args.junction_offset)
    junction_offset = shape.parameters['junction_offset']
    pdb = args.input if args.input else shape.save_pdb(args.output / 'geometry.pdb')
    with (args.output / 'build.log').open('w') as log, contextlib.redirect_stdout(log):
        system = System(file=str(pdb), make_dir=False, print_actions=False)
        selected = [int(row['num']) for _, row in system.balls.iterrows() if str(row['name']).strip() == 'INT']
        if not selected:
            raise ValueError('No INT selected cells found in PDB.')
        group = Group(system, name='topology_diagnostic', atoms=selected,
                      settings=sett('mv', ['30']), make_net=True)
        group.build()
        group.get_info()
        # Recover membership from loaded coordinates, not target topology.
        xyz = np.array([system.balls.iloc[i]['loc'] for i in selected])
        loop_cells = []
        for loop in range(args.torus_count):
            center = [(2 * loop - args.torus_count + 1) * args.major_radius, 0, 0]
            relative = xyz - center
            plane, perpendicular = (1, 2) if loop % 2 == 0 else (2, 1)
            distances = np.hypot(np.hypot(relative[:, 0], relative[:, plane]) - args.major_radius,
                                 relative[:, perpendicular])
            tolerance = max(0.002, 1.5 * np.hypot(junction_offset, 1.37 * junction_offset))
            loop_cells.append([selected[i] for i in np.flatnonzero(distances < tolerance)])
        positions = [[(2 * j - args.torus_count + 2) * args.major_radius, 0, 0]
                     for j in range(args.torus_count - 1)]
        report = inspect_group(group, loop_cells, positions)
    report['input'] = str(pdb.resolve())
    report['settings'] = {'net_type': 'aw', 'max_vert': 30, 'surf_res': 0.2}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2))
    for key, value in report.items():
        if key not in {'selected_selected_faces', 'selected_cell_component_members', 'junctions'}:
            print(f'{key}: {value}')
    for junction in report['junctions']:
        print(json.dumps(junction, indent=2))
    if args.expect_genus is not None:
        errors = acceptance_errors(report, args.expect_genus)
        print('Acceptance: ' + ('FAIL\n' + '\n'.join(errors) if errors else 'PASS'))
        if errors:
            raise SystemExit(1)


def cavity_acceptance_errors(report, n_cavities=1, atol=1e-8):
    """Analytic targets are used only here, after independent measurement."""
    expected = {
        'selected_cell_components': 1, 'boundary_components': n_cavities + 1,
        'euler_characteristic': 2 * (n_cavities + 1),
        'complete': True, 'closed': True, 'manifold': True, 'orientable': True,
    }
    errors = [f'{k}: measured {report[k]!r}, expected {v!r}'
              for k, v in expected.items() if report[k] != v]
    if report['boundary_face_symmetric_difference']:
        errors.append('Boundary differs from selected/excluded interfaces.')
    if not report['independent_boundary_checks']['orientable_by_face_winding']:
        errors.append('Independent winding check failed.')
    measured = report['integrated_gaussian_curvature']
    if not np.isclose(measured, 4 * np.pi * (n_cavities + 1), rtol=0, atol=atol):
        errors.append('Total Gaussian curvature does not match cavity target.')
    parts = report['component_reports']
    if len(parts) != n_cavities + 1:
        errors.append('Independent component count failed.')
    if not np.isclose(sum(p['integrated_gaussian_curvature'] for p in parts), measured,
                      rtol=0, atol=atol):
        errors.append('Component curvature sum differs from group measurement.')
    for p in parts:
        if p['chi'] != 2 or not all(p[k] for k in (
                'closed', 'manifold', 'orientable', 'solid_outward_winding_consistent')):
            errors.append(f'Component {p["component"]} topology/orientation failed.')
        if not np.isclose(p['integrated_gaussian_curvature'], 4 * np.pi, rtol=0, atol=atol):
            errors.append(f'Component {p["component"]} curvature failed.')
        # Radial checks use each component's own vertex centroid.
        if p['component'] == 0:
            oriented = p['signed_volume'] > 0 and p['normal_dot_position_min'] > 0
        else:
            oriented = p['signed_volume'] < 0 and p['normal_dot_position_max'] < 0
        if not oriented:
            errors.append(f'Component {p["component"]} solid-outward normal failed.')
    return errors


def run_cavity_diagnostic(args):
    from vorpy.src.geometry.validation_shapes import sphere_with_cavities
    from vorpy.src.command.set import sett
    from vorpy.src.group.group import Group
    from vorpy.src.system.system import System

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    shape = sphere_with_cavities(args.n_cavities, args.outer_radius,
                                 args.cavity_radius, args.resolution)
    max_vert = max(30., 2 * shape.parameters['outer_radius'])
    if max_vert > 5000:
        raise ValueError('Required vertex radius exceeds VorPy settings limit; reduce cavity_radius.')
    pdb = shape.save_pdb(output / 'geometry.pdb')
    shape.save_xyzr(output / 'geometry.xyzr')
    with (output / 'build.log').open('w') as log, contextlib.redirect_stdout(log):
        system = System(file=str(pdb), make_dir=False, print_actions=False)
        selected = [int(row['num']) for _, row in system.balls.iterrows()
                    if str(row['name']).strip() == 'INT']
        group = Group(system, name='cavity_diagnostic', atoms=selected,
                      settings=sett('mv', [str(max_vert)]), make_net=True)
        group.build()
        group.get_info()
        report = inspect_group(group)
        report['component_reports'] = boundary_component_reports(group)
        export_planar_boundary(group, report['component_reports'], output)
    target = shape.expected_int_gaussian_curvature
    report.update({
        'input': str(pdb), 'n_cavities': args.n_cavities,
        'outer_radius': shape.parameters['outer_radius'],
        'cavity_centers': shape.parameters['cavity_centers'].tolist(),
        'settings': {'net_type': 'aw', 'max_vert': max_vert, 'surf_res': 0.2},
        'expected_gaussian_curvature': target,
        'gaussian_error': report['integrated_gaussian_curvature'] - target,
        'gaussian_relative_error': abs(report['integrated_gaussian_curvature'] - target) / target,
    })
    errors = cavity_acceptance_errors(report, args.n_cavities)
    report['acceptance_errors'] = errors
    (output / 'report.json').write_text(json.dumps(report, indent=2))
    for key in ('selected_cells', 'complete_selected_cells', 'selected_cell_components',
                'boundary_components', 'complete', 'closed', 'manifold', 'orientable',
                'euler_characteristic', 'integrated_gaussian_curvature', 'gaussian_error'):
        print(f'{key}: {report[key]}')
    for part in report['component_reports']:
        print(json.dumps({k: v for k, v in part.items() if k != 'faces'}, indent=2))
    print('Acceptance: ' + ('FAIL\n' + '\n'.join(errors) if errors else 'PASS'))
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
