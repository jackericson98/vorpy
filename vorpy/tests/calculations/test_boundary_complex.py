from collections import Counter
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

from vorpy.src.calculations.boundary_complex import boundary_complex


def triangular_complex(faces):
    edge_ids, edges, surfaces = {}, [], []
    for face in faces:
        ids = []
        for first, second in zip(face, (*face[1:], face[0])):
            key = tuple(sorted((first, second)))
            if key not in edge_ids:
                edge_ids[key] = len(edges)
                edges.append({'verts': key})
            ids.append(edge_ids[key])
        surfaces.append({'edges': ids, 'tris': np.array([[0, 1, 2]])})
    return SimpleNamespace(edges=pd.DataFrame(edges), surfs=pd.DataFrame(surfaces))


def tetrahedron(offset=0):
    return [tuple(v + offset for v in face) for face in
            [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]]


def test_disconnected_disks_in_one_pairwise_surface_do_not_merge_closed_components():
    net = triangular_complex(tetrahedron() + tetrahedron(4))
    merged = {'edges': net.surfs.iloc[0]['edges'] + net.surfs.iloc[4]['edges'],
              'tris': np.array([[0, 1, 2], [3, 4, 5]])}
    net.surfs = pd.DataFrame([merged] + [net.surfs.iloc[i].to_dict() for i in [1, 2, 3, 5, 6, 7]])
    result = boundary_complex(net, range(7))
    assert result['face_euler'][0] == 2
    assert result['components'] == 2
    assert result['orientable']
    assert result['unresolved_faces'] == {}
    assert 8 - 12 + sum(result['face_euler'].values()) == 4


def test_annular_face_is_not_silently_filled_as_two_disks():
    net = SimpleNamespace(
        edges=pd.DataFrame({'verts': [(0, 1), (1, 2), (2, 3), (3, 0),
                                      (4, 5), (5, 6), (6, 7), (7, 4)]}),
        surfs=pd.DataFrame([{'edges': list(range(8)), 'tris': np.array([
            [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
            [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])}]))
    result = boundary_complex(net, [0])
    assert result['unresolved_faces']
    assert not result['orientable']


def test_orientation_constraints_reject_projective_plane():
    phi = (1 + np.sqrt(5)) / 2
    points = np.array([(0, a, b * phi) for a in (-1, 1) for b in (-1, 1)]
                      + [(a, b * phi, 0) for a in (-1, 1) for b in (-1, 1)]
                      + [(a * phi, 0, b) for a in (-1, 1) for b in (-1, 1)])
    representatives, quotient = {}, []
    for point in points:
        canonical = point if point[np.flatnonzero(point)[0]] > 0 else -point
        key = tuple(canonical)
        representatives.setdefault(key, len(representatives))
        quotient.append(representatives[key])
    faces = sorted({tuple(sorted(quotient[v] for v in face)) for face in ConvexHull(points).simplices})
    net = triangular_complex(faces)
    assert len(faces) == 10
    assert set(Counter(e for ids in net.surfs['edges'] for e in ids).values()) == {2}
    assert not boundary_complex(net, range(len(faces)))['orientable']


def test_aggregate_euler_cannot_hide_a_sphere_and_an_annulus():
    net = SimpleNamespace(
        edges=pd.DataFrame({'verts': [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]}),
        surfs=pd.DataFrame([{'edges': list(range(6)), 'tris': np.array(
            tetrahedron() + [(4, 5, 9), (4, 9, 8), (5, 6, 10), (5, 10, 9),
                             (6, 7, 11), (6, 11, 10), (7, 4, 8), (7, 8, 11)])}]))
    # Two mesh components and total chi=2 are insufficient to prove two disks.
    result = boundary_complex(net, [0])
    assert result['unresolved_faces']
    assert not result['orientable']


def test_negative_edge_reference_is_not_a_python_row_selection():
    net = triangular_complex(tetrahedron())
    net.surfs.at[0, 'edges'] = [-1, *net.surfs.at[0, 'edges'][1:]]
    assert boundary_complex(net, range(4))['unresolved_faces']


def test_even_euler_characteristic_does_not_certify_orientability():
    def vertex(i, j):
        if j == 4:
            i, j = -i, 0  # Twisted periodic seam of a Klein bottle.
        return (i % 4) * 4 + j

    faces = []
    for i in range(4):
        for j in range(4):
            a, b, c, d = vertex(i, j), vertex(i + 1, j), vertex(i + 1, j + 1), vertex(i, j + 1)
            faces.extend([(a, b, c), (a, c, d)])
    net = triangular_complex(faces)
    assert 16 - len(net.edges) + len(faces) == 0
    assert set(Counter(e for ids in net.surfs['edges'] for e in ids).values()) == {2}
    assert not boundary_complex(net, range(len(faces)))['orientable']
