"""Mesh equivalence checks for triangulation preparation optimizations."""
import importlib

import numpy as np
import pytest
from shapely import Point, Polygon

triangulation = importlib.import_module('vorpy.src.network.triangulate')


def reference_sort(perimeter, triangles, polygon, points, timing=None):
    inside, outside = [], []
    for triangle in triangles:
        if all(index < len(perimeter) for index in triangle):
            a, b, c = [points[index] for index in triangle]
            centroid = ((a[0] + b[0] + c[0]) / 3., (a[1] + b[1] + c[1]) / 3.)
            destination = inside if polygon.contains(Point(centroid)) else outside
        else:
            destination = inside
        destination.append(triangle)
    designations = {index: ('e' if index < len(perimeter) else 'i') for index in range(len(points))}
    return inside, outside, [], designations


@pytest.mark.parametrize('scalar', [False, True])
@pytest.mark.parametrize('perimeter', [
    [[0., 0.], [2., 0.], [2., 2.], [0., 2.]],
    [[0., 0.], [3., 0.], [3., 1.], [1., 1.], [1., 3.], [0., 3.]],
    [[0., 0.], [2., 0.], [1., 0.12]],
])
def test_complete_mesh_matches_reference(monkeypatch, perimeter, scalar):
    if scalar:
        monkeypatch.setattr(triangulation, 'contains_xy', None)
    perimeter = [np.asarray(point) for point in perimeter]
    timing = {}
    points, triangles = triangulation.triangulate_2D_Surface(perimeter, res=.2, timing=timing)
    monkeypatch.setattr(triangulation, 'sort_tris', reference_sort)
    monkeypatch.setattr(triangulation, '_cached_ring_directions', triangulation._ring_directions)
    reference_points, reference_triangles = triangulation.triangulate_2D_Surface(perimeter, res=.2)
    np.testing.assert_array_equal(points, reference_points)
    np.testing.assert_array_equal(triangles, reference_triangles)
    assert timing['tri_raw_triangles'] == timing['tri_inside_triangles'] + timing['tri_outside_triangles']


def test_triangle_sort_preserves_boundary_rejection_order_and_counters():
    perimeter = [[0., 0.], [3., 0.], [3., 1.], [1., 1.], [1., 3.], [0., 3.]]
    points = perimeter + [[.2, .2], [.5, .5], [.7, .7]]
    triangles = np.asarray([[0, 1, 2], [2, 4, 5], [0, 6, 7], [6, 7, 8]])
    timing = {}
    result = triangulation.sort_tris(perimeter, triangles, Polygon(perimeter), points, timing)
    reference = reference_sort(perimeter, triangles, Polygon(perimeter), points)
    for actual, expected in zip(result[:3], reference[:3]):
        np.testing.assert_array_equal(actual, expected)
    assert result[3] == reference[3]
    assert timing['tri_all_perimeter'] == 2
    assert timing['tri_all_interior'] == timing['tri_mixed'] == 1


def test_ring_cache_is_bounded_and_read_only():
    triangulation._cached_ring_directions.cache_clear()
    directions = triangulation._cached_ring_directions(12)
    with pytest.raises(ValueError):
        directions[0][0] = 3.
    for count in range(1, 270):
        triangulation._cached_ring_directions(count)
    assert triangulation._cached_ring_directions.cache_info().currsize == 256
