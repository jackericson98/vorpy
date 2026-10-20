import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.calculations import calc_surf_sa, calc_tetra_vol
from vorpy.src.network.build_surfs import build_surfs, _surface_tetra_volumes


build_surfs_module = importlib.import_module('vorpy.src.network.build_surfs')
STRICT_TOLERANCE = 1.0e-14


def _expected_volumes(references, points, triangles):
    return np.array([
        sum(
            calc_tetra_vol(reference, points[triangle[0]], points[triangle[1]], points[triangle[2]])
            for triangle in triangles
        )
        for reference in references
    ])


def _single_surface_net(points, triangles, curvature_area):
    references = np.array([[0.25, -0.5, 1.75], [2.0, 1.5, -0.75]])
    net = SimpleNamespace(
        surfs=pd.DataFrame([{'balls': [0, 1], 'edges': []}]),
        balls=pd.DataFrame({
            'rad': [1.0, 1.0],
            'loc': list(references),
            'num': [10, 20],
        }),
        edges=pd.DataFrame({'points': []}),
        verts=pd.DataFrame(),
        settings={
            'surf_res': 0.1,
            'net_type': 'pow',
            'surf_scheme': 'mean',
            'verbose': False,
        },
        metrics={'start': 0.0, 'vert': 0.0, 'con': 0.0},
        update_progress=lambda *args: None,
    )
    surface = (
        points, triangles,
        [], 0.0, 0.0,
        [], 0.0, 0.0,
        0.0, 0.0, 0.0,
        curvature_area,
        None, None, curvature_area is None, None,
    )
    return net, references, surface


def test_flat_surface_uses_fused_area_and_preserves_volumes(monkeypatch):
    references = np.array([
        [0.25, -0.5, 1.75],
        [2.0, 1.5, -0.75],
    ])
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ])
    triangles = np.array([[0, 1, 2], [1, 3, 2]])
    net, _, surface = _single_surface_net(points, triangles, curvature_area=None)
    monkeypatch.setattr(build_surfs_module, 'build_surf', lambda **kwargs: surface)

    build_surfs(net)

    assert net.surfs.at[0, 'sa'] == pytest.approx(
        calc_surf_sa(triangles, points), rel=STRICT_TOLERANCE, abs=STRICT_TOLERANCE
    )
    assert list(net.surfs.at[0, 'vols'].values()) == pytest.approx(
        _expected_volumes(references, points, triangles),
        rel=STRICT_TOLERANCE,
        abs=STRICT_TOLERANCE,
    )


def test_curved_surface_keeps_curvature_area_and_preserves_volumes(monkeypatch):
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 0.5],
    ])
    triangles = np.array([[0, 1, 2], [0, 3, 1], [1, 3, 4], [1, 4, 2]])
    curvature_area = calc_surf_sa(triangles, points) + 0.125
    net, references, surface = _single_surface_net(points, triangles, curvature_area)
    monkeypatch.setattr(build_surfs_module, 'build_surf', lambda **kwargs: surface)

    build_surfs(net)

    assert net.surfs.at[0, 'sa'] == pytest.approx(
        curvature_area, rel=STRICT_TOLERANCE, abs=STRICT_TOLERANCE
    )
    assert list(net.surfs.at[0, 'vols'].values()) == pytest.approx(
        _expected_volumes(references, points, triangles),
        rel=STRICT_TOLERANCE,
        abs=STRICT_TOLERANCE,
    )


def test_surface_counts_record_geometry_work(monkeypatch):
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    triangles = np.array([[0, 1, 2]])
    net, _, surface = _single_surface_net(points, triangles, curvature_area=None)
    monkeypatch.setattr(build_surfs_module, 'build_surf', lambda **kwargs: surface)

    build_surfs(net)

    assert net.surface_counts == {
        'candidate_facets': 1,
        'accepted_facets': 1,
        'discarded_facets': 0,
        'generated_points': 3,
        'generated_triangles': 1,
        'surface_integration_calls': 1,
        'analytic_bisector_calls': 0,
        'curvature_triangle_evaluations': 0,
        'projection_evaluations': 0,
        'triangulation_calls': 0,
        'containment_queries': 0,
    }


def test_batched_surface_volume_preserves_empty_surface_zero():
    actual_volumes, actual_area = _surface_tetra_volumes(
        [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        np.empty((0, 3)),
        np.empty((0, 3), dtype=int),
    )

    assert np.array_equal(actual_volumes, np.zeros(2))
    assert actual_area == 0.0


def test_batched_surface_volume_is_invariant_to_triangle_orientation():
    references = np.array([[0.25, -0.5, 1.75]])
    points = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])

    forward_volumes, forward_area = _surface_tetra_volumes(references, points, [[0, 1, 2]])
    reverse_volumes, reverse_area = _surface_tetra_volumes(references, points, [[0, 2, 1]])

    assert reverse_volumes == pytest.approx(forward_volumes, rel=0.0, abs=0.0)
    assert reverse_area == pytest.approx(forward_area, rel=0.0, abs=0.0)
    assert forward_volumes == pytest.approx(
        _expected_volumes(references, points, [[0, 1, 2]]),
        rel=STRICT_TOLERANCE,
        abs=STRICT_TOLERANCE,
    )
    assert forward_area == pytest.approx(
        calc_surf_sa([[0, 1, 2]], points),
        rel=STRICT_TOLERANCE,
        abs=STRICT_TOLERANCE,
    )
