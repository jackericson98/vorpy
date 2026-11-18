import numpy as np

from vorpy.src.geometry.periodic_box import (
    point_in_vertex_box,
    scaled_vertex_region,
    vertex_region_from_cryst1,
)
from vorpy.src.network import Network


CRYST1 = 'CRYST1   39.854   39.854   39.854  60.00  60.00  90.00 P 1           1'


def test_rhombic_dodecahedron_rejects_aabb_corner():
    region = vertex_region_from_cryst1(CRYST1)
    center = np.asarray(region['center'])
    box = [[-100.0] * 3, [100.0] * 3, scaled_vertex_region(region, 1.0)]

    assert point_in_vertex_box(center, box)
    assert not point_in_vertex_box(center + [19.0, 19.0, 19.0], box)


def test_cryst1_region_replaces_atom_derived_aabb_limit():
    region = vertex_region_from_cryst1(CRYST1)
    center = np.asarray(region['center'])
    box = [center.tolist(), center.tolist(), scaled_vertex_region(region, 1.0)]

    assert point_in_vertex_box(center + [1.0, 0.0, 0.0], box)


def test_box_multiplier_scales_cryst1_region_about_center():
    region = vertex_region_from_cryst1(CRYST1)
    center = np.asarray(region['center'])
    point = center + [0.0, 0.0, 30.0]

    assert not point_in_vertex_box(
        point, [[-100.0] * 3, [100.0] * 3, scaled_vertex_region(region, 1.0)]
    )
    assert point_in_vertex_box(
        point, [[-100.0] * 3, [100.0] * 3, scaled_vertex_region(region, 1.25)]
    )


def test_orthogonal_cryst1_uses_six_face_cell():
    region = vertex_region_from_cryst1(
        'CRYST1   20.000   30.000   40.000  90.00  90.00  90.00 P 1           1'
    )

    assert region['face_count'] == 6


def test_legacy_aabb_is_unchanged_without_cryst1_region():
    box = [[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]]

    assert point_in_vertex_box([1.0, 1.0, 1.0], box)
    assert not point_in_vertex_box([3.0, 1.0, 1.0], box)


def test_spatial_grid_keeps_cryst1_region_separate_from_corner_array():
    region = vertex_region_from_cryst1(CRYST1)
    network = Network(
        locs=[[0.0, 0.0, 0.0], [40.0, 0.0, 0.0],
              [0.0, 40.0, 0.0], [0.0, 0.0, 40.0]],
        rads=[1.0] * 4,
        settings={'box_size': 0.95, 'vertex_region': region},
    )

    network.sort_balls(num_boxes=8)

    assert np.asarray(network.box['verts']).shape == (2, 3)
    assert network.box['vertex_region']['face_count'] == 12
    assert network.box['vertex_region']['multiplier'] == 0.95
    assert all(
        0 <= index < network.settings['num_splits']
        for indices in network.balls['box'] for index in indices
    )
