import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.calculations.polyhedral_gaussian_curvature import (
    POLYHEDRAL_STATUS_CALCULATED,
    POLYHEDRAL_STATUS_UNAVAILABLE,
    _face_angles,
    calculate_polyhedral_gaussian_curvature,
)
from vorpy.src.group.sort import get_info


_CUBE_FACES = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (3, 7, 6, 2),
    (0, 4, 7, 3),
    (1, 2, 6, 5),
)


def _cube_network(scheme="pow", offset=(0.0, 0.0, 0.0), vertex_offset=0, edge_offset=0):
    base = np.asarray([
        (-1.0, -1.0, -1.0), (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0), (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0),
    ]) + np.asarray(offset, dtype=float)
    vertices = [
        {"loc": point, "num": index + vertex_offset}
        for index, point in enumerate(base)
    ]
    edges = []
    edge_ids = {}
    faces = []
    for face_vertices in _CUBE_FACES:
        face_edges = []
        for index, start in enumerate(face_vertices):
            end = face_vertices[(index + 1) % len(face_vertices)]
            key = frozenset((start + vertex_offset, end + vertex_offset))
            if key not in edge_ids:
                edge_ids[key] = edge_offset + len(edges)
                edges.append({"verts": sorted(key)})
            face_edges.append(edge_ids[key])
        faces.append({
            "verts": sorted(vertex + vertex_offset for vertex in face_vertices),
            "edges": face_edges,
            "flat": True,
        })
    return SimpleNamespace(
        settings={"net_type": scheme},
        verts=pd.DataFrame(vertices),
        edges=pd.DataFrame(edges),
        surfs=pd.DataFrame(faces),
    )


def _merge_cubes(first, second):
    return SimpleNamespace(
        settings={"net_type": "pow"},
        verts=pd.concat([first.verts, second.verts], ignore_index=True),
        edges=pd.concat([first.edges, second.edges], ignore_index=True),
        surfs=pd.concat([first.surfs, second.surfs], ignore_index=True),
    )


def _torus_network():
    major_radius, minor_radius = 3.0, 1.0
    u_count, v_count = 8, 6
    vertices = []
    for u_index in range(u_count):
        u = 2.0 * math.pi * u_index / u_count
        for v_index in range(v_count):
            v = 2.0 * math.pi * v_index / v_count
            vertices.append({
                "loc": np.asarray((
                    (major_radius + minor_radius * math.cos(v)) * math.cos(u),
                    (major_radius + minor_radius * math.cos(v)) * math.sin(u),
                    minor_radius * math.sin(v),
                )),
                "num": len(vertices),
            })

    edges = []
    edge_ids = {}
    faces = []

    def vertex_id(u_index, v_index):
        return (u_index % u_count) * v_count + (v_index % v_count)

    def edge_id(first, second):
        key = frozenset((first, second))
        if key not in edge_ids:
            edge_ids[key] = len(edges)
            edges.append({"verts": sorted(key)})
        return edge_ids[key]

    for u_index in range(u_count):
        for v_index in range(v_count):
            a = vertex_id(u_index, v_index)
            b = vertex_id(u_index + 1, v_index)
            c = vertex_id(u_index + 1, v_index + 1)
            d = vertex_id(u_index, v_index + 1)
            for loop in ((a, b, c), (a, c, d)):
                faces.append({
                    "verts": sorted(loop),
                    "edges": [
                        edge_id(loop[index], loop[(index + 1) % 3])
                        for index in range(3)
                    ],
                    "flat": True,
                })
    return SimpleNamespace(
        settings={"net_type": "pow"},
        verts=pd.DataFrame(vertices),
        edges=pd.DataFrame(edges),
        surfs=pd.DataFrame(faces),
    )


@pytest.mark.parametrize("scheme", ["pow", "prm"])
def test_cube_polyhedral_gaussian_curvature_matches_independent_gauss_bonnet(scheme):
    result = calculate_polyhedral_gaussian_curvature(_cube_network(scheme), range(6))

    assert result.status == POLYHEDRAL_STATUS_CALCULATED
    assert result.total == pytest.approx(4.0 * math.pi, abs=1e-12)
    assert result.euler_characteristic == 2
    assert result.gauss_bonnet_residual == pytest.approx(0.0, abs=1e-12)
    assert all(value == pytest.approx(math.pi / 2.0) for value in result.vertex_curvatures.values())
    assert all(len(loop) == 4 for loop in result.face_loops.values())


def test_concave_polygon_uses_reflex_interior_angle():
    points = np.asarray([
        (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0),
        (1.0, 1.0, 0.0), (0.0, 2.0, 0.0),
    ])
    angles = _face_angles(points, np.asarray((0.0, 0.0, 1.0)), 0, 1e-8)

    assert max(angles) > math.pi
    assert angles[3] == pytest.approx(1.5 * math.pi)


def test_multiple_closed_components_are_summed_without_chi_substitution():
    first = _cube_network(offset=(0.0, 0.0, 0.0))
    second = _cube_network(offset=(5.0, 0.0, 0.0), vertex_offset=8, edge_offset=12)
    result = calculate_polyhedral_gaussian_curvature(
        _merge_cubes(first, second), range(12)
    )

    assert result.status == POLYHEDRAL_STATUS_CALCULATED
    assert result.component_count == 2
    assert result.euler_characteristic == 4
    assert result.total == pytest.approx(8.0 * math.pi, abs=1e-12)


def test_closed_higher_genus_boundary_is_supported():
    result = calculate_polyhedral_gaussian_curvature(_torus_network(), range(96))

    assert result.status == POLYHEDRAL_STATUS_CALCULATED
    assert result.component_count == 1
    assert result.euler_characteristic == 0
    assert result.total == pytest.approx(0.0, abs=1e-12)
    assert result.gauss_bonnet_residual == pytest.approx(0.0, abs=1e-12)


def test_nonplanar_face_is_explicitly_unavailable():
    network = _cube_network("prm")
    network.verts.at[2, "loc"] = np.asarray((1.0, 1.0, -0.8))
    result = calculate_polyhedral_gaussian_curvature(network, range(6))

    assert result.status == POLYHEDRAL_STATUS_UNAVAILABLE
    assert result.total is None
    assert any("not planar" in reason for reason in result.reasons)


def test_open_boundary_is_not_reported_as_zero():
    result = calculate_polyhedral_gaussian_curvature(_cube_network(), range(5))

    assert result.status == POLYHEDRAL_STATUS_UNAVAILABLE
    assert result.total is None
    assert any("exactly two faces" in reason for reason in result.reasons)


def test_aw_network_is_rejected_by_polyhedral_calculator():
    with pytest.raises(ValueError, match="POW or PRM"):
        calculate_polyhedral_gaussian_curvature(_cube_network("aw"), range(6))


@pytest.mark.parametrize("scheme", ["pow", "prm"])
def test_group_analysis_dispatches_polyhedral_gaussian_curvature(scheme):
    network = _cube_network(scheme)
    locations = [
        (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
        (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1), (0, 0, 0),
    ]
    network.balls = pd.DataFrame({
        "loc": [np.asarray(point, dtype=float) for point in locations],
        "num": list(range(9)), "rad": [1.0] * 9,
        "complete": [True] * 9, "vol": [1.0] * 9,
        "vdw_vol": [1.0] * 9, "mass": [1.0] * 9,
        "com": [np.asarray(point, dtype=float) for point in locations],
    })
    network.surfs["balls"] = [[0, 8]] * 6
    network.surfs["sa"] = [1.0] * 6
    network.surfs["points"] = [[]] * 6
    network.surfs["tris"] = [[]] * 6
    network.surfs["int_mean_curv"] = [0.0] * 6
    network.surfs["int_mean_curv_sq"] = [0.0] * 6
    network.surfs["int_gauss_curv"] = [0.0] * 6
    network.surfs["flat"] = [True] * 6
    network.surfs["com"] = [np.zeros(3)] * 6
    network.edges["length"] = [1.0] * len(network.edges)
    group = SimpleNamespace(
        net=network,
        ball_ndxs=list(range(8)),
        layer_surfs=[list(range(6))],
        layer_net_atoms=[list(range(8))],
        get_layers=lambda max_layers=1: None,
    )

    get_info(group)

    assert group.gaussian_curvature_status == POLYHEDRAL_STATUS_CALCULATED
    assert group.int_gauss_curv_face == 0.0
    assert group.int_gauss_curv_edge == 0.0
    assert group.int_gauss_curv_vertex == pytest.approx(4.0 * math.pi)
    assert group.int_gauss_curv == pytest.approx(4.0 * math.pi)
