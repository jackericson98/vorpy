import numpy as np
import pytest

from vorpy.src.calculations.aw_interface_orientation import (
    aw_interface_edge_orientation,
    induced_boundary_tangent,
    signed_dihedral,
)
from vorpy.src.calculations.surface_geometry import QuadraticSurfaceGeometry


def _aw_case(group_a={0}, group_b={1, 2}, permutation=None, reverse_rows=False):
    xyz = {
        0: np.array([0., 0., 0.]),
        1: np.array([2., 0., 0.]),
        2: np.array([0., 2., 0.]),
    }
    if permutation is None:
        permutation = {0: 0, 1: 1, 2: 2}
    remapped_a = {permutation[x] for x in group_a}
    remapped_b = {permutation[x] for x in group_b}
    point = np.array([1., 1., 0.])
    pairs = [(permutation[0], permutation[1]), (permutation[0], permutation[2])]
    if reverse_rows:
        pairs.reverse()
    locations = {permutation[old]: xyz[old] for old in xyz}
    return point, pairs, tuple(permutation[x] for x in (0, 1, 2)), remapped_a, remapped_b, locations


def test_aw_edge_orientation_is_independent_of_generator_and_surface_order():
    reference = aw_interface_edge_orientation(*_aw_case())
    permuted = aw_interface_edge_orientation(
        *_aw_case(permutation={0: 2, 1: 0, 2: 1}, reverse_rows=True)
    )
    assert reference["beta"] == pytest.approx(np.pi / 2)
    assert permuted["beta"] == pytest.approx(reference["beta"])
    assert permuted["tangent"] == pytest.approx(reference["tangent"])


def test_reversing_groups_reverses_smooth_orientation_and_edge_beta():
    forward = aw_interface_edge_orientation(*_aw_case())
    reverse = aw_interface_edge_orientation(*_aw_case({1, 2}, {0}))
    assert reverse["beta"] == pytest.approx(-forward["beta"])
    assert abs(reverse["beta"]) == pytest.approx(abs(forward["beta"]))


def test_cube_has_positive_convex_edge_measure_and_zero_smooth_term():
    # Unit cube: twelve unit edges, each with a 90-degree outward turn.
    contributions = []
    for axis_a in range(3):
        for axis_b in range(axis_a + 1, 3):
            axis_c = 3 - axis_a - axis_b
            for sign_a in (-1, 1):
                for sign_b in (-1, 1):
                    n1 = np.eye(3)[axis_a] * sign_a
                    n2 = np.eye(3)[axis_b] * sign_b
                    # The first face's outward co-normal is the second face
                    # normal at their common cube edge.
                    tangent = induced_boundary_tangent(n1, n2)
                    beta = signed_dihedral(tangent, n1, n2)
                    contributions.append(beta)
    c_smooth = 0.0
    c_edge_signed = float(sum(contributions))
    c_edge_unsigned = float(sum(abs(value) for value in contributions))
    c_h = c_smooth + 0.5 * c_edge_signed
    assert len(contributions) == 12
    assert c_edge_signed == pytest.approx(6 * np.pi)
    assert c_edge_unsigned == pytest.approx(6 * np.pi)
    assert c_h == pytest.approx(3 * np.pi)


def test_concave_ridge_has_negative_signed_edge_measure():
    n1 = np.array([1., 0., 0.])
    n2 = np.array([0., -1., 0.])
    m1 = np.array([0., 1., 0.])
    tangent = induced_boundary_tangent(n1, m1)
    beta = signed_dihedral(tangent, n1, n2)
    assert beta == pytest.approx(-np.pi / 2)
    assert abs(beta) == pytest.approx(np.pi / 2)


def test_sphere_quadratic_geometry_has_known_smooth_integral():
    radius = 2.5
    geometry = QuadraticSurfaceGeometry([
        1, 1, 1, 0, 0, 0, 0, 0, 0, -radius**2, 0, 0, 0, 0,
    ])
    nodes, weights = np.polynomial.legendre.leggauss(32)
    phi_count = 128
    phi_weight = 2 * np.pi / phi_count
    c_smooth = 0.0
    for mu, weight in zip(nodes, weights, strict=True):
        rho = np.sqrt(1 - mu * mu)
        for phi_index in range(phi_count):
            phi = (phi_index + 0.5) * phi_weight
            point = radius * np.array([rho * np.cos(phi), rho * np.sin(phi), mu])
            c_smooth += geometry.mean_curvature(point) * radius**2 * float(weight) * phi_weight
    c_edge_signed = 0.0
    c_edge_unsigned = 0.0
    c_h = c_smooth + 0.5 * c_edge_signed
    assert c_smooth == pytest.approx(4 * np.pi * radius, abs=1e-11)
    assert c_edge_signed == c_edge_unsigned == 0.0
    assert c_h == pytest.approx(10 * np.pi, abs=1e-11)


def test_edge_orientation_is_rigid_motion_invariant():
    reference = aw_interface_edge_orientation(*_aw_case())
    angle = 0.73
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.],
        [np.sin(angle), np.cos(angle), 0.],
        [0., 0., 1.],
    ])
    translation = np.array([13., -7., 4.])
    point, pairs, ids, ga, gb, locations = _aw_case()
    transformed = aw_interface_edge_orientation(
        (rotation @ point) + translation,
        pairs,
        ids,
        ga,
        gb,
        {key: rotation @ value + translation for key, value in locations.items()},
    )
    assert transformed["beta"] == pytest.approx(reference["beta"], abs=1e-12)
    assert transformed["tangent"] == pytest.approx(rotation @ reference["tangent"])


def test_reversing_edge_parameter_flips_signed_dihedral_only():
    """Parameter reversal changes the signed tangent convention, not geometry."""
    tangent = np.array([0.0, 0.0, 1.0])
    normal_1 = np.array([1.0, 0.0, 0.0])
    normal_2 = np.array([0.0, 1.0, 0.0])

    forward = signed_dihedral(tangent, normal_1, normal_2)
    reversed_parameter = signed_dihedral(-tangent, normal_1, normal_2)

    assert forward == pytest.approx(np.pi / 2.0)
    assert reversed_parameter == pytest.approx(-np.pi / 2.0)
    assert reversed_parameter == pytest.approx(-forward)
    assert abs(reversed_parameter) == pytest.approx(abs(forward))


def test_surface_boundary_tangent_uses_right_hand_rule():
    normal = np.array([0.0, 0.0, 1.0])
    outward_conormal = np.array([1.0, 0.0, 0.0])

    assert induced_boundary_tangent(normal, outward_conormal) == pytest.approx(
        [0.0, 1.0, 0.0]
    )
    assert induced_boundary_tangent(-normal, outward_conormal) == pytest.approx(
        [0.0, -1.0, 0.0]
    )
