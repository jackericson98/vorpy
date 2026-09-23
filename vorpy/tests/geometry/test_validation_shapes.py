import numpy as np
import pytest

from vorpy.src.geometry.validation_shapes import (
    box, cube, sphere, spherocylinder, spherocylinder_multicell,
    torus_multicell,
)


def test_sphere_radius_is_target_power_cell_radius_for_equal_weights():
    shape = sphere(10.0, resolution=20, ball_radius=1.0, center_radius=1.0)
    distances = np.linalg.norm(shape.xyzr[1:, :3], axis=1)
    assert distances == pytest.approx(np.full(20, 20.0))
    assert shape.parameters["target_power_cell_radius"] == 10.0
    assert shape.parameters["surrounding_generator_distance"] == pytest.approx(20.0)


def test_sphere_unequal_weights_reproduce_requested_power_radius():
    target, center_radius, ball_radius = 10.0, 2.0, 3.0
    shape = sphere(
        target, resolution=20,
        ball_radius=ball_radius,
        center_radius=center_radius,
    )
    distance = shape.parameters["surrounding_generator_distance"]
    recovered = (distance**2 + center_radius**2 - ball_radius**2) / (2.0 * distance)
    assert recovered == pytest.approx(target)
    assert np.linalg.norm(shape.xyzr[1:, :3], axis=1) == pytest.approx(
        np.full(20, distance)
    )
    assert shape.expected_int_mean_curvature == pytest.approx(40.0 * np.pi)
    assert shape.expected_int_gaussian_curvature == pytest.approx(4.0 * np.pi)


def test_sphere_rejects_impossible_power_bisector():
    with pytest.raises(ValueError, match="power-cell face"):
        sphere(1.0, resolution=20, center_radius=2.0, ball_radius=1.0)


def test_cube_is_seven_weighted_generators_with_expected_planes():
    shape = cube(20.0, resolution=10, sphere_radius=1.0)
    assert shape.n_atoms == 7
    assert np.allclose(shape.xyzr[1:, 3], 1.0)
    assert np.allclose(np.abs(shape.xyzr[1:, :3]),
                       [[20, 0, 0], [20, 0, 0], [0, 20, 0], [0, 20, 0], [0, 0, 20], [0, 0, 20]])
    assert shape.expected_int_mean_curvature == pytest.approx(60.0 * np.pi)
    assert shape.expected_int_gaussian_curvature == pytest.approx(4.0 * np.pi)


def test_box_uses_independent_power_bisector_distances_and_radii():
    shape = box(20.0, 10.0, 6.0, sphere_radius=2.0, center_radius=1.0)
    assert shape.n_atoms == 7
    centers = shape.xyzr[1:, :3]
    distances = sorted(np.linalg.norm(centers, axis=1))
    assert distances == pytest.approx(sorted([
        10.0 + np.sqrt(10.0**2 + 3.0),
        10.0 + np.sqrt(10.0**2 + 3.0),
        5.0 + np.sqrt(5.0**2 + 3.0),
        5.0 + np.sqrt(5.0**2 + 3.0),
        3.0 + np.sqrt(3.0**2 + 3.0),
        3.0 + np.sqrt(3.0**2 + 3.0),
    ]))
    assert shape.xyzr[0, 3] == 1.0
    assert np.allclose(shape.xyzr[1:, 3], 2.0)
    assert shape.expected_int_mean_curvature == pytest.approx(36.0 * np.pi)


def test_spherocylinder_has_one_center_and_linear_generator_count():
    shape = spherocylinder(10.0, 20.0, resolution=50)
    assert shape.n_atoms == 51
    assert np.allclose(shape.xyzr[0], [0.0, 0.0, 0.0, 1.0])
    assert spherocylinder(10.0, 20.0, resolution=20).n_atoms == 21
    assert spherocylinder(10.0, 20.0, resolution=100).n_atoms == 101


def test_spherocylinder_support_normals_and_power_planes():
    radius, length = 10.0, 20.0
    shape = spherocylinder(radius, length, resolution=50)
    centers = shape.xyzr[1:, :3]
    distances = np.linalg.norm(centers, axis=1)
    normals = centers / distances[:, None]
    r0 = shape.xyzr[0, 3]
    r1 = shape.xyzr[1, 3]

    assert np.linalg.norm(normals, axis=1) == pytest.approx(np.ones(50))
    assert np.unique(np.round(normals, 10), axis=0).shape[0] == 50
    assert np.allclose(centers / distances[:, None], normals)

    support = radius + 0.5 * length * np.abs(normals[:, 2])
    recovered = (distances**2 + r0**2 - r1**2) / (2.0 * distances)
    assert recovered == pytest.approx(support)
    assert np.min(support) == pytest.approx(radius)
    assert np.max(support) == pytest.approx(radius + length / 2.0)
    assert shape.expected_int_mean_curvature == pytest.approx(np.pi * length + 4 * np.pi * radius)
    assert shape.expected_int_gaussian_curvature == pytest.approx(4 * np.pi)


def test_spherocylinder_unequal_radii_preserve_support_planes():
    shape = spherocylinder(
        10.0, 20.0, resolution=20,
        center_radius=2.0, neighbor_radius=3.0,
    )
    centers = shape.xyzr[1:, :3]
    distances = np.linalg.norm(centers, axis=1)
    normals = centers / distances[:, None]
    recovered = (distances**2 + 2.0**2 - 3.0**2) / (2.0 * distances)
    assert recovered == pytest.approx(10.0 + 10.0 * np.abs(normals[:, 2]))


def test_multicell_spherocylinder_group_and_constraint_geometry():
    shape = spherocylinder_multicell(10.0, 20.0, interior_count=5,
                                     angular_resolution=8, cap_resolution=8)
    assert shape.n_atoms == 5 + 5 * (8 * 3 + 8)
    assert shape.parameters["interior_indices"] == [0, 1, 2, 3, 4]
    assert np.asarray(shape.parameters["interior_coordinates"])[:, 2] == pytest.approx(
        [-10.0, -5.0, 0.0, 5.0, 10.0]
    )

    origins = np.asarray(shape.parameters["interior_coordinates"])
    assignments = shape.parameters["constraint_assignments"]
    exterior = shape.xyzr[5:, :3]
    local_vectors = exterior - origins[np.asarray(assignments)]
    distances = np.linalg.norm(local_vectors, axis=1)
    normals = local_vectors / distances[:, None]
    support = np.einsum("ij,ij->i", normals, local_vectors)

    assert distances == pytest.approx(np.full(len(distances), 20.0))
    assert support == pytest.approx(np.full(len(support), 20.0))
    assert shape.parameters["target_volume"] == pytest.approx(10000.0 * np.pi / 3.0)
    assert shape.parameters["target_surface_area"] == pytest.approx(800.0 * np.pi)
    assert shape.expected_int_mean_curvature == pytest.approx(60.0 * np.pi)
    assert shape.expected_int_gaussian_curvature == pytest.approx(4.0 * np.pi)


def test_torus_multicell_has_ring_interior_and_local_constraints():
    shape = torus_multicell(10.0, 3.0, interior_count=12,
                            cross_section_resolution=12)
    assert shape.n_atoms == 12 + 12 * 12
    assert shape.parameters["interior_indices"] == list(range(12))
    assert shape.parameters["exterior_count"] == 144
    assert np.allclose(np.linalg.norm(shape.xyzr[:12, :3], axis=1), 10.0)
    assert np.allclose(shape.xyzr[:12, 2], 0.0)
    # The selected ring never occupies the central hole.
    assert np.min(np.linalg.norm(shape.xyzr[:12, :2], axis=1)) == pytest.approx(10.0)
    assert shape.expected_int_mean_curvature == pytest.approx(20.0 * np.pi**2)
    assert shape.expected_int_gaussian_curvature == pytest.approx(0.0)


def test_torus_multicell_constraints_reproduce_tube_power_radius():
    shape = torus_multicell(10.0, 3.0, interior_count=12,
                            cross_section_resolution=12,
                            center_radius=2.0, neighbor_radius=3.0)
    d = shape.parameters["constraint_distance"]
    recovered = (d * d + 2.0**2 - 3.0**2) / (2.0 * d)
    assert recovered == pytest.approx(3.0)
    centers = shape.xyzr[:12, :3]
    constraints = shape.xyzr[12:, :3]
    assignment = shape.parameters["constraint_assignments"]
    distances = np.linalg.norm(constraints - centers[np.asarray(assignment)], axis=1)
    assert distances == pytest.approx(np.full(144, d))


def test_multicell_default_is_three_cell_short_spherocylinder():
    shape = spherocylinder_multicell(10.0, 10.0)
    assert shape.parameters["interior_count"] == 3
    assert np.asarray(shape.parameters["interior_coordinates"])[:, 2] == pytest.approx(
        [-5.0, 0.0, 5.0]
    )
    assert shape.n_atoms == 3 + 3 * (12 * 3 + 12)


def test_multicell_pdb_marks_interior_chain_separately(tmp_path):
    shape = spherocylinder_multicell(10.0, 20.0, interior_count=5)
    path = shape.save_pdb(tmp_path / "multi.pdb")
    atom_lines = [line for line in path.read_text().splitlines() if line.startswith("HETATM")]
    assert sum(" INT " in line for line in atom_lines) == 5
    assert not any(" EXT " in line for line in atom_lines)
