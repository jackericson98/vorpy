import numpy as np
import pytest

from vorpy.src.geometry.validation_shapes import (
    box, cube, sphere, spherocylinder, spherocylinder_multicell,
    torus, torus_multicell, double_torus_multicell,
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


def test_double_torus_multicell_has_genus_two_figure_eight_centerline():
    loop_resolution, cross_resolution = 8, 10
    shape = double_torus_multicell(
        10.0, 3.0, interior_count=loop_resolution,
        cross_section_resolution=cross_resolution,
    )

    interior_count = 2 * loop_resolution - 1
    assert shape.parameters["genus"] == 2
    assert shape.parameters["euler_characteristic"] == -2
    assert shape.expected_int_gaussian_curvature == pytest.approx(-4.0 * np.pi)
    assert shape.parameters["interior_indices"] == list(range(interior_count))
    omitted = shape.parameters["overlap_constraints_omitted"]
    assert omitted >= 0
    assert shape.n_atoms == interior_count + interior_count * cross_resolution - omitted
    assert shape.parameters["junction_constraints_omitted"] == cross_resolution
    assert np.all(np.bincount(
        shape.parameters["constraint_assignments"], minlength=interior_count
    ) <= cross_resolution)

    centers = shape.xyzr[:interior_count, :3]
    junction = np.asarray([0.0, shape.parameters["junction_offset"],
                           1.37 * shape.parameters["junction_offset"]])
    assert np.count_nonzero(np.linalg.norm(centers - junction, axis=1) < 1e-12) == 1
    assert np.unique(np.round(centers, 10), axis=0).shape[0] == interior_count

    assignments = np.asarray(shape.parameters["constraint_assignments"])
    constraints = shape.xyzr[interior_count:, :3]
    distances = np.linalg.norm(constraints - centers[assignments], axis=1)
    assert distances == pytest.approx(
        np.full(len(constraints),
                shape.parameters["constraint_distance"])
    )


def test_dense_double_torus_uses_an_unconstrained_but_bounded_junction():
    loop_resolution, cross_resolution = 15, 12
    shape = double_torus_multicell(
        10.0, 3.0, interior_count=loop_resolution,
        cross_section_resolution=cross_resolution,
    )
    interior_count = 2 * loop_resolution - 1
    assignments = np.asarray(shape.parameters["constraint_assignments"])

    assert shape.parameters["dense_junction"]
    assert shape.parameters["junction_constraints_omitted"] == 2 * cross_resolution
    assert shape.parameters["junction_constraint_loop_retained"] is None
    assert not np.any(assignments == 0)
    assert np.all(np.bincount(assignments, minlength=interior_count)[1:] <= cross_resolution)
    assert shape.n_atoms == (
        interior_count + (interior_count - 1) * cross_resolution
        - shape.parameters["overlap_constraints_omitted"]
    )


@pytest.mark.parametrize("resolution", [8, 15, 40])
@pytest.mark.parametrize("include_center", [True, False])
@pytest.mark.parametrize("radii", [(1.0, 1.0), (2.0, 3.0)])
def test_double_torus_removes_buried_shell_spheres(resolution, include_center, radii):
    shape = double_torus_multicell(
        10.0, 3.0, interior_count=resolution, include_center=include_center,
        center_radius=radii[0], neighbor_radius=radii[1],
    )
    params = shape.parameters
    exterior = shape.xyzr[params["exterior_start"]:, :3]
    left_distance = np.hypot(np.hypot(exterior[:, 0] + 10, exterior[:, 1]) - 10,
                             exterior[:, 2])
    right_distance = np.hypot(np.hypot(exterior[:, 0] - 10, exterior[:, 2]) - 10,
                              exterior[:, 1])
    # The shared junction center is locally offset to remove multi-cell
    # co-spherical vertices; allow its small corresponding shell displacement.
    tolerance = 0.01 * params["junction_offset"]
    assert np.all(left_distance >= params["constraint_distance"] - tolerance)
    assert np.all(right_distance >= params["constraint_distance"] - tolerance)
    if resolution >= 15:
        assert params["overlap_constraints_omitted"] > 0
    assert len(exterior) == params["exterior_count"] == len(params["constraint_assignments"])
    origins = params["interior_coordinates"][params["constraint_assignments"]]
    assert np.linalg.norm(exterior - origins, axis=1) == pytest.approx(
        np.full(len(exterior), params["constraint_distance"])
    )


def test_double_torus_multicell_is_available_from_cli():
    from vorpy.src.geometry.validation_shapes import _build_parser, _generate_from_args

    args = _build_parser().parse_args([
        "--shape", "double_torus_multicell",
        "--major-radius", "10", "--minor-radius", "3",
        "--interior-count", "8", "--cross-section-resolution", "10",
    ])
    shape = _generate_from_args(args)
    assert shape.name == "double_torus_multicell"
    assert shape.parameters["genus"] == 2

    default_args = _build_parser().parse_args(["--shape", "double_torus"])
    default_shape = _generate_from_args(default_args)
    assert default_shape.parameters["interior_count_per_loop"] == 10
    assert default_shape.parameters["interior_count"] == 19


@pytest.mark.parametrize("loop_resolution", [10, 15])
def test_double_torus_group_boundary_has_genus_two(
        tmp_path, loop_resolution,
):
    from vorpy.src.command.set import sett
    from vorpy.src.group.group import Group
    from vorpy.src.system.system import System

    shape = double_torus_multicell(10.0, 3.0, interior_count=loop_resolution,
                                   cross_section_resolution=12)
    pdb = shape.save_pdb(tmp_path / "double_torus.pdb")
    system = System(file=str(pdb), make_dir=False)
    selected_count = shape.parameters["interior_count"]
    group = Group(
        system,
        name="double_torus_interior",
        atoms=list(range(selected_count)),
        settings=sett("mv", ["30"]),
        make_net=True,
    )
    group.build()
    group.get_info()

    assert all(bool(value) for value in group.net.balls.iloc[:selected_count]["complete"])
    assert group.boundary_is_complete
    assert group.boundary_is_closed
    assert group.boundary_is_manifold
    assert group.boundary_is_orientable
    assert group.boundary_component_count == 1
    assert (group.boundary_vertex_count - group.boundary_edge_count
            + group.boundary_face_count) == -2
    assert group.euler_characteristic == -2
    assert group.boundary_genus == pytest.approx(2.0)
    assert group.int_gauss_curv == pytest.approx(-4.0 * np.pi)
    # Interfaces between selected cells are excluded from the external boundary.
    assert len(group.net.surfs) > group.boundary_face_count


def test_torus_constraints_reproduce_tube_power_radius_and_alternate_phase():
    major_count, cross_count = 8, 12
    shape = torus(
        10.0, 3.0, major_count, cross_count,
        sphere_radius=1.0, center_radius=1.0,
    )
    assert shape.n_atoms == major_count + major_count * cross_count
    assert shape.parameters["constraint_distance"] == pytest.approx(6.0)
    centers = shape.xyzr[:major_count, :3]
    constraints = shape.xyzr[major_count:, :3]
    assignments = np.asarray(shape.parameters["constraint_assignments"])
    vectors = constraints - centers[assignments]
    distances = np.linalg.norm(vectors, axis=1)
    recovered = (distances**2 + 1.0**2 - 1.0**2) / (2.0 * distances)
    assert distances == pytest.approx(np.full(len(distances), 6.0))
    assert recovered == pytest.approx(np.full(len(recovered), 3.0))

    phase = 0.35 * np.pi / cross_count
    first_by_center = vectors[::cross_count] / distances[::cross_count, None]
    cross_section_angles = []
    for index, normal in enumerate(first_by_center):
        angle = 2.0 * np.pi * index / major_count
        radial = np.array([np.cos(angle), np.sin(angle), 0.0])
        cross_section_angles.append(np.arctan2(normal[2], np.dot(normal, radial)))
    assert cross_section_angles == pytest.approx(
        [0.0 if i % 2 == 0 else phase for i in range(major_count)]
    )


def test_torus_unequal_radii_keep_requested_power_support():
    shape = torus(
        10.0, 3.0, 6, 8, sphere_radius=3.0, center_radius=2.0,
    )
    distance = shape.parameters["constraint_distance"]
    recovered = (distance**2 + 2.0**2 - 3.0**2) / (2.0 * distance)
    assert recovered == pytest.approx(3.0)


def test_torus_rejects_incompatible_power_support_radii():
    with pytest.raises(ValueError, match="power-support radius"):
        torus(10.0, 3.0, 6, 8, sphere_radius=1.0, center_radius=5.0)


def test_torus_analytic_summary_reports_counts_instead_of_points():
    from vorpy.src.geometry.validation_shapes import analytic_summary

    shape = torus(10.0, 3.0, 20, 24)
    summary = analytic_summary(shape)

    assert "interior_coordinates: 20 points" in summary
    assert "constraint_assignments: 480 entries" in summary
    assert "interior_indices: 20 entries" in summary
    assert "array([[" not in summary


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


@pytest.mark.parametrize("serial, expected", [
    (99999, (99999, "A")), (100000, (0, "B")),
    (100001, (1, "B")), (1000000, (0, "L")),
])
def test_validation_pdb_serial_rollover(serial, expected):
    from vorpy.src.geometry.validation_shapes import _validation_pdb_identity
    assert _validation_pdb_identity(serial, "A") == expected


def test_validation_pdb_rollover_preserves_columns_and_reader_indices(tmp_path):
    from vorpy.src.geometry.validation_shapes import ValidationShape
    from vorpy.src.system.system import System

    shape = ValidationShape("torus", np.tile([12.5, -3.25, 1.5, 1.0], (100001, 1)), 0, 0)
    path = shape.save_pdb(tmp_path / "rollover.pdb")
    records = [line for line in path.read_text().splitlines() if line.startswith("HETATM")]
    selected = [records[0], *records[99998:]]
    assert [(int(line[6:11]), line[21]) for line in selected] == [
        (1, "A"), (99999, "A"), (0, "B"), (1, "B"),
    ]
    for line in selected:
        assert line[17:20] == "VAL"
        assert int(line[22:26]) == 2
        assert [float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]] == [12.5, -3.25, 1.5]
    sample = tmp_path / "sample.pdb"
    sample.write_text("HEADER test\n" + "\n".join(selected) + "\nEND\n")
    system = System(file=str(sample), make_dir=False, print_actions=False)
    assert system.balls["num"].tolist() == [0, 1, 2, 3]
    assert [chain.name for chain in system.chains] == ["A", "B"]


def test_validation_pdb_rejects_exhausted_chains():
    from vorpy.src.geometry.validation_shapes import _validation_pdb_identity
    with pytest.raises(ValueError, match="exhausted"):
        _validation_pdb_identity(6100000, "A")
    with pytest.raises(ValueError, match="99,999"):
        _validation_pdb_identity(100000, "A", interior=True)


def test_pdb_formatter_rejects_serial_overflow():
    from vorpy.src.output import make_pdb_line
    with pytest.raises(ValueError, match="serial"):
        make_pdb_line(ser_num=1000000)
