import numpy as np
import pytest

from vorpy.src.geometry.validation_shapes import (
    _build_parser, _default_basename, _generate_from_args,
    double_torus_multicell, multitorus, torus_multicell,
)
from vorpy.src.geometry.topology_diagnostic import acceptance_errors, inspect_group


@pytest.mark.parametrize("count", [1, 2, 3, 7, 20])
@pytest.mark.parametrize("resolution", [15, 16])
def test_multitorus_shared_junctions_and_shell_union(count, resolution):
    shape = multitorus(10, 3, torus_count=count, interior_count=resolution)
    p = shape.parameters
    samples = resolution + int(count > 2 and resolution % 2 != 0)
    assert p["interior_count_per_loop"] == samples
    assert p["interior_count"] == count * samples - (count - 1)
    assert p["genus"] == count
    assert p["euler_characteristic"] == 2 - 2 * count
    assert shape.expected_int_gaussian_curvature == pytest.approx(4 * np.pi * (1 - count))
    interior = p["interior_coordinates"]
    assert len(np.unique(np.round(interior, 10), axis=0)) == len(interior)
    for join in range(count - 1):
        junction = [(2 * join - count + 2) * 10,
                    p["junction_offset"], 1.37 * p["junction_offset"]]
        assert np.count_nonzero(np.linalg.norm(interior - junction, axis=1) < 1e-9) == 1
    exterior = shape.xyzr[p["exterior_start"]:, :3]
    for loop in range(count):
        relative = exterior - [(2 * loop - count + 1) * 10, 0, 0]
        in_plane, perpendicular = (1, 2) if loop % 2 == 0 else (2, 1)
        distance = np.hypot(np.hypot(relative[:, 0], relative[:, in_plane]) - 10,
                            relative[:, perpendicular])
        assert np.all(distance >= p["constraint_distance"] - 1e-9)
    assert len(exterior) == p["exterior_count"] == len(p["constraint_assignments"])
    origins = interior[p["constraint_assignments"]]
    assert np.linalg.norm(exterior - origins, axis=1) == pytest.approx(
        np.full(len(exterior), p["constraint_distance"])
    )


def test_multitorus_one_and_two_match_existing_generators():
    assert multitorus(10, 3, 1, 16).xyzr == pytest.approx(torus_multicell(10, 3, 16).xyzr)
    assert multitorus(10, 3, 2, 15).xyzr == pytest.approx(double_torus_multicell(10, 3, 15).xyzr)


@pytest.mark.parametrize("count", [0, -1, 1.5, True, float("inf"), None])
def test_multitorus_rejects_invalid_counts(count):
    with pytest.raises(ValueError, match="torus_count"):
        multitorus(10, 3, count)


@pytest.mark.parametrize("include_center", [True, False])
def test_multitorus_cli_and_exports(tmp_path, include_center):
    args = ["--shape", "multitorus", "--torus-count", "4", "--interior-count", "15",
            "--center-radius", "2", "--neighbor-radius", "3"]
    if not include_center:
        args.append("--no-center")
    shape = _generate_from_args(_build_parser().parse_args(args))
    assert shape.name == "multitorus"
    assert _default_basename(shape) == "multitorus_n4_R10_r3_i16_x12"
    p = shape.parameters
    assert p["interior_indices"] == (list(range(61)) if include_center else [])
    assert p["exterior_start"] == (61 if include_center else 0)
    xyzr, pdb, pml = shape.save(tmp_path, _default_basename(shape))
    assert np.loadtxt(xyzr).shape == shape.xyzr.shape
    records = [line for line in pdb.read_text().splitlines() if line.startswith("HETATM")]
    assert len(records) == shape.n_atoms
    assert sum(" INT " in line for line in records) == (61 if include_center else 0)
    assert "load " + pdb.name + ", multitorus" in pml.read_text()


@pytest.mark.parametrize("count", [1, 2, 3, 5])
def test_multitorus_boundary_topology(count):
    # Equal weights reduce to ordinary Voronoi cells. Inspect the union's
    # boundary directly, independently of the generator's target metadata.
    from scipy.spatial import Voronoi

    shape = multitorus(10, 3, torus_count=count, interior_count=16)
    vor = Voronoi(shape.xyzr[:, :3])
    selected = shape.parameters["interior_count"]
    vertices, edge_uses, faces = set(), {}, 0
    for pair, ridge in zip(vor.ridge_points, vor.ridge_vertices):
        if (pair[0] < selected) == (pair[1] < selected):
            continue
        assert -1 not in ridge
        faces += 1
        vertices.update(ridge)
        for a, b in zip(ridge, ridge[1:] + ridge[:1]):
            edge = tuple(sorted((a, b)))
            edge_uses[edge] = edge_uses.get(edge, 0) + 1
    assert set(edge_uses.values()) == {2}
    assert len(vertices) - len(edge_uses) + faces == 2 - 2 * count
    neighbors = {vertex: set() for vertex in vertices}
    for a, b in edge_uses:
        neighbors[a].add(b)
        neighbors[b].add(a)
    visited, pending = set(), [next(iter(vertices))]
    while pending:
        vertex = pending.pop()
        if vertex not in visited:
            visited.add(vertex)
            pending.extend(neighbors[vertex] - visited)
    assert visited == vertices


def test_vorpy_recovers_genus_two_from_assembled_boundary(tmp_path):
    from vorpy.src.command.set import sett
    from vorpy.src.group.group import Group
    from vorpy.src.system.system import System

    shape = multitorus(10, 3, torus_count=2, interior_count=16)
    pdb = shape.save_pdb(tmp_path / "multitorus_g2.pdb")
    system = System(file=str(pdb), make_dir=False, print_actions=False)
    selected = list(range(shape.parameters["interior_count"]))
    group = Group(system, name="multitorus_g2", atoms=selected,
                  settings=sett("mv", ["30"]), make_net=True,
                  print_metrics=False)
    group.build()
    group.get_info()
    coords = shape.parameters["interior_coordinates"]
    offset = shape.parameters["junction_offset"]
    loops = []
    for loop in range(2):
        center = [(2 * loop - 1) * 10, 0, 0]
        relative = coords - center
        plane, perpendicular = (1, 2) if loop % 2 == 0 else (2, 1)
        distance = np.hypot(np.hypot(relative[:, 0], relative[:, plane]) - 10,
                            relative[:, perpendicular])
        loops.append([i for i in range(len(coords)) if distance[i] <
                      1.5 * np.hypot(offset, 1.37 * offset)])
    report = inspect_group(group, loops, [[0, 0, 0]])
    assert acceptance_errors(report, genus=2) == []
    assert report["selected_cells"] == report["complete_selected_cells"] == 31
    assert (report["boundary_vertices"], report["boundary_edges"],
            report["boundary_faces"]) == (1354, 2369, 1013)
    assert report["euler_characteristic"] == -2
    assert report["integrated_gaussian_curvature"] == pytest.approx(-4 * np.pi, abs=1e-8)
