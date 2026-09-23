import numpy as np
import pytest

from vorpy.src.calculations import boundary_turning_angle
from vorpy.src.geometry.validation_shapes import spherocylinder_multicell
from vorpy.src.group.group import Group
from vorpy.src.system.system import System


@pytest.fixture(scope="module")
def multicell_group(tmp_path_factory):
    path = tmp_path_factory.mktemp("multicell") / "validation.pdb"
    spherocylinder_multicell(
        10.0, 10.0, interior_count=3, angular_resolution=12,
        cap_resolution=12, axial_ring_count=3,
    ).save_pdb(path)
    system = System(file=str(path), make_dir=False, print_actions=False)
    group = Group(
        sys=system, name="multicell_test", chains=[system.chains[0]],
        make_net=True, print_metrics=False,
    )
    group.build()
    group.get_info()
    return group


def _external_boundary_imc(group):
    """Reconstruct group IMC from retained faces, independently of reduction."""
    net = group.net
    body = {int(value) for value in group.layer_net_atoms[0]}
    retained = {int(value) for value in group.layer_surfs[0]}
    owner = {}
    for surf_id in retained:
        balls = {int(value) for value in net.surfs.iloc[surf_id]["balls"]}
        cells = balls.intersection(body)
        assert len(cells) == 1
        owner[surf_id] = next(iter(cells))

    edge_faces = {}
    for surf_id in retained:
        for edge_id in net.surfs.iloc[surf_id].get("edges", []):
            edge_faces.setdefault(int(edge_id), []).append(surf_id)

    total = 0.0
    for edge_id, faces in edge_faces.items():
        assert len(faces) == 2
        normals = []
        for surf_id in faces:
            balls = tuple(int(value) for value in net.surfs.iloc[surf_id]["balls"])
            cell = owner[surf_id]
            other = next(value for value in balls if value != cell)
            vector = (
                np.asarray(net.balls.iloc[other]["loc"], dtype=float)
                - np.asarray(net.balls.iloc[cell]["loc"], dtype=float)
            )
            normals.append(vector / np.linalg.norm(vector))
        theta = boundary_turning_angle(normals[0], normals[1])
        total += 0.5 * float(net.edges.iloc[edge_id]["length"]) * theta
    return total


def test_multicell_group_uses_external_boundary_mean_curvature(multicell_group):
    group = multicell_group
    expected = _external_boundary_imc(group)
    cell_sum = float(group.net.balls.iloc[:3]["int_mean_curv"].sum())

    assert group.int_mean_curv == pytest.approx(expected, abs=1e-8)
    assert group.int_mean_curv != pytest.approx(cell_sum, abs=1e-4)
    assert group.int_mean_curv == pytest.approx(160.0410289477, abs=1e-7)
    assert group.int_gauss_curv == pytest.approx(4.0 * np.pi, abs=1e-8)
    assert group.vol == pytest.approx(3749.455756, abs=1e-6)
    assert group.sa == pytest.approx(1241.017365, abs=1e-6)


def test_single_cell_group_mean_curvature_is_not_halved(multicell_group):
    system = multicell_group.sys
    group = Group(
        sys=system, name="single_cell_test", atoms=[0],
        make_net=True, print_metrics=False,
    )
    group.build()
    group.get_info()
    assert group.int_mean_curv == pytest.approx(
        float(group.net.balls.iloc[0]["int_mean_curv"]), abs=1e-8
    )
