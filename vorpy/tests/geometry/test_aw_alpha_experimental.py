from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.apollonius import ApolloniusComplex
from vorpy.src.geometry.aw_alpha.births import (
    edge_birth,
    generator_birth,
    surface_birth,
    vertex_birth,
)
from vorpy.src.geometry.aw_alpha.cli import extract_aw_alpha_options
from vorpy.src.geometry.aw_alpha.complex import AWAlphaFiltration


def _network(locations, radii, *, surfaces=(), edges=(), vertices=()):
    return SimpleNamespace(
        settings={"net_type": "aw"},
        group_name="synthetic-aw-alpha",
        balls=pd.DataFrame(
            {
                "loc": list(locations),
                "rad": list(radii),
                "complete": [True] * len(radii),
                "name": [f"A{i}" for i in range(len(radii))],
                "chain_name": ["A"] * len(radii),
                "res_name": ["GLY"] * len(radii),
                "res_seq": list(range(1, len(radii) + 1)),
                "pdb_ins_code": [""] * len(radii),
            }
        ),
        surfs=pd.DataFrame(surfaces, columns=["balls", "edges", "verts"]),
        edges=pd.DataFrame(edges, columns=["balls", "verts", "points"]),
        verts=pd.DataFrame(vertices, columns=["balls", "loc"]),
    )


def _pair_network(distance, radii):
    net = _network(
        [(0.0, 0.0, 0.0), (distance, 0.0, 0.0)], radii, surfaces=[([0, 1], [], [])]
    )
    simplex = ApolloniusComplex(net).simplices[1][(0, 1)]
    return net, simplex


@pytest.mark.parametrize(
    "distance,radii",
    [
        (8.0, (2.0, 2.0)),
        (9.0, (2.0, 4.0)),
        (5.0, (2.0, 2.0)),
    ],
)
def test_analytic_pair_surface_birth_for_equal_unequal_and_overlapping(distance, radii):
    net, simplex = _pair_network(distance, radii)
    value, diagnostic = surface_birth(net, simplex)
    assert value == pytest.approx((distance - sum(radii)) / 2.0, abs=1e-12)
    assert diagnostic["birth_source"] == "interior_analytic"
    assert diagnostic["residual_A"] <= 1e-12


def test_single_generator_birth_is_formal_negative_radius():
    assert generator_birth(2.75) == -2.75


def test_aw_edge_birth_uses_existing_analytic_edge_solver():
    from vorpy.tests.network.test_build_net import _analytic_diagnostic_network

    net = _analytic_diagnostic_network()
    value, diagnostic = edge_birth(net, 0)
    assert value == pytest.approx(2.0, abs=1e-10)
    assert diagnostic["birth_source"] in {"endpoint", "stationary_analytic"}


def test_four_generator_vertex_birth_verifies_equal_clearances():
    locations = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    radius = 1.0
    birth = np.sqrt(3.0) - radius
    net = _network(locations, [radius] * 4, vertices=[([0, 1, 2, 3], (0.0, 0.0, 0.0))])
    value, diag = vertex_birth(net, 0, (0, 1, 2, 3))
    assert value == pytest.approx(birth)
    assert diag["residual_A"] < 1e-12
    assert len(diag["individual_clearances_A"]) == 4


def test_pair_birth_rigid_transform_radius_shift_and_permutation_invariance():
    base, simplex = _pair_network(8.0, (2.0, 4.0))
    expected, _ = surface_birth(base, simplex)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([11.0, -5.0, 2.0])
    points = np.asarray(base.balls["loc"].tolist()) @ rotation.T + shift
    moved, moved_simplex = _pair_network(
        np.linalg.norm(points[1] - points[0]), (2.0, 4.0)
    )
    moved.balls["loc"] = list(points)
    value, _ = surface_birth(moved, moved_simplex)
    assert value == pytest.approx(expected, abs=1e-12)

    shifted, shifted_simplex = _pair_network(8.0, (3.0, 5.0))
    shifted_value, _ = surface_birth(shifted, shifted_simplex)
    assert shifted_value == pytest.approx(expected - 1.0, abs=1e-12)

    permuted = _network(
        [(8.0, 0.0, 0.0), (0.0, 0.0, 0.0)], (4.0, 2.0), surfaces=[([0, 1], [], [])]
    )
    permuted_value, _ = surface_birth(
        permuted, ApolloniusComplex(permuted).simplices[1][(0, 1)]
    )
    assert permuted_value == pytest.approx(expected, abs=1e-12)


def test_infeasible_pair_contact_without_known_boundary_is_unresolved():
    net = _network(
        [(0.0, 0.0, 0.0), (8.0, 0.0, 0.0), (4.0, 0.0, 0.0)],
        (1.0, 1.0, 10.0),
        surfaces=[([0, 1], [], [])],
    )
    value, diagnostic = surface_birth(net, ApolloniusComplex(net).simplices[1][(0, 1)])
    assert value is None
    assert diagnostic["birth_source"] == "unresolved"


def test_filtration_detects_geometric_order_violation_without_hiding_it():
    net, _ = _pair_network(8.0, (2.0, 2.0))
    filtration = AWAlphaFiltration(net)
    filtration._record(0, (0,), -2.0, supported=True, source="test")
    filtration._record(0, (1,), -2.0, supported=True, source="test")
    filtration._record(1, (0, 1), -3.0, supported=True, source="test")
    filtration._close_filtration()
    filtration._built = True
    assert filtration.validate_filtration()[0]["kind"] == "monotonicity_violation"
    assert filtration.records[1][(0, 1)].filtration_birth == -2.0
    assert filtration.records[1][(0, 1)].closure_adjustment == 1.0


def test_cli_option_is_separate_and_retains_multiple_alpha_values():
    remaining, alphas = extract_aw_alpha_options(
        [
            "molecule.pdb",
            "--alpha-complex",
            "-2",
            "--aw-alpha-experimental",
            "0",
            "-s",
            "mv",
            "5",
            "--aw-alpha-experimental",
            "1.4",
        ]
    )
    assert remaining == ["molecule.pdb", "--alpha-complex", "-2", "-s", "mv", "5"]
    assert alphas == [0.0, 1.4]
    from vorpy.src.geometry.aw_alpha.cli import extract_aw_alpha_pair_only

    without_pair_flag, pair_only = extract_aw_alpha_pair_only(
        ["-g", "c", "A", "--aw-alpha-pairs-only-experimental"]
    )
    assert without_pair_flag == ["-g", "c", "A"]
    assert pair_only


def test_pair_overlap_can_build_only_generator_and_surface_births():
    net, _ = _pair_network(8.0, (2.0, 2.0))
    filtration = AWAlphaFiltration(net).calculate_births(max_dimension=1)
    assert filtration._max_dimension_built == 1
    assert len(filtration.records[0]) == 2
    assert len(filtration.records[1]) == 1
    assert not filtration.records[2]
    assert not filtration.records[3]


def test_truncated_filtration_rejects_negative_alpha_query():
    net, _ = _pair_network(8.0, (2.0, 2.0))
    filtration = AWAlphaFiltration(net).calculate_births(max_dimension=1)
    with pytest.raises(ValueError, match="truncated to alpha>=0"):
        filtration.alpha_complex(-0.1, truncate_negative=True, max_dimension=1)


def test_bicolor_export_uses_parent_system_stable_atom_metadata():
    from vorpy.src.geometry.aw_alpha.complex import AWAlphaFiltration

    net, _ = _pair_network(2.0, (2.0, 2.0))
    net.surfs.at[0, "edges"] = [0]
    net.surfs.at[0, "verts"] = [0, 1]
    net.verts = pd.DataFrame(
        {
            "balls": [[0, 1], [0, 1]],
            "loc": [(1.0, 1.0, 0.0), (1.0, -1.0, 0.0)],
        }
    )
    net.sys = SimpleNamespace(
        balls=pd.DataFrame(
            {
                "loc": [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
                "rad": [2.0, 2.0],
                "num": [101, 202],
                "chain_name": ["A", "I"],
                "res_seq": [15, 3],
                "pdb_ins_code": ["", "B"],
                "res_name": ["LYS", "TYR"],
                "name": ["NZ", "CA"],
            }
        )
    )
    filtration = AWAlphaFiltration(net).calculate_births(max_dimension=1)
    rows, _ = filtration.bicolor_surface_records(0.0, {0}, {1}, max_dimension=1)
    assert len(rows) == 1
    assert rows[0]["atom_i_stable_id"] == "A|15|LYS|NZ"
    assert rows[0]["atom_i_metadata_status"] == "available"
    assert rows[0]["atom_i_atom_index"] == 101
    assert rows[0]["atom_j_stable_id"] == "I|3B|TYR|CA"
    assert rows[0]["atom_i_radius_A"] == 2.0


def test_parser_safe_water_chain_normalization_preserves_atoms_and_coordinates(
    tmp_path,
):
    from vorpy.src.geometry.aw_alpha.cli import parser_safe_aw_alpha_structure

    source = tmp_path / "water.pdb"
    protein = "ATOM      1  CA  GLY A   1       1.000   2.000   3.000  1.00 10.00           C  \n"
    water = "HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00 10.00           O  \n"
    source.write_text(protein + water, encoding="utf-8")
    temporary, count = parser_safe_aw_alpha_structure(source)
    try:
        original_lines = source.read_text(encoding="utf-8").splitlines()
        transformed_lines = temporary.read_text(encoding="utf-8").splitlines()
        assert count == 1
        assert len(original_lines) == len(transformed_lines) == 2
        assert transformed_lines[0] == original_lines[0]
        assert transformed_lines[1][21] == " "
        assert transformed_lines[1][30:54] == original_lines[1][30:54]
        assert water in source.read_text(encoding="utf-8")
    finally:
        temporary.unlink()
