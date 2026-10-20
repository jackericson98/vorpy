from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
from vorpy.src.calculations.aw_interface_orientation import aw_interface_edge_orientation
from vorpy.src.calculations.edge_geometry import (
    AdditivelyWeightedTrisectorBranch,
    AdditivelyWeightedTrisectorConic,
    match_aw_trisector_conic_to_samples,
)
from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge
from vorpy.src.geometry.interfaces.curvature import aw_curvature_summary


def _interface_network(unequal_pair=False, missing_third_face=False, incomplete_cell=False):
    radii = [1.0, 1.0, 1.0]
    if unequal_pair:
        radii[1] = 0.5
    balls = pd.DataFrame({
        "loc": [(0., 0., 0.), (2., 0., 0.), (0., 2., 0.)],
        "rad": radii,
        "num": [10, 20, 30],
        "complete": [True, not incomplete_cell, True],
    })
    surfaces = pd.DataFrame({
        "balls": [[0, 1], [0, 2], [1, 2]],
        "edges": [[0], [0], [0]],
        "verts": [[0, 1], [0, 1], [0, 1]],
        "tris": [[], [], []],
        "sa": [4., 4., 4.],
        "func": [
            [0., 0., 0., 0., 0., 0., 1., 0., 0., -1.25, 0., 0., 0., 0.],
            [0.] * 14,
            [0.] * 14,
        ],
        "com": [(1.25, 0., 0.), (0., 1., 0.), (1., 1., 0.)],
        "flat": [not unequal_pair, True, True],
        "int_mean_curv": [2.0 if unequal_pair else 0.0, 0.0, 0.0],
    })
    balls_on_edge = [0, 1] if missing_third_face else [0, 1, 2]
    edges = pd.DataFrame({
        "balls": [[0, 1, 2]],
        "verts": [[0, 1]],
        "surfs": [balls_on_edge],
        "points": [[(1., 1., -1.), (1., 1., 0.), (1., 1., 1.)]],
        "length": [2.],
    })
    vertices = pd.DataFrame({
        "balls": [[0, 1, 2], [0, 1, 2]],
        "loc": [(1., 1., -1.), (1., 1., 1.)],
    })
    return SimpleNamespace(
        settings={"net_type": "aw"}, balls=balls, surfs=surfaces,
        edges=edges, verts=vertices, group_name="synthetic-interface",
    )


def test_equal_radius_bicolor_interface_is_planar_and_edge_integral_is_analytic():
    result = analyze_aw_interface_curvature(
        _interface_network(), {0}, {1, 2},
        quadrature_order=8, max_quadrature_order=64,
    )

    assert len(result.selected_surfaces) == 2
    assert result.surface_curvature == pytest.approx(0.0)
    assert len(result.selected_edges) == 1
    edge = result.selected_edges[0]
    assert edge.integration_status == "converged"
    assert edge.length == pytest.approx(2.0, abs=1e-12)
    assert edge.signed_contribution == pytest.approx(np.pi, abs=1e-12)
    assert edge.unsigned_contribution == pytest.approx(np.pi, abs=1e-12)
    assert edge.beta_mean == pytest.approx(np.pi / 2., abs=1e-12)
    assert edge.signed_mean_curvature == pytest.approx(np.pi / 2, abs=1e-12)
    assert result.edge_mean_curvature == pytest.approx(np.pi / 2, abs=1e-12)
    assert result.edge_contribution_complete is True
    assert result.total_curvature_certified is False
    summary = aw_curvature_summary(result)
    assert summary.conventional_edge_H_A == pytest.approx(np.pi / 2, abs=1e-12)
    assert summary.edge_contribution_complete is True
    assert summary.total_curvature_certified is False
    accounting = result.integrated_mean_curvature
    assert accounting["smooth_patch"] == pytest.approx(0.0)
    assert accounting["signed_seam"] == pytest.approx(np.pi / 2, abs=1e-12)
    assert accounting["signed_seam_raw_beta_ds"] == pytest.approx(np.pi, abs=1e-12)
    assert accounting["partial"] == pytest.approx(np.pi / 2, abs=1e-12)
    assert accounting["total"] is None
    assert accounting["certified"] is False
    assert accounting["boundary"] is None
    assert accounting["corner"] is None
    assert accounting["units"] == "angstrom"
    assert result.edge_angle_statistics["signed_degrees"] == pytest.approx(90.0)
    assert result.combined_curvature == pytest.approx(np.pi / 2)


def test_unequal_radius_pair_reuses_nonzero_oriented_smooth_surface_curvature():
    network = _interface_network(unequal_pair=True)
    forward = analyze_aw_interface_curvature(network, {0}, {1, 2})
    reverse = analyze_aw_interface_curvature(network, {1, 2}, {0})

    assert forward.surface_curvature == pytest.approx(2.0)
    assert reverse.surface_curvature == pytest.approx(-2.0)
    assert forward.edge_curvature["signed"] == pytest.approx(-reverse.edge_curvature["signed"])
    assert forward.edge_curvature["unsigned"] == pytest.approx(reverse.edge_curvature["unsigned"])
    assert forward.edge_mean_curvature is None
    assert reverse.edge_mean_curvature is None
    assert forward.edge_contribution_complete is False
    assert reverse.edge_contribution_complete is False
    assert forward.integrated_mean_curvature["smooth_patch"] == pytest.approx(2.0)
    assert forward.integrated_mean_curvature["signed_seam"] is None
    assert forward.integrated_mean_curvature["partial"] is None
    assert forward.integrated_mean_curvature["unresolved"]


def test_apollonius_triangle_requires_actual_three_surface_edge_incidence():
    network = _interface_network(missing_third_face=True)
    result = analyze_aw_interface_curvature(network, {0}, {1, 2}, selection_mode="supported")

    assert len(result.edges) == 1
    assert result.edges[0].status == "excluded"
    assert not result.edges[0].supported
    assert result.edges[0].reason
    assert result.edges[0].signed_mean_curvature is None
    assert result.edge_mean_curvature is None
    assert result.edge_contribution_complete is False
    summary = aw_curvature_summary(result)
    assert summary.edge_contribution_complete is False
    assert summary.conventional_edge_H_A is None
    assert summary.combined_H_A is None
    assert not result.selected_edges
    assert len(result.selected_surfaces) == 2
    accounting = result.integrated_mean_curvature
    assert accounting["signed_seam"] is None
    assert accounting["partial"] is None
    assert accounting["total"] is None
    assert accounting["certified"] is False
    assert accounting["unresolved"]


def test_edge_without_exactly_two_retained_bicolor_surfaces_is_unresolved():
    network = _interface_network()
    network.edges.at[0, "surfs"] = [0]
    result = analyze_aw_interface_curvature(network, {0}, {1, 2})

    assert len(result.edges) == 1
    assert result.edges[0].surface_ids == (0,)
    assert result.edges[0].status == "excluded"
    assert "exactly two bicolor" in result.edges[0].reason
    assert result.edge_contribution_complete is False
    assert result.edge_mean_curvature is None
    assert result.exclusion_diagnostics["incomplete_bicolor_surface_incidence"] == 1
    accounting = result.integrated_mean_curvature
    assert accounting["signed_seam"] is None
    assert accounting["signed_seam_raw_beta_ds"] is None
    assert accounting["partial"] is None
    assert accounting["total"] is None


def test_surface_referenced_edge_missing_from_network_is_unresolved():
    network = _interface_network()
    network.edges = network.edges.iloc[0:0].copy()
    result = analyze_aw_interface_curvature(network, {0}, {1, 2})

    assert result.missing_edge_ids == (0,)
    assert result.edge_contribution_complete is False
    assert result.edge_mean_curvature is None
    assert result.exclusion_diagnostics["missing_interface_edge_record"] == 1
    accounting = result.integrated_mean_curvature
    assert accounting["signed_seam"] is None
    assert accounting["signed_seam_raw_beta_ds"] is None
    assert accounting["partial"] is None
    assert accounting["total"] is None
    assert "missing interface edge record 0" in accounting["unresolved"]


def test_alpha_selection_fails_clearly_while_pair_surface_births_are_unresolved():
    with pytest.raises(ValueError, match="pair-surface alpha births are unresolved"):
        analyze_aw_interface_curvature(_interface_network(), {0}, {1, 2}, selection_mode="alpha")


def test_incomplete_generator_cells_are_classified_and_excluded():
    result = analyze_aw_interface_curvature(
        _interface_network(incomplete_cell=True), {0}, {1, 2}, selection_mode="raw"
    )
    assert len(result.selected_surfaces) == 1
    assert not result.selected_edges
    assert result.exclusion_diagnostics["incomplete_surface_generator_cells"] == 1
    assert result.exclusion_diagnostics["incomplete_edge_generator_cells"] == 1


def test_edge_at_maximum_quadrature_order_is_unresolved_not_counted():
    result = analyze_aw_interface_curvature(
        _interface_network(), {0}, {1, 2}, quadrature_order=8, max_quadrature_order=8
    )
    assert not result.selected_edges
    assert result.edges[0].integration_status == "maximum_order_reached"
    assert result.edges[0].signed_mean_curvature is None
    assert result.edge_mean_curvature is None
    assert result.edge_contribution_complete is False
    assert result.exclusion_diagnostics["quadrature_not_converged"] == 1
    assert result.integrated_mean_curvature["signed_seam"] is None
    assert result.integrated_mean_curvature["partial"] is None


def test_degenerate_analytic_edge_is_unresolved_not_zero():
    network = _interface_network()
    network.verts.loc[1, "loc"] = network.verts.loc[0, "loc"]
    network.edges.at[0, "points"] = [(1., 1., 0.), (1., 1., 0.)]
    result = analyze_aw_interface_curvature(network, {0}, {1, 2})

    edge = result.edges[0]
    assert edge.status == "excluded"
    assert edge.integration_status == "unresolved"
    assert "distinct" in edge.reason
    assert edge.signed_contribution is None
    assert result.integrated_mean_curvature["signed_seam"] is None
    assert result.integrated_mean_curvature["total"] is None


def test_orientation_keeps_signed_and_absolute_quantities_separate():
    forward = analyze_aw_interface_curvature(_interface_network(), {0}, {1, 2})
    reverse = analyze_aw_interface_curvature(_interface_network(), {1, 2}, {0})
    assert forward.edge_curvature["positive"] > 0.0
    assert forward.edge_curvature["negative"] == 0.0
    assert forward.edge_curvature["signed"] == pytest.approx(-reverse.edge_curvature["signed"])
    assert forward.edge_curvature["unsigned"] == pytest.approx(reverse.edge_curvature["unsigned"])
    assert forward.edge_mean_curvature == pytest.approx(-reverse.edge_mean_curvature)
    assert forward.integrated_mean_curvature["signed_seam"] == pytest.approx(
        -reverse.integrated_mean_curvature["signed_seam"]
    )
    assert forward.integrated_mean_curvature["partial"] == pytest.approx(
        -reverse.integrated_mean_curvature["partial"]
    )
    assert forward.combined_curvature == pytest.approx(
        forward.surface_curvature + 0.5*forward.edge_curvature["signed"]
    )
    assert reverse.combined_curvature == pytest.approx(-forward.combined_curvature)
    assert forward.selected_edges[0].estimated_error < 1e-12


def _analytic_edge_network(source, parameters):
    locations = np.asarray(source.generator_locations, dtype=float)
    radii = np.asarray(source.generator_radii, dtype=float)
    parameters = tuple(float(value) for value in parameters)
    points = [source.point(value) for value in parameters]
    balls = pd.DataFrame({
        "loc": list(locations), "rad": radii, "num": [10, 20, 30],
        "complete": [True, True, True],
    })
    surfaces = pd.DataFrame({
        "balls": [[0, 1], [0, 2], [1, 2]],
        "edges": [[0], [0], [0]],
        "verts": [[0, 1], [0, 1], [0, 1]],
        "sa": [4., 4., 4.],
        "func": [[0.] * 14, [0.] * 14, [0.] * 14],
        "com": [(0., 0., 0.)] * 3,
        "flat": [False, False, False],
        "int_mean_curv": [0., 0., 0.],
    })
    edges = pd.DataFrame({
        "balls": [[0, 1, 2]], "verts": [[0, 1]], "surfs": [[0, 1, 2]],
        "points": [points], "length": [source.arc_length()],
    })
    vertices = pd.DataFrame({
        "balls": [[0, 1, 2], [0, 1, 2]],
        "loc": [points[0], points[-1]],
    })
    return SimpleNamespace(
        settings={"net_type": "aw"}, balls=balls, surfs=surfaces,
        edges=edges, verts=vertices, group_name="analytic-interface",
    )


def _independent_signed_integral(network, resolved, order=256):
    geometry = resolved.geometry
    nodes, weights = np.polynomial.legendre.leggauss(order)
    half = 0.5 * (geometry.t_max - geometry.t_min)
    center = 0.5 * (geometry.t_max + geometry.t_min)
    locations = {
        index: np.asarray(network.balls.loc[index, "loc"], dtype=float)
        for index in (0, 1, 2)
    }
    signed = length = 0.0
    for node, weight in zip(nodes, weights, strict=True):
        parameter = center + half * float(node)
        point = geometry.point(parameter)
        tangent = geometry.tangent(parameter)
        oriented = aw_interface_edge_orientation(
            point, [(0, 1), (0, 2)], (0, 1, 2), {0}, {1, 2}, locations,
        )
        ds_weight = half * float(weight) * np.linalg.norm(tangent)
        signed += oriented["beta"] * ds_weight
        length += ds_weight
    return signed, length


@pytest.mark.parametrize(
    ("source", "parameters", "resolver_status"),
    [
        (
            AdditivelyWeightedTrisectorBranch(
                [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
                [1., 1.5, 2.], 2., 4., branch=-1,
            ),
            (4., 2.),
            "matched_curve",
        ),
        (
            AdditivelyWeightedTrisectorConic(
                [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
                [1., 1.5, 2.], -.7, 1.1, component=1,
            ),
            (-.7, 1.1),
            "matched_nonsingular",
        ),
    ],
)
def test_regular_analytic_edge_integral_matches_independent_quadrature(
        source, parameters, resolver_status):
    network = _analytic_edge_network(source, parameters)
    resolved = resolve_aw_network_edge(network, 0)
    assert resolved.status == resolver_status
    expected_signed, expected_length = _independent_signed_integral(network, resolved)

    result = analyze_aw_interface_curvature(
        network, {0}, {1, 2}, quadrature_order=8, max_quadrature_order=256,
    )
    edge = result.edges[0]
    assert edge.status == "included"
    assert edge.integration_status == "converged"
    assert edge.length == pytest.approx(expected_length, rel=1e-11, abs=1e-11)
    assert edge.signed_contribution == pytest.approx(
        expected_signed, rel=1e-10, abs=1e-10,
    )
    assert edge.signed_mean_curvature == pytest.approx(
        expected_signed / 2., rel=1e-10, abs=1e-10,
    )
    accounting = result.integrated_mean_curvature
    assert accounting["signed_seam"] == pytest.approx(
        expected_signed / 2., rel=1e-10, abs=1e-10,
    )
    assert accounting["partial"] == accounting["signed_seam"]
    assert accounting["units"] == "angstrom"
    assert accounting["convention"] == "H=(k1+k2)/2; C_e=0.5*integral(beta ds)"
    assert accounting["certified"] is False


def test_near_turning_branch_uses_validated_nonsingular_reparameterization():
    source = AdditivelyWeightedTrisectorBranch(
        [[0., 0., 0.], [4., 0., 0.], [0., 5., 0.]],
        [1., 1.5, 2.], 1.517, 1.54, branch=1,
    )
    assert min(source._q(value) for value in (source.t_min, source.t_max)) < 3e-3
    network = _analytic_edge_network(source, (source.t_min, source.t_max))
    resolved = resolve_aw_network_edge(network, 0)
    assert resolved.status == "matched_curve"

    conic = match_aw_trisector_conic_to_samples(
        [network.balls.loc[index, "loc"] for index in (0, 1, 2)],
        [network.balls.loc[index, "rad"] for index in (0, 1, 2)],
        network.edges.loc[0, "verts"],
        [network.verts.loc[index, "loc"] for index in network.edges.loc[0, "verts"]],
        network.edges.loc[0, "points"],
    )
    assert conic.status == "matched_nonsingular"
    expected_signed, expected_length = _independent_signed_integral(
        network, SimpleNamespace(geometry=conic.geometry), order=256,
    )

    result = analyze_aw_interface_curvature(
        network, {0}, {1, 2}, quadrature_order=8, max_quadrature_order=64,
        relative_tolerance=1e-9, absolute_tolerance=1e-11,
    )
    edge = result.edges[0]
    assert edge.status == "included"
    assert edge.integration_status == "converged"
    assert edge.length == pytest.approx(expected_length, rel=1e-10, abs=1e-10)
    assert edge.signed_contribution == pytest.approx(
        expected_signed, rel=1e-10, abs=1e-10,
    )
    assert edge.estimated_error < 1e-10


def test_generator_id_and_surface_storage_order_do_not_change_physical_curvature():
    original = _interface_network()
    reference = analyze_aw_interface_curvature(original, {0}, {1, 2})

    # Permute atom/generator IDs while preserving their physical coordinates.
    mapping = {0: 2, 1: 0, 2: 1}
    reordered = _interface_network()
    reordered.balls = reordered.balls.rename(index=mapping).sort_index()
    reordered.balls.index.name = original.balls.index.name
    for table in (reordered.surfs, reordered.edges, reordered.verts):
        table["balls"] = table["balls"].map(lambda values: [mapping[int(v)] for v in values])
    # Reverse surface row order and the order in which each edge stores faces.
    reordered.surfs = reordered.surfs.iloc[::-1].copy()
    reordered.edges["surfs"] = reordered.edges["surfs"].map(lambda values: list(reversed(values)))
    permuted = analyze_aw_interface_curvature(reordered, {mapping[0]}, {mapping[1], mapping[2]})

    assert permuted.surface_curvature == pytest.approx(reference.surface_curvature)
    assert permuted.edge_curvature["signed"] == pytest.approx(reference.edge_curvature["signed"])
    assert permuted.edge_curvature["unsigned"] == pytest.approx(reference.edge_curvature["unsigned"])
    assert permuted.combined_curvature == pytest.approx(reference.combined_curvature)
