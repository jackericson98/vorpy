from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature


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


def test_apollonius_triangle_requires_actual_three_surface_edge_incidence():
    network = _interface_network(missing_third_face=True)
    result = analyze_aw_interface_curvature(network, {0}, {1, 2}, selection_mode="supported")

    assert len(result.edges) == 1
    assert result.edges[0].status == "excluded"
    assert not result.edges[0].supported
    assert result.edges[0].reason
    assert not result.selected_edges
    assert len(result.selected_surfaces) == 2


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
    assert result.exclusion_diagnostics["quadrature_not_converged"] == 1


def test_orientation_keeps_signed_and_absolute_quantities_separate():
    forward = analyze_aw_interface_curvature(_interface_network(), {0}, {1, 2})
    reverse = analyze_aw_interface_curvature(_interface_network(), {1, 2}, {0})
    assert forward.edge_curvature["positive"] > 0.0
    assert forward.edge_curvature["negative"] == 0.0
    assert forward.edge_curvature["signed"] == pytest.approx(-reverse.edge_curvature["signed"])
    assert forward.edge_curvature["unsigned"] == pytest.approx(reverse.edge_curvature["unsigned"])
    assert forward.combined_curvature == pytest.approx(
        forward.surface_curvature + 0.5*forward.edge_curvature["signed"]
    )
    assert reverse.combined_curvature == pytest.approx(-forward.combined_curvature)
    assert forward.selected_edges[0].estimated_error < 1e-12


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
