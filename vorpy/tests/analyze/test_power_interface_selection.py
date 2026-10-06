import numpy as np
import pytest

from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
from vorpy.src.analyze.power_interface_selection import power_selection_views


def test_full_power_uses_bounded_geometry_without_alpha_or_m_filters():
    pytest.importorskip("gudhi")
    points = np.random.default_rng(43).uniform(-3., 3., (25, 3))
    a = set(np.flatnonzero(points[:, 0] < 0))
    b = set(range(len(points))) - a
    frozen = analyze_power_interface_curvature(points, np.full(len(points), .4), a, b,
                                               condition_beta_m=.01)
    original = [(f.final_selected, f.status, f.reason) for f in frozen.facets]
    views, polygons = power_selection_views(frozen)
    from vorpy.src.analyze.cazals_validation import power_facet_polygon_area
    assert polygons == {f.generator_ids: power_facet_polygon_area(frozen, f) for f in frozen.facets}
    full = views["full_power"]
    assert len(full.final_facets) > 0
    assert {f.generator_ids for f in full.final_facets} == {
        ids for ids, polygon in polygons.items() if polygon["status"] == "finite_bounded"}
    assert any(not f.alpha_zero_selected for f in full.final_facets)
    assert len(full.final_facets) > len(frozen.final_facets)
    assert views["cazals_alpha0_power"] is frozen
    assert original == [(f.final_selected, f.status, f.reason) for f in frozen.facets]
    selected = {f.generator_ids for f in full.final_facets}
    assert {e.regular_triangle for e in full.included_edges} == {
        e.regular_triangle for e in frozen.edges
        if e.geometry_status == "finite" and all(f in selected for f in e.bicolor_facets)}
    originals = {e.regular_triangle: e for e in frozen.edges}
    for edge in full.included_edges:
        assert edge.length_beta == originals[edge.regular_triangle].length_beta
        assert edge.length_abs_beta == originals[edge.regular_triangle].length_abs_beta
    alpha = views["alpha0_power_before_M"]
    assert {f.generator_ids for f in alpha.final_facets} == {
        f.generator_ids for f in full.final_facets if f.alpha_zero_selected}


def test_full_power_reversal_preserves_geometry_and_reverses_signed_measure():
    pytest.importorskip("gudhi")
    points = np.random.default_rng(8).normal(size=(20, 3))
    a, b = set(range(10)), set(range(10, 20))
    forward, fp = power_selection_views(analyze_power_interface_curvature(points, np.ones(20), a, b))
    reverse, rp = power_selection_views(analyze_power_interface_curvature(points, np.ones(20), b, a))
    for method in forward:
        f, r = forward[method], reverse[method]
        assert {s.generator_ids for s in f.final_facets} == {s.generator_ids for s in r.final_facets}
        assert f.edge_totals["signed"] == pytest.approx(-r.edge_totals["signed"])
        assert f.edge_totals["unsigned"] == pytest.approx(r.edge_totals["unsigned"])
        assert sum(fp[s.generator_ids]["area_A2"] for s in f.final_facets) == pytest.approx(
            sum(rp[s.generator_ids]["area_A2"] for s in r.final_facets))
