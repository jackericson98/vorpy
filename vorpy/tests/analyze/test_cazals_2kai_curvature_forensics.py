from __future__ import annotations

from pathlib import Path

import pytest

from vorpy.src.analyze.power_interface_cli import run_power_interface_pdb
from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_2kai_curvature_forensics import (
    EXPECTED,
    export_2kai_forensics,
)
from vorpy.src.analyze.power_interface_curvature import cazals_signed_beta


@pytest.fixture(scope="module")
def forensic_run(tmp_path_factory):
    pdb = Path(__file__).resolve().parents[2] / "data" / "2KAI.pdb"
    out = tmp_path_factory.mktemp("2kai_cazals_forensic")
    run = run_power_interface_pdb(pdb, {"A", "B"}, {"I"}, out / "production_source",
                                  probe_radius=1.4, alpha=0.0, condition_beta_m=5.0)
    result, prepared = run["result"], run["prepared"]
    # Snapshot result state that the exporter must only read.
    snapshot = {
        "points": result.points.copy(), "radii": result.expanded_radii.copy(),
        "groups": (frozenset(result.group_a), frozenset(result.group_b)),
        "facets": tuple((f.generator_ids, f.alpha_zero_selected, f.condition_beta_accepted,
                         f.final_selected, f.m_over_r) for f in result.facets),
        "edges": tuple((e.regular_triangle, e.length, e.beta_radians, e.length_beta,
                        e.length_abs_beta, e.status) for e in result.edges),
    }
    exported = export_2kai_forensics(result, prepared, out / "forensics")
    yield result, prepared, exported, snapshot


def test_forensic_export_reconstructs_frozen_2kai_totals(forensic_run):
    _, _, exported, _ = forensic_run
    total = exported["total"]
    assert total["edge_count"] == EXPECTED["edge_count"] == 864
    assert total["total_length_A"] == pytest.approx(EXPECTED["total_length_A"], abs=2e-6)
    assert total["signed_beta_l_A_rad"] == pytest.approx(EXPECTED["signed_beta_l_A_rad"], abs=2e-6)
    assert total["unsigned_beta_l_A_rad"] == pytest.approx(EXPECTED["unsigned_beta_l_A_rad"], abs=2e-6)
    assert total["signed_sH_deg"] == pytest.approx(EXPECTED["signed_sH_deg"], abs=2e-6)
    assert exported["m5_counts"]["alpha0_bicolor_facets"] == 362
    assert exported["m5_counts"]["real_orthogonal_ball_candidate_q_ge_0"] == 88
    assert exported["m5_counts"]["q_lt_0_no_real_candidate"] == 274
    assert exported["m5_counts"]["M5_rejected"] == 0
    assert {"facet_1_m_if_real_A", "facet_1_m_over_r_if_real", "facet_1_q_classification"} <= set(exported["edges"][0])


def test_edge_lengths_beta_and_pattern_sums_match_production(forensic_run):
    result, _, exported, _ = forensic_run
    rows = exported["edges"]
    assert len(rows) == len(result.included_edges) == 864
    for row, edge in zip(rows, result.included_edges):
        assert row["edge_length_A"] == pytest.approx(edge.length, abs=0.0)
        assert row["beta_abs_rad"] == pytest.approx(abs(edge.beta_radians), abs=0.0)
        assert row["beta_signed_rad"] == pytest.approx(edge.beta_radians, abs=0.0)
        assert row["signed_beta_times_length_A_rad"] == pytest.approx(edge.length_beta, abs=0.0)
    patterns = {row["pattern"]: row for row in exported["patterns"]}
    assert patterns["AAB"]["signed_beta_l_A_rad"] + patterns["ABB"]["signed_beta_l_A_rad"] == pytest.approx(
        patterns["TOTAL"]["signed_beta_l_A_rad"], abs=1e-10)
    total = patterns["TOTAL"]
    assert total["positive_beta_l_A_rad"] + total["negative_beta_l_A_rad"] == pytest.approx(
        total["signed_beta_l_A_rad"], abs=1e-10)
    assert total["unsigned_beta_l_A_rad"] == pytest.approx(
        sum(abs(edge.beta_radians) * edge.length for edge in result.included_edges), abs=1e-10)


def test_group_reversal_flips_signed_curvature_only(forensic_run):
    result, _, exported, _ = forensic_run
    reversed_betas = [cazals_signed_beta(edge.regular_triangle, result.points,
                                         result.group_b, result.group_a)[2]
                      for edge in result.included_edges]
    assert sum(float(edge.length) * beta for edge, beta in zip(result.included_edges, reversed_betas)) == pytest.approx(
        -exported["total"]["signed_beta_l_A_rad"], abs=1e-10)
    assert [abs(beta) for beta in reversed_betas] == pytest.approx(
        [abs(edge.beta_radians) for edge in result.included_edges], abs=1e-12)
    assert [edge.length for edge in result.included_edges] == pytest.approx(
        [row["edge_length_A"] for row in exported["edges"]], abs=0.0)


def test_diagnostic_export_does_not_mutate_production_result(forensic_run):
    result, _, _, snapshot = forensic_run
    assert (result.points == snapshot["points"]).all()
    assert (result.expanded_radii == snapshot["radii"]).all()
    assert (frozenset(result.group_a), frozenset(result.group_b)) == snapshot["groups"]
    assert tuple((f.generator_ids, f.alpha_zero_selected, f.condition_beta_accepted,
                  f.final_selected, f.m_over_r) for f in result.facets) == snapshot["facets"]
    assert tuple((e.regular_triangle, e.length, e.beta_radians, e.length_beta,
                  e.length_abs_beta, e.status) for e in result.edges) == snapshot["edges"]
