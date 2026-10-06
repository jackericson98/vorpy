from pathlib import Path
from types import SimpleNamespace
import csv

import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.power_interface_cli import (
    aw_comparison_structure_without_waters,
    chain_groups_from_legacy_args,
    extract_power_interface_options,
    power_output_directory,
    run_power_interface_pdb,
)
from vorpy.src.analyze.power_vs_aw import compare_power_and_aw
from vorpy.src.analyze.aw_interface_curvature import (
    AWInterfaceCurvatureResult,
    InterfaceEdgeRecord,
    InterfaceSurfaceRecord,
)


def test_power_cli_options_and_existing_chain_group_syntax(tmp_path):
    args, settings = extract_power_interface_options([
        "2KAI.pdb", "--power-interface", "cazals", "--power-alpha", "0",
        "--power-probe-radius", "1.4", "--power-M", "5",
        "-g", "c", "A", "and", "c", "B", "-g", "c", "I",
    ])
    assert args[1:3] == ["-g", "c"]
    assert settings["preset"] == "cazals"
    assert settings["probe_radius"] == 1.4
    assert settings["alpha"] == 0
    assert settings["condition_beta_m"] == 5
    assert chain_groups_from_legacy_args(args[1:]) == ({"A", "B"}, {"I"})
    assert power_output_directory(["-e", "dir", str(tmp_path), "-g", "c", "A"], "2KAI.pdb") == tmp_path / "2KAI"


def test_power_option_errors_are_actionable():
    with pytest.raises(SystemExit, match="require --power-interface cazals"):
        extract_power_interface_options(["2KAI.pdb", "--power-alpha", "0"])


def test_power_runner_writes_standard_summary_and_csv_exports(tmp_path, monkeypatch):
    pdb = tmp_path / "mini.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    atom = SimpleNamespace(index=0, pdb_serial=1, group="A", chain="A",
                           residue_name="ALA", residue_number=1, insertion_code="",
                           atom_name="CA", element="C", x=0., y=0., z=0.,
                           base_radius=1.7, expanded_radius=3.1, radius_class="C",
                           group_unused=None)
    prepared = SimpleNamespace(points=np.array([[0., 0., 0.], [1., 0., 0.]]),
                               expanded_radii=np.array([3.1, 3.1]),
                               atoms=[atom], group_a={0}, group_b={1},
                               fallback_atoms=[], probe_radius=1.4)
    totals = {"length": 0., "positive_edges": 0, "negative_edges": 0,
              "signed": 0., "unsigned": 0., "signed_mean_rad": 0.,
              "unsigned_mean_rad": 0., "signed_mean_deg": 0., "unsigned_mean_deg": 0.}
    result = SimpleNamespace(final_facets=[], group_a={0}, group_b={1}, edges=[],
                             all_bicolor_facets=[], included_edges=[], edge_totals=totals,
                             alpha=0., alpha_tolerance=1e-12, condition_beta_m=5.)
    monkeypatch.setattr("vorpy.src.analyze.power_interface_cli.prepare_cazals_pdb",
                        lambda *a, **k: prepared)
    monkeypatch.setattr("vorpy.src.analyze.power_interface_cli.analyze_power_interface_curvature",
                        lambda *a, **k: result)
    run = run_power_interface_pdb(pdb, {"A"}, {"I"}, tmp_path / "exports")
    expected = {"power_interface_summary.txt", "power_interface_atoms.csv",
                "power_interface_facets.csv", "power_interface_edges.csv",
                "power_interface_components.csv"}
    assert expected <= {path.name for path in (tmp_path / "exports").iterdir()}
    assert "Power/Laguerre Interface Analysis" in run["summary"]
    with (tmp_path / "exports" / "power_interface_atoms.csv").open(newline="", encoding="utf-8") as handle:
        assert list(csv.DictReader(handle)) == []


def test_aw_compare_water_filter_is_temporary_and_preserves_source(tmp_path):
    original = ("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
                "HETATM    2  O   HOH A   2       1.000   0.000   0.000  1.00  0.00           O\n")
    pdb = tmp_path / "input.pdb"
    pdb.write_text(original, encoding="utf-8")
    filtered, count = aw_comparison_structure_without_waters(pdb)
    try:
        assert count == 1
        assert "HOH" not in filtered.read_text(encoding="utf-8")
        assert pdb.read_text(encoding="utf-8") == original
    finally:
        filtered.unlink(missing_ok=True)


def test_power_aw_comparison_passes_through_standalone_aw_curvature(tmp_path):
    power = SimpleNamespace(
        final_facets=[], included_edges=[], group_a=frozenset({0}),
        group_b=frozenset({1}), alpha=0., condition_beta_m=5.,
        edge_totals={"signed": 0., "unsigned": 0., "positive_edges": 0,
                     "negative_edges": 0, "length": 0.},
    )
    prepared = SimpleNamespace(points=np.array([[0., 0., 0.], [1., 0., 0.]]),
                               atoms=[], group_a={0}, group_b={1}, probe_radius=1.4)
    power_run = {"result": power, "prepared": prepared,
                 "facet_records": [], "facet_area_by_facet": {}}
    network = SimpleNamespace(balls=pd.DataFrame(columns=["loc"]))
    surface = InterfaceSurfaceRecord(
        surface_id=3, generator_ids=(0, 1), atom_numbers=(1, 2),
        group_orientation="A->B", area=4., bounded=True,
        feature_complete=True, generator_cells_complete=True, supported=True,
        included=True, dual_edge=(0, 1), integrated_mean_curvature=2.5,
        status="included",
    )
    edge = InterfaceEdgeRecord(
        edge_id=9, generator_ids=(0, 1, 2), atom_numbers=(1, 2, 3),
        dual_triangle=(0, 1, 2), surface_ids=(3, 4), geometry_type="analytic",
        supported=True, feature_complete=True, generator_cells_complete=True,
        length=3., signed_contribution=-4., unsigned_contribution=7.,
        positive_contribution=0., negative_contribution=-4., beta_min=-2.,
        beta_max=2., beta_mean=0., beta_statistics_method="test",
        integration_status="converged", estimated_error=0., quadrature_order=16,
        status="included",
    )
    aw = AWInterfaceCurvatureResult(
        network=network, complex=None, group_a=frozenset({0}),
        group_b=frozenset({1, 2}), selection_mode="supported",
        surfaces=[surface], edges=[edge],
    )
    result = compare_power_and_aw(power_run, aw, tmp_path,
                                  group_labels=({"A", "B"}, {"I"}))
    report = (tmp_path / "power_vs_aw_summary.txt").read_text(encoding="utf-8")
    assert "smooth integral(H dA)=2.5 A" in report
    assert "raw signed edge integral(integral beta ds)=-4 A rad" in report
    assert "raw unsigned edge integral(integral |beta| ds)=7 A rad" in report
    assert "conventional edge contribution(1/2 integral beta ds)=-2 A" in report
    assert "combined conventional AW integrated mean curvature C_H=0.5 A" in report
    assert "interface normal A->B" in report
    assert "reversing A/B reverses signed smooth, edge, and combined quantities" in report
    assert "group A chains=['A', 'B']; group B chains=['I']" in report
    metrics = {row["metric"]: row for row in result["metrics"] if row["model"] == "AW"}
    assert metrics["C_surface_H"]["value"] == aw.surface_curvature
    assert metrics["C_edge_raw_signed"]["value"] == aw.edge_curvature["signed"]
    assert metrics["C_edge_raw_unsigned"]["value"] == aw.edge_curvature["unsigned"]
    assert metrics["C_edge_H_signed"]["value"] == aw.edge_curvature["signed"] / 2
    assert metrics["C_total_H"]["value"] == aw.combined_curvature
    assert metrics["orientation_convention"]["value"].startswith("interface normal A->B")
    assert metrics["C_total_H"]["value"] == metrics["C_surface_H"]["value"] + metrics["C_edge_H_signed"]["value"]
    assert (tmp_path / "power_vs_aw_metrics.csv").is_file()


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("unequal_pair", [False, True])
def test_comparison_matches_analyzed_aw_and_exported_totals(tmp_path, reverse, unequal_pair):
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    from vorpy.tests.analyze.test_aw_interface_curvature import _interface_network

    groups = ({1, 2}, {0}) if reverse else ({0}, {1, 2})
    aw = analyze_aw_interface_curvature(_interface_network(unequal_pair=unequal_pair), *groups)
    aw.export_csv(tmp_path / "standalone")
    power = SimpleNamespace(
        final_facets=[], included_edges=[], group_a={0}, group_b={1},
        alpha=0., condition_beta_m=5., edge_totals={"signed": 0.},
    )
    prepared = SimpleNamespace(points=np.empty((0, 3)), atoms=[],
                               group_a={0}, group_b={1}, probe_radius=1.4)
    report = compare_power_and_aw({"result": power, "prepared": prepared}, aw, tmp_path)
    with (tmp_path / "power_vs_aw_metrics.csv").open(newline="", encoding="utf-8") as handle:
        metrics = {r["metric"]: r for r in csv.DictReader(handle) if r["model"] == "AW"}
    expected = {
        "C_surface_H": aw.surface_curvature,
        "C_edge_raw_signed": aw.edge_curvature["signed"],
        "C_edge_raw_unsigned": aw.edge_curvature["unsigned"],
        "C_edge_H_signed": 0.5 * aw.edge_curvature["signed"],
        "C_total_H": aw.combined_curvature,
    }
    for metric, value in expected.items():
        assert float(metrics[metric]["value"]) == pytest.approx(value)
        assert f"{value:.9g}" in report["summary"]
    assert metrics["C_edge_raw_signed"]["unit"] == "A rad"
    assert metrics["C_edge_H_signed"]["unit"] == "A"
    assert "reversing A/B reverses signed smooth, edge, and combined" in metrics["orientation_convention"]["value"]
    assert "withheld" not in report["summary"].lower()
    assert "pending orientation" not in report["summary"].lower()
    standalone = (tmp_path / "standalone" / "aw_interface_summary.txt").read_text(encoding="utf-8")
    assert f"C_H = C_smooth + 1/2 C_edge,signed: {aw.combined_curvature:.9g} A" in standalone


def test_comparison_consumes_combined_analyzer_property(tmp_path, monkeypatch):
    # A sentinel catches accidental reimplementation of the combined formula.
    monkeypatch.setattr(AWInterfaceCurvatureResult, "combined_curvature", property(lambda self: 123.))
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    from vorpy.tests.analyze.test_aw_interface_curvature import _interface_network
    aw = analyze_aw_interface_curvature(_interface_network(), {0}, {1, 2})
    power = SimpleNamespace(final_facets=[], included_edges=[], group_a={0}, group_b={1},
                            alpha=0., condition_beta_m=5., edge_totals={"signed": 0.})
    prepared = SimpleNamespace(points=np.empty((0, 3)), atoms=[],
                               group_a={0}, group_b={1}, probe_radius=1.4)
    report = compare_power_and_aw({"result": power, "prepared": prepared}, aw, tmp_path)
    assert next(r["value"] for r in report["metrics"]
                if r["model"] == "AW" and r["metric"] == "C_total_H") == 123.
    assert "C_H=123 A" in report["summary"]
