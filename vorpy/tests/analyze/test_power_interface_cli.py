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


def test_power_aw_comparison_writes_separate_convention_report(tmp_path):
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
    aw = SimpleNamespace(selected_surfaces=[], selected_edges=[], network=network,
                         interface_area=0., edge_curvature={"signed": 0.},
                         surface_curvature=0., selection_mode="supported",
                         group_a={0}, group_b={1})
    result = compare_power_and_aw(power_run, aw, tmp_path,
                                  group_labels=({"A", "B"}, {"I"}))
    report = (tmp_path / "power_vs_aw_summary.txt").read_text(encoding="utf-8")
    assert "C_total_H is withheld" in report
    assert "group A chains=['A', 'B']; group B chains=['I']" in report
    assert not any(row["metric"] == "C_total_H" and row["value"] is not None
                   for row in result["metrics"])
    assert "AW curvature: unavailable" in report
    assert next(row for row in result["metrics"] if row["metric"] == "C_surface_H")["value"] is None
    assert (tmp_path / "power_vs_aw_metrics.csv").is_file()
