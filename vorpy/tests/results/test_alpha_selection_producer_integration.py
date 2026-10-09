"""Cached InterfaceGeometryAnalysis alpha selection -> typed Results regression."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.output import visualization
from vorpy.src.results import (
    AlphaSelectionResult,
    INTERFACE_LOG_COLUMNS,
    InterfaceResult,
    Provenance,
    ResultIdentity,
    Status,
)


FIXTURES = Path(__file__).parent / "fixtures"


def _identity(representation):
    return ResultIdentity(
        result_kind="interface",
        system="2KAI",
        system_id="2KAI:fixture",
        frame="static",
        network_id="2KAI:aw:fixture",
        network_scope="system",
        environment="dry",
        partition="aw",
        representation=representation,
        radius_configuration_id="2KAI:fixture",
        interface_id="2KAI_matched_probe",
        group_A="A",
        group_B="B",
        side="shared",
        orientation="group1 -> group2",
        alpha_method="additive",
        alpha_value=1.4,
        alpha_units="A",
        probe_radius=1.4,
        probe_units="Å",
    )


def _cached_analysis():
    state = json.loads((FIXTURES / "2kai_aw_cached_analysis.json").read_text())
    analysis = InterfaceGeometryAnalysis(**state)
    analysis.archive_provenance = {"fixture": "2kai_aw_cached_analysis"}
    return analysis


def test_cached_alpha_selection_becomes_typed_result_without_a_geometry_query(monkeypatch):
    analysis = _cached_analysis()
    iface = SimpleNamespace(geometry_analysis=analysis)
    monkeypatch.setattr(
        visualization,
        "_cached_result_identity",
        lambda _iface, **_kwargs: _identity("alpha_selection"),
    )

    result = visualization._cached_alpha_selection_result(iface)

    assert isinstance(result, AlphaSelectionResult)
    assert result.identity.representation == "alpha_selection"
    assert result.identity.alpha_method == "additive"
    assert result.identity.alpha_value == 1.4
    assert result.identity.alpha_units == "A"
    assert result.status is Status.PARTIAL
    assert result.selection["selected_features"].value == 339
    assert result.selection["coverage"].value == 1.0
    assert result.selection["eligible_pairs"].value is None
    assert result.selection["eligible_pairs"].status is Status.NOT_CALCULATED
    assert result.provenance.details["fixture"] == "2kai_aw_cached_analysis"
    assert result.selected_generator_pairs == ()


def test_cached_alpha_result_uses_existing_eight_column_log_rows(monkeypatch):
    analysis = _cached_analysis()
    alpha_result = visualization.adapt_alpha_selection(
        analysis,
        _identity("alpha_selection"),
        provenance=Provenance("fixture alpha cache"),
    )
    physical_result = InterfaceResult(
        _identity("voronoi"),
        Provenance("fixture physical cache"),
        Status.NOT_CALCULATED,
    )
    iface = SimpleNamespace(
        geometry_analysis=analysis,
        water_topology=None,
        representation_caches={},
    )
    monkeypatch.setattr(visualization, "_cached_interface_result", lambda _iface: physical_result)
    monkeypatch.setattr(
        visualization,
        "_cached_result_identity",
        lambda _iface, representation="voronoi", **_kwargs: _identity(representation),
    )
    monkeypatch.setattr(
        visualization, "_cached_alpha_selection_result", lambda _iface: alpha_result
    )

    rows = visualization._compact_result_rows(iface)
    logged = {(row["representation"], row["quantity"]): row for row in rows}

    selected = logged[("alpha_selection", "alpha_selected_count")]
    coverage = logged[("alpha_selection", "alpha_coverage")]
    assert selected["value"] == "339"
    assert selected["status"] == "PARTIAL"
    assert coverage["value"] == "1.0"
    assert all(tuple(row) == INTERFACE_LOG_COLUMNS for row in rows)


def test_cached_alpha_selection_writes_canonical_log_without_rebuilding(monkeypatch, tmp_path):
    """The normal log writer consumes the frozen cache, not geometry builders."""
    analysis = _cached_analysis()
    physical_result = InterfaceResult(
        _identity("voronoi"), Provenance("fixture physical cache"), Status.NOT_CALCULATED,
    )
    iface = SimpleNamespace(
        geometry_analysis=analysis,
        water_topology=None,
        representation_caches={},
        analyze_geometry=lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("canonical log rebuilt geometry")
        ),
    )
    monkeypatch.setattr(
        visualization,
        "_cached_result_identity",
        lambda _iface, representation="voronoi", **_kwargs: _identity(representation),
    )
    monkeypatch.setattr(visualization, "_cached_interface_result", lambda _iface: physical_result)

    path = tmp_path / "logs.csv"
    visualization._compact_log(iface, path)
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert tuple(rows[0]) == INTERFACE_LOG_COLUMNS
    logged = {(row["representation"], row["quantity"]): row for row in rows}
    selected = logged[("alpha_selection", "alpha_selected_count")]
    coverage = logged[("alpha_selection", "alpha_coverage")]
    unavailable = logged[("alpha_selection", "alpha_eligible_count")]
    alpha_value = logged[("alpha_selection", "alpha_value")]
    assert (selected["value"], selected["status"]) == ("339", "PARTIAL")
    assert (coverage["value"], coverage["status"]) == ("1.0", "PARTIAL")
    assert (unavailable["value"], unavailable["status"]) == ("", "NOT_CALCULATED")
    assert (alpha_value["value"], alpha_value["units"]) == ("1.4", "A")
    assert "2kai_aw_cached_analysis" in alpha_value["provenance"]
