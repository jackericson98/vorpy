"""Focused regressions for cache-only Results and canonical interface logs."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.output import visualization
from vorpy.src.results import (
    AlphaSelectionResult,
    INTERFACE_LOG_COLUMNS,
    InterfaceResult,
    Provenance,
    ResultIdentity,
    Status,
    adapt_alpha_selection,
)


FIXTURE = Path(__file__).parent / "fixtures" / "2kai_aw_cached_analysis.json"


def _identity(representation, side="shared"):
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
        side=side,
        orientation="group1 -> group2",
        alpha_method="additive",
        alpha_value=1.4,
        alpha_units="A",
        probe_radius=1.4,
        probe_units="Å",
    )


def _dual(network):
    edge = SimpleNamespace(
        simplex_id="1:3,8",
        dimension=1,
        generator_ids=(3, 8),
        primal_feature_ids=(("surface", 17),),
        bounded=True,
        complete=True,
        supported=True,
        status="bounded_complete_supported",
    )
    audit = SimpleNamespace(
        valid=True,
        generator_count=2,
        feature_counts={"surface": 1},
        unmapped_feature_count=0,
    )
    return SimpleNamespace(
        network=network,
        simplices={1: {(3, 8): edge}},
        audit=lambda: audit,
    )


def _dimensioned_dual(network):
    def simplex(dimension, generator_ids, feature):
        return SimpleNamespace(
            simplex_id=f"{dimension}:" + ",".join(map(str, generator_ids)),
            dimension=dimension,
            generator_ids=generator_ids,
            primal_feature_ids=((feature, 1),),
            bounded=True,
            complete=True,
            supported=True,
            status="bounded_complete_supported",
        )

    audit = SimpleNamespace(
        valid=True,
        generator_count=3,
        feature_counts={"cell": 3},
        unmapped_feature_count=0,
    )
    return SimpleNamespace(
        network=network,
        simplices={
            0: {(3,): simplex(0, (3,), "cell")},
            1: {
                (3, 8): simplex(1, (3, 8), "surf"),
                (8, 9): simplex(1, (8, 9), "surf"),
            },
            2: {(3, 8, 9): simplex(2, (3, 8, 9), "edge")},
        },
        metadata={
            "explicit_group_membership": {
                "generator_ids": {"A": [3], "B": [8, 9]},
                "simplex_ids": {
                    "A": ["0:3"],
                    "B": ["1:8,9"],
                    "AB/contact": ["1:3,8", "2:3,8,9"],
                },
            }
        },
        archive_provenance={
            "network_id": "2KAI:aw:fixture",
            "interface_id": "2KAI_matched_probe",
            "cache_revision": "ATLAS-A5",
        },
        audit=lambda: audit,
    )


def _stub_identity(monkeypatch):
    monkeypatch.setattr(
        visualization,
        "_cached_result_identity",
        lambda _iface, representation="voronoi", side="shared", **_kwargs: _identity(
            representation, side
        ),
    )


def _mismatched_alpha_source():
    return {
        "alpha_selection": {
            "status": "partial",
            "certified": False,
            "eligible_bicolor_pairs": 2,
            "candidate_bicolor_pairs": 2,
            "selected_pairs": 2,
            "certified_rejected_pairs": 0,
            "unresolved_pairs": 1,
            "selected_physical_pairs": 1,
            "excluded_pairs": 0,
            "coverage": 0.5,
            "selected_generator_pairs": [[3, 8], [1320, 4061]],
            "selected_system_pairs": [["system:A", "system:B"], ["system:X", "system:Y"]],
            "selected_physical_system_pairs": [["system:A", "system:B"]],
        }
    }


def test_alpha_selection_preserves_missing_ids_and_fails_closed(monkeypatch, tmp_path):
    network = object()
    source = _mismatched_alpha_source()
    result = adapt_alpha_selection(
        source,
        _identity("alpha_selection"),
        dual=_dual(network),
        provenance=Provenance("synthetic alpha selection"),
    )

    assert isinstance(result, AlphaSelectionResult)
    assert result.status is Status.UNRESOLVED
    assert result.selected_dual_feature_ids == ("1:3,8", "1:1320,4061")
    assert [feature.simplex_id for feature in result.selected_dual_features] == ["1:3,8"]
    assert result.provenance.details["missing_dual_simplex_ids"] == ("1:1320,4061",)
    assert "1:1320,4061" in json.dumps(result.to_dict())

    analysis = InterfaceGeometryAnalysis(
        metadata={"alpha_value": 1.4, "alpha_units": "A", "alpha_method": "additive"},
        alpha_selection=source["alpha_selection"],
    )
    iface = SimpleNamespace(
        net=network,
        _geometry_dual=_dual(network),
        geometry_analysis=analysis,
        water_topology=None,
        representation_caches={},
    )
    _stub_identity(monkeypatch)
    monkeypatch.setattr(
        visualization,
        "_cached_interface_result",
        lambda _iface: SimpleNamespace(
            identity=_identity("voronoi"), provenance=Provenance("cached physical state"),
            metrics={}, selection={}, topology=None, mean_curvature=None,
            gaussian_curvature=None, contacts=(),
        ),
    )

    cached = visualization._cached_alpha_selection_result(iface)
    assert cached.status is Status.UNRESOLVED
    assert cached.selected_dual_feature_ids == result.selected_dual_feature_ids

    path = tmp_path / "logs.csv"
    visualization._compact_log(iface, path)
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    alpha_rows = [row for row in rows if row["representation"] == "alpha_selection"]
    assert alpha_rows
    assert all(row["status"] == Status.UNRESOLVED.value for row in alpha_rows)
    assert all(row["value"] == "" for row in alpha_rows)
    assert all("1:1320,4061" in row["provenance"] for row in alpha_rows)


@pytest.mark.parametrize("cache_state", ["missing", "stale", "valid"])
def test_cached_dual_results_keep_cache_state_explicit(monkeypatch, cache_state):
    network = object()
    dual = None if cache_state == "missing" else _dual(
        object() if cache_state == "stale" else network
    )
    iface = SimpleNamespace(net=network, _geometry_dual=dual)
    _stub_identity(monkeypatch)

    contact, generators = visualization._cached_dual_results(iface)

    expected = {
        "missing": Status.NOT_CALCULATED,
        "stale": Status.STALE,
        "valid": Status.CERTIFIED,
    }[cache_state]
    assert contact.status is generators.status is expected
    if expected is Status.CERTIFIED:
        assert contact.metrics["simplex_count_dim_1"].value == 1
        assert generators.generator_complex.simplex_counts[1].value == 1
    else:
        assert contact.metrics == generators.metrics == {}
        assert generators.generator_complex is None
    a_result, b_result = visualization._cached_dual_result_bundle(iface)[2:]
    assert a_result.status is b_result.status is {
        "missing": Status.NOT_CALCULATED,
        "stale": Status.STALE,
        "valid": Status.NOT_SUPPORTED,
    }[cache_state]


def test_atlas_membership_preserves_dimensions_and_filters_contact_edges(monkeypatch):
    network = object()
    iface = SimpleNamespace(net=network, _geometry_dual=_dimensioned_dual(network))
    _stub_identity(monkeypatch)

    contact, combined, a_result, b_result = visualization._cached_dual_result_bundle(iface)

    assert [record.simplex_id for record in contact.contact_features] == ["1:3,8"]
    assert [record.simplex_id for record in combined.simplex_records] == ["1:3,8", "2:3,8,9"]
    assert [record.simplex_id for record in a_result.simplex_records] == ["0:3"]
    assert [record.simplex_id for record in b_result.simplex_records] == ["1:8,9"]
    assert all(result.status is Status.CERTIFIED for result in (contact, combined, a_result, b_result))
    assert all(result.identity.network_id == "2KAI:aw:fixture" for result in (contact, combined, a_result, b_result))
    assert all(result.identity.partition == "aw" for result in (contact, combined, a_result, b_result))
    assert all(result.provenance.details["cache_revision"] == "ATLAS-A5"
               for result in (contact, combined, a_result, b_result))


def test_incomplete_dual_is_partial_without_inventing_a_or_b_membership(monkeypatch):
    network = object()
    dual = _dual(network)
    edge = dual.simplices[1][(3, 8)]
    edge.complete = False
    edge.status = "partial"
    iface = SimpleNamespace(net=network, _geometry_dual=dual)
    _stub_identity(monkeypatch)

    contact, combined, a_result, b_result = visualization._cached_dual_result_bundle(iface)

    assert contact.status is combined.status is Status.PARTIAL
    assert contact.contact_features[0].status is Status.PARTIAL
    assert a_result.status is b_result.status is Status.NOT_SUPPORTED
    assert a_result.generator_complex is b_result.generator_complex is None


def test_canonical_log_emits_separate_a_and_b_generator_rows(monkeypatch, tmp_path):
    network = object()
    iface = SimpleNamespace(
        net=network,
        _geometry_dual=_dimensioned_dual(network),
        geometry_analysis=None,
        water_topology=None,
        representation_caches={},
    )
    _stub_identity(monkeypatch)
    monkeypatch.setattr(visualization, "_cached_interface_result", lambda _iface: SimpleNamespace(
        identity=_identity("voronoi"), provenance=Provenance("cached physical state"),
        metrics={}, selection={}, topology=None, mean_curvature=None,
        gaussian_curvature=None, contacts=(),
    ))
    monkeypatch.setattr(visualization, "_cached_alpha_selection_result", lambda _iface: None)

    path = tmp_path / "logs.csv"
    visualization._compact_log(iface, path)
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    logged = {(row["representation"], row["side"], row["quantity"]): row for row in rows}
    assert logged[("dual_generator_complex", "A", "simplex_count_dim_0")]["value"] == "1"
    assert logged[("dual_generator_complex", "B", "simplex_count_dim_1")]["value"] == "1"
    assert logged[("dual_generator_complex", "AB/contact", "simplex_count_dim_2")]["value"] == "1"
    assert logged[("dual_contact_complex", "shared", "simplex_count_dim_1")]["value"] == "1"
    assert all(
        row["status"] == Status.CERTIFIED.value
        for row in rows
        if row["representation"].startswith("dual_")
    )


def test_partition_mismatched_dual_is_stale_without_traversal(monkeypatch):
    network = object()
    dual = _dual(network)
    dual.scheme = "pow"
    dual.audit = lambda: (_ for _ in ()).throw(AssertionError("stale dual was traversed"))
    iface = SimpleNamespace(net=network, _geometry_dual=dual)
    _stub_identity(monkeypatch)

    results = visualization._cached_dual_result_bundle(iface)

    assert all(result.status is Status.STALE for result in results)


@pytest.mark.parametrize("field", ["network_id", "interface_id"])
def test_cached_dual_dependency_identity_mismatch_is_stale(monkeypatch, field):
    network = object()
    dual = _dimensioned_dual(network)
    dual.archive_provenance[field] = "different-cache-owner"
    iface = SimpleNamespace(net=network, _geometry_dual=dual)
    _stub_identity(monkeypatch)

    results = visualization._cached_dual_result_bundle(iface)

    assert all(result.status is Status.STALE for result in results)
    assert all(
        "cached dual references another" in " ".join(result.provenance.diagnostics)
        for result in results
    )


def test_canonical_log_rows_are_deterministic_provenance_bearing_and_null_when_unavailable(
    monkeypatch, tmp_path
):
    analysis = InterfaceGeometryAnalysis(**json.loads(FIXTURE.read_text(encoding="utf-8")))
    analysis.archive_provenance = {"fixture": "2kai_aw_cached_analysis"}
    iface = SimpleNamespace(
        net=object(),
        geometry_analysis=analysis,
        _geometry_dual=None,
        water_topology=None,
        representation_caches={},
    )
    _stub_identity(monkeypatch)
    monkeypatch.setattr(
        visualization,
        "_cached_interface_result",
        lambda _iface: InterfaceResult(
            _identity("voronoi"), Provenance("cached physical state"), Status.NOT_CALCULATED
        ),
    )

    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    visualization._compact_log(iface, first)
    visualization._compact_log(iface, second)

    assert first.read_bytes() == second.read_bytes()
    with first.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert tuple(rows[0]) == INTERFACE_LOG_COLUMNS
    identities = [
        tuple(row[name] for name in ("interface_id", "representation", "side", "quantity"))
        for row in rows
    ]
    assert identities == sorted(identities)
    assert len(identities) == len(set(identities))
    assert {row["status"] for row in rows} >= {Status.PARTIAL.value, Status.NOT_CALCULATED.value}
    assert all(
        row["value"] == ""
        for row in rows
        if row["status"] in {Status.NOT_CALCULATED.value, Status.NOT_SUPPORTED.value, Status.UNRESOLVED.value}
    )
    selected = next(
        row
        for row in rows
        if row["representation"] == "alpha_selection"
        and row["quantity"] == "alpha_selected_count"
    )
    assert selected["status"] == Status.PARTIAL.value
    assert "2kai_aw_cached_analysis" in selected["provenance"]
