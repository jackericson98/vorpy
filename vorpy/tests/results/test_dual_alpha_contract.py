"""Focused contracts for dual incidence and alpha-selection Results adapters."""

import json
from types import SimpleNamespace

from vorpy.src.results import (
    AlphaSelectionResult,
    DualContactComplexResult,
    DualGeneratorComplexResult,
    GeometryLineage,
    InterfaceLogRow,
    INTERFACE_LOG_COLUMNS,
    Provenance,
    ResultIdentity,
    Status,
    adapt_alpha_selection,
    adapt_dual_contact_complex,
    adapt_dual_generator_complex,
)


def identity(representation, side):
    return ResultIdentity(
        result_kind="interface",
        system="synthetic",
        system_id="synthetic:atoms",
        frame="static",
        network_id="synthetic:power",
        network_scope="system",
        environment="dry",
        partition="power",
        representation=representation,
        radius_configuration_id="synthetic:base",
        interface_id="A_B",
        group_A="A",
        group_B="B",
        side=side,
        orientation="A -> B",
        alpha_method="power_distance",
        alpha_value=0.0,
        alpha_units="A^2",
        probe_radius=1.4,
        probe_units="Å",
    )


def dual_cache():
    vertex = SimpleNamespace(
        simplex_id="0:3",
        dimension=0,
        generator_ids=(3,),
        primal_feature_ids=(("cell", 3),),
        bounded=True,
        complete=True,
        supported=True,
        status="bounded_complete_supported",
    )
    contact = SimpleNamespace(
        simplex_id="1:3,8",
        dimension=1,
        generator_ids=(3, 8),
        primal_feature_ids=(("surf", 17),),
        bounded=True,
        complete=True,
        supported=True,
        status="bounded_complete_supported",
    )
    audit = SimpleNamespace(
        valid=True,
        generator_count=2,
        feature_counts={"cell": 1, "surf": 1},
        unmapped_feature_count=0,
    )
    return SimpleNamespace(simplices={0: {(3,): vertex}, 1: {(3, 8): contact}}, audit=lambda: audit)


def provenance():
    return Provenance(
        "synthetic dual cache",
        details={
            "geometry_lineage": {
                "origin": "UNRESOLVED",
                "reason": "Synthetic incidence fixture has no geometric solve.",
            }
        },
    )


def test_dual_contact_adapter_serializes_incidences_without_primal_measurements():
    result = adapt_dual_contact_complex(
        dual_cache(), identity("dual", "shared"), provenance=provenance()
    )

    assert isinstance(result, DualContactComplexResult)
    assert result.status is Status.CERTIFIED
    assert result.contact_features[0].simplex_id == "1:3,8"
    assert result.contact_features[0].primal_feature_ids == ("surf:17",)
    assert result.metrics["simplex_count_dim_1"].value == 1
    assert result.metrics["dual_incidence_coverage"].value == 1.0
    assert result.metrics["dual_incidence_coverage"].support_count == 2
    assert "area" not in result.metrics
    assert result.provenance.geometry_lineage == GeometryLineage(
        "UNRESOLVED", reason="Synthetic incidence fixture has no geometric solve."
    )
    json.dumps(result.to_dict())


def test_dual_generator_adapter_traverses_dimensioned_simplices_without_topology_inference():
    result = adapt_dual_generator_complex(
        dual_cache(), identity("dual_a", "A"), provenance=provenance()
    )

    assert isinstance(result, DualGeneratorComplexResult)
    assert result.generator_complex.side == "A"
    assert result.generator_complex.simplex_counts[0].value == 1
    assert result.generator_complex.simplex_counts[1].value == 1
    assert result.generator_complex.topology == {}
    assert [record.simplex_id for record in result.simplex_records] == ["0:3", "1:3,8"]
    json.dumps(result.to_dict())


def test_alpha_selection_has_its_own_identity_and_selected_dual_correspondence():
    payload = {
        "status": "available",
        "certified": True,
        "eligible_bicolor_pairs": 1,
        "candidate_bicolor_pairs": 1,
        "selected_pairs": 1,
        "selected_physical_pairs": 1,
        "excluded_pairs": 0,
        "unresolved_pairs": 0,
        "coverage": 1.0,
        "selected_generator_pairs": [[3, 8]],
        "selected_system_pairs": [["chain:A", "chain:B"]],
        "selected_physical_system_pairs": [["system:A", "system:B"]],
    }
    result = adapt_alpha_selection(
        payload,
        identity("alpha_selection", "shared"),
        dual=dual_cache(),
        provenance=provenance(),
    )

    assert isinstance(result, AlphaSelectionResult)
    assert result.status is Status.CERTIFIED
    assert result.selected_generator_pairs == ((3, 8),)
    assert result.selected_dual_feature_ids == ("1:3,8",)
    assert result.selected_dual_features[0].simplex_id == "1:3,8"
    assert result.selection["coverage"].value == 1.0
    json.dumps(result.to_dict())


def test_unavailable_alpha_selection_does_not_turn_unknown_counts_into_zero():
    result = adapt_alpha_selection(None, identity("alpha_selection", "shared"))

    assert result.status is Status.NOT_CALCULATED
    assert result.selection == {}
    assert result.selected_generator_pairs == ()
    assert result.selected_dual_feature_ids == ()


def test_unresolved_alpha_cache_coverage_sentinel_remains_unavailable():
    result = adapt_alpha_selection(
        {"status": "unresolved", "certified": False, "coverage": 0.0},
        identity("alpha_selection", "shared"),
    )

    assert result.status is Status.UNRESOLVED
    assert result.selection["coverage"].value is None
    assert result.selection["coverage"].status is Status.UNRESOLVED


def test_alpha_selection_uses_existing_eight_column_interface_log_schema():
    row = InterfaceLogRow(
        "A_B",
        "alpha_selection",
        "shared",
        "alpha_selected_count",
        1,
        "1",
        Status.CERTIFIED,
    )

    assert tuple(row.as_csv_row()) == INTERFACE_LOG_COLUMNS
    assert len(INTERFACE_LOG_COLUMNS) == 8
