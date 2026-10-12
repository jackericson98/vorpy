from dataclasses import replace
from types import SimpleNamespace

import pytest

from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.results import (
    GeometryLineage,
    GeometryLineageOrigin,
    Provenance,
    Quantity,
    ResultsDependency,
    ResultIdentity,
    SharedNetworkIdentity,
    Status,
    adapt_cached_result,
    adapt_group_view_analysis,
    adapt_interface_analysis,
)
from vorpy.src.results.model import dictionary


def _lineage(origin, *, complete=True, exact=False, references=False):
    values = {
        "origin": origin,
        "solve_universe_id": "solve:aw:123",
        "geometry_source_id": "geometry:sha256:456",
        "geometry_version": "aw-geometry-reuse-v1",
        "solver_settings_digest": "settings:sha256:789",
        "environment_digest": "environment:sha256:abc",
        "source_network_id": "network:source",
        "selection_scope": {"groups": {"A": [0], "B": [1]}},
        "orientation": {"convention": "physical_interface_A_to_B"},
        "completeness": {
            "completeness_proven": complete,
            "vertex_completeness_proven": complete,
        },
    }
    if references:
        values.update(
            solver_settings=None,
            solver_settings_ref="settings-ref:sha256:789",
            environment=None,
            environment_ref="environment-ref:sha256:abc",
        )
    else:
        values.update(
            solver_settings={"net_type": "aw", "surface_resolution": 32},
            environment={"network_environment": "dry"},
        )
    if exact:
        values.update(
            source_view_id="view:sha256:def",
            compatibility={
                "decision": "EXACT_REUSE",
                "generator_compatibility": {"match": True, "digest": "g"},
                "competitor_compatibility": {"match": True, "digest": "c"},
                "solve_context": {"match": True, "digest": "s"},
            },
        )
    return GeometryLineage(**values)


def _interface_identity():
    return ResultIdentity(
        result_kind="interface",
        system="synthetic",
        frame="static",
        network_id="archive:network",
        network_scope="system",
        environment="dry",
        partition="aw",
        representation="voronoi",
        radius_configuration_id="radii:synthetic",
        interface_id="A_B",
        group_A="A",
        group_B="B",
        side="shared",
        orientation="group1 -> group2",
        alpha_method="additive",
        alpha_value=1.0,
        alpha_units="A",
        probe_radius=1.4,
        probe_units="Å",
    )


def _cached_interface_analysis():
    return InterfaceGeometryAnalysis(metadata={}, voronoi_side_1={})


def _shared_identity():
    return SharedNetworkIdentity(
        source_network_id="producer:network",
        source_identity_status="AUTHORITATIVE",
        solve_universe_id="solve:context",
        frame_id="static",
        analysis_context_kind="GLOBAL_NETWORK",
        geometry_scheme="aw",
        partition_scheme="aw",
        scope_mode="full",
        solvent_mode="MOLECULE_ONLY",
        solvent_cells_complete=False,
        selection_id="view:group",
    )


@pytest.mark.parametrize(
    "origin, exact",
    [
        (GeometryLineageOrigin.DIRECT_SOLVE, False),
        (GeometryLineageOrigin.EXACT_REUSE, True),
        (GeometryLineageOrigin.MATERIALIZED_VIEW, False),
    ],
)
def test_resolved_lineage_origins_are_valid(origin, exact):
    lineage = _lineage(origin, exact=exact)
    assert lineage.origin is origin
    assert lineage.solve_universe_id == "solve:aw:123"


@pytest.mark.parametrize(
    "lineage",
    [
        _lineage(GeometryLineageOrigin.DIRECT_SOLVE),
        _lineage(GeometryLineageOrigin.MATERIALIZED_VIEW),
        _lineage(GeometryLineageOrigin.EXACT_REUSE, exact=True),
        GeometryLineage(
            GeometryLineageOrigin.UNRESOLVED,
            reason="producer did not preserve geometry provenance",
        ),
    ],
)
def test_typed_interface_adapter_preserves_producer_lineage(lineage):
    result = adapt_interface_analysis(
        _cached_interface_analysis(),
        _interface_identity(),
        geometry_lineage=lineage,
    )

    assert result.provenance.geometry_lineage == lineage
    assert result.provenance.details["geometry_lineage"]["origin"] == (
        lineage.origin.value
    )


def test_adapters_do_not_invent_lineage_and_merge_existing_details():
    identity = _interface_identity()
    result = adapt_interface_analysis(_cached_interface_analysis(), identity)
    assert result.provenance.geometry_lineage is None

    lineage = _lineage(GeometryLineageOrigin.DIRECT_SOLVE)
    result = adapt_interface_analysis(
        _cached_interface_analysis(),
        identity,
        provenance=Provenance("producer", details={"cache_key": "cache:1"}),
        geometry_lineage=lineage,
    )
    assert result.provenance.details["cache_key"] == "cache:1"
    assert result.provenance.geometry_lineage == lineage


def test_cached_result_rejects_conflicting_producer_lineage():
    identity = _interface_identity()
    first = _lineage(GeometryLineageOrigin.DIRECT_SOLVE)
    second = _lineage(GeometryLineageOrigin.MATERIALIZED_VIEW)

    with pytest.raises(ValueError, match="Conflicting geometry lineage"):
        adapt_cached_result(
            identity,
            provenance=Provenance(
                "producer", details={"geometry_lineage": first}
            ),
            geometry_lineage=second,
        )


def test_unresolved_handoff_preserves_existing_archived_lineage():
    existing = _lineage(GeometryLineageOrigin.MATERIALIZED_VIEW)
    fallback = GeometryLineage(
        GeometryLineageOrigin.UNRESOLVED,
        reason="solver evidence unavailable after archive reload",
    )

    result = adapt_cached_result(
        _interface_identity(),
        provenance=Provenance(
            "archive cache",
            details={"geometry_lineage": existing, "cache_key": "archive:1"},
        ),
        status=Status.PARTIAL,
        geometry_lineage=fallback,
    )

    assert result.provenance.geometry_lineage == existing
    assert result.provenance.details["cache_key"] == "archive:1"
    assert result.status is Status.PARTIAL


def test_cached_adapter_preserves_identity_status_and_provenance_details():
    identity = replace(
        _interface_identity(),
        source_network_id="producer:network",
        source_identity_status="PROVISIONAL",
    )
    lineage = replace(
        _lineage(GeometryLineageOrigin.EXACT_REUSE, exact=True),
        source_network_id="producer:network",
    )
    provenance = Provenance(
        "producer cache",
        source_revision="producer-v2",
        scientific_state_id="state:cached",
        details={"cache_key": "cache:1"},
        diagnostics=("producer diagnostic",),
    )

    result = adapt_cached_result(
        identity,
        provenance=provenance,
        status=Status.UNRESOLVED,
        geometry_lineage=lineage,
    )

    assert result.identity.network_id == "archive:network"
    assert result.identity.source_network_id == "producer:network"
    assert result.identity.source_identity_status == "PROVISIONAL"
    assert result.status is Status.UNRESOLVED
    assert result.provenance.details["cache_key"] == "cache:1"
    assert result.provenance.source_revision == "producer-v2"
    assert result.provenance.scientific_state_id == "state:cached"
    assert result.provenance.diagnostics == ("producer diagnostic",)
    assert result.provenance.geometry_lineage == lineage


def test_adapter_rejects_unvalidated_producer_lineage():
    with pytest.raises(ValueError, match="Resolved geometry lineage is incomplete"):
        adapt_cached_result(
            _interface_identity(),
            provenance=Provenance("producer cache"),
            geometry_lineage={"origin": "DIRECT_SOLVE"},
        )


def test_nexus_group_adapter_preserves_materialized_lineage():
    identity = _shared_identity()
    analysis = SimpleNamespace(
        source_network_id=identity.source_network_id,
        solve_universe_id=identity.solve_universe_id,
        selection_id=identity.selection_id,
        boundary_surface_ids=("surface:1",),
        area=2.0,
        volume=1.0,
        complete=True,
        boundary_dependent_area=0.0,
        physical_interface_area=2.0,
    )
    dependency = ResultsDependency(
        identity=identity,
        source_id="view-cache:group",
        source_revision="producer-v1",
    )
    lineage = _lineage(GeometryLineageOrigin.MATERIALIZED_VIEW)

    result = adapt_group_view_analysis(
        analysis,
        identity,
        dependency,
        group_id="group-A",
        system="synthetic",
        radius_configuration_id="radii:synthetic",
        source_solve_context=SimpleNamespace(
            solve_universe_id=identity.solve_universe_id,
            frame_id=identity.frame_id,
            scheme=identity.geometry_scheme,
            partition_scheme=identity.partition_scheme,
            scope_mode=identity.scope_mode,
        ),
        source_solvent_provenance={
            "solvent_mode": identity.solvent_mode,
            "solvent_generator_digest": identity.solvent_generator_digest,
            "solvent_cells_complete": identity.solvent_cells_complete,
        },
        source_retained_generator_provenance={
            "policy": "explicit_system_indices",
            "retained_system_indices": (),
            "retained_generator_ids": (),
            "retained_generator_digest": None,
            "boundary_independent": True,
        },
        geometry_lineage=lineage,
    )

    assert result.provenance.geometry_lineage == lineage


def test_resolved_lineage_accepts_immutable_metadata_references():
    lineage = _lineage(GeometryLineageOrigin.DIRECT_SOLVE, references=True)
    assert lineage.solver_settings is None
    assert lineage.solver_settings_ref == "settings-ref:sha256:789"
    assert lineage.environment is None
    assert lineage.environment_ref == "environment-ref:sha256:abc"


def test_unresolved_lineage_allows_unknown_identifiers_without_fabrication():
    lineage = GeometryLineage(
        GeometryLineageOrigin.UNRESOLVED,
        reason="legacy geometry provenance unavailable",
    )
    payload = lineage.as_dict()
    assert payload["origin"] == "UNRESOLVED"
    assert payload["solve_universe_id"] is None
    assert payload["geometry_source_id"] is None


def test_invalid_origin_is_rejected():
    with pytest.raises(ValueError, match="Unknown geometry lineage origin"):
        GeometryLineage("REUSED_SOMEHOW", reason="unsupported")


def test_resolved_lineage_requires_solver_and_environment_identity():
    with pytest.raises(ValueError, match="solver settings"):
        replace(
            _lineage(GeometryLineageOrigin.DIRECT_SOLVE),
            solver_settings=None,
            solver_settings_ref=None,
        )

    with pytest.raises(ValueError, match="environment metadata"):
        GeometryLineage(
            GeometryLineageOrigin.DIRECT_SOLVE,
            solve_universe_id="solve:aw:123",
            geometry_source_id="geometry:sha256:456",
            geometry_version="aw-geometry-reuse-v1",
            solver_settings_digest="settings:sha256:789",
            solver_settings={"net_type": "aw"},
            environment_digest="environment:sha256:abc",
            source_network_id="network:source",
            selection_scope={"scope": "full"},
            orientation={"convention": "physical_interface_A_to_B"},
            completeness={"completeness_proven": True},
        )


def test_exact_reuse_requires_positive_compatibility_and_completeness():
    with pytest.raises(ValueError, match="compatibility evidence"):
        GeometryLineage(
            GeometryLineageOrigin.EXACT_REUSE,
            solve_universe_id="solve:aw:123",
            geometry_source_id="geometry:sha256:456",
            geometry_version="aw-geometry-reuse-v1",
            solver_settings_digest="settings:sha256:789",
            solver_settings={"net_type": "aw"},
            environment_digest="environment:sha256:abc",
            environment={"network_environment": "dry"},
            source_network_id="network:source",
            source_view_id="view:sha256:def",
            selection_scope={"scope": "interface"},
            orientation={"convention": "physical_interface_A_to_B"},
            completeness={"completeness_proven": True},
        )

    with pytest.raises(ValueError, match="proven completeness"):
        _lineage(GeometryLineageOrigin.EXACT_REUSE, complete=False, exact=True)


def test_materialized_view_does_not_require_certification_or_promote_status():
    lineage = _lineage(GeometryLineageOrigin.MATERIALIZED_VIEW, complete=False)
    provenance = Provenance("materialized source view", details={"geometry_lineage": lineage})
    quantity = Quantity(
        "area",
        1.0,
        "Å²",
        Status.PARTIAL,
        "stored primitive view",
        "materialized_view",
        provenance=provenance,
    )
    assert provenance.geometry_lineage.origin is GeometryLineageOrigin.MATERIALIZED_VIEW
    assert provenance.geometry_lineage.completeness["completeness_proven"] is False
    assert quantity.status is Status.PARTIAL


def test_geometry_origin_does_not_change_scientific_status():
    provenance = Provenance(
        "exact reused geometry",
        details={"geometry_lineage": _lineage(GeometryLineageOrigin.EXACT_REUSE, exact=True)},
    )
    quantity = Quantity(
        "area",
        12.5,
        "Å²",
        Status.CERTIFIED,
        "independently certified cached measurement",
        "full_interface",
        provenance=provenance,
    )
    assert provenance.geometry_lineage.origin is GeometryLineageOrigin.EXACT_REUSE
    assert quantity.status is Status.CERTIFIED


def test_legacy_missing_lineage_is_not_fabricated():
    provenance = Provenance("legacy archive", details={"partition": "aw"})
    assert provenance.geometry_lineage is None
    assert "geometry_lineage" not in provenance.details


def test_provenance_mapping_round_trip_preserves_nulls_and_evidence():
    lineage = _lineage(GeometryLineageOrigin.EXACT_REUSE, exact=True)
    provenance = Provenance(
        "typed geometry result",
        source_revision="producer-v2",
        scientific_state_id="state:123",
        details={"geometry_lineage": lineage.as_dict(), "partition": "aw"},
    )

    serialized = dictionary(provenance)
    restored = Provenance(
        serialized["source"],
        source_revision=serialized["source_revision"],
        scientific_state_id=serialized["scientific_state_id"],
        details=serialized["details"],
    )

    assert restored.geometry_lineage == lineage
    assert serialized["details"]["geometry_lineage"]["source_interface_id"] is None
    assert serialized["details"]["geometry_lineage"]["compatibility"][
        "competitor_compatibility"
    ]["match"] is True
