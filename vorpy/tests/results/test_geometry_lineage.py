from dataclasses import replace

import pytest

from vorpy.src.results import (
    GeometryLineage,
    GeometryLineageOrigin,
    Provenance,
    Quantity,
    Status,
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
