from vorpy.src.results import (
    GeometryLineageOrigin,
    make_geometry_lineage,
)


def _metadata(**overrides):
    values = {
        "geometry_source_id": "geometry:synthetic",
        "geometry_version": "aw-geometry-reuse-v1",
        "solver_settings_digest": "settings:synthetic",
        "solver_settings": {"net_type": "aw"},
        "environment_digest": "environment:synthetic",
        "environment": {"mode": "dry"},
        "source_network_id": "network:synthetic",
        "selection_scope": {"scope": "interface"},
        "orientation": {"convention": "physical_interface_A_to_B"},
        "completeness": {"completeness_proven": False},
    }
    values.update(overrides)
    return values


def test_direct_solve_requires_enabled_solver_evidence():
    unresolved = make_geometry_lineage(
        _metadata(solve_universe_id="solve:synthetic"),
        solver_evidence={
            "schema_version": 1,
            "enabled": False,
            "solve_universe_id": "solve:synthetic",
        },
    )
    assert unresolved.origin is GeometryLineageOrigin.UNRESOLVED

    lineage = make_geometry_lineage(
        _metadata(),
        solver_evidence={
            "schema_version": 1,
            "enabled": True,
            "solve_key": {"digest": "synthetic"},
        },
    )
    assert lineage.origin is GeometryLineageOrigin.DIRECT_SOLVE
    assert lineage.solve_universe_id == "solve:synthetic"


def test_materialized_view_requires_validated_view_evidence():
    evidence = {
        "selection_id": "view:synthetic",
        "source_network_id": "network:view-source",
        "solve_universe_id": "solve:view-source",
        "selection_kind": "interface",
        "group_a_ids": (1,),
        "group_b_ids": (2,),
        "completeness": {"proven": True},
    }
    lineage = make_geometry_lineage(
        _metadata(
            completeness=None,
            source_network_id="network:view-source",
            selection_scope=None,
        ),
        view_evidence=evidence,
    )
    assert lineage.origin is GeometryLineageOrigin.MATERIALIZED_VIEW
    assert lineage.source_view_id == "view:synthetic"
    assert lineage.selection_scope["selection_kind"] == "interface"
    assert lineage.completeness["completeness_proven"] is True

    evidence["completeness"] = {"proven": False}
    unresolved = make_geometry_lineage(
        _metadata(completeness=None, source_network_id="network:view-source"),
        view_evidence=evidence,
    )
    assert unresolved.origin is GeometryLineageOrigin.UNRESOLVED


def test_exact_reuse_requires_positive_certificate():
    certificate = {
        "decision": "EXACT_REUSE",
        "certificate": {"proven": True},
        "comparison": {
            "equivalent": True,
            "same_competitor_universe": True,
            "compatible_geometry_context": True,
            "mismatches": {},
        },
    }
    lineage = make_geometry_lineage(
        _metadata(
            solve_universe_id="solve:reuse-source",
            source_view_id="view:reuse-source",
            completeness=None,
        ),
        reuse_certificate=certificate,
    )
    assert lineage.origin is GeometryLineageOrigin.EXACT_REUSE
    assert lineage.completeness["completeness_proven"] is True

    certificate["comparison"]["equivalent"] = False
    unresolved = make_geometry_lineage(
        _metadata(
            solve_universe_id="solve:reuse-source",
            source_view_id="view:reuse-source",
        ),
        reuse_certificate=certificate,
    )
    assert unresolved.origin is GeometryLineageOrigin.UNRESOLVED


def test_factory_rejects_ambiguous_or_conflicting_evidence():
    unresolved = make_geometry_lineage(
        _metadata(),
        solver_evidence={"schema_version": 1, "enabled": True, "solve_universe_id": "solve:1"},
        view_evidence={
            "selection_id": "view:1",
            "source_network_id": "network:1",
            "solve_universe_id": "solve:1",
            "selection_kind": "group",
            "completeness": {"proven": True},
        },
    )
    assert unresolved.origin is GeometryLineageOrigin.UNRESOLVED

    unresolved = make_geometry_lineage(
        _metadata(solve_universe_id="solve:wrong"),
        solver_evidence={
            "schema_version": 1,
            "enabled": True,
            "solve_universe_id": "solve:right",
        },
    )
    assert unresolved.origin is GeometryLineageOrigin.UNRESOLVED
