from dataclasses import replace
from types import SimpleNamespace

import pytest

from vorpy.src.results import (
    ResultsDependency,
    ResultIdentity,
    SharedNetworkIdentity,
    Status,
    canonical_retained_generator_digest,
    is_provisional_source_network_id,
    validate_retained_generator_provenance,
    validate_results_dependency,
    validate_solve_context_identity,
    validate_solvent_provenance,
)
from vorpy.src.results.model import dictionary
from vorpy.tests.network.test_retained_generators import _edta_group


def _identity(**updates):
    values = {
        "source_network_id": "archive-network-7",
        "source_identity_status": "AUTHORITATIVE",
        "solve_universe_id": "solve:context-1",
        "frame_id": "7",
        "analysis_context_kind": "GLOBAL_NETWORK",
        "geometry_scheme": "prm",
        "partition_scheme": "prm",
        "scope_mode": "full",
        "selection_id": "view:selection-1",
        "group_A_selection_id": "view:group-a",
        "group_B_selection_id": "view:group-b",
        "orientation": "A -> B",
    }
    values.update(updates)
    return SharedNetworkIdentity(**values)


def _dependency(identity, **updates):
    values = {
        "identity": identity,
        "source_id": "dual-cache-1",
        "source_revision": "producer-1",
    }
    values.update(updates)
    return ResultsDependency(**values)


def _result_identity(**updates):
    values = {
        "result_kind": "interface",
        "system": "synthetic",
        "frame": "7",
        "network_id": "legacy-network",
        "network_scope": "system",
        "environment": "dry",
        "partition": "power",
        "representation": "dual",
        "radius_configuration_id": "radii-1",
        "interface_id": "interface-1",
        "group_A": "A",
        "group_B": "B",
        "side": "shared",
        "orientation": "A -> B",
    }
    values.update(updates)
    return ResultIdentity(**values)


def test_view_selection_changes_dependency_but_not_solve_context_identity():
    first = _identity(selection_id="view:first")
    second = replace(first, selection_id="view:second")

    assert first.solve_universe_id == second.solve_universe_id
    assert first.context_payload() == second.context_payload()
    assert _dependency(first).dependency_digest != _dependency(second).dependency_digest


def test_geometry_affecting_solve_change_changes_dependency_identity():
    baseline = _dependency(_identity())
    changed_solve = _dependency(
        replace(baseline.identity, solve_universe_id="solve:context-2")
    )

    assert baseline.dependency_digest != changed_solve.dependency_digest
    assert validate_results_dependency(baseline, changed_solve) is Status.STALE


def test_explicit_solve_context_is_required_to_match_canonical_identity():
    identity = _identity(
        group_A_selection_id=None,
        group_B_selection_id=None,
        orientation=None,
    )
    context = SimpleNamespace(
        solve_universe_id="solve:context-1",
        frame_id="7",
        scheme="prm",
        partition_scheme="prm",
        scope_mode="full",
    )

    assert validate_solve_context_identity(identity, context) is Status.CERTIFIED
    assert validate_solve_context_identity(
        identity,
        SimpleNamespace(**{**vars(context), "solve_universe_id": "solve:context-2"}),
    ) is Status.UNRESOLVED


def test_authoritative_solvent_context_fields_are_checked_when_present():
    identity = _identity(
        group_A_selection_id=None,
        group_B_selection_id=None,
        orientation=None,
        solvent_mode="SOLVENT_CONTEXT",
        solvent_generator_digest="solvent-generators-1",
        solvent_cells_complete=False,
    )
    context = SimpleNamespace(
        solve_universe_id="solve:context-1",
        frame_id="7",
        scheme="prm",
        partition_scheme="prm",
        scope_mode="full",
        solvent_mode="SOLVENT_CONTEXT",
        solvent_generator_digest="solvent-generators-1",
        solvent_cells_complete=False,
    )

    assert validate_solve_context_identity(identity, context) is Status.CERTIFIED
    assert validate_solve_context_identity(
        identity,
        SimpleNamespace(**{**vars(context), "solvent_mode": "FULL_SOLVENT"}),
    ) is Status.UNRESOLVED
    del context.solvent_cells_complete
    assert validate_solve_context_identity(identity, context) is Status.UNRESOLVED


def test_f10_solvent_provenance_mapping_is_explicit_and_fail_closed():
    identity = _identity(
        group_A_selection_id=None,
        group_B_selection_id=None,
        orientation=None,
        solvent_mode="SOLVENT_CONTEXT",
        solvent_generator_digest="solvent-generators-1",
        solvent_cells_complete=False,
    )
    provenance = {
        "solvent_mode": "SOLVENT_CONTEXT",
        "solvent_generator_digest": "solvent-generators-1",
        "solvent_cells_complete": False,
        "solvent_generator_ids": ("solvent-1",),
    }

    assert validate_solvent_provenance(identity, provenance) is Status.CERTIFIED
    assert validate_solvent_provenance(
        identity,
        {**provenance, "solvent_mode": "MOLECULE_ONLY"},
    ) is Status.UNRESOLVED
    assert validate_solvent_provenance(identity, {"solvent_mode": "SOLVENT_CONTEXT"}) is Status.UNRESOLVED


def test_solvent_mode_changes_dependency_without_changing_solve_id():
    molecule_only = _identity()
    solvent_context = replace(
        molecule_only,
        solvent_mode="SOLVENT_CONTEXT",
        solvent_generator_digest="solvent-generators-1",
        solvent_cells_complete=False,
    )
    full_solvent = replace(
        solvent_context,
        solvent_mode="FULL_SOLVENT",
        solvent_cells_complete=True,
    )

    assert molecule_only.solve_universe_id == solvent_context.solve_universe_id
    assert solvent_context.solve_universe_id == full_solvent.solve_universe_id
    assert len({
        _dependency(molecule_only).dependency_digest,
        _dependency(solvent_context).dependency_digest,
        _dependency(full_solvent).dependency_digest,
    }) == 3
    assert validate_results_dependency(
        _dependency(solvent_context), _dependency(molecule_only)
    ) is Status.STALE


def test_invalid_solvent_provenance_is_rejected():
    try:
        replace(_identity(), solvent_mode="SOLVENT_CONTEXT")
    except ValueError as error:
        assert "generator digest" in str(error)
    else:
        raise AssertionError("missing solvent generator evidence must fail closed")

    try:
        replace(
            _identity(
                solvent_mode="FULL_SOLVENT",
                solvent_generator_digest="solvent-generators-1",
            )
        )
    except ValueError as error:
        assert "complete solvent cells" in str(error)
    else:
        raise AssertionError("incomplete full-solvent evidence must fail closed")


def test_reversing_ordered_interface_sides_changes_selection_identity():
    forward = _identity(
        selection_id="view:forward",
        group_A_selection_id="view:group-a",
        group_B_selection_id="view:group-b",
    )
    reverse = _identity(
        selection_id="view:reverse",
        group_A_selection_id="view:group-b",
        group_B_selection_id="view:group-a",
        orientation="B -> A",
    )

    assert forward.solve_universe_id == reverse.solve_universe_id
    assert forward.payload() != reverse.payload()
    assert _dependency(forward).dependency_digest != _dependency(reverse).dependency_digest


def test_dependency_revision_mismatch_is_stale_without_rebuild():
    expected = _dependency(_identity())
    cached = _dependency(_identity(), source_revision="producer-0")

    assert validate_results_dependency(expected, cached) is Status.STALE
    assert validate_results_dependency(expected, None) is Status.NOT_CALCULATED


def test_dependency_status_is_evidence_and_stale_is_never_promoted():
    certified = _dependency(_identity())
    stale = replace(certified, status=Status.STALE)
    unresolved = replace(certified, status=Status.UNRESOLVED)

    assert stale.dependency_digest == certified.dependency_digest
    assert validate_results_dependency(certified, stale) is Status.STALE
    assert validate_results_dependency(stale, certified) is Status.STALE
    assert validate_results_dependency(unresolved, certified) is Status.UNRESOLVED


def test_solve_derived_network_namespace_is_unresolved_not_certified():
    provisional = "network:context-1"
    assert is_provisional_source_network_id(provisional, "solve:context-1")

    dependency = _dependency(
        _identity(
            source_network_id=provisional,
            source_identity_status="PROVISIONAL",
        )
    )
    assert validate_results_dependency(dependency, dependency) is Status.UNRESOLVED


def test_archive_namespace_without_authority_stays_unresolved_even_with_digest():
    identity = _identity(
        source_network_id="archive-only-identity",
        source_identity_status="UNRESOLVED",
    )
    dependency = _dependency(
        identity,
        source_revision="archive-only-identity",
        source_digest="digest-only",
    )

    assert validate_results_dependency(dependency, dependency) is Status.UNRESOLVED


def test_authoritative_status_is_required_and_provisional_status_cannot_certify():
    authoritative = _dependency(_identity())
    assert validate_results_dependency(authoritative, authoritative) is Status.CERTIFIED

    provisional = _dependency(
        _identity(
            source_network_id="archive-only-identity",
            source_identity_status="PROVISIONAL",
        )
    )
    assert validate_results_dependency(provisional, provisional) is Status.UNRESOLVED

    with pytest.raises(ValueError, match="cannot be authoritative"):
        _identity(
            source_network_id="network:context-1",
            source_identity_status="AUTHORITATIVE",
        )


def test_conflicting_source_identity_status_is_unresolved_not_stale():
    expected = _dependency(_identity())
    cached = _dependency(
        replace(expected.identity, source_identity_status="UNRESOLVED")
    )

    assert validate_results_dependency(expected, cached) is Status.UNRESOLVED


def test_result_identity_projects_canonical_solve_and_view_ids():
    legacy = _result_identity()
    shared = _result_identity(
        source_network_id="archive-network-7",
        solve_universe_id="solve:context-1",
        analysis_context_kind="GLOBAL_NETWORK",
        selection_id="view:selection-1",
    )

    assert legacy.result_id != shared.result_id
    assert shared.solve_universe_id == "solve:context-1"
    assert shared.selection_id == "view:selection-1"


def test_result_identity_projects_solvent_mode():
    result = _result_identity(
        source_network_id="archive-network-7",
        solve_universe_id="solve:context-1",
        analysis_context_kind="GLOBAL_NETWORK",
        selection_id="view:selection-1",
        solvent_mode="SOLVENT_CONTEXT",
    )

    assert result.solvent_mode == "SOLVENT_CONTEXT"


def test_legacy_result_identity_omits_optional_solvent_provenance():
    assert "solvent_mode" not in dictionary(_result_identity())


def _retention_provenance(indices=(32,), generator_ids=("system:32",)):
    return {
        "policy": "explicit_system_indices",
        "retained_system_indices": indices,
        "retained_generator_ids": generator_ids,
        "retained_generator_digest": canonical_retained_generator_digest(
            indices, generator_ids
        ),
        "boundary_independent": True,
    }


def test_f12_default_edta_retention_evidence_is_explicit_and_legacy_safe():
    identity = _identity()
    default_provenance = {
        "policy": "explicit_system_indices",
        "retained_system_indices": (),
        "retained_generator_ids": (),
        "retained_generator_digest": None,
        "boundary_independent": True,
    }

    assert validate_retained_generator_provenance(identity, default_provenance) is Status.CERTIFIED
    assert "retained_generator_digest" not in identity.context_payload()


def test_f12_retained_mg_changes_results_dependency_identity():
    default = _identity()
    retained = replace(
        default,
        retained_system_indices=(32,),
        retained_generator_ids=("system:32",),
        retained_generator_digest=canonical_retained_generator_digest(
            (32,), ("system:32",)
        ),
    )

    assert _dependency(default).dependency_digest != _dependency(retained).dependency_digest
    assert validate_retained_generator_provenance(
        retained, _retention_provenance()
    ) is Status.CERTIFIED


def test_real_edta_default_and_retained_mg_have_distinct_solve_and_dependency_ids():
    _, default_group = _edta_group()
    _, retained_group = _edta_group((32,))

    def project(network, selection_id):
        context = network.get_solve_context_key(frame_id="static")
        solvent = network.solvent_provenance
        retention = network.retained_generator_provenance
        return SharedNetworkIdentity(
            source_network_id="edta-source-network",
            source_identity_status="AUTHORITATIVE",
            solve_universe_id=context.solve_universe_id,
            frame_id=context.frame_id,
            analysis_context_kind="GLOBAL_NETWORK",
            geometry_scheme=context.scheme,
            partition_scheme=context.partition_scheme,
            scope_mode=context.scope_mode,
            solvent_mode=solvent["solvent_mode"],
            solvent_generator_digest=solvent["solvent_generator_digest"],
            solvent_cells_complete=solvent["solvent_cells_complete"],
            retained_system_indices=retention["retained_system_indices"],
            retained_generator_ids=retention["retained_generator_ids"],
            retained_generator_digest=retention["retained_generator_digest"],
            selection_id=selection_id,
        )

    default_identity = project(default_group.net, "view:edta-default")
    retained_identity = project(retained_group.net, "view:edta-retained-mg")
    default_dependency = _dependency(default_identity)
    retained_dependency = _dependency(retained_identity)

    assert default_identity.solve_universe_id != retained_identity.solve_universe_id
    assert default_dependency.dependency_digest != retained_dependency.dependency_digest
    assert retained_identity.retained_generator_ids == ("system:32",)


def test_f12_retention_evidence_is_checked_when_solve_context_exposes_it():
    retained = replace(
        _identity(),
        retained_system_indices=(32,),
        retained_generator_ids=("system:32",),
        retained_generator_digest=canonical_retained_generator_digest(
            (32,), ("system:32",)
        ),
    )
    context = SimpleNamespace(
        solve_universe_id=retained.solve_universe_id,
        frame_id=retained.frame_id,
        scheme=retained.geometry_scheme,
        partition_scheme=retained.partition_scheme,
        scope_mode=retained.scope_mode,
        retained_generator_provenance=_retention_provenance(),
    )

    assert validate_solve_context_identity(retained, context) is Status.CERTIFIED
    context.retained_generator_provenance = {
        **_retention_provenance(),
        "retained_generator_ids": ("system:31",),
    }
    assert validate_solve_context_identity(retained, context) is Status.UNRESOLVED


def test_f12_retention_metadata_is_fail_closed_when_partial_or_contradictory():
    retained = replace(
        _identity(),
        retained_system_indices=(32,),
        retained_generator_ids=("system:32",),
        retained_generator_digest=canonical_retained_generator_digest(
            (32,), ("system:32",)
        ),
    )

    assert validate_retained_generator_provenance(
        retained, {"retained_system_indices": (32,)}
    ) is Status.UNRESOLVED
    assert validate_retained_generator_provenance(
        retained, {**_retention_provenance(), "retained_generator_ids": ("system:31",)}
    ) is Status.UNRESOLVED
    with pytest.raises(ValueError, match="inconsistent retained-generator"):
        replace(
            _identity(),
            retained_system_indices=(32,),
            retained_generator_ids=(),
            retained_generator_digest="digest",
        )
