from dataclasses import replace
from types import SimpleNamespace

import pytest

from vorpy.src.boundary import SolventMode
from vorpy.src.group.shared_view import GroupViewAnalysis
from vorpy.src.io.network_archive import load_network, save_network
from vorpy.src.network import Network
from vorpy.src.results import (
    ResultsDependency,
    SharedNetworkIdentity,
    Status,
    adapt_group_view_analysis,
    canonical_retained_generator_digest,
    collect_shared_results,
)
from vorpy.src.system import System
from vorpy.tests.network.test_views import _network


def _identity(**updates):
    values = {
        "source_network_id": "archive-network-7",
        "source_identity_status": "AUTHORITATIVE",
        "solve_universe_id": "solve:context-1",
        "frame_id": "7",
        "analysis_context_kind": "GLOBAL_NETWORK",
        "geometry_scheme": "aw",
        "partition_scheme": "aw",
        "scope_mode": "full",
        "solvent_cells_complete": False,
        "selection_id": "view:group-a",
    }
    values.update(updates)
    return SharedNetworkIdentity(**values)


def _dependency(identity, **updates):
    values = {
        "identity": identity,
        "source_id": "group-view-cache-1",
        "source_revision": "producer-1",
    }
    values.update(updates)
    return ResultsDependency(**values)


def _analysis(identity, **updates):
    values = {
        "source_network_id": identity.source_network_id,
        "solve_universe_id": identity.solve_universe_id,
        "selection_id": identity.selection_id,
        "system_ids": (10,),
        "boundary_surface_ids": (101, 102),
        "competitor_system_ids": (11,),
        "area": 4.5,
        "volume": 2.0,
        "complete": True,
        "boundary_dependent_area": 0.0,
        "physical_interface_area": 4.5,
    }
    values.update(updates)
    return GroupViewAnalysis(**values)


def _adapt(identity=None, dependency=None, analysis=None, **kwargs):
    identity = identity or _identity()
    dependency = dependency or _dependency(identity)
    analysis = analysis or _analysis(identity)
    kwargs.setdefault(
        "source_solve_context",
        SimpleNamespace(
            solve_universe_id=identity.solve_universe_id,
            frame_id=identity.frame_id,
            scheme=identity.geometry_scheme,
            partition_scheme=identity.partition_scheme,
            scope_mode=identity.scope_mode,
        ),
    )
    kwargs.setdefault(
        "source_solvent_provenance",
        {
            "solvent_mode": identity.solvent_mode,
            "solvent_generator_digest": identity.solvent_generator_digest,
            "solvent_cells_complete": identity.solvent_cells_complete,
        },
    )
    kwargs.setdefault(
        "source_retained_generator_provenance",
        {
            "policy": "explicit_system_indices",
            "retained_system_indices": (),
            "retained_generator_ids": (),
            "retained_generator_digest": None,
            "boundary_independent": True,
        },
    )
    return adapt_group_view_analysis(
        analysis,
        identity,
        dependency,
        group_id="group-A",
        system="synthetic",
        radius_configuration_id="radii-1",
        **kwargs,
    )


def test_group_view_adapter_returns_typed_identity_metrics_and_provenance():
    result = _adapt()

    assert result.status is Status.CERTIFIED
    assert result.identity.group_id == "group-A"
    assert result.identity.network_id == "archive-network-7"
    assert result.identity.source_network_id == "archive-network-7"
    assert result.identity.source_identity_status == "AUTHORITATIVE"
    assert result.identity.solve_universe_id == "solve:context-1"
    assert result.identity.selection_id == "view:group-a"
    assert result.identity.solvent_mode == "MOLECULE_ONLY"
    assert result.identity.frame == "7"
    assert result.metrics["area"].value == pytest.approx(4.5)
    assert result.metrics["volume"].value == pytest.approx(2.0)
    assert result.provenance.details["dependency_digest"].startswith("dependency-v1:")
    assert result.provenance.details["cache_dependencies"][0]["source_id"] == (
        "group-view-cache-1"
    )
    assert result.provenance.details["completeness"]["proven"] is True


def test_boundary_area_components_remain_group_scoped_and_auditable():
    analysis = _analysis(
        _identity(),
        boundary_dependent_area=1.25,
        physical_interface_area=3.25,
    )

    result = _adapt(analysis=analysis)

    assert result.status is Status.CERTIFIED
    assert result.identity.result_kind == "group"
    assert result.identity.interface_id is None
    assert result.metrics["area"].scope == "group_view"
    assert result.metrics["area"].value == pytest.approx(4.5)
    assert result.provenance.details["area_components"] == {
        "total": 4.5,
        "boundary_dependent": 1.25,
        "physical_interface": 3.25,
        "components_proven": True,
        "units": "\u00c5\u00b2",
    }


@pytest.mark.parametrize(
    "field", ("boundary_dependent_area", "physical_interface_area")
)
@pytest.mark.parametrize("value", (None, -1.0, float("nan")))
def test_invalid_area_component_evidence_cannot_certify(field, value):
    result = _adapt(analysis=_analysis(_identity(), **{field: value}))

    assert result.status is Status.UNRESOLVED
    assert result.metrics["area"].value is None
    assert result.provenance.details["area_components"]["components_proven"] is False


@pytest.mark.parametrize(
    "evidence",
    (
        "source_solve_context",
        "source_solvent_provenance",
        "source_retained_generator_provenance",
    ),
)
def test_missing_certification_evidence_is_unresolved(evidence):
    result = _adapt(**{evidence: None})

    assert result.status is Status.UNRESOLVED
    assert result.metrics["area"].value is None
    assert result.metrics["volume"].value is None


@pytest.mark.parametrize(
    "complete, expected",
    [
        (False, Status.PARTIAL),
        ("False", Status.UNRESOLVED),
        (None, Status.UNRESOLVED),
    ],
)
def test_geometry_completeness_requires_an_explicit_boolean(complete, expected):
    result = _adapt(analysis=_analysis(_identity(), complete=complete))

    assert result.status is expected
    if expected is Status.PARTIAL:
        assert result.metrics["area"].value == pytest.approx(4.5)
    else:
        assert result.metrics["area"].value is None


def test_selection_changes_result_dependency_but_not_solve_identity():
    first_identity = _identity(selection_id="view:group-a")
    second_identity = replace(first_identity, selection_id="view:group-b")
    first = _adapt(first_identity, _dependency(first_identity), _analysis(first_identity))
    second = _adapt(
        second_identity,
        _dependency(second_identity),
        _analysis(second_identity, system_ids=(11,)),
    )

    assert first.identity.solve_universe_id == second.identity.solve_universe_id
    assert first.result_id != second.result_id
    assert first.provenance.details["dependency_digest"] != second.provenance.details[
        "dependency_digest"
    ]


def test_group_adapter_validates_explicit_source_solve_context():
    context = SimpleNamespace(
        solve_universe_id="solve:context-1",
        frame_id="7",
        scheme="aw",
        partition_scheme="aw",
        scope_mode="full",
    )
    assert _adapt(source_solve_context=context).status is Status.CERTIFIED

    stale_context = SimpleNamespace(**{**vars(context), "frame_id": "8"})
    unresolved = _adapt(source_solve_context=stale_context)
    assert unresolved.status is Status.UNRESOLVED
    assert unresolved.metrics["area"].value is None


def test_group_adapter_validates_f10_solvent_provenance():
    identity = replace(_identity(), solvent_cells_complete=False)
    provenance = {
        "solvent_mode": "MOLECULE_ONLY",
        "solvent_generator_digest": None,
        "solvent_cells_complete": False,
        "solvent_generator_ids": (),
    }
    result = _adapt(identity, source_solvent_provenance=provenance)
    assert result.status is Status.CERTIFIED
    assert result.provenance.details["source_solvent_provenance"] == provenance

    unresolved = _adapt(
        identity,
        source_solvent_provenance={**provenance, "solvent_mode": "FULL_SOLVENT"},
    )
    assert unresolved.status is Status.UNRESOLVED


@pytest.mark.parametrize(
    "mode, cells_complete",
    [("SOLVENT_CONTEXT", False), ("FULL_SOLVENT", True)],
)
def test_group_adapter_accepts_authoritative_participating_solvent_modes(
    mode, cells_complete
):
    identity = replace(
        _identity(),
        solvent_mode=mode,
        solvent_generator_digest="solvent-generators-1",
        solvent_cells_complete=cells_complete,
    )
    provenance = {
        "solvent_mode": mode,
        "solvent_generator_digest": "solvent-generators-1",
        "solvent_cells_complete": cells_complete,
        "solvent_generator_ids": ("solvent-1",),
    }

    result = _adapt(
        identity,
        _dependency(identity),
        _analysis(identity),
        source_solvent_provenance=provenance,
    )

    assert result.status is Status.CERTIFIED
    assert result.identity.solvent_mode == mode


@pytest.mark.parametrize(
    "identity_updates, dependency_updates",
    [
        ({"solve_universe_id": "solve:context-2"}, {}),
        ({}, {"source_revision": "producer-2"}),
        ({"frame_id": "8"}, {}),
        (
            {
                "solvent_mode": "SOLVENT_CONTEXT",
                "solvent_generator_digest": "solvent-generators-2",
                "solvent_cells_complete": False,
            },
            {},
        ),
    ],
)
def test_geometry_revision_or_frame_change_stales_cached_group_result(
    identity_updates, dependency_updates
):
    cached_identity = _identity()
    cached_dependency = _dependency(cached_identity)
    current_identity = replace(cached_identity, **identity_updates)
    current_dependency = _dependency(current_identity, **dependency_updates)

    result = _adapt(
        current_identity,
        current_dependency,
        _analysis(cached_identity),
        cached_dependency=cached_dependency,
    )

    assert result.status is Status.STALE
    assert result.metrics["area"].status is Status.STALE
    assert result.metrics["volume"].status is Status.STALE


def test_stale_cached_dependency_remains_stale_with_values_preserved():
    identity = _identity()
    stale_cache = _dependency(identity, status=Status.STALE)

    result = _adapt(
        identity,
        _dependency(identity),
        _analysis(identity),
        cached_dependency=stale_cache,
    )

    assert result.status is Status.STALE
    assert result.metrics["area"].status is Status.STALE
    assert result.metrics["area"].value == pytest.approx(4.5)
    assert result.provenance.details["cached_dependency_status"] == "STALE"


def test_missing_view_lineage_and_provisional_source_remain_unresolved():
    identity = _identity(
        source_network_id="network:context-1",
        source_identity_status="PROVISIONAL",
    )
    dependency = _dependency(identity)
    provisional = _adapt(identity, dependency)

    assert provisional.status is Status.UNRESOLVED
    assert provisional.metrics["area"].value is None
    assert provisional.metrics["volume"].value is None

    missing_view_lineage = _adapt(
        analysis=_analysis(identity, source_network_id=None),
    )
    assert missing_view_lineage.status is Status.UNRESOLVED
    assert missing_view_lineage.metrics["area"].value is None


def test_archive_namespace_without_authoritative_source_identity_is_unresolved():
    identity = _identity(
        source_network_id="archive-only-identity",
        source_identity_status="UNRESOLVED",
    )
    result = _adapt(identity, _dependency(identity), _analysis(identity))

    assert result.status is Status.UNRESOLVED
    assert result.identity.source_identity_status == "UNRESOLVED"
    assert result.metrics["area"].value is None
    assert result.metrics["volume"].value is None


def _retained_provenance():
    indices = (32,)
    generator_ids = ("system:32",)
    return {
        "policy": "explicit_system_indices",
        "retained_system_indices": indices,
        "retained_generator_ids": generator_ids,
        "retained_generator_digest": canonical_retained_generator_digest(
            indices, generator_ids
        ),
        "boundary_independent": True,
    }


def test_group_adapter_preserves_explicit_edta_mg_retention_and_fails_closed_without_it():
    provenance = _retained_provenance()
    identity = replace(
        _identity(),
        retained_system_indices=(32,),
        retained_generator_ids=("system:32",),
        retained_generator_digest=provenance["retained_generator_digest"],
    )

    result = _adapt(identity, source_retained_generator_provenance=provenance)
    assert result.status is Status.CERTIFIED
    assert result.provenance.details["source_retained_generator_provenance"] == provenance

    unresolved = _adapt(identity)
    assert unresolved.status is Status.UNRESOLVED
    assert unresolved.metrics["area"].value is None

    contradictory = _adapt(
        identity,
        source_retained_generator_provenance={
            **provenance,
            "retained_generator_digest": "not-the-f12-digest",
        },
    )
    assert contradictory.status is Status.UNRESOLVED
    assert contradictory.metrics["volume"].value is None


def test_group_adapter_accepts_default_edta_f12_retention_mapping():
    provenance = {
        "policy": "explicit_system_indices",
        "retained_system_indices": (),
        "retained_generator_ids": (),
        "retained_generator_digest": None,
        "boundary_independent": True,
    }

    result = _adapt(source_retained_generator_provenance=provenance)

    assert result.status is Status.CERTIFIED
    assert result.provenance.details["source_retained_generator_provenance"] == provenance


def test_retained_membership_change_marks_cached_group_result_stale():
    cached_identity = _identity()
    retained_provenance = _retained_provenance()
    current_identity = replace(
        cached_identity,
        retained_system_indices=(32,),
        retained_generator_ids=("system:32",),
        retained_generator_digest=retained_provenance["retained_generator_digest"],
    )
    current_dependency = _dependency(current_identity)
    cached_dependency = _dependency(cached_identity)

    result = _adapt(
        current_identity,
        current_dependency,
        _analysis(current_identity),
        cached_dependency=cached_dependency,
        source_retained_generator_provenance=retained_provenance,
    )

    assert result.status is Status.STALE
    assert result.metrics["area"].status is Status.STALE


def test_retained_result_validates_after_corrected_archive_reload(tmp_path, monkeypatch):
    network = _network(scheme="aw", source_network_id="retained-archive-source")
    network.group = [0, 1, 2, 3]
    network.sys = System(
        simple=True,
        make_dir=False,
        retained_system_indices=(32,),
    )
    network.sys.balls = network.balls.copy(deep=True)
    network.sys.balls["res"] = None
    network.sys.balls["chn"] = None
    network.solvent_mode = SolventMode.MOLECULE_ONLY.value
    network.retained_system_indices = (32,)
    network.retained_generator_ids = ("system:32",)
    network.refresh_solvent_provenance()
    archive = save_network(network, tmp_path / "retained-results.vpy")

    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("archive reload rebuilt geometry"),
    )
    restored = load_network(archive)
    context = restored.get_solve_context_key()
    view = restored.make_view(group_ids=[10])
    analysis = GroupViewAnalysis(
        source_network_id=view.source_network_id,
        solve_universe_id=view.solve_universe_id,
        selection_id=view.selection_id,
        system_ids=(10,),
        boundary_surface_ids=(0,),
        competitor_system_ids=(11,),
        area=1.0,
        volume=1.0,
        complete=True,
    )
    retention = restored.retained_generator_provenance
    identity = replace(
        _identity(
            source_network_id=view.source_network_id,
            solve_universe_id=context.solve_universe_id,
            frame_id=context.frame_id,
            geometry_scheme=context.scheme,
            partition_scheme=context.partition_scheme,
            scope_mode=context.scope_mode,
            solvent_mode=restored.solvent_provenance["solvent_mode"],
            solvent_generator_digest=restored.solvent_provenance[
                "solvent_generator_digest"
            ],
            solvent_cells_complete=restored.solvent_provenance[
                "solvent_cells_complete"
            ],
            selection_id=view.selection_id,
        ),
        retained_system_indices=retention["retained_system_indices"],
        retained_generator_ids=retention["retained_generator_ids"],
        retained_generator_digest=retention["retained_generator_digest"],
    )
    dependency = _dependency(identity, source_revision=restored.archive_id)

    result = adapt_group_view_analysis(
        analysis,
        identity,
        dependency,
        group_id="retained-group",
        system="synthetic-reloaded",
        radius_configuration_id="radii-1",
        results_partition="aw",
        source_solve_context=context,
        source_solvent_provenance=restored.solvent_provenance,
        source_retained_generator_provenance=retention,
    )

    assert result.status is Status.CERTIFIED
    assert result.provenance.details["source_retained_generator_provenance"] == retention
    assert result.provenance.details["shared_network_identity"][
        "retained_generator_digest"
    ] == retention["retained_generator_digest"]


def test_shared_results_collection_keeps_distinct_selections_on_one_network():
    first_identity = _identity(selection_id="view:group-a")
    second_identity = replace(first_identity, selection_id="view:group-b")
    first = _adapt(first_identity, _dependency(first_identity), _analysis(first_identity))
    second = _adapt(
        second_identity,
        _dependency(second_identity),
        _analysis(second_identity, system_ids=(11,)),
    )

    collection = collect_shared_results(
        (first, second),
        source_network_id="archive-network-7",
        solve_universe_id="solve:context-1",
    )

    assert len(collection) == 2
    assert collection.selection_ids == ("view:group-a", "view:group-b")
    assert collection.by_selection("view:group-b") is second
    with pytest.raises(ValueError, match="distinct selection"):
        collect_shared_results(
            (first, first),
            source_network_id="archive-network-7",
            solve_universe_id="solve:context-1",
        )
