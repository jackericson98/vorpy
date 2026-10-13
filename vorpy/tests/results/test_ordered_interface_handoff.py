"""Synthetic FORGE-to-ATLAS ordered-interface handoff regressions."""

from dataclasses import replace

import pandas as pd
import pytest
import numpy as np

from vorpy.src.boundary import SolventMode
from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.io import load_network, save_network
from vorpy.src.results import (
    GeometryLineage,
    GeometryLineageOrigin,
    InterfaceResult,
    Provenance,
    ResultIdentity,
    ResultsDependency,
    SharedNetworkIdentity,
    Status,
    adapt_cached_result,
    adapt_interface_analysis,
    validate_results_dependency,
)
from vorpy.src.system import System
from vorpy.workbench.services.vorpy_backend import InterfaceSolveSelection
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


def _selection():
    return InterfaceSolveSelection(
        interface_id="chains-A-B__chain-I",
        group_a_id="chains-A-B",
        group_a_system_indices=(0, 1),
        group_b_id="chain-I",
        group_b_system_indices=(2,),
    )


def _shared_identity(*, source_network_id="producer:interface-network-1",
                     source_identity_status="AUTHORITATIVE", solvent_mode="MOLECULE_ONLY",
                     solvent_generator_digest=None, solvent_cells_complete=False):
    return SharedNetworkIdentity(
        source_network_id=source_network_id,
        source_identity_status=source_identity_status,
        solve_universe_id="solve:interface-context-1",
        frame_id="static",
        analysis_context_kind="GLOBAL_NETWORK",
        geometry_scheme="aw",
        partition_scheme="aw",
        scope_mode="interface",
        solvent_mode=solvent_mode,
        solvent_generator_digest=solvent_generator_digest,
        solvent_cells_complete=solvent_cells_complete,
        selection_id="view:chains-A-B__chain-I",
        group_A_selection_id="view:chains-A-B",
        group_B_selection_id="view:chain-I",
        orientation="A -> B",
    )


def _result_identity(selection, shared, *, network_id="archive:interface-network-1",
                     source_network_id=None, source_identity_status=None):
    return ResultIdentity(
        result_kind="interface",
        system="synthetic-interface",
        frame=shared.frame_id,
        # Archive identity remains distinct from the producer namespace.
        network_id=network_id,
        network_scope="interface",
        environment="dry",
        partition="aw",
        representation="voronoi",
        radius_configuration_id="radii:synthetic-v1",
        interface_id=selection.interface_id,
        group_A=selection.group_a_id,
        group_B=selection.group_b_id,
        side="shared",
        orientation="group1 -> group2",
        source_network_id=(
            shared.source_network_id
            if source_network_id is None
            else source_network_id
        ),
        source_identity_status=(
            shared.source_identity_status
            if source_identity_status is None
            else source_identity_status
        ),
        source_interface_id=selection.interface_id,
        solve_universe_id=shared.solve_universe_id,
        analysis_context_kind=shared.analysis_context_kind,
        selection_id=shared.selection_id,
        solvent_mode=shared.solvent_mode,
    )


def _lineage(selection, shared):
    return GeometryLineage(
        origin=GeometryLineageOrigin.DIRECT_SOLVE,
        solve_universe_id=shared.solve_universe_id,
        geometry_source_id="geometry:interface-1",
        geometry_version="aw-interface-v1",
        solver_settings_digest="settings:interface-1",
        solver_settings={"net_type": "aw"},
        environment_digest="environment:dry-1",
        environment={"network_environment": "dry"},
        source_network_id=shared.source_network_id,
        source_interface_id=selection.interface_id,
        selection_scope={
            "group_A_system_indices": selection.group_a_system_indices,
            "group_B_system_indices": selection.group_b_system_indices,
        },
        orientation={"convention": "physical_interface_A_to_B"},
        completeness={"completeness_proven": False},
    )


def test_ordered_forger_interface_maps_to_typed_result_and_dependency():
    selection = _selection()
    shared = _shared_identity()
    dependency = ResultsDependency(
        shared,
        source_id=selection.interface_id,
        source_revision="producer-interface-v1",
    )
    identity = _result_identity(selection, shared)
    lineage = _lineage(selection, shared)
    analysis = InterfaceGeometryAnalysis(
        metadata={"result_status": "PARTIAL"},
        voronoi_side_1={},
    )
    provenance = Provenance(
        "FORGE ordered interface producer",
        source_revision=dependency.source_revision,
        details={
            "shared_network_identity": shared.payload(),
            "dependency_digest": dependency.dependency_digest,
        },
    )

    result = adapt_interface_analysis(
        analysis,
        identity,
        provenance=provenance,
        geometry_lineage=lineage,
    )

    assert isinstance(result, InterfaceResult)
    assert result.identity.network_id == "archive:interface-network-1"
    assert result.identity.source_network_id == shared.source_network_id
    assert result.identity.source_identity_status == "AUTHORITATIVE"
    assert result.identity.group_A == "chains-A-B"
    assert result.identity.group_B == "chain-I"
    assert result.identity.source_interface_id == selection.interface_id
    assert result.provenance.geometry_lineage == lineage
    assert result.provenance.details["dependency_digest"] == dependency.dependency_digest
    assert result.status is Status.PARTIAL
    assert validate_results_dependency(dependency, dependency) is Status.CERTIFIED


def test_ordered_interface_selection_rejects_overlap():
    with pytest.raises(ValueError, match="overlap"):
        InterfaceSolveSelection(
            "interface-overlap",
            "chains-A-B",
            (0, 1),
            "chain-I",
            (1, 2),
        )


def test_interface_solvent_mode_changes_dependency_without_changing_solve_id():
    molecule_only = _shared_identity()
    solvent_context = replace(
        molecule_only,
        solvent_mode="SOLVENT_CONTEXT",
        solvent_generator_digest="solvent:generators-1",
        solvent_cells_complete=False,
    )
    full_solvent = replace(
        solvent_context,
        solvent_mode="FULL_SOLVENT",
        solvent_cells_complete=True,
    )
    dependencies = tuple(
        ResultsDependency(identity, "chains-A-B__chain-I", source_revision="v1")
        for identity in (molecule_only, solvent_context, full_solvent)
    )

    assert len({item.dependency_digest for item in dependencies}) == 3
    assert {
        item.identity.solve_universe_id for item in dependencies
    } == {molecule_only.solve_universe_id}
    assert validate_results_dependency(dependencies[1], dependencies[0]) is Status.STALE


def test_unresolved_producer_evidence_stays_unresolved():
    selection = _selection()
    shared = _shared_identity(source_identity_status="UNRESOLVED")
    identity = _result_identity(selection, shared)
    lineage = GeometryLineage(
        GeometryLineageOrigin.UNRESOLVED,
        reason="FORGE producer identity was not supplied",
    )
    result = adapt_cached_result(
        identity,
        provenance=Provenance(
            "FORGE ordered interface producer",
            details={"shared_network_identity": shared.payload()},
        ),
        status=Status.UNRESOLVED,
        geometry_lineage=lineage,
    )
    dependency = ResultsDependency(
        shared,
        source_id=selection.interface_id,
        source_revision="producer-interface-v1",
    )

    assert result.status is Status.UNRESOLVED
    assert result.provenance.geometry_lineage == lineage
    assert validate_results_dependency(dependency, dependency) is Status.UNRESOLVED


def _live_forger_interface():
    system = System(
        spheres=pd.DataFrame({
            "loc": [
                (0.0, 0.0, 0.0),
                (2.0, 0.2, 0.1),
                (0.2, 2.1, 0.3),
                (0.1, 0.4, 2.2),
            ],
            "rad": [1.0, 1.1, 0.9, 1.05],
            "mass": [12.0, 16.0, 18.0, 14.0],
            "num": [100, 101, 102, 103],
            "res_name": ["ALA", "GLY", "SER", "THR"],
            "name": ["CA", "O", "CA", "CA"],
        }),
        name="live-interface-handoff",
        make_dir=False,
        solvent_mode=SolventMode.MOLECULE_ONLY,
        retained_system_indices=(),
    )
    selection = InterfaceSolveSelection(
        interface_id="chains-A-B__chain-I",
        group_a_id="chains-A-B",
        group_a_system_indices=(0, 1),
        group_b_id="chain-I",
        group_b_system_indices=(2, 3),
    )
    backend = VorPyBackend(
        VorPySolveSettings(source_network_id="producer:live-interface")
    )
    interface = backend.make_interface_network(system, selection)
    return backend, selection, interface


def _prepare_archive_system(network):
    """Complete only the molecular columns required by the archive codec."""
    network.sys.balls["res"] = [None] * len(network.sys.balls)
    network.sys.balls["chn"] = [None] * len(network.sys.balls)


def test_live_forger_producer_adapts_to_typed_interface_result():
    backend, selection, interface = _live_forger_interface()
    network = interface.net
    assert network.network_mode == "interface"
    assert network.iface_grps == ({0, 1}, {2, 3})
    assert network.source_network_id == "producer:live-interface"

    interface.build(analyze_waters=False, calculate_curvature=False)
    handoff = backend.make_interface_result_handoff(
        interface,
        selection,
        InterfaceGeometryAnalysis(
            metadata={"result_status": "PARTIAL"},
            voronoi_side_1={},
        ),
        system_id="synthetic-interface",
        result_network_id="archive:live-interface",
        radius_configuration_id="radii:synthetic-v1",
        selection_id="view:chains-A-B__chain-I",
        group_a_selection_id="view:chains-A-B",
        group_b_selection_id="view:chain-I",
        orientation="A -> B",
        source_identity_status="AUTHORITATIVE",
        source_revision="producer:live-interface-v1",
        environment="dry",
        completeness_proven=False,
    )
    result = handoff.result
    dependency = handoff.dependency
    lineage = handoff.geometry_lineage

    assert result.identity.source_network_id == network.source_network_id
    assert result.identity.network_id == "archive:live-interface"
    assert result.identity.source_interface_id == selection.interface_id
    assert result.identity.group_A == selection.group_a_id
    assert result.identity.group_B == selection.group_b_id
    assert result.identity.selection_id == "view:chains-A-B__chain-I"
    assert result.provenance.geometry_lineage == lineage
    assert result.status is Status.PARTIAL
    assert validate_results_dependency(dependency, dependency) is Status.CERTIFIED


def test_archived_forger_interface_re_adaptation_preserves_lineage_and_status(tmp_path):
    backend, selection, interface = _live_forger_interface()
    interface.build(analyze_waters=False, calculate_curvature=False)
    _prepare_archive_system(interface.net)
    archive = save_network(interface.net, tmp_path / "live-interface.vpy")
    restored = load_network(archive)
    restored_interface = next(
        item for item in restored.sys.ifaces
        if item.interface_id == selection.interface_id
    )
    assert restored.source_network_id == "producer:live-interface"
    assert restored_interface.net is restored
    assert restored.iface_grps == ({0, 1}, {2, 3})

    handoff = backend.make_interface_result_handoff(
        restored_interface,
        selection,
        InterfaceGeometryAnalysis(
            metadata={"result_status": "UNRESOLVED"},
            voronoi_side_1={},
        ),
        system_id="synthetic-interface",
        result_network_id=restored.archive_id,
        radius_configuration_id="radii:synthetic-v1",
        selection_id="view:chains-A-B__chain-I",
        group_a_selection_id="view:chains-A-B",
        group_b_selection_id="view:chain-I",
        orientation="A -> B",
        source_identity_status="AUTHORITATIVE",
        source_revision=restored.archive_id,
        environment="dry",
        completeness_proven=False,
    )
    result = handoff.result
    lineage = handoff.geometry_lineage

    assert result.identity.network_id == restored.archive_id
    assert result.identity.source_network_id == restored.source_network_id
    assert result.provenance.source_revision == restored.archive_id
    assert result.provenance.geometry_lineage == lineage
    assert result.status is Status.UNRESOLVED
