"""Synthetic FORGE-to-ATLAS ordered-interface handoff regressions."""

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import numpy as np

from vorpy.src.boundary import SolventMode
from vorpy.src.interface import Interface
from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.io import load_network, save_network
from vorpy.tests.network.test_views import _network
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
from vorpy.src.results.model import stable_id
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
        metadata={"result_status": "PARTIAL", "curvature_requested": True},
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


def _live_forger_interface(
    *, network_type="aw", source_network_id="producer:live-interface"
):
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
        VorPySolveSettings(
            network_type=network_type,
            source_network_id=source_network_id,
        )
    )
    interface = backend.make_interface_network(system, selection)
    return backend, selection, interface


def test_dry_pow_interface_result_acceptance_contract():
    """The dry POW pilot accepts only explicit producer and evidence fields."""
    backend, selection, interface = _live_forger_interface(
        network_type="pow",
        source_network_id="producer:dry-pow-interface",
    )
    network = interface.net
    context = network.get_solve_context_key()
    analysis = InterfaceGeometryAnalysis(
        metadata={
            "result_status": "PARTIAL",
            "system_id": "dry-pow-interface",
            "network_id": "archive:dry-pow-interface",
            "network_scope": "interface",
            "frame": context.frame_id,
            "scheme": "pow",
            "radius_configuration_id": "radii:dry-pow-v1",
            "solvent_context": False,
            "curvature_requested": True,
        },
        alpha_selection={"selected_pairs": 3},
        voronoi_side_1={
            "area": 12.5,
            "area_selection_complete": True,
            "topology_vertices": 4,
            "topology_edges": 5,
            "topology_faces": 3,
            "gaussian_curvature_certified": False,
            "gaussian_accounting_scope": "full_selected_interface",
            "total_K": None,
            "gauss_bonnet_lhs": None,
            "gauss_bonnet_expected": None,
            "gauss_bonnet_residual": None,
            "euler_characteristic": 2,
            "partial_integrated_mean_curvature": 1.25,
            "mean_curvature_certified": False,
            "mean_supported_faces": 2,
            "mean_unsupported_faces": 1,
            "smooth_K_partial": 0.5,
            "gaussian_curvature_certified": False,
            "gaussian_accounting_scope": "supported_subcomplex",
            "gaussian_supported_faces": 2,
            "gaussian_unsupported_faces": 1,
            # Present-but-null fields are unresolved evidence, not zeroes.
            "total_K": None,
            "gauss_bonnet_lhs": None,
            "gauss_bonnet_expected": None,
            "gauss_bonnet_residual": None,
        },
        unresolved=["Gaussian coverage is incomplete for the selected interface"],
    )
    handoff = backend.make_interface_result_handoff(
        interface,
        selection,
        analysis,
        system_id="dry-pow-interface",
        result_network_id="archive:dry-pow-interface",
        radius_configuration_id="radii:dry-pow-v1",
        selection_id="view:chains-A-B__chain-I",
        group_a_selection_id="view:chains-A-B",
        group_b_selection_id="view:chain-I",
        orientation="A -> B",
        source_identity_status="AUTHORITATIVE",
        source_revision="producer:dry-pow-interface-v1",
        environment="dry",
        completeness_proven=False,
    )
    result = handoff.result
    identity = handoff.identity
    dependency = handoff.dependency
    lineage = handoff.geometry_lineage

    # Producer identity and archive identity are separate namespaces.
    assert network.source_network_id == "producer:dry-pow-interface"
    assert identity.source_network_id == "producer:dry-pow-interface"
    assert identity.source_identity_status == "AUTHORITATIVE"
    assert identity.authoritative_source_network
    assert result.identity.network_id == "archive:dry-pow-interface"
    assert result.identity.network_id != result.identity.source_network_id

    # A/B order is preserved in both the dependency and typed Result identity.
    assert network.iface_grps == ({0, 1}, {2, 3})
    assert identity.group_A_selection_id == "view:chains-A-B"
    assert identity.group_B_selection_id == "view:chain-I"
    assert identity.orientation == "A -> B"
    assert result.identity.group_A == selection.group_a_id
    assert result.identity.group_B == selection.group_b_id
    assert result.identity.orientation == "group1 -> group2"
    assert result.identity.source_interface_id == selection.interface_id

    # The dependency is tied to the POW solve context, not merely the view.
    assert identity.geometry_scheme == "pow"
    assert identity.partition_scheme == "power"
    assert identity.scope_mode == "interface"
    assert result.identity.frame == context.frame_id
    assert result.identity.solve_universe_id == identity.solve_universe_id
    assert result.identity.partition == "power"
    assert identity.solvent_mode == "MOLECULE_ONLY"
    assert dependency.identity == identity
    assert dependency.source_id == selection.interface_id
    assert dependency.source_revision == "producer:dry-pow-interface-v1"
    assert dependency.dependency_digest
    assert validate_results_dependency(dependency, dependency) is Status.CERTIFIED
    assert result.provenance.details["dependency_digest"] == dependency.dependency_digest

    # Lineage is direct producer evidence with explicit ordered scope.
    assert lineage.origin is GeometryLineageOrigin.DIRECT_SOLVE
    assert lineage.solve_universe_id == identity.solve_universe_id
    assert lineage.source_network_id == identity.source_network_id
    assert lineage.source_interface_id == selection.interface_id
    assert lineage.selection_scope["group_A_system_indices"] == (0, 1)
    assert lineage.selection_scope["group_B_system_indices"] == (2, 3)
    assert lineage.selection_scope["selection_id"] == "view:chains-A-B__chain-I"
    assert lineage.orientation["convention"] == "physical_interface_A_to_B"
    assert lineage.completeness["completeness_proven"] is False
    assert result.provenance.geometry_lineage == lineage

    # Surface and geometry quantities remain typed and retain their evidence status.
    assert result.status is Status.PARTIAL
    assert result.metrics["area"].value == pytest.approx(12.5)
    assert result.metrics["area"].units == "Å²"
    assert result.metrics["area"].status is Status.CERTIFIED
    assert result.topology.quantities["vertices"].value == 4
    assert result.topology.quantities["edges"].value == 5
    assert result.topology.quantities["faces"].value == 3
    assert result.mean_curvature.integrated_mean_curvature.value == pytest.approx(1.25)
    assert result.mean_curvature.integrated_mean_curvature.status is Status.PARTIAL
    assert result.metrics["mean_unsupported_faces"].value == 1
    assert result.metrics["mean_unsupported_faces"].status is Status.PARTIAL
    assert result.gaussian_curvature.smooth_face_contribution.value == pytest.approx(0.5)
    assert result.gaussian_curvature.smooth_face_contribution.status is Status.PARTIAL

    # Requested but incomplete curvature cannot be promoted or replaced by zero.
    gaussian = result.gaussian_curvature
    assert gaussian.integrated_gaussian_curvature.value is None
    assert gaussian.integrated_gaussian_curvature.status is Status.UNRESOLVED
    assert gaussian.integrated_gaussian_curvature.scope == "full_selected_interface"
    assert gaussian.gauss_bonnet_completion.value is None
    assert gaussian.gauss_bonnet_completion.status is Status.UNRESOLVED
    assert gaussian.gauss_bonnet_residual.value is None
    assert gaussian.gauss_bonnet_residual.status is Status.UNRESOLVED
    assert result.extensions[0].payload["unresolved"] == (
        "Gaussian coverage is incomplete for the selected interface",
    )


def test_incomplete_source_view_keeps_ordered_partial_geometry_typed():
    """A 191-facet-style incomplete view is not an A/B membership failure."""
    source = _network(scheme="pow", source_network_id="producer:1bna-pow-dry")
    source.network_mode = "complete"
    source.solvent_mode = SolventMode.MOLECULE_ONLY.value
    source.settings["solvent_mode"] = source.solvent_mode
    source.balls.loc[2, "complete"] = False
    extra = source.surfs.iloc[1].copy()
    extra["balls"] = [0, 1, 2]
    extra["sa"] = 1.0
    source.surfs.loc[2] = extra
    view = source.make_view(iface_grps=([10], [11]))
    assert view.group_a_ids == (10,)
    assert view.group_b_ids == (11,)
    assert view.completeness.proven is False

    group_a = SimpleNamespace(
        name="group-A", group_id="group-A", ball_ndxs=[10],
        settings=dict(source.settings), interfaces={},
    )
    group_b = SimpleNamespace(
        name="group-B", group_id="group-B", ball_ndxs=[11],
        settings=dict(source.settings), interfaces={},
    )
    def register_for(group):
        return lambda **kwargs: group.interfaces.__setitem__(
            kwargs["interface_id"], kwargs["interface"]
        )

    group_a.register_interface = register_for(group_a)
    group_b.register_interface = register_for(group_b)
    system = SimpleNamespace(
        balls=source.balls, groups=[group_a, group_b], ifaces=[],
        files={"dir": None}, frame_index=7,
    )
    interface = Interface(system, group_a, group_b)
    selection = InterfaceSolveSelection(
        interface.interface_id, "group-A", (10,), "group-B", (11,)
    )
    interface.net = source
    interface.source_network = source
    interface.source_view = view
    interface.geometry_analysis = InterfaceGeometryAnalysis(
        metadata={"result_status": "PARTIAL", "curvature_requested": False},
        voronoi_side_1={
            "area": 179.944,
            "area_selection_complete": False,
            "area_scope": "supported source facets; unresolved facets omitted",
            "topology_faces": 191,
            "topology_edges": 512,
            "topology_vertices": 322,
            "partial_integrated_mean_curvature": 4.0,
            "mean_curvature_certified": False,
            "gaussian_curvature_certified": False,
        },
        unresolved=["source generator 12 cell completeness is unproven"],
    )
    handoff = VorPyBackend().make_interface_result_handoff(
        interface, selection, interface.geometry_analysis,
        system_id="1bna-pow-dry",
        result_network_id="archive:1bna-pow-dry-interface",
        radius_configuration_id="radii:synthetic-v1",
        selection_id=view.selection_id,
        group_a_selection_id="view:group-A",
        group_b_selection_id="view:group-B",
        orientation="A -> B",
        source_identity_status="AUTHORITATIVE",
        source_revision="producer:1bna-pow-dry:v1",
        environment="dry",
        completeness_proven=False,
    )

    result = handoff.result
    assert result.status is Status.PARTIAL
    assert result.provenance.geometry_lineage.origin is GeometryLineageOrigin.MATERIALIZED_VIEW
    assert result.provenance.geometry_lineage.completeness["completeness_proven"] is False
    assert result.identity.group_A == "group-A"
    assert result.identity.group_B == "group-B"
    assert result.metrics["area"].value == pytest.approx(179.944)
    assert result.metrics["area"].status is Status.PARTIAL
    assert result.extensions[0].payload["physical_accounting"]["topology_faces"] == 191
    assert result.extensions[0].payload["physical_accounting"]["topology_edges"] == 512
    assert result.extensions[0].payload["physical_accounting"]["topology_vertices"] == 322
    assert result.mean_curvature.integrated_mean_curvature.value == pytest.approx(4.0)
    assert result.mean_curvature.integrated_mean_curvature.status is Status.PARTIAL
    assert result.gaussian_curvature.integrated_gaussian_curvature.value is None
    assert result.gaussian_curvature.integrated_gaussian_curvature.status is Status.NOT_CALCULATED


def test_2kai_pow_cli_adapter_uses_minimum_live_handoff_without_solving(monkeypatch):
    """The CLI-shaped producer path adapts one ordered 2KAI POW result."""
    source = Path("2KAI_MD_centered_no_sol_corrected.pdb")
    backend, selection, interface = _live_forger_interface(
        network_type="pow",
        source_network_id="producer:2KAI:pow:dry:A+B-I",
    )
    context = interface.net.get_solve_context_key()
    interface.geometry_analysis = InterfaceGeometryAnalysis(
        metadata={
            "result_status": "PARTIAL",
            "scheme": "pow",
            "frame": context.frame_id,
            "solvent_context": False,
            "curvature_requested": True,
        },
        voronoi_side_1={
            "area": 123.0,
            "area_selection_complete": True,
            "topology_vertices": 11,
            "topology_edges": 22,
            "topology_faces": 12,
            "partial_integrated_mean_curvature": 3.0,
            "total_integrated_mean_curvature": None,
            "mean_curvature_certified": False,
            "gaussian_curvature_certified": False,
            "gaussian_accounting_scope": "full_selected_interface",
            "total_K": None,
            "gauss_bonnet_lhs": None,
            "gauss_bonnet_expected": None,
            "gauss_bonnet_residual": None,
        },
        unresolved=["POW interface edge/vertex coverage is incomplete"],
    )

    # This is the live producer-to-adapter call used by the Workbench worker;
    # the view is stubbed only because view construction requires solved
    # topology, which this adapter-contract test deliberately avoids.
    monkeypatch.setattr(
        interface.net,
        "make_view",
        lambda **_kwargs: SimpleNamespace(selection_id="view:2kai-pow-dry-interface"),
    )
    typed = backend._interface_typed_results(interface, selection, source, 1)

    assert len(typed) == 1
    result = typed[0]
    expected_archive_id = stable_id(
        "archive", {"source": str(source), "interface_id": selection.interface_id}
    )
    expected_group_a_selection_id = stable_id(
        "view", {"group_id": selection.group_a_id, "system_ids": (0, 1)}
    )
    expected_group_b_selection_id = stable_id(
        "view", {"group_id": selection.group_b_id, "system_ids": (2, 3)}
    )
    assert result.identity.network_id == expected_archive_id
    assert result.identity.network_id.startswith("archive:")
    assert result.identity.source_network_id == "producer:2KAI:pow:dry:A+B-I"
    assert result.identity.source_identity_status == "AUTHORITATIVE"
    assert result.identity.network_id != result.identity.source_network_id

    assert result.identity.group_A == "chains-A-B"
    assert result.identity.group_B == "chain-I"
    assert result.identity.orientation == "group1 -> group2"
    assert result.identity.selection_id == "view:2kai-pow-dry-interface"
    assert result.identity.partition == "power"
    assert result.identity.solve_universe_id == context.solve_universe_id
    assert result.identity.solvent_mode == "MOLECULE_ONLY"

    details = result.provenance.details
    shared = details["shared_network_identity"]
    assert shared["source_network_id"] == result.identity.source_network_id
    assert shared["geometry_scheme"] == "pow"
    assert shared["partition_scheme"] == "power"
    assert shared["solve_universe_id"] == context.solve_universe_id
    assert shared["scope_mode"] == "interface"
    assert shared["solvent_mode"] == "MOLECULE_ONLY"
    assert shared["selection_id"] == "view:2kai-pow-dry-interface"
    assert shared["group_A_selection_id"] == expected_group_a_selection_id
    assert shared["group_B_selection_id"] == expected_group_b_selection_id
    assert shared["orientation"] == "A -> B"
    assert details["dependency_digest"].startswith("dependency-v1:")
    assert result.provenance.source_revision == "interface:chains-A-B__chain-I"

    lineage = result.provenance.geometry_lineage
    assert lineage.origin is GeometryLineageOrigin.DIRECT_SOLVE
    assert lineage.source_network_id == result.identity.source_network_id
    assert lineage.source_interface_id == selection.interface_id
    assert lineage.selection_scope["group_A_system_indices"] == (0, 1)
    assert lineage.selection_scope["group_B_system_indices"] == (2, 3)
    assert lineage.selection_scope["selection_id"] == "view:2kai-pow-dry-interface"
    assert lineage.orientation["convention"] == "physical_interface_A_to_B"
    assert lineage.completeness["completeness_proven"] is False

    assert result.status is Status.PARTIAL
    assert result.metrics["area"].value == pytest.approx(123.0)
    assert result.metrics["area"].units == "Å²"
    assert result.metrics["area"].status is Status.CERTIFIED
    assert result.topology.quantities["vertices"].value == 11
    assert result.topology.quantities["edges"].value == 22
    assert result.topology.quantities["faces"].value == 12
    assert result.mean_curvature.integrated_mean_curvature.value == pytest.approx(3.0)
    assert result.mean_curvature.integrated_mean_curvature.status is Status.PARTIAL
    assert result.mean_curvature.complete_total.value is None
    assert result.mean_curvature.complete_total.status is Status.UNRESOLVED
    assert result.gaussian_curvature.integrated_gaussian_curvature.value is None
    assert result.gaussian_curvature.integrated_gaussian_curvature.status is Status.UNRESOLVED
    assert result.extensions[0].payload["unresolved"] == (
        "POW interface edge/vertex coverage is incomplete",
    )


def _prepare_archive_system(network):
    """Complete only the molecular columns required by the archive codec."""
    network.sys.balls["res"] = [None] * len(network.sys.balls)
    network.sys.balls["chn"] = [None] * len(network.sys.balls)


def _archiveable_typed_interface():
    """Create a complete archive graph without invoking a molecular solve."""
    backend, selection, interface = _live_forger_interface(
        network_type="pow",
        source_network_id="producer:typed-interface-roundtrip",
    )
    network = interface.net
    network.archive_id = "archive:typed-interface-roundtrip"
    analysis = InterfaceGeometryAnalysis(
        metadata={"result_status": "PARTIAL", "curvature_requested": True},
        voronoi_side_1={
            "area": 12.5,
            "area_selection_complete": True,
            "topology_vertices": 4,
            "topology_edges": 5,
            "topology_faces": 3,
            "gaussian_curvature_certified": False,
            "gaussian_accounting_scope": "full_selected_interface",
            "total_K": None,
            "gauss_bonnet_lhs": None,
            "gauss_bonnet_expected": None,
            "gauss_bonnet_residual": None,
        },
        unresolved=["curvature coverage is incomplete"],
    )
    handoff = backend.make_interface_result_handoff(
        interface,
        selection,
        analysis,
        system_id="synthetic-interface",
        result_network_id=network.archive_id,
        radius_configuration_id="radii:synthetic-v1",
        selection_id="view:chains-A-B__chain-I",
        group_a_selection_id="view:chains-A-B",
        group_b_selection_id="view:chain-I",
        orientation="A -> B",
        source_identity_status="AUTHORITATIVE",
        source_revision="producer:typed-interface-roundtrip-v1",
        environment="dry",
        completeness_proven=False,
    )
    interface.typed_results = (handoff.result,)
    interface.results_dependency = handoff.dependency

    # The archive validator needs solved-shaped tables; empty tables are a
    # bounded fixture and do not claim any geometric measurements.
    network.verts = pd.DataFrame({
        "loc": pd.Series(dtype=object),
        "rad": pd.Series(dtype=float),
        "balls": pd.Series(dtype=object),
    })
    network.edges = pd.DataFrame({
        "balls": pd.Series(dtype=object),
        "verts": pd.Series(dtype=object),
        "points": pd.Series(dtype=object),
    })
    network.surfs = pd.DataFrame({
        "balls": pd.Series(dtype=object),
        "verts": pd.Series(dtype=object),
        "edges": pd.Series(dtype=object),
        "points": pd.Series(dtype=object),
        "tris": pd.Series(dtype=object),
    })
    _prepare_archive_system(network)
    return interface, handoff


def test_typed_interface_result_vpy_round_trip_preserves_identity_evidence(tmp_path):
    interface, handoff = _archiveable_typed_interface()
    archive = save_network(interface.net, tmp_path / "typed-interface.vpy")

    restored = load_network(archive)
    restored_interface = next(
        item for item in restored.sys.ifaces
        if item.interface_id == handoff.result.identity.interface_id
    )
    restored_result = restored_interface.typed_results[0]
    restored_dependency = restored_interface.results_dependency

    assert restored_result.to_dict() == handoff.result.to_dict()
    assert restored_result.identity.network_id == restored.archive_id
    assert restored_result.identity.network_id != restored_result.identity.source_network_id
    assert restored_result.identity.source_network_id == restored.source_network_id
    assert restored_result.identity.group_A == "chains-A-B"
    assert restored_result.identity.group_B == "chain-I"
    assert restored_result.identity.orientation == "group1 -> group2"
    assert restored_dependency.identity.group_A_selection_id == "view:chains-A-B"
    assert restored_dependency.identity.group_B_selection_id == "view:chain-I"
    assert restored_dependency.identity.orientation == "A -> B"
    assert restored_dependency.status is Status.CERTIFIED
    assert restored_result.provenance.geometry_lineage == handoff.geometry_lineage
    assert restored_result.provenance.geometry_lineage.completeness == {
        "completeness_proven": False,
    }
    assert restored_result.metrics["area"].units == "Å²"
    assert restored_result.metrics["area"].value == pytest.approx(12.5)
    assert restored_result.status is Status.PARTIAL
    assert handoff.result.gaussian_curvature.integrated_gaussian_curvature.status is Status.UNRESOLVED
    assert restored_result.gaussian_curvature.integrated_gaussian_curvature.status is Status.UNRESOLVED


def test_typed_interface_result_vpy_archive_missing_dependency_fails_closed(tmp_path):
    interface, handoff = _archiveable_typed_interface()
    interface.results_dependency = None
    archive = save_network(interface.net, tmp_path / "typed-interface-no-dependency.vpy")

    restored = load_network(archive)
    result = restored.sys.ifaces[0].typed_results[0]

    assert result.status is Status.UNRESOLVED
    assert result.metrics["area"].status is Status.CERTIFIED
    assert result.metrics["area"].value == pytest.approx(12.5)
    assert result.provenance.geometry_lineage == handoff.geometry_lineage


def test_typed_interface_result_vpy_archive_missing_lineage_fails_closed(tmp_path):
    interface, handoff = _archiveable_typed_interface()
    details = dict(handoff.result.provenance.details)
    details.pop("geometry_lineage")
    interface.typed_results = (replace(
        handoff.result,
        provenance=Provenance(
            handoff.result.provenance.source,
            source_revision=handoff.result.provenance.source_revision,
            scientific_state_id=handoff.result.provenance.scientific_state_id,
            details=details,
            diagnostics=handoff.result.provenance.diagnostics,
        ),
    ),)
    archive = save_network(interface.net, tmp_path / "typed-interface-no-lineage.vpy")

    restored = load_network(archive)
    result = restored.sys.ifaces[0].typed_results[0]

    assert result.status is Status.UNRESOLVED
    assert result.provenance.geometry_lineage is None


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
