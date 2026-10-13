import os
from dataclasses import replace
from types import SimpleNamespace

import pandas as pd

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from vorpy.src.boundary import BoundaryConfig
from vorpy.src.network import Network
from vorpy.src.group.shared_view import GroupViewAnalysis
from vorpy.src.results import (
    InterfaceResult,
    Provenance,
    ResultIdentity,
    ResultsDependency,
    SharedNetworkIdentity,
    Status,
    adapt_group_view_analysis,
)
from vorpy.workbench.ui.panels import (
    AnalysisCatalogEntryPresentation,
    AnalysisCatalogPanel,
    AnalysisCatalogPresentation,
    catalog_presentation_from_results,
    context_presentation_from_owner,
    context_presentation_from_results,
    NetworkContextPanel,
    NetworkContextPresentation,
    validated_group_system_ids,
)


def app():
    return QApplication.instance() or QApplication([])


def test_network_context_empty_state_does_not_claim_group_analysis():
    app()
    panel = NetworkContextPanel()

    assert panel.identity_badge.text() == "No network loaded"
    assert all(value.text() == "—" for value in panel.fields.values())
    assert "not analyzed network contexts" in panel.note.text()


def test_network_context_renders_validated_global_identity():
    app()
    panel = NetworkContextPanel()
    panel.set_context(NetworkContextPresentation(
        source_network_id="network-7",
        solve_universe_id="solve-9",
        archive_identity="archive-7",
        geometry_lineage="DIRECT_SOLVE · geometry-7",
        geometry_scheme="AW",
        frame="12",
        context_scope="global_network",
        solvent_mode="molecule_only",
        solvent_mode_state="validated",
        retained_generators="32 retained (explicit policy)",
        completeness="CERTIFIED",
        status="CERTIFIED",
        identity_state="validated",
    ))

    assert panel.identity_badge.text() == "Validated identity"
    assert panel.fields["source_network_id"].text() == "network-7"
    assert panel.fields["solve_universe_id"].text() == "solve-9"
    assert panel.fields["archive_identity"].text() == "archive-7"
    assert panel.fields["geometry_lineage"].text() == "DIRECT_SOLVE · geometry-7"
    assert panel.fields["geometry_scheme"].text() == "AW"
    assert panel.fields["frame"].text() == "12"
    assert panel.fields["context_scope"].text() == "Global network"
    assert panel.fields["solvent_mode"].text() == "Molecule-only"
    assert panel.fields["retained_generators"].text() == "32 retained (explicit policy)"
    assert panel.fields["completeness"].text() == "CERTIFIED"
    assert panel.fields["status"].text() == "CERTIFIED"


def test_network_context_marks_provisional_and_unresolved_identity():
    app()
    panel = NetworkContextPanel()
    panel.set_context(NetworkContextPresentation(
        geometry_scheme="Power",
        frame="3",
        context_scope="unknown",
        solvent_mode="full_solvent",
        solvent_mode_state="legacy",
        completeness="UNKNOWN",
        status="UNRESOLVED",
        identity_state="provisional",
        note="Identity was supplied without a validated namespace.",
    ))
    assert panel.identity_badge.text() == "Provisional identity"
    assert panel.fields["source_network_id"].text() == "—"
    assert panel.fields["context_scope"].text() == "Unknown"
    assert panel.fields["solvent_mode"].text() == "Legacy / unresolved (Full solvent)"
    assert panel.fields["status"].text() == "UNRESOLVED"

    panel.set_context(NetworkContextPresentation(
        status="UNRESOLVED",
        completeness="UNRESOLVED",
        identity_state="unresolved",
    ))
    assert panel.identity_badge.text() == "Identity unresolved"
    assert panel.fields["solvent_mode"].text() == "Unknown / unresolved"
    assert panel.fields["status"].text() == "UNRESOLVED"


def test_analysis_catalog_keeps_two_selections_on_one_network_distinct():
    app()
    panel = AnalysisCatalogPanel()
    panel.set_catalog(AnalysisCatalogPresentation(
        source_network_id="network-7",
        identity_state="validated",
        entries=(
            AnalysisCatalogEntryPresentation(
                selection_label="Group: protein",
                selection_id="view:group-protein",
                analysis_id="analysis:protein",
                context_scope="global_network",
                completeness="COMPLETE",
                status="CERTIFIED",
                area="12.5 Å²",
                volume="31.0 Å³",
            ),
            AnalysisCatalogEntryPresentation(
                selection_label="Interface: ligand/contact",
                selection_id="view:interface-contact",
                analysis_id="analysis:contact",
                context_scope="global_network",
                completeness="COMPLETE",
                status="CERTIFIED",
                area="7.5 Å²",
                volume="9.0 Å³",
            ),
        ),
    ))

    assert panel.network_label.text() == "Network: network-7"
    assert panel.table.rowCount() == 2
    assert panel.table.item(0, 1).text() == "view:group-protein"
    assert panel.table.item(1, 1).text() == "view:interface-contact"
    assert panel.table.item(0, 2).text() != panel.table.item(1, 2).text()
    assert panel.table.item(0, 3).text() == "Global network"
    assert "boundary-dependent" in panel.table.horizontalHeaderItem(6).text()
    panel.table.cellClicked.emit(1, 0)
    assert panel.active_entry.selection_id == "view:interface-contact"


def test_analysis_catalog_empty_state_does_not_claim_analysis():
    app()
    panel = AnalysisCatalogPanel()
    panel.set_catalog(AnalysisCatalogPresentation(
        source_network_id="network-7",
        identity_state="provisional",
    ))

    assert panel.network_label.text() == "Provisional network: network-7"
    assert panel.table.rowCount() == 0
    assert "No analyses associated" in panel.empty.text()
    assert "producer-attested" in panel.empty.text()


def test_context_adapter_consumes_validated_f10_solvent_evidence():
    network = SimpleNamespace(
        source_network_id="producer-network-7",
        iface_grps=None,
        solvent_provenance={
            "solvent_mode": "SOLVENT_CONTEXT",
            "solvent_generator_ids": ("system:1",),
            "solvent_generator_digest": "digest-1",
            "solvent_cells_complete": False,
            "generator_membership": ((0, 0, "system:0", False, False),
                                      (1, 1, "system:1", True, False)),
        },
        retained_generator_provenance={
            "policy": "explicit_system_indices",
            "retained_system_indices": (0,),
            "retained_generator_ids": ("system:0",),
            "retained_generator_digest": "retained-1",
            "boundary_independent": True,
        },
        get_solve_context_key=lambda **_: SimpleNamespace(
            digest="context-7",
            solve_universe_id="solve:context-7",
            scheme="aw",
            frame_id="static",
            scope_mode="full",
        ),
    )

    context = context_presentation_from_owner(SimpleNamespace(net=network))

    assert context.identity_state == "validated"
    assert context.source_network_id == "producer-network-7"
    assert context.solve_universe_id == "solve:context-7"
    assert context.context_scope == "global_network"
    assert context.solvent_mode == "SOLVENT_CONTEXT"
    assert context.solvent_mode_state == "validated"
    assert context.retained_generators == "1 retained (explicit policy)"


def test_context_adapter_reads_actual_network_producer_metadata_without_building():
    network = Network(
        locs=[(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        rads=[1.5, 1.5],
        group=[0, 1],
        settings={
            "net_type": "aw",
            "surf_res": 0.2,
            "solvent_mode": "SOLVENT_CONTEXT",
        },
        balls=pd.DataFrame({
            "loc": [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
            "rad": [1.5, 1.5],
            "num": [0, 1],
        }),
        system=SimpleNamespace(frame_index=4),
        source_network_id="producer-network-7",
        generator_system_indices=(0, 1),
        generator_stable_ids=("system:0", "system:1"),
        solvent_local_indices=(1,),
        solvent_generator_ids=("system:1",),
        solvent_mode="SOLVENT_CONTEXT",
    )
    network.boundary_config = BoundaryConfig("none")

    context = context_presentation_from_owner(SimpleNamespace(net=network))

    assert context.identity_state == "validated"
    assert context.solve_universe_id.startswith("solve:")
    assert context.geometry_scheme == "aw"
    assert context.frame == "4"
    assert context.solvent_mode_state == "validated"


def test_context_adapter_does_not_promote_archive_id_to_source_identity():
    network = SimpleNamespace(
        source_network_id=None,
        archive_id="archive-object-7",
        iface_grps=None,
        solvent_provenance=None,
        solvent_mode=None,
        settings={},
        get_solve_context_key=lambda **_: SimpleNamespace(
            digest="context-archive-only",
            solve_universe_id="solve:context-archive-only",
            scheme="aw",
            frame_id="static",
            scope_mode="full",
        ),
    )

    context = context_presentation_from_owner(SimpleNamespace(net=network))

    assert context.source_network_id is None
    assert context.archive_identity == "archive-object-7"
    assert context.identity_state == "unresolved"
    assert context.status == "UNRESOLVED"
    assert "Archive object identity is not used" in context.note


def test_context_adapter_fails_closed_for_legacy_solvent_metadata():
    network = SimpleNamespace(
        source_network_id=None,
        iface_grps=None,
        solvent_provenance={"solvent_mode": "FULL_SOLVENT"},
        get_solve_context_key=None,
    )

    context = context_presentation_from_owner(SimpleNamespace(net=network))

    assert context.identity_state == "unresolved"
    assert context.solvent_mode == "FULL_SOLVENT"
    assert context.solvent_mode_state == "unresolved"
    assert "incomplete" in context.note


def test_context_adapter_labels_raw_legacy_solvent_mode():
    network = SimpleNamespace(
        source_network_id=None,
        iface_grps=None,
        solvent_mode="FULL_SOLVENT",
        settings={},
        solvent_provenance=None,
        get_solve_context_key=None,
    )

    context = context_presentation_from_owner(SimpleNamespace(net=network))
    panel = NetworkContextPanel()
    panel.set_context(context)

    assert context.solvent_mode_state == "legacy"
    assert panel.fields["solvent_mode"].text() == "Legacy / unresolved (Full solvent)"


def _typed_group_result(selection_id, group_id, system_id, area, volume):
    identity = SharedNetworkIdentity(
        source_network_id="producer-network-7",
        source_identity_status="AUTHORITATIVE",
        solve_universe_id="solve:context-1",
        frame_id="static",
        analysis_context_kind="GLOBAL_NETWORK",
        geometry_scheme="aw",
        partition_scheme="aw",
        scope_mode="full",
        selection_id=selection_id,
        solvent_mode="MOLECULE_ONLY",
        solvent_generator_digest=None,
        solvent_cells_complete=False,
    )
    dependency = ResultsDependency(
        identity=identity,
        source_id=f"view-cache-{group_id}",
        source_revision="producer-1",
    )
    analysis = GroupViewAnalysis(
        source_network_id=identity.source_network_id,
        solve_universe_id=identity.solve_universe_id,
        selection_id=selection_id,
        system_ids=(system_id,),
        boundary_surface_ids=(f"surface-{group_id}",),
        competitor_system_ids=(99,),
        area=area,
        volume=volume,
        complete=True,
    )
    return adapt_group_view_analysis(
        analysis,
        identity,
        dependency,
        group_id=group_id,
        system="synthetic",
        radius_configuration_id="radii-1",
        source_solve_context=SimpleNamespace(
            solve_universe_id="solve:context-1",
            frame_id="static",
            scheme="aw",
            partition_scheme="aw",
            scope_mode="full",
            solvent_mode="MOLECULE_ONLY",
            solvent_generator_digest=None,
            solvent_cells_complete=False,
        ),
        source_solvent_provenance={
            "solvent_mode": "MOLECULE_ONLY",
            "solvent_generator_digest": None,
            "solvent_cells_complete": False,
        },
        source_retained_generator_provenance={
            "policy": "explicit_system_indices",
            "retained_system_indices": (),
            "retained_generator_ids": (),
            "retained_generator_digest": None,
            "boundary_independent": True,
        },
    )


def test_catalog_and_context_adapters_consume_validated_typed_results():
    results = (
        _typed_group_result("view:group-a", "group-A", 10, 4.5, 2.0),
        _typed_group_result("view:group-b", "group-B", 11, 6.5, 3.0),
    )

    catalog = catalog_presentation_from_results(results)
    context = context_presentation_from_results(results)

    assert catalog.identity_state == "validated"
    assert catalog.source_network_id == "producer-network-7"
    assert [entry.selection_id for entry in catalog.entries] == [
        "view:group-a", "view:group-b"
    ]
    assert catalog.entries[0].area.startswith("4.5")
    assert catalog.entries[1].volume.startswith("3")
    assert context.identity_state == "validated"
    assert context.context_scope == "global_network"
    assert context.solve_universe_id == "solve:context-1"
    assert context.solvent_mode == "MOLECULE_ONLY"
    assert context.solvent_mode_state == "validated"
    assert context.completeness == "PROVEN"
    assert context.status == "CERTIFIED"


def test_catalog_preserves_stale_result_status_and_boundary_area_label():
    app()
    stale = replace(
        _typed_group_result("view:group-stale", "group-stale", 10, 4.5, 2.0),
        status=Status.STALE,
    )

    catalog = catalog_presentation_from_results((stale,))
    panel = AnalysisCatalogPanel()
    panel.set_catalog(catalog)

    assert catalog.identity_state == "validated"
    assert catalog.entries[0].status == "STALE"
    assert "boundary-dependent" in panel.table.horizontalHeaderItem(6).text()


def test_catalog_presents_group_and_interface_results_from_one_network():
    group = _typed_group_result("view:group-a", "group-A", 10, 4.5, 2.0)
    shared_identity = {
        "source_network_id": "producer-network-7",
        "solve_universe_id": "solve:context-1",
        "frame_id": "static",
        "analysis_context_kind": "GLOBAL_NETWORK",
        "geometry_scheme": "aw",
        "partition_scheme": "aw",
        "scope_mode": "full",
        "selection_id": "view:interface-ab",
        "solvent_mode": "MOLECULE_ONLY",
        "solvent_generator_digest": None,
        "solvent_cells_complete": False,
    }
    interface = InterfaceResult(
        identity=ResultIdentity(
            result_kind="interface",
            system="synthetic",
            frame="static",
            network_id="producer-network-7",
            network_scope="interface",
            environment="dry",
            partition="aw",
            representation="voronoi",
            radius_configuration_id="radii-1",
            source_network_id="producer-network-7",
            solve_universe_id="solve:context-1",
            analysis_context_kind="GLOBAL_NETWORK",
            selection_id="view:interface-ab",
            solvent_mode="MOLECULE_ONLY",
            interface_id="interface-ab",
            group_A="group-A",
            group_B="group-B",
            side="shared",
            orientation="A->B",
        ),
        provenance=Provenance(
            source="nexus",
            details={
                "shared_network_identity": shared_identity,
                "source_solvent_provenance": {
                    "solvent_mode": "MOLECULE_ONLY",
                    "solvent_generator_digest": None,
                    "solvent_cells_complete": False,
                },
                "completeness": {"completeness_proven": True},
            },
        ),
        status=Status.CERTIFIED,
    )

    catalog = catalog_presentation_from_results((group, interface))

    assert [entry.selection_kind for entry in catalog.entries] == ["group", "interface"]
    assert [entry.selection_id for entry in catalog.entries] == [
        "view:group-a", "view:interface-ab"
    ]
    assert catalog.entries[1].selection_target == "interface-ab"


def test_catalog_labels_unresolved_interface_status_without_certifying_it():
    panel = AnalysisCatalogPanel()
    panel.set_catalog(AnalysisCatalogPresentation(
        source_network_id=None,
        identity_state="unresolved",
        entries=(
            AnalysisCatalogEntryPresentation(
                selection_label="Interface: A ↔ B",
                selection_id="view:interface-ab",
                selection_kind="interface",
                selection_target="interface-ab",
                completeness="UNKNOWN",
                status="UNRESOLVED",
            ),
        ),
    ))

    assert panel.network_label.text() == "Unresolved network: —"
    assert panel.table.item(0, 0).text() == "Interface: A ↔ B"
    assert panel.table.item(0, 4).text() == "UNKNOWN"
    assert panel.table.item(0, 5).text() == "UNRESOLVED"


def test_context_keeps_solvent_mode_unresolved_without_source_evidence():
    result = _typed_group_result("view:group-a", "group-A", 10, 4.5, 2.0)
    incomplete = replace(
        result,
        provenance=Provenance(
            source="legacy",
            details={"shared_network_identity": {"geometry_scheme": "aw"}},
        ),
    )

    context = context_presentation_from_results((incomplete,))

    assert context.identity_state == "validated"
    assert context.solvent_mode == "MOLECULE_ONLY"
    assert context.solvent_mode_state == "unresolved"
    assert "omit" in context.note


def test_group_overlay_mapping_uses_network_view_without_building():
    view = SimpleNamespace(
        selection_kind="group",
        system_ids=(1, 3),
        completeness=SimpleNamespace(proven=True),
    )

    class NetworkStub:
        def __init__(self):
            self.sys = SimpleNamespace()

        def get_solve_context_key(self):
            return SimpleNamespace(local_to_system=((0, 1), (1, 3)))

        def make_view(self, **_kwargs):
            return view

        def build(self):
            raise AssertionError("catalog display must not build the network")

    network = NetworkStub()
    group = SimpleNamespace(group_id="group-A", name="Group A", ball_ndxs=[1, 3])
    group.net = network
    network.sys.groups = (group,)
    result = SimpleNamespace(
        export_group=group,
        atoms=(SimpleNamespace(index=1), SimpleNamespace(index=3)),
    )
    entry = AnalysisCatalogEntryPresentation(
        selection_label="Group: group-A",
        selection_kind="group",
        selection_target="group-A",
    )

    assert validated_group_system_ids(result, entry) == (1, 3)


def test_catalog_rejects_provisional_typed_results_without_presenting_rows():
    result = _typed_group_result("view:group-a", "group-A", 10, 4.5, 2.0)
    identity = replace(
        result.identity,
        source_network_id="network:context-1",
        network_id="network:context-1",
    )
    provisional = replace(result, identity=identity)

    catalog = catalog_presentation_from_results((provisional,))

    assert catalog.identity_state == "provisional"
    assert catalog.entries == ()
    assert "provisional" in catalog.note
