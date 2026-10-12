from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from vorpy.workbench.services.vorpy_backend import _layers_from_network, _ProgressBridge


def test_progress_bridge_forwards_network_context_and_cancellation():
    updates = []
    bridge = _ProgressBridge(lambda label, value: updates.append((label, value)), lambda: False)
    bridge.update_progress("Building edges", 48.6, "protein")
    assert updates == [("protein: Building edges", 49)]

    cancelled = _ProgressBridge(lambda *_: None, lambda: True)
    with pytest.raises(RuntimeError, match="cancelled"):
        cancelled.update_progress("Building", 10)


def test_network_geometry_is_converted_to_viewer_layers():
    network = SimpleNamespace(
        edges=pd.DataFrame({"points": [[np.array([0, 0, 0]), np.array([1, 0, 0])]]}),
        verts=pd.DataFrame({"loc": [np.array([0, 0, 0]), np.array([1, 0, 0])]}),
        surfs=pd.DataFrame({
            "points": [np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])],
            "tris": [np.array([[0, 1, 2]])],
            "mean_tri_curvs": [np.array([0.25])],
            "gauss_tri_curvs": [np.array([0.125])],
        }),
    )

    layers = _layers_from_network(network)

    assert [layer.kind for layer in layers] == ["edges", "vertices", "surfaces"]
    assert layers[0].points.shape == (2, 3)
    assert layers[0].lines.tolist() == [[0, 1]]
    assert layers[1].points.shape == (2, 3)
    assert layers[2].faces.tolist() == [[0, 1, 2]]
    assert not layers[2].visible
    assert layers[2].cell_scalars["surface_energy"].tolist() == [0.125]
    assert set(layers[2].cell_scalars) == {"gaussian_curvature", "mean_curvature", "surface_energy", "distance", "inside_outside", "integrated_mean_curvature", "integrated_gaussian_curvature"}


def test_boundary_facing_geometry_is_labeled_and_boundary_internal_surfaces_hidden():
    network = SimpleNamespace(
        group=[0],
        boundary_indices=(3, 4),
        balls=pd.DataFrame({
            "num": [0, 1, 3, 4],
            "loc": [np.array([0., 0., 0.]), np.array([1., 0., 0.]),
                    np.array([2., 0., 0.]), np.array([3., 0., 0.])],
            "rad": [1., 1., 1.5, 1.5],
        }),
        edges=None,
        verts=None,
        surfs=pd.DataFrame({
            "balls": [[0, 1], [0, 3], [3, 4]],
            "points": [np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])] * 3,
            "tris": [np.array([[0, 1, 2]])] * 3,
            "mean_tri_curvs": [np.array([0.25])] * 3,
            "gauss_tri_curvs": [np.array([0.125])] * 3,
        }),
    )
    layers = _layers_from_network(network)
    assert [layer.name for layer in layers] == ["Voronoi surfaces", "Virtual boundary surfaces"]
    assert len(layers[0].faces) == len(layers[1].faces) == 1
    assert layers[1].visible is False


def test_worker_forwards_selection_snapshot_to_backend():
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.workers.solve_worker import SolveWorker

    seen = {}

    class BackendStub:
        def solve(self, source, progress, is_cancelled, selected_indices=None):
            seen["indices"] = selected_indices
            return AnalysisResult(source=source, name="selected")

    worker = SolveWorker(BackendStub(), None, (4, 9))
    worker.run()

    assert seen["indices"] == (4, 9)


def test_worker_forwards_selected_frame_to_backend():
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.workers.solve_worker import SolveWorker

    seen = {}

    class BackendStub:
        def solve(self, source, progress, is_cancelled, selected_indices=None,
                  frame_index=1, frame_ranges=None):
            seen.update(frame_index=frame_index, frame_ranges=frame_ranges)
            return AnalysisResult(source=source, name="frame")

    ranges = ((0, 10), (10, 20))
    SolveWorker(BackendStub(), None, frame_index=2, frame_ranges=ranges).run()
    assert seen == {"frame_index": 2, "frame_ranges": ranges}


def test_worker_preserves_live_typed_results_handoff():
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.workers.solve_worker import SolveWorker

    typed = object()

    class BackendStub:
        def solve(self, source, progress, is_cancelled, selected_indices=None):
            return AnalysisResult(
                source=source,
                name="typed",
                typed_results=(typed,),
            )

    completed = []
    worker = SolveWorker(BackendStub(), None)
    worker.completed.connect(completed.append)
    worker.run()

    assert completed[0].typed_results == (typed,)


def test_n8_edta_backend_handoff_populates_two_catalog_results(tmp_path, monkeypatch):
    from pathlib import Path
    from types import SimpleNamespace

    from vorpy.src.io import load_network, save_network
    from vorpy.src.network import Network
    from vorpy.src.results.model import Status
    from vorpy.src.results.shared_network import is_provisional_source_network_id
    from vorpy.src.group.shared_view import attach_group_view_results
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings
    from vorpy.workbench.ui.panels import (
        catalog_presentation_from_results,
        context_presentation_from_owner,
        context_presentation_from_results,
    )

    build_calls = []
    original_build = Network.build

    def counted_build(*args, **kwargs):
        build_calls.append(True)
        return original_build(*args, **kwargs)

    monkeypatch.setattr(Network, "build", counted_build)
    source = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"
    result = VorPyBackend(
        VorPySolveSettings(
            network_type="aw",
            boundary_mode="shell",
            typed_group_selections={
                "edta-group-0": (0,),
                "edta-group-1": (1,),
            },
        )
    ).solve(source, lambda *_args: None, lambda: False)
    network = result.export_group.net
    assert len(build_calls) == 1
    assert getattr(network, "archive_id", None) is None

    # A fresh solve may carry only a solve-scoped provisional namespace.
    assert len(result.typed_results) == 2
    assert result.typed_selection_system_ids == {
        "edta-group-0": (0,), "edta-group-1": (1,)
    }
    fresh_context = context_presentation_from_owner(SimpleNamespace(net=network))
    assert fresh_context.identity_state == "unresolved"
    assert fresh_context.status == "UNRESOLVED"
    assert all(item.status is Status.UNRESOLVED for item in result.typed_results)
    assert all(
        item.identity.source_identity_status == "PROVISIONAL"
        for item in result.typed_results
    )
    fresh_catalog = catalog_presentation_from_results(result.typed_results)
    assert fresh_catalog.identity_state == "provisional"
    assert fresh_catalog.entries == ()
    fresh_results_context = context_presentation_from_results(result.typed_results)
    assert fresh_results_context.identity_state == "provisional"
    assert fresh_results_context.status == "UNRESOLVED"

    unresolved_archive = tmp_path / "edta-n8-unresolved.vpy"
    save_network(network, unresolved_archive)
    unresolved_restored = load_network(unresolved_archive)
    unresolved_context = context_presentation_from_owner(
        SimpleNamespace(net=unresolved_restored)
    )
    restored_source_id = getattr(unresolved_restored, "source_network_id", None)
    if restored_source_id is not None:
        restored_key = unresolved_restored.get_solve_context_key()
        assert is_provisional_source_network_id(
            restored_source_id, restored_key.solve_universe_id
        )
    assert unresolved_context.identity_state in {"unresolved", "provisional"}
    assert unresolved_context.status == "UNRESOLVED"
    unresolved_result = AnalysisResult(source=source, name="EDTA unresolved restored")
    unresolved_typed = attach_group_view_results(
        unresolved_result,
        unresolved_restored,
        (
            SimpleNamespace(group_id="edta-group-0", ball_ndxs=(0,)),
            SimpleNamespace(group_id="edta-group-1", ball_ndxs=(1,)),
        ),
        system="EDTA",
        radius_configuration_id="workbench:atomic-defaults-v1",
        source_revision=unresolved_restored.archive_id,
    )
    unresolved_catalog = catalog_presentation_from_results(unresolved_typed)
    assert all(item.status is Status.UNRESOLVED for item in unresolved_typed)
    assert all(
        item.identity.source_identity_status == "PROVISIONAL"
        for item in unresolved_typed
    )
    assert unresolved_catalog.identity_state == "provisional"
    assert unresolved_catalog.entries == ()
    unresolved_results_context = context_presentation_from_results(unresolved_typed)
    assert unresolved_results_context.identity_state == "provisional"
    assert unresolved_results_context.status == "UNRESOLVED"

    producer_source_id = "edta-producer-n8"
    network.source_network_id = producer_source_id
    archive = tmp_path / "edta-n8-workbench.vpy"
    save_network(network, archive)
    restored = load_network(archive)
    assert restored.source_network_id == producer_source_id
    restored_result = AnalysisResult(source=source, name="EDTA restored")

    typed = attach_group_view_results(
        restored_result,
        restored,
        (
            SimpleNamespace(group_id="edta-group-0", ball_ndxs=(0,)),
            SimpleNamespace(group_id="edta-group-1", ball_ndxs=(1,)),
        ),
        system="EDTA",
        radius_configuration_id="workbench:atomic-defaults-v1",
        source_revision=restored.archive_id,
    )
    catalog = catalog_presentation_from_results(typed)

    assert len(typed) == 2
    assert len(catalog.entries) == 2
    assert catalog.identity_state == "validated"
    assert all(
        typed_result.identity.source_network_id == producer_source_id
        for typed_result in typed
    )
    assert all(entry.status == "CERTIFIED" for entry in catalog.entries)
    assert all(entry.completeness == "PROVEN" for entry in catalog.entries)
    assert all(entry.area and entry.volume for entry in catalog.entries)
    assert typed[0].identity.solve_universe_id == typed[1].identity.solve_universe_id
    assert typed[0].identity.selection_id != typed[1].identity.selection_id
    assert len(build_calls) == 1


def test_p11_edta_explicit_source_identity_is_certified_after_reload(tmp_path, monkeypatch):
    from pathlib import Path
    from types import SimpleNamespace

    from vorpy.src.group.shared_view import attach_group_view_results
    from vorpy.src.io import load_network, save_network
    from vorpy.src.network import Network
    from vorpy.src.results.model import Status
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings
    from vorpy.workbench.ui.panels import catalog_presentation_from_results

    build_calls = []
    original_build = Network.build

    def counted_build(*args, **kwargs):
        build_calls.append(True)
        return original_build(*args, **kwargs)

    monkeypatch.setattr(Network, "build", counted_build)
    source = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"
    source_network_id = "edta-source-p11"
    settings = VorPySolveSettings(
        network_type="aw",
        boundary_mode="shell",
        source_network_id=source_network_id,
        typed_group_selections={
            "edta-group-0": (0,),
            "edta-group-1": (1,),
        },
    )
    result = VorPyBackend(settings).solve(source, lambda *_args: None, lambda: False)
    network = result.export_group.net
    assert len(build_calls) == 1
    assert getattr(network, "archive_id", None) is None
    assert result.typed_selection_system_ids == {
        "edta-group-0": (0,), "edta-group-1": (1,)
    }

    typed = result.typed_results
    catalog = catalog_presentation_from_results(typed)
    assert len(typed) == len(catalog.entries) == 2
    assert catalog.identity_state == "validated"
    assert all(item.status is Status.CERTIFIED for item in typed)
    assert all(entry.status == "CERTIFIED" for entry in catalog.entries)
    assert all(entry.completeness == "PROVEN" for entry in catalog.entries)
    assert all(entry.area and entry.volume for entry in catalog.entries)
    assert all(item.identity.source_network_id == source_network_id for item in typed)
    assert typed[0].identity.solve_universe_id == typed[1].identity.solve_universe_id
    assert typed[0].identity.selection_id != typed[1].identity.selection_id

    archive = tmp_path / "edta-p11-explicit-source.vpy"
    save_network(network, archive)
    restored = load_network(archive)
    assert restored.source_network_id == source_network_id
    restored_result = AnalysisResult(source=source, name="EDTA restored")
    restored_typed = attach_group_view_results(
        restored_result,
        restored,
        (
            SimpleNamespace(group_id="edta-group-0", ball_ndxs=(0,)),
            SimpleNamespace(group_id="edta-group-1", ball_ndxs=(1,)),
        ),
        system="EDTA",
        radius_configuration_id="workbench:atomic-defaults-v1",
        source_revision=restored.archive_id,
    )
    restored_catalog = catalog_presentation_from_results(restored_typed)

    assert restored.archive_id
    assert len(restored_typed) == len(restored_catalog.entries) == 2
    assert restored_catalog.identity_state == "validated"
    assert all(item.status is Status.CERTIFIED for item in restored_typed)
    assert all(item.identity.source_network_id == source_network_id for item in restored_typed)
    assert [item.identity.selection_id for item in restored_typed] == [
        item.identity.selection_id for item in typed
    ]
    assert restored_typed[0].identity.solve_universe_id == typed[0].identity.solve_universe_id
    assert len(build_calls) == 1

    from vorpy.src.io import load_session, save_session
    from vorpy.workbench.project import Project
    from vorpy.workbench.session import session_from_workbench, workbench_from_session

    session_archive = save_session(
        session_from_workbench(Project(name="EDTA"), {source: result}, source),
        tmp_path / "edta-p11-workbench-session.vpy",
    )
    loaded_session = load_session(session_archive)
    _, session_results, active_source, _ = workbench_from_session(
        loaded_session, session_archive
    )
    session_typed = session_results[active_source].typed_results
    session_catalog = catalog_presentation_from_results(session_typed)

    assert len(session_typed) == len(session_catalog.entries) == 2
    assert all(item.status is Status.CERTIFIED for item in session_typed)
    assert all(item.identity.source_network_id == source_network_id for item in session_typed)
    assert [item.identity.selection_id for item in session_typed] == [
        item.identity.selection_id for item in typed
    ]
    assert session_typed[0].identity.solve_universe_id == typed[0].identity.solve_universe_id
    assert session_typed[0].metrics["area"].value == typed[0].metrics["area"].value
    assert session_typed[0].metrics["volume"].value == typed[0].metrics["volume"].value
    assert len(build_calls) == 1


def test_n10_retained_mg_archive_catalog_preserves_context_and_metrics(tmp_path, monkeypatch):
    from pathlib import Path
    from types import SimpleNamespace

    from vorpy.src.group.shared_view import attach_group_view_results, make_group_view
    from vorpy.src.io import load_network, save_network
    from vorpy.src.network import Network
    from vorpy.src.network.views import NetworkViewError
    from vorpy.src.results.model import Status
    from vorpy.workbench.domain import AnalysisResult
    from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings
    from vorpy.workbench.ui.panels import (
        AnalysisCatalogPanel,
        catalog_presentation_from_results,
        context_presentation_from_owner,
    )

    build_calls = []
    original_build = Network.build

    def counted_build(*args, **kwargs):
        build_calls.append(True)
        return original_build(*args, **kwargs)

    monkeypatch.setattr(Network, "build", counted_build)
    source = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"
    result = VorPyBackend(
        VorPySolveSettings(network_type="aw", boundary_mode="shell", retained_system_indices=(32,))
    ).solve(source, lambda *_args: None, lambda: False)
    network = result.export_group.net
    producer_source_id = (
        getattr(network, "source_network_id", None)
        or "edta-producer-retained-mg-p11"
    )
    network.source_network_id = producer_source_id
    archive = tmp_path / "edta-retained-mg-workbench.vpy"
    save_network(network, archive)
    restored = load_network(archive)
    assert restored.source_network_id == producer_source_id

    groups = []
    for system_id in range(32):
        candidate = SimpleNamespace(
            group_id=f"edta-atom-{system_id}",
            ball_ndxs=(system_id,),
        )
        try:
            make_group_view(candidate, restored)
        except NetworkViewError:
            continue
        groups.append(candidate)
        if len(groups) == 2:
            break
    assert len(groups) == 2

    restored_result = AnalysisResult(source=source, name="EDTA retained Mg restored")
    typed = attach_group_view_results(
        restored_result,
        restored,
        tuple(groups),
        system="EDTA",
        radius_configuration_id="workbench:atomic-defaults-v1",
        source_revision=restored.archive_id,
    )
    catalog = catalog_presentation_from_results(typed)
    context = context_presentation_from_owner(SimpleNamespace(net=restored))

    assert len(build_calls) == 1
    assert catalog.identity_state == "validated"
    assert catalog.source_network_id == producer_source_id
    assert len(catalog.entries) == 2
    assert all(item.status is Status.CERTIFIED for item in typed)
    assert all(item.identity.source_network_id == producer_source_id for item in typed)
    assert typed[0].identity.solve_universe_id == typed[1].identity.solve_universe_id
    assert typed[0].identity.selection_id != typed[1].identity.selection_id
    assert all(entry.status == "CERTIFIED" for entry in catalog.entries)
    assert all(entry.completeness == "PROVEN" for entry in catalog.entries)
    assert all(entry.area and entry.volume for entry in catalog.entries)
    assert context.solvent_mode == "MOLECULE_ONLY"
    assert context.retained_generators == "1 retained (explicit policy)"
    assert "boundary-dependent" in AnalysisCatalogPanel._HEADERS[6].lower()
