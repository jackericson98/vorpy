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
    )

    layers = _layers_from_network(network)

    assert [layer.kind for layer in layers] == ["edges", "vertices"]
    assert layers[0].points.shape == (2, 3)
    assert layers[0].lines.tolist() == [[0, 1]]
    assert layers[1].points.shape == (2, 3)


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
