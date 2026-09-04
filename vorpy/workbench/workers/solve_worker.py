"""Run an analysis backend without blocking the Qt event loop."""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from threading import Event

from PySide6.QtCore import QObject, Signal, Slot

from vorpy.workbench.services.backend import AnalysisBackend


class SolveWorker(QObject):
    progress = Signal(str, int)
    completed = Signal(object)
    failed = Signal(str)
    cancelled = Signal()
    finished = Signal()

    def __init__(
        self,
        backend: AnalysisBackend,
        source: Path | None,
        selected_indices: Sequence[int] | None = None,
    ):
        super().__init__()
        self._backend = backend
        self._source = source
        self._selected_indices = tuple(selected_indices or ())
        self._cancelled = Event()

    @Slot()
    def run(self) -> None:
        try:
            result = self._backend.solve(
                self._source,
                lambda label, value: self.progress.emit(label, value),
                self._cancelled.is_set,
                self._selected_indices or None,
            )
            self.completed.emit(result)
        except Exception as error:  # noqa: BLE001 - backends report failures through this boundary.
            if self._cancelled.is_set():
                self.cancelled.emit()
            else:
                self.failed.emit(str(error))
        finally:
            self.finished.emit()

    @Slot()
    def cancel(self) -> None:
        self._cancelled.set()


