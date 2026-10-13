"""Stable boundary between the workbench and a solver implementation."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Protocol

from vorpy.workbench.domain import AnalysisResult

ProgressCallback = Callable[[str, int], None]
CancellationCheck = Callable[[], bool]

class AnalysisBackend(Protocol):
    def solve(
        self,
        source: Path | None,
        progress: ProgressCallback,
        is_cancelled: CancellationCheck,
        selected_indices: tuple[int, ...] | None = None,
        interface_selection=None,
        frame_index: int = 1,
        frame_ranges: tuple[tuple[int, int], ...] | None = None,
    ) -> AnalysisResult:
        """Analyze one structure and return display-ready, format-neutral data.

        A producer that has ATLAS/NEXUS Results may attach one authoritative
        collection to ``AnalysisResult.typed_results``. The Workbench consumes
        that live handoff without deriving Results from defined selections;
        legacy backends may leave the field unset.
        """


