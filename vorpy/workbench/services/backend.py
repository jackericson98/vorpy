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
    ) -> AnalysisResult:
        """Analyze one structure and return display-ready, format-neutral data."""



