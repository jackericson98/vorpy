"""Format-neutral data passed between an analysis backend and the GUI."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Atom:
    index: int
    serial: int
    name: str
    element: str
    position: tuple[float, float, float]
    residue_name: str = ""
    residue_sequence: str = ""
    chain: str = ""
    radius: float = 0.35


@dataclass(frozen=True)
class Bond:
    atom_a: int
    atom_b: int


@dataclass
class GeometryLayer:
    name: str
    kind: str
    points: np.ndarray = field(default_factory=lambda: np.empty((0, 3), dtype=float))
    lines: np.ndarray | None = None
    source_path: Path | None = None
    color: str = "#55a9d9"
    opacity: float = 1.0
    visible: bool = True


@dataclass
class AnalysisResult:
    source: Path | None
    name: str
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)
    layers: list[GeometryLayer] = field(default_factory=list)
    complete_cells: int = 0
    surface_count: int = 0
    elapsed_seconds: float = 0.0


