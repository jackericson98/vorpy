"""Deterministic mock solver used while the GUI remains separate from VorPy."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from vorpy.workbench.domain import AnalysisResult, Atom, Bond, GeometryLayer
from vorpy.workbench.services.backend import CancellationCheck, ProgressCallback


class MockBackend:
    def solve(
        self,
        source: Path | None,
        progress: ProgressCallback,
        is_cancelled: CancellationCheck,
        selected_indices: tuple[int, ...] | None = None,
    ) -> AnalysisResult:
        start = time.perf_counter()
        stages = [
            ("Reading structure", 8),
            ("Sorting atoms", 20),
            ("Finding vertices", 48),
            ("Building surfaces", 75),
            ("Analyzing network", 92),
            ("Preparing viewer", 100),
        ]
        for label, percent in stages:
            if is_cancelled():
                raise RuntimeError("Analysis cancelled")
            progress(label, percent)
            time.sleep(0.16)

        positions = [
            (-2.4, 0.0, 0.1), (-1.35, 1.05, 0.15), (-0.05, 1.45, -0.1),
            (1.15, 0.75, 0.2), (2.25, 1.45, -0.2), (-1.1, -0.85, 0.05),
            (0.15, -1.25, -0.15), (1.45, -0.85, 0.1), (2.55, -0.05, 0.0),
        ]
        elements = ["O", "N", "C", "N", "O", "O", "N", "C", "O"]
        atoms = [
            Atom(
                index=i,
                serial=i + 1,
                name=f"{element}{i + 1}",
                element=element,
                position=position,
                residue_name="EDT",
                residue_sequence="1",
                chain="A",
                radius=0.34 if element == "C" else 0.38,
            )
            for i, (element, position) in enumerate(zip(elements, positions))
        ]
        bonds = [Bond(i, i + 1) for i in range(len(atoms) - 1)] + [Bond(1, 5), Bond(3, 7)]

        vertex_points = np.array([
            (-1.85, .55, .65), (-.72, 1.22, .55), (.55, 1.1, .7),
            (-.48, -.65, .62), (.85, -.75, .58), (1.92, .12, .65),
        ])
        edge_points = np.array([
            vertex_points[0], vertex_points[1], vertex_points[1], vertex_points[2],
            vertex_points[0], vertex_points[3], vertex_points[3], vertex_points[4],
            vertex_points[4], vertex_points[5], vertex_points[2], vertex_points[5],
        ])
        edge_lines = np.arange(len(edge_points), dtype=np.int64).reshape(-1, 2)
        layers = [
            GeometryLayer("Voronoi edges", "edges", edge_points, edge_lines, color="#55a9d9"),
            GeometryLayer("Voronoi vertices", "vertices", vertex_points, color="#efb84f"),
        ]
        name = source.stem if source else "EDTA demo"
        return AnalysisResult(
            source=source,
            name=name,
            atoms=atoms,
            bonds=bonds,
            layers=layers,
            complete_cells=len(atoms),
            surface_count=24,
            elapsed_seconds=time.perf_counter() - start,
        )


