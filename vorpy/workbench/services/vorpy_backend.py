"""Production adapter from the Analysis Studio to VorPy's solver objects."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from vorpy.src.group import Group
from vorpy.src.system import System
from vorpy.workbench.domain import AnalysisResult, Atom, GeometryLayer
from vorpy.workbench.services.backend import CancellationCheck, ProgressCallback
from vorpy.workbench.services.structure_loader import load_pdb


@dataclass(frozen=True)
class VorPySolveSettings:
    network_type: str = "aw"
    max_vertices: int = 40
    box_size: float = 1.25
    surface_resolution: float = 0.2
    build_surfaces: bool = True
    build_vertices: bool = True
    build_edges: bool = True


class _ProgressBridge:
    def __init__(self, progress: ProgressCallback, is_cancelled: CancellationCheck):
        self._progress = progress
        self._is_cancelled = is_cancelled

    def update_progress(
        self, process: str, value: float, network: str | None = None
    ) -> None:
        if self._is_cancelled():
            raise RuntimeError("Analysis cancelled")
        label = f"{network}: {process}" if network else process
        self._progress(label, round(value))


class VorPyBackend:
    """Run VorPy and convert its output to the viewer's stable data model."""

    def __init__(self, settings: VorPySolveSettings | None = None):
        self.settings = settings or VorPySolveSettings()

    def solve(
        self,
        source: Path | None,
        progress: ProgressCallback,
        is_cancelled: CancellationCheck,
        selected_indices: tuple[int, ...] | None = None,
    ) -> AnalysisResult:
        if source is None:
            raise ValueError("Load a structure before running analysis")
        if source.suffix.lower() != ".pdb":
            raise ValueError("The integrated solver currently accepts PDB structures")

        started = time.perf_counter()
        progress("Reading structure", 0)
        system = System(
            file=str(source),
            gui=_ProgressBridge(progress, is_cancelled),
            print_actions=False,
        )
        if is_cancelled():
            raise RuntimeError("Analysis cancelled")

        group = Group(
            system,
            name=source.stem,
            atoms=list(selected_indices or ()),
            net_type=self.settings.network_type,
            max_vert=self.settings.max_vertices,
            box_size=self.settings.box_size,
            surf_res=self.settings.surface_resolution,
            print_metrics=False,
        )
        if group not in system.groups:
            system.groups.append(group)

        system.start_run()
        group.build()
        system.finish_run()
        if is_cancelled():
            raise RuntimeError("Analysis cancelled")

        result = load_pdb(source)
        result.atoms = _merge_atoms(result.atoms, _atoms_from_system(system))
        result.layers = _layers_from_network(group.net)
        requested_layers = {"vertices", "edges", "surfaces"}
        if not self.settings.build_vertices:
            requested_layers.discard("vertices")
        if not self.settings.build_edges:
            requested_layers.discard("edges")
        if not self.settings.build_surfaces:
            requested_layers.discard("surfaces")
        result.layers = [
            layer for layer in result.layers
            if layer.kind.lower() in requested_layers
        ]
        complete = group.net.balls.get("complete")
        result.complete_cells = (
            int(sum(bool(value) for value in complete)) if complete is not None else 0
        )
        result.surface_count = (
            len(group.net.surfs) if group.net.surfs is not None else 0
        )
        result.elapsed_seconds = time.perf_counter() - started
        progress("Preparing viewer", 100)
        return result


def _atoms_from_system(system: System) -> list[Atom]:
    atoms: list[Atom] = []
    for index, row in system.balls.reset_index(drop=True).iterrows():
        atoms.append(
            Atom(
                index=index,
                serial=int(row.get("num", index)) + 1,
                name=str(row.get("name", "")),
                element=str(row.get("element", "")),
                position=tuple(float(value) for value in row["loc"]),
                residue_name=str(row.get("res_name", row.get("residue", ""))),
                residue_sequence=str(row.get("res_seq", "")),
                chain=str(row.get("chain_name", row.get("chain", ""))),
                radius=float(row["rad"]),
            )
        )
    return atoms


def _merge_atoms(display_atoms: list[Atom], solved_atoms: list[Atom]) -> list[Atom]:
    """Keep reliable PDB identities while adopting solved coordinates."""
    if len(display_atoms) != len(solved_atoms):
        return display_atoms
    return [
        Atom(
            index=display.index,
            serial=display.serial,
            name=display.name,
            element=display.element,
            position=solved.position,
            residue_name=display.residue_name,
            residue_sequence=display.residue_sequence,
            chain=display.chain,
            radius=display.radius,
        )
        for display, solved in zip(display_atoms, solved_atoms, strict=True)
    ]


def _layers_from_network(network) -> list[GeometryLayer]:
    layers: list[GeometryLayer] = []
    if network.edges is not None and "points" in network.edges:
        points: list[np.ndarray] = []
        lines: list[tuple[int, int]] = []
        for edge_points in network.edges["points"]:
            edge_points = np.asarray(edge_points, dtype=float)
            start = len(points)
            points.extend(edge_points)
            lines.extend(
                (start + index, start + index + 1)
                for index in range(max(len(edge_points) - 1, 0))
            )
        layers.append(
            GeometryLayer(
                "Voronoi edges",
                "edges",
                np.asarray(points, dtype=float).reshape((-1, 3)),
                np.asarray(lines, dtype=np.int64).reshape((-1, 2)),
                color="#55a9d9",
            )
        )
    if network.verts is not None and "loc" in network.verts:
        layers.append(
            GeometryLayer(
                "Voronoi vertices",
                "vertices",
                np.asarray(list(network.verts["loc"]), dtype=float).reshape((-1, 3)),
                color="#efb84f",
            )
        )
    return layers
