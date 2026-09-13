"""Production adapter from the Analysis Studio to VorPy's solver objects."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from vorpy.src.group import Group
from vorpy.src.system import System
from vorpy.workbench.domain import AnalysisResult, Atom, GeometryLayer
from vorpy.workbench.atomic_defaults import apply_system_defaults
from vorpy.workbench.services.backend import CancellationCheck, ProgressCallback
from vorpy.workbench.services.structure_loader import load_pdb
from vorpy.src.output.curvature_colors import (
    component_value,
    mean_vertex_display_value,
)


@dataclass(frozen=True)
class VorPySolveSettings:
    network_type: str = "aw"
    max_vertices: int = 40
    box_size: float = 1.25
    surface_resolution: float = 0.2
    build_surfaces: bool = True
    build_vertices: bool = True
    build_edges: bool = True
    atomic_defaults: dict[str, dict[str, float]] = field(default_factory=dict)


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

        apply_system_defaults(system, self.settings.atomic_defaults)
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
        shell_surfs = shell_edges = shell_verts = ()
        if all(
            getattr(group.net, name, None) is not None
            for name in ("surfs", "edges", "verts")
        ):
            group.get_layers(max_layers=1, group_resids=False)
            if group.layer_surfs:
                shell_surfs = group.layer_surfs[0]
            if group.layer_edges:
                shell_edges = group.layer_edges[0]
            if group.layer_verts:
                shell_verts = group.layer_verts[0]
        system.finish_run()
        if is_cancelled():
            raise RuntimeError("Analysis cancelled")

        result = load_pdb(source)
        result.atoms = _merge_atoms(result.atoms, _atoms_from_system(system))
        result.layers = _layers_from_network(
            group.net,
            shell_surfs=shell_surfs,
            shell_edges=shell_edges,
            shell_verts=shell_verts,
        )
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


def _numeric_charge(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


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
                mass=float(row["mass"]) if row.get("mass") is not None else None,
                charge=_numeric_charge(row.get("charge")),
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
            # The solved System is the source of truth for geometry:
            # preserve PDB identity fields, but keep the exact radius that
            # VorPy used to build the network.
            radius=solved.radius,
            mass=solved.mass,
            charge=solved.charge if solved.charge is not None else display.charge,
            source_properties={"radius": display.radius, "mass": display.mass, "charge": display.charge},
        )
        for display, solved in zip(display_atoms, solved_atoms, strict=True)
    ]


def _surface_geometry(
    surfaces, balls, target_cells=None, interpretation="magnitude"
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    points: list[np.ndarray] = []
    faces: list[np.ndarray] = []
    scalar_values = {
        "gaussian_curvature": [],
        "mean_curvature": [],
        "integrated_mean_curvature": [],
        "integrated_gaussian_curvature": [],
        "surface_energy": [],
        "distance": [],
        "inside_outside": [],
    }
    offset = 0
    for _, surface in surfaces.iterrows():
        surface_points = np.asarray(surface["points"], dtype=float).reshape((-1, 3))
        surface_faces = np.asarray(surface["tris"], dtype=np.int64).reshape((-1, 3))
        face_count = len(surface_faces)
        points.extend(surface_points)
        faces.extend(surface_faces + offset)
        offset += len(surface_points)

        for scheme, column in (
            ("gaussian_curvature", "gauss_tri_curvs"),
            ("mean_curvature", "mean_tri_curvs"),
        ):
            values = np.asarray(surface.get(column, []), dtype=float).ravel()
            if len(values) != face_count:
                fallback = "gauss_curv" if scheme == "gaussian_curvature" else "mean_curv"
                values = np.full(face_count, float(surface.get(fallback, 0.0) or 0.0))
            scalar_values[scheme].extend(values)

        mean_values = np.asarray(scalar_values["mean_curvature"][-face_count:], dtype=float)
        scalar_values["surface_energy"].extend(2.0 * np.square(mean_values))

        for display_name, core_name in (
            ("integrated_mean_curvature", "int_mean_curv"),
            ("integrated_gaussian_curvature", "int_gauss_curv"),
        ):
            value = component_value(
                surface, "surface", core_name, target_cells, mode=interpretation
            )
            scalar_values[display_name].extend(
                np.full(face_count, 0.0 if value is None else float(value))
            )

        center = np.asarray(surface.get("loc", surface_points.mean(axis=0)), dtype=float)
        point_distances = np.linalg.norm(surface_points - center, axis=1)
        scalar_values["distance"].extend(np.max(point_distances[surface_faces], axis=1))

        inside = np.zeros(len(surface_points), dtype=bool)
        defining = surface.get("balls", [])
        if balls is not None and len(defining):
            matches = balls
            if "num" in balls:
                matches = balls.loc[balls["num"] == int(defining[0])]
            elif int(defining[0]) in balls.index:
                matches = balls.loc[[int(defining[0])]]
            if len(matches):
                atom = matches.iloc[0]
                inside = np.linalg.norm(
                    surface_points - np.asarray(atom["loc"], dtype=float), axis=1
                ) < float(atom["rad"])
        scalar_values["inside_outside"].extend(
            np.all(inside[surface_faces], axis=1).astype(float)
        )

    return (
        np.asarray(points, dtype=float).reshape((-1, 3)),
        np.asarray(faces, dtype=np.int64).reshape((-1, 3)),
        {key: np.asarray(values, dtype=float) for key, values in scalar_values.items()},
    )


def _edge_geometry(
    edges, target_cells=None, interpretation="magnitude"
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    points: list[np.ndarray] = []
    lines: list[tuple[int, int]] = []
    scalars = {
        "integrated_mean_curvature": [],
        "integrated_gaussian_curvature": [],
    }
    for _, edge in edges.iterrows():
        edge_points = np.asarray(edge["points"], dtype=float)
        start = len(points)
        points.extend(edge_points)
        segment_count = max(len(edge_points) - 1, 0)
        lines.extend((start + index, start + index + 1) for index in range(segment_count))
        for display_name, core_name in (
            ("integrated_mean_curvature", "int_mean_curv"),
            ("integrated_gaussian_curvature", "int_gauss_curv"),
        ):
            value = component_value(
                edge, "edge", core_name, target_cells, mode=interpretation
            )
            scalars[display_name].extend(
                np.full(segment_count, 0.0 if value is None else float(value))
            )
    line_array = np.asarray(lines, dtype=np.int64).reshape((-1, 2))
    scalar_arrays = {
        key: np.asarray(values, dtype=float)
        for key, values in scalars.items()
    }
    for name, values in scalar_arrays.items():
        if len(values) != len(line_array):
            raise ValueError(
                f"Edge scalar {name!r} has {len(values)} values for "
                f"{len(line_array)} rendered line segments"
            )

    return (
        np.asarray(points, dtype=float).reshape((-1, 3)),
        line_array,
        scalar_arrays,
    )


def _vertex_geometry(
    network, vertices, target_cells=None, interpretation="magnitude"
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    points = np.asarray(list(vertices["loc"]), dtype=float).reshape((-1, 3))
    mean_values, gauss_values = [], []
    for _, vertex in vertices.iterrows():
        mean_value = mean_vertex_display_value(
            network, vertex, target_cells, mode=interpretation
        )
        gauss_value = component_value(
            vertex, "vertex", "int_gauss_curv", target_cells, mode=interpretation
        )
        mean_values.append(0.0 if mean_value is None else float(mean_value))
        gauss_values.append(0.0 if gauss_value is None else float(gauss_value))
    return points, {
        "integrated_mean_curvature": np.asarray(mean_values, dtype=float),
        "integrated_gaussian_curvature": np.asarray(gauss_values, dtype=float),
    }


def _layers_from_network(
    network,
    shell_surfs=(),
    shell_edges=(),
    shell_verts=(),
) -> list[GeometryLayer]:
    layers: list[GeometryLayer] = []
    raw_group = getattr(network, "group", None)
    target_cells = () if raw_group is None else tuple(int(value) for value in raw_group)

    if network.edges is not None and "points" in network.edges:
        points, lines, scalars = _edge_geometry(
            network.edges, target_cells, interpretation="magnitude"
        )
        layers.append(GeometryLayer(
            "Voronoi edges", "edges", points, lines, color="#55a9d9",
            cell_scalars=scalars, interpretation="magnitude",
        ))
        shell_edge_rows = network.edges.iloc[list(shell_edges)]
        if not shell_edge_rows.empty:
            points, lines, scalars = _edge_geometry(
                shell_edge_rows, target_cells, interpretation="boundary"
            )
            layers.append(GeometryLayer(
                "Voronoi shell edges", "edges", points, lines, color="#42d6c7",
                cell_scalars=scalars, interpretation="boundary",
            ))

    if network.verts is not None and "loc" in network.verts:
        points, scalars = _vertex_geometry(
            network, network.verts, target_cells, interpretation="magnitude"
        )
        layers.append(GeometryLayer(
            "Voronoi vertices", "vertices", points, color="#efb84f",
            cell_scalars=scalars, interpretation="magnitude",
        ))
        shell_vertex_rows = network.verts.iloc[list(shell_verts)]
        if not shell_vertex_rows.empty:
            points, scalars = _vertex_geometry(
                network, shell_vertex_rows, target_cells, interpretation="boundary"
            )
            layers.append(GeometryLayer(
                "Voronoi shell vertices", "vertices", points, color="#f29f67",
                cell_scalars=scalars, interpretation="boundary",
            ))

    if network.surfs is not None and "points" in network.surfs and "tris" in network.surfs:
        points, faces, scalars = _surface_geometry(
            network.surfs, getattr(network, "balls", None), target_cells,
            interpretation="magnitude",
        )
        layers.append(GeometryLayer(
            "Voronoi surfaces", "surfaces", points=points, faces=faces,
            color="#4f9fcf", cell_scalars=scalars, opacity=0.45,
            visible=False, interpretation="magnitude",
        ))
        shell_surface_rows = network.surfs.iloc[list(shell_surfs)]
        if not shell_surface_rows.empty:
            points, faces, scalars = _surface_geometry(
                shell_surface_rows, getattr(network, "balls", None), target_cells,
                interpretation="boundary",
            )
            layers.append(GeometryLayer(
                "Voronoi shell surfaces", "surfaces", points=points, faces=faces,
                color="#806df0", cell_scalars=scalars, opacity=0.45,
                interpretation="boundary",
            ))
    return layers

