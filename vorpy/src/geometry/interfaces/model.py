"""Shared interface/topology result records; physical geometry remains primal."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class InterfaceSurface:
    feature_id: int
    generator_ids: tuple[int, int]
    dual_simplex_id: str | None
    area_A2: float | None
    bounded: bool
    complete: bool
    generator_cells_complete: bool
    supported: bool
    selected: bool
    status: str
    reason: str = ""
    boundary_edge_ids: tuple[int, ...] = ()
    boundary_vertex_ids: tuple[int, ...] = ()
    smooth_H_A: float | None = None
    stable_atom_i: str = ""
    stable_atom_j: str = ""
    group_i: str = ""
    group_j: str = ""


@dataclass(frozen=True)
class InterfaceEdge:
    feature_id: int
    generator_ids: tuple[int, ...]
    dual_simplex_id: str | None
    incident_candidate_surface_ids: tuple[int, ...]
    incident_selected_surface_ids: tuple[int, ...]
    bounded: bool
    complete: bool
    supported: bool
    incidence_status: str
    reason: str = ""
    endpoint_vertex_ids: tuple[int, ...] = ()
    classification: str = "unresolved"
    length_A: float | None = None
    signed_turn_or_beta: float | None = None
    signed_integral_A_rad: float | None = None
    unsigned_integral_A_rad: float | None = None
    dual_triangle_present_at_alpha: bool | None = None
    dual_triangle_birth: float | None = None
    dual_triangle_status: str = "not_applicable"
    curvature_status: str = "not_computed"


@dataclass(frozen=True)
class InterfaceVertex:
    feature_id: int
    generator_ids: tuple[int, ...]
    incident_candidate_surface_ids: tuple[int, ...]
    incident_selected_surface_ids: tuple[int, ...]
    bounded: bool
    supported: bool
    status: str
    xyz: tuple[float, float, float] | None = None


@dataclass
class InterfaceRepresentation:
    scheme: str
    selection_mode: str
    alpha: float | None
    group_a: frozenset[int]
    group_b: frozenset[int]
    surfaces: list[InterfaceSurface] = field(default_factory=list)
    edges: list[InterfaceEdge] = field(default_factory=list)
    vertices: list[InterfaceVertex] = field(default_factory=list)
    curvature: object | None = None
    metadata: dict = field(default_factory=dict)
    selection_mappings: list[dict] = field(default_factory=list)
    alpha_incidence_audit: list[dict] = field(default_factory=list)

    @property
    def candidate_surfaces(self):
        return [record for record in self.surfaces if record.status != "not_bicolor"]

    @property
    def selected_surfaces(self):
        return [record for record in self.surfaces if record.selected]

    @property
    def excluded_surfaces(self):
        return [record for record in self.candidate_surfaces if not record.selected]

    @property
    def interface_edges(self):
        return [record for record in self.edges
                if record.incident_selected_surface_ids]

    @property
    def interface_vertices(self):
        return [record for record in self.vertices
                if record.incident_selected_surface_ids]

    @property
    def interface_atoms(self):
        return frozenset(generator for surface in self.selected_surfaces
                         for generator in surface.generator_ids)

    @property
    def area_A2(self):
        return float(sum(surface.area_A2 or 0.0 for surface in self.selected_surfaces))

    @property
    def component_count(self):
        selected = self.selected_surfaces
        if not selected:
            return 0
        parent = {surface.feature_id: surface.feature_id for surface in selected}

        def find(value):
            while parent[value] != value:
                parent[value] = parent[parent[value]]
                value = parent[value]
            return value

        def union(first, second):
            a, b = find(first), find(second)
            if a != b:
                parent[b] = a

        for edge in self.interface_edges:
            ids = edge.incident_selected_surface_ids
            for other in ids[1:]:
                union(ids[0], other)
        return len({find(surface.feature_id) for surface in selected})

    @property
    def euler_characteristic(self):
        # This is the Euler characteristic of the retained incidence complex;
        # it does not imply a closed surface or justify a genus assignment.
        return (len(self.interface_vertices) - len(self.interface_edges)
                + len(self.selected_surfaces))

    @property
    def unresolved_feature_count(self):
        return sum(surface.status == "unresolved" for surface in self.surfaces) + sum(
            edge.incidence_status in {"unresolved", "nonmanifold"} for edge in self.edges
        )

    def summary(self):
        status_counts = {}
        curvature_status_counts = {}
        for edge in self.edges:
            status_counts[edge.incidence_status] = status_counts.get(edge.incidence_status, 0) + 1
            curvature_status_counts[edge.curvature_status] = (
                curvature_status_counts.get(edge.curvature_status, 0) + 1
            )
        return {
            "scheme": self.scheme,
            "selection_mode": self.selection_mode,
            "alpha": self.alpha,
            "group_a_generator_count": len(self.group_a),
            "group_b_generator_count": len(self.group_b),
            "candidate_surfaces": len(self.candidate_surfaces),
            "selected_surfaces": len(self.selected_surfaces),
            "excluded_surfaces": len(self.excluded_surfaces),
            "interface_atoms": len(self.interface_atoms),
            "interface_edges": len(self.interface_edges),
            "interface_vertices": len(self.interface_vertices),
            "connected_components": self.component_count,
            "area_A2": self.area_A2,
            "euler_characteristic_incidence_complex": self.euler_characteristic,
            "edge_incidence_status_counts": status_counts,
            "edge_curvature_status_counts": curvature_status_counts,
            "unresolved_features": self.unresolved_feature_count,
            "curvature": None if self.curvature is None else asdict(self.curvature),
            "selection_mapping_count": len(self.selection_mappings),
            "selection_mapping_failures": sum(not row.get("mapped", False) for row in self.selection_mappings),
            "alpha_incidence_audit_count": len(self.alpha_incidence_audit),
            "metadata": self.metadata,
        }

    def export(self, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        for name, rows in (("interface_surfaces.csv", self.surfaces),
                           ("interface_edges.csv", self.edges),
                           ("interface_vertices.csv", self.vertices)):
            payload = [asdict(row) for row in rows]
            fields = list(payload[0]) if payload else []
            with (directory / name).open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                if fields:
                    writer.writeheader()
                    writer.writerows({key: _csv_value(value) for key, value in row.items()}
                                     for row in payload)
        _write_dict_rows(directory / "alpha_incidence_audit.csv", self.alpha_incidence_audit)
        _write_dict_rows(directory / "alpha_surface_mapping.csv", self.selection_mappings)
        if self.curvature is not None:
            curvature = asdict(self.curvature)
            _write_dict_rows(directory / "curvature_summary.csv", [curvature])
        (directory / "interface_summary.json").write_text(
            json.dumps(self.summary(), indent=2, sort_keys=True, default=str, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return directory


def _csv_value(value):
    if isinstance(value, (tuple, list, set, frozenset, dict)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def _write_dict_rows(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if fields:
            writer.writeheader()
            writer.writerows({key: _csv_value(value) for key, value in row.items()} for row in rows)
