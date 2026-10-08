"""Curved selected partner-restricted union-of-balls contact surfaces.

The canonical geometry is a spherical arrangement.  Tessellated rendering is
deliberately not part of this module's scientific representation.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import Iterable, Mapping, Optional, Sequence


Vec3 = tuple[float, float, float]


def _add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _scale(a: Vec3, value: float) -> Vec3:
    return (a[0] * value, a[1] * value, a[2] * value)


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vec3, tolerance: float) -> Optional[Vec3]:
    length = _norm(a)
    if length <= tolerance:
        return None
    return _scale(a, 1.0 / length)


def _clamp(value: float, low: float = -1.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def _token(prefix: str, *values: object) -> str:
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"), default=str)
    return prefix + "-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


@dataclass(frozen=True)
class GeometryTolerance:
    distance: float = 1e-9
    angular: float = 1e-9
    area: float = 1e-10
    quadrature: float = 1e-8


@dataclass(frozen=True)
class AtomBall:
    """One weighted atom ball used by the expanded union-of-balls surface."""

    generator_index: int
    stable_id: str
    group: str
    center: Vec3
    radius: float

    def __post_init__(self):
        if self.group not in {"A", "B", "S"}:
            raise ValueError("AtomBall group must be A, B, or S")
        if len(self.center) != 3 or not all(math.isfinite(value) for value in self.center):
            raise ValueError("AtomBall center must contain finite xyz coordinates")
        if not math.isfinite(self.radius) or self.radius <= 0:
            raise ValueError("AtomBall radius must be positive and finite")


@dataclass(frozen=True)
class SelectedContact:
    """A selected Power dual contact, retained as input provenance."""

    contact_id: str
    atom_a_index: int
    atom_b_index: int
    source_network_id: str = ""
    source_interface_id: str = ""
    status: str = "SELECTED"

    def __post_init__(self):
        if not self.contact_id:
            raise ValueError("Selected contacts require contact_id")
        if self.status not in {"SELECTED", "UNRESOLVED", "REJECTED"}:
            raise ValueError("Unknown selected-contact status")


@dataclass(frozen=True)
class JunctionVertex:
    vertex_id: str
    carrier_atom_id: str
    xyz: Vec3
    unit_direction: Vec3
    source_circle_ids: tuple[str, ...]
    incident_arc_ids: tuple[str, ...] = ()
    status: str = "CERTIFIED"
    source_carrier_atom_ids: tuple[str, ...] = ()
    incident_edge_ids: tuple[str, ...] = ()
    local_link: str = "UNRESOLVED"


@dataclass(frozen=True)
class CircularArc:
    arc_id: str
    carrier_atom_id: str
    circle_id: str
    start_vertex_id: Optional[str]
    end_vertex_id: Optional[str]
    start_xyz: Vec3
    end_xyz: Vec3
    circle_normal: Vec3
    circle_offset: float
    angular_span: float
    source_atom_ids: tuple[str, ...]
    source_contact_ids: tuple[str, ...]
    retained_side: int
    status: str = "CERTIFIED"


@dataclass(frozen=True)
class SphericalPatch:
    patch_id: str
    carrier_atom_id: str
    carrier_center: Vec3
    carrier_radius: float
    side: str
    partner_atom_ids: tuple[str, ...]
    self_occluding_atom_ids: tuple[str, ...]
    source_contact_ids: tuple[str, ...]
    source_network_id: str
    source_interface_id: str
    definition_version: str
    boundary_arc_ids: tuple[str, ...]
    junction_ids: tuple[str, ...]
    component_id: str
    area_A2: float
    area_error_A2: float
    region_kind: str
    status: str = "CERTIFIED"
    source_contact_selection_id: str = ""


@dataclass(frozen=True)
class RefinedSeamPiece:
    """One canonical interval of a same-group carrier intersection circle."""

    seam_id: str
    carrier_atom_ids: tuple[str, str]
    circle_id: str
    circle_center: Vec3
    circle_radius: float
    circle_normal: Vec3
    parameter_interval: tuple[float, float]
    endpoint_junction_ids: tuple[str, str]
    start_xyz: Vec3
    end_xyz: Vec3
    incident_patch_ids: tuple[str, ...]
    source_arc_ids: tuple[str, ...]
    source_contact_ids: tuple[str, ...]
    orientation_on_patches: tuple[tuple[str, int], ...]
    classification: str
    status: str = "CERTIFIED"


@dataclass(frozen=True)
class TopologicalCut:
    """Auxiliary cut used only to regularize a multiply-boundary face cell."""

    cut_id: str
    patch_id: str
    carrier_atom_id: str
    cycle_indices: tuple[int, int]
    endpoint_junction_ids: tuple[str, str]
    start_xyz: Vec3
    end_xyz: Vec3
    status: str = "TOPOLOGY_ONLY"


@dataclass(frozen=True)
class MolecularContactSurfaceCurvature:
    """Curvature audit derived from canonical spherical geometry.

    Topology-only cuts are intentionally absent from every quantity here.
    """

    smooth_mean_curvature: Optional[float]
    seam_mean_curvature: Optional[float]
    total_integrated_mean_curvature: Optional[float]
    smooth_gaussian_curvature: Optional[float]
    intrinsic_discrete_gaussian_curvature: Optional[float]
    intrinsic_seam_gaussian_curvature: Optional[float]
    interior_vertex_defects: Optional[float]
    boundary_geodesic_curvature: Optional[float]
    boundary_corner_turning: Optional[float]
    completed_gaussian_curvature: Optional[float]
    gauss_bonnet_expected: Optional[float]
    gauss_bonnet_residual: Optional[float]
    mean_status: str
    gaussian_status: str
    status: str
    method: str
    diagnostics: tuple[str, ...] = ()
    physical_seam_count: int = 0
    exterior_boundary_piece_count: int = 0
    audited_junction_count: int = 0

    def summary(self, area_A2: float) -> dict:
        return {
            "smooth_mean_curvature": self.smooth_mean_curvature,
            "seam_mean_curvature": self.seam_mean_curvature,
            "total_integrated_mean_curvature": self.total_integrated_mean_curvature,
            "smooth_gaussian_curvature": self.smooth_gaussian_curvature,
            "intrinsic_discrete_gaussian_curvature": self.intrinsic_discrete_gaussian_curvature,
            "intrinsic_seam_gaussian_curvature": self.intrinsic_seam_gaussian_curvature,
            "interior_vertex_defects": self.interior_vertex_defects,
            "boundary_geodesic_curvature": self.boundary_geodesic_curvature,
            "boundary_corner_turning": self.boundary_corner_turning,
            "completed_gaussian_curvature": self.completed_gaussian_curvature,
            "gauss_bonnet_expected": self.gauss_bonnet_expected,
            "gauss_bonnet_residual": self.gauss_bonnet_residual,
            "smooth_mean_curvature_per_area": (
                self.smooth_mean_curvature / area_A2 if area_A2 else None
            ),
            "smooth_gaussian_curvature_per_area": (
                self.smooth_gaussian_curvature / area_A2 if area_A2 else None
            ),
            "mean_status": self.mean_status,
            "gaussian_status": self.gaussian_status,
            "status": self.status,
            "method": self.method,
            "diagnostics": list(self.diagnostics),
            "physical_seam_count": self.physical_seam_count,
            "exterior_boundary_piece_count": self.exterior_boundary_piece_count,
            "audited_junction_count": self.audited_junction_count,
        }


@dataclass(frozen=True)
class MolecularContactSurface:
    """Canonical 2D molecular-contact boundary surface.

    Its Euler characteristic and Betti numbers describe this surface, not
    the 3D solid enclosed by it.
    """

    side: str
    convention: str
    definition_version: str
    patches: tuple[SphericalPatch, ...]
    arcs: tuple[CircularArc, ...]
    junctions: tuple[JunctionVertex, ...]
    carrier_atom_ids: tuple[str, ...]
    selected_contact_ids: tuple[str, ...]
    total_area_A2: float
    naive_pairwise_area_A2: float
    components: int
    boundary_loops: int
    euler_characteristic: Optional[int]
    nonmanifold_vertices: int
    unresolved_contact_ids: tuple[str, ...]
    area_method: str
    area_error_A2: float
    status: str
    diagnostics: tuple[str, ...] = ()
    contact_selection_id: str = ""
    refined_seams: tuple[RefinedSeamPiece, ...] = ()
    global_vertex_count: int = 0
    global_edge_count: int = 0
    global_face_count: int = 0
    manifold_vertices: int = 0
    singular_vertices: int = 0
    junction_link_status: tuple[tuple[str, str], ...] = ()
    orientability: Optional[bool] = None
    genus_by_component: Optional[tuple[int, ...]] = None
    area_status: str = "CERTIFIED"
    topological_cuts: tuple[TopologicalCut, ...] = ()
    topological_face_count: int = 0
    spherical_patch_record_count: int = 0
    face_boundary_cycle_counts: tuple[tuple[str, int], ...] = ()
    geometric_edge_count: int = 0
    betti_numbers: Optional[tuple[int, int, int]] = None
    incidence_status: str = "UNRESOLVED"
    edge_incidence_distribution: tuple[tuple[str, int], ...] = ()
    curvature: Optional[MolecularContactSurfaceCurvature] = None
    contact_selection_source: Optional[str] = None

    @property
    def boundary_arc_count(self) -> int:
        return len(self.arcs)

    @property
    def junction_count(self) -> int:
        return len(self.junctions)

    def summary(self) -> dict:
        return {
            "side": self.side,
            "convention": self.convention,
            "definition_version": self.definition_version,
            "contact_selection_source": self.contact_selection_source,
            "contact_selection_id": self.contact_selection_id,
            "carrier_atoms": len(self.carrier_atom_ids),
            "patches": len(self.patches),
            "arcs": len(self.arcs),
            "junction_vertices": len(self.junctions),
            "area_A2": self.total_area_A2,
            "naive_pairwise_area_A2": self.naive_pairwise_area_A2,
            "area_reduction_A2": self.naive_pairwise_area_A2 - self.total_area_A2,
            "area_error_A2": self.area_error_A2,
            "components": self.components,
            "boundary_loops": self.boundary_loops,
            "euler_characteristic": self.euler_characteristic,
            "nonmanifold_vertices": self.nonmanifold_vertices,
            "global_vertices": self.global_vertex_count,
            "global_edges": self.global_edge_count,
            "global_faces": self.global_face_count,
            "refined_seams": len(self.refined_seams),
            "manifold_vertices": self.manifold_vertices,
            "singular_vertices": self.singular_vertices,
            "junction_link_status": dict(self.junction_link_status),
            "orientability": self.orientability,
            "genus_by_component": self.genus_by_component,
            "area_status": self.area_status,
            "topological_cuts": len(self.topological_cuts),
            "topological_faces": self.topological_face_count,
            "spherical_patch_records": self.spherical_patch_record_count,
            "face_boundary_cycle_counts": dict(self.face_boundary_cycle_counts),
            "geometric_edges": self.geometric_edge_count,
            "betti_numbers": self.betti_numbers,
            "incidence_status": self.incidence_status,
            "edge_incidence_distribution": dict(self.edge_incidence_distribution),
            "curvature": self.curvature.summary(self.total_area_A2) if self.curvature else None,
            "patchwise_gauss_bonnet": audit_molecular_contact_surface_gauss_bonnet(self),
            "topological_cut_ids": [cut.cut_id for cut in self.topological_cuts],
            "selected_contacts": len(self.selected_contact_ids),
            "unresolved_contacts": list(self.unresolved_contact_ids),
            "area_method": self.area_method,
            "status": self.status,
            "diagnostics": list(self.diagnostics),
        }

    def to_contract(self, identity, provenance, interaction_links=None):
        """Adapt measurements/provenance to Agent #1's immutable result contract."""
        from vorpy.src.results import (
            GaussianCurvature,
            MeanCurvature,
            MolecularContactSelection,
            MolecularContactSurfaceGeometry,
            MolecularContactSurfaceResult,
            MolecularGeometryCounts,
            MolecularPatchProvenance,
            Quantity,
            Status,
            TopologicalCellComplex,
        )

        def canonical_status(value):
            try:
                return value if isinstance(value, Status) else Status(value)
            except (TypeError, ValueError):
                return Status.UNRESOLVED

        status = canonical_status(self.status)
        area_status = canonical_status(self.area_status)

        def quantity(name, value, units="1", scope="surface", quantity_status=None):
            effective_status = quantity_status or status
            if effective_status in {
                Status.UNRESOLVED,
                Status.NOT_CALCULATED,
                Status.NOT_SUPPORTED,
            }:
                value = None
            elif value is None and effective_status in {Status.CERTIFIED, Status.PARTIAL}:
                effective_status = Status.NOT_CALCULATED
            return Quantity(
                name, value, units,
                effective_status,
                self.area_method, scope, provenance=provenance,
            )

        geometry = MolecularContactSurfaceGeometry(
            quantities={
                "carrier_atom_count": quantity("carrier_atom_count", len(self.carrier_atom_ids)),
                "geometric_edge_count": quantity("geometric_edge_count", self.geometric_edge_count),
                "spherical_patch_count": quantity("spherical_patch_count", len(self.patches)),
                "circular_arc_count": quantity("circular_arc_count", len(self.arcs)),
                "junction_record_count": quantity("junction_record_count", len(self.junctions)),
                "refined_seam_count": quantity("refined_seam_count", len(self.refined_seams)),
                "boundary_arc_count": quantity("boundary_arc_count", self.boundary_arc_count),
            },
            status=area_status,
        )
        geometric_counts = MolecularGeometryCounts(
            quantities={
                "carrier_atom_count": quantity("carrier_atom_count", len(self.carrier_atom_ids)),
                "geometric_edge_count": quantity("geometric_edge_count", self.geometric_edge_count),
                "spherical_patch_count": quantity("spherical_patch_count", len(self.patches)),
                "circular_arc_count": quantity("circular_arc_count", len(self.arcs)),
                "junction_record_count": quantity("junction_record_count", len(self.junctions)),
                "refined_seam_count": quantity("refined_seam_count", len(self.refined_seams)),
                "boundary_arc_count": quantity("boundary_arc_count", self.boundary_arc_count),
            },
            status=area_status,
        )
        selection = MolecularContactSelection(
            quantities={
                "selected_contacts": quantity("selected_contacts", len(self.selected_contact_ids), scope="selection"),
                "unresolved_contacts": quantity("unresolved_contacts", len(self.unresolved_contact_ids), scope="selection"),
            },
            selected_ab_contact_ids=self.selected_contact_ids,
            direct_contact_atom_ids=self.carrier_atom_ids,
            contributing_surface_atom_ids=tuple(sorted({patch.carrier_atom_id for patch in self.patches})),
            occluding_same_group_atom_ids=tuple(sorted({atom for patch in self.patches for atom in patch.self_occluding_atom_ids})),
            partner_clipping_atom_ids=tuple(sorted({atom for patch in self.patches for atom in patch.partner_atom_ids})),
            unresolved_contact_ids=self.unresolved_contact_ids,
            status=status,
        )
        topology_status = (
            Status.CERTIFIED
            if status is Status.CERTIFIED and self.incidence_status == "CERTIFIED"
            else Status.PARTIAL
            if status is Status.CERTIFIED or status is Status.PARTIAL
            else status
        )

        def topology_quantity(name, value, scope="topology"):
            return quantity(name, value, scope=scope, quantity_status=topology_status)

        def certification(prefix=None):
            if self.incidence_status == "CERTIFIED":
                return True
            if prefix is None and self.diagnostics:
                return False
            if prefix is not None and any(diagnostic.startswith(prefix) for diagnostic in self.diagnostics):
                return False
            return None

        genus = (self.genus_by_component[0]
                 if self.genus_by_component is not None and len(self.genus_by_component) == 1
                 else None)
        nonorientable_genus = genus if self.orientability is False else None
        topology_quantities = {
            "vertex_count": topology_quantity("vertex_count", self.global_vertex_count),
            "edge_count": topology_quantity("edge_count", self.global_edge_count),
            "topological_edge_count": topology_quantity("topological_edge_count", self.global_edge_count),
            "topological_cut_count": topology_quantity("topological_cut_count", len(self.topological_cuts)),
            "patch_attachment_record_count": topology_quantity(
                "patch_attachment_record_count", self.spherical_patch_record_count),
            "face_count": topology_quantity("face_count", self.global_face_count),
            "connected_component_count": topology_quantity("connected_component_count", self.components),
            "boundary_component_count": topology_quantity("boundary_component_count", self.boundary_loops),
            "euler_characteristic": topology_quantity("euler_characteristic", self.euler_characteristic),
            "orientable": topology_quantity("orientable", self.orientability),
            "genus": topology_quantity("genus", genus),
            "nonorientable_genus": topology_quantity("nonorientable_genus", nonorientable_genus),
            "manifold_vertex_count": topology_quantity("manifold_vertex_count", self.manifold_vertices),
            "singular_vertex_count": topology_quantity("singular_vertex_count", self.singular_vertices),
            "chain_complex_valid": topology_quantity("chain_complex_valid", certification()),
            "boundary_of_boundary_zero": topology_quantity(
                "boundary_of_boundary_zero", certification("boundary_of_boundary:")),
            "beta0": topology_quantity("beta0", self.betti_numbers[0] if self.betti_numbers else None),
            "beta1": topology_quantity("beta1", self.betti_numbers[1] if self.betti_numbers else None),
            "beta2": topology_quantity("beta2", self.betti_numbers[2] if self.betti_numbers else None),
        }
        topology = TopologicalCellComplex(
            quantities=topology_quantities,
            status=topology_status,
            provenance=provenance,
        )
        mean_curvature = None
        gaussian_curvature = None
        if self.curvature is not None:
            curvature = self.curvature

            def curvature_quantity(name, value, units, scope, quantity_status):
                quantity_status = Status(quantity_status)
                if value is None and quantity_status in {Status.CERTIFIED, Status.PARTIAL}:
                    quantity_status = Status.UNRESOLVED
                return Quantity(
                    name, value, units, quantity_status, curvature.method, scope,
                    provenance=provenance,
                )

            mean_status = Status(curvature.mean_status)
            gaussian_status = Status(curvature.gaussian_status)
            mean_curvature = MeanCurvature(
                curvature_quantity("integrated_mean_curvature", curvature.total_integrated_mean_curvature, "Å", "canonical_surface", mean_status),
                curvature_quantity("complete_mean_curvature", curvature.total_integrated_mean_curvature, "Å", "canonical_surface", mean_status),
                curvature_quantity("smooth_face_contribution", curvature.smooth_mean_curvature, "Å", "spherical_patches", Status.CERTIFIED),
                curvature_quantity("edge_contribution", curvature.seam_mean_curvature, "Å", "physical_refined_seams", mean_status),
                curvature_quantity("boundary_contribution", None, "Å", "free_boundary_convention", Status.NOT_SUPPORTED),
                "outward molecular normal; H=(k1+k2)/2",
            )
            interior_gaussian = (
                curvature.interior_vertex_defects
                if curvature.interior_vertex_defects is not None
                else curvature.intrinsic_discrete_gaussian_curvature
            )
            intrinsic_gaussian = None
            if all(
                value is not None
                for value in (
                    curvature.smooth_gaussian_curvature,
                    curvature.intrinsic_seam_gaussian_curvature,
                    interior_gaussian,
                )
            ):
                intrinsic_gaussian = (
                    curvature.smooth_gaussian_curvature
                    + curvature.intrinsic_seam_gaussian_curvature
                    + interior_gaussian
                )
            gaussian_curvature = GaussianCurvature(
                curvature_quantity("integrated_gaussian_curvature", intrinsic_gaussian, "1", "intrinsic_surface", gaussian_status),
                curvature_quantity("smooth_face_contribution", curvature.smooth_gaussian_curvature, "1", "spherical_patches", Status.CERTIFIED),
                curvature_quantity("intrinsic_discrete_contribution", curvature.intrinsic_discrete_gaussian_curvature, "1", "interior_junctions", gaussian_status),
                curvature_quantity("intrinsic_seam_contribution", curvature.intrinsic_seam_gaussian_curvature, "1", "physical_refined_seams", Status.CERTIFIED if curvature.intrinsic_seam_gaussian_curvature is not None else gaussian_status),
                curvature_quantity("interior_vertex_defects", curvature.interior_vertex_defects, "1", "interior_junctions", gaussian_status),
                curvature_quantity("boundary_geodesic_contribution", curvature.boundary_geodesic_curvature, "1", "physical_boundary_arcs", Status.CERTIFIED if curvature.boundary_geodesic_curvature is not None else gaussian_status),
                curvature_quantity("boundary_corner_contribution", curvature.boundary_corner_turning, "1", "physical_boundary_junctions", gaussian_status),
                curvature_quantity("gauss_bonnet_completion", curvature.completed_gaussian_curvature, "1", "gauss_bonnet_completed_surface", gaussian_status),
                curvature_quantity("gauss_bonnet_expected", curvature.gauss_bonnet_expected, "1", "certified_topology", Status.CERTIFIED if curvature.gauss_bonnet_expected is not None else Status.UNRESOLVED),
                curvature_quantity("gauss_bonnet_residual", curvature.gauss_bonnet_residual, "1", "gauss_bonnet_audit", gaussian_status),
                "outward molecular normal; intrinsic spherical curvature; physical boundary completion; topology-only cuts excluded",
            )
        patch_provenance = tuple(
            MolecularPatchProvenance(
                patch_id=patch.patch_id,
                source_atom_stable_id=patch.carrier_atom_id,
                source_sphere_center=patch.carrier_center,
                source_radius=patch.carrier_radius,
                source_radius_units="Å",
                partner_contact_atom_ids=patch.partner_atom_ids,
                parent_selected_contact_ids=patch.source_contact_ids,
                clipping_occlusion_atom_ids=patch.self_occluding_atom_ids,
                component_id=patch.component_id,
                boundary_arc_ids=patch.boundary_arc_ids,
                junction_ids=patch.junction_ids,
                status=status,
            )
            for patch in self.patches
        )
        return MolecularContactSurfaceResult(
            identity=identity,
            provenance=provenance,
            status=status,
            metrics={
                "area": quantity(
                    "area", self.total_area_A2, "Å²", "molecular_contact_surface", area_status
                ),
            },
            molecular_selection=selection,
            geometry=geometry,
            geometric_counts=geometric_counts,
            topology=topology,
            mean_curvature=mean_curvature,
            gaussian_curvature=gaussian_curvature,
            patch_provenance=patch_provenance,
            interaction_links=interaction_links,
        )


@dataclass
class _Constraint:
    kind: str
    source_index: int
    source_atom_id: str
    normal: Vec3
    threshold: float
    contact_ids: tuple[str, ...] = ()
    circle_id: Optional[str] = None


@dataclass
class _Circle:
    circle_id: str
    normal: Vec3
    offset: float
    constraints: list[_Constraint] = field(default_factory=list)
    nodes: list[int] = field(default_factory=list)
    basis_1: Vec3 = (1.0, 0.0, 0.0)
    basis_2: Vec3 = (0.0, 1.0, 0.0)
    radial: float = 0.0

    @property
    def source_atom_ids(self) -> tuple[str, ...]:
        return tuple(sorted({constraint.source_atom_id for constraint in self.constraints}))

    @property
    def source_contact_ids(self) -> tuple[str, ...]:
        return tuple(sorted({contact for constraint in self.constraints for contact in constraint.contact_ids}))


@dataclass
class _Node:
    direction: Vec3
    source_circle_ids: set[str] = field(default_factory=set)
    synthetic: bool = False


@dataclass
class _Segment:
    segment_id: int
    circle_index: int
    start_node: int
    end_node: int
    theta_start: float
    delta: float


@dataclass
class _HalfEdge:
    halfedge_id: int
    segment_id: int
    reverse: bool
    start_node: int
    end_node: int
    theta_start: float
    delta: float
    tangent_start: Vec3
    tangent_end: Vec3
    twin: int
    next_halfedge: Optional[int] = None
    face_id: Optional[int] = None


def _sphere_relation(radius: float, other_radius: float, distance: float, tolerance: GeometryTolerance) -> str:
    if distance <= tolerance.distance:
        if abs(radius - other_radius) <= tolerance.distance:
            return "coincident"
        return "other_contains" if other_radius > radius else "carrier_contains"
    if distance > radius + other_radius + tolerance.distance:
        return "disjoint"
    if abs(distance - (radius + other_radius)) <= tolerance.distance:
        return "external_tangent"
    if distance < abs(radius - other_radius) - tolerance.distance:
        return "other_contains" if other_radius > radius else "carrier_contains"
    if abs(distance - abs(radius - other_radius)) <= tolerance.distance:
        return "internal_tangent"
    return "partial_overlap"


def _circle_basis(normal: Vec3, tolerance: GeometryTolerance) -> tuple[Vec3, Vec3]:
    references = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    reference = min(references, key=lambda value: abs(_dot(value, normal)))
    first = _unit(_cross(reference, normal), tolerance.angular)
    if first is None:
        raise ValueError("Could not construct a spherical-circle basis")
    second = _cross(normal, first)
    return first, second


def _circle_point(circle: _Circle, theta: float) -> Vec3:
    tangent = _add(
        _scale(circle.basis_1, math.cos(theta)),
        _scale(circle.basis_2, math.sin(theta)),
    )
    return _add(_scale(circle.normal, circle.offset), _scale(tangent, circle.radial))


def _circle_tangent(circle: _Circle, theta: float) -> Vec3:
    return _scale(
        _add(
            _scale(circle.basis_1, -math.sin(theta)),
            _scale(circle.basis_2, math.cos(theta)),
        ),
        circle.radial,
    )


def _circle_intersections(first: _Circle, second: _Circle, tolerance: GeometryTolerance) -> list[Vec3]:
    cross = _cross(first.normal, second.normal)
    cross_norm = _norm(cross)
    if cross_norm <= tolerance.angular:
        return []
    denominator = cross_norm * cross_norm
    base = _scale(
        _add(
            _scale(_cross(second.normal, cross), first.offset),
            _scale(_cross(cross, first.normal), second.offset),
        ),
        1.0 / denominator,
    )
    height_squared = 1.0 - _dot(base, base)
    if height_squared < -tolerance.angular:
        return []
    if height_squared <= tolerance.angular:
        return [base]
    direction = _scale(cross, 1.0 / cross_norm)
    height = math.sqrt(max(0.0, height_squared))
    return [_add(base, _scale(direction, height)), _sub(base, _scale(direction, height))]


def _same_direction(first: Vec3, second: Vec3, tolerance: GeometryTolerance) -> bool:
    return _norm(_sub(first, second)) <= tolerance.angular


def _canonical_shared_circle(
    first: AtomBall,
    second: AtomBall,
    tolerance: GeometryTolerance,
) -> Optional[tuple[str, Vec3, float, Vec3, Vec3, Vec3]]:
    """Return the order-independent circle geometry for two intersecting balls."""
    ordered = tuple(sorted((first, second), key=lambda atom: atom.stable_id))
    first, second = ordered
    displacement = _sub(second.center, first.center)
    distance = _norm(displacement)
    if distance <= tolerance.distance:
        return None
    normal = _unit(displacement, tolerance.angular)
    if normal is None:
        return None
    offset = (
        first.radius * first.radius
        - second.radius * second.radius
        + distance * distance
    ) / (2.0 * distance)
    radial_squared = first.radius * first.radius - offset * offset
    if radial_squared <= tolerance.distance * tolerance.distance:
        return None
    center = _add(first.center, _scale(normal, offset))
    radial = math.sqrt(max(0.0, radial_squared))
    basis_1, basis_2 = _circle_basis(normal, tolerance)
    circle_id = _token(
        "seam-circle",
        (first.stable_id, second.stable_id),
        tuple(round(value, 10) for value in center),
        round(radial, 10),
        tuple(round(value, 10) for value in normal),
    )
    return circle_id, center, radial, normal, basis_1, basis_2


def _circle_theta(
    point: Vec3,
    center: Vec3,
    radial: float,
    basis_1: Vec3,
    basis_2: Vec3,
) -> float:
    direction = _scale(_sub(point, center), 1.0 / radial)
    return math.atan2(_dot(direction, basis_2), _dot(direction, basis_1)) % math.tau


def _split_common_arc(
    arc: CircularArc,
    center: Vec3,
    radial: float,
    basis_1: Vec3,
    basis_2: Vec3,
    tolerance: GeometryTolerance,
) -> list[tuple[float, float, int]]:
    """Map one stored carrier-local arc into [0, 2π), retaining direction."""
    if arc.start_vertex_id is None or arc.end_vertex_id is None:
        return [(0.0, math.tau, 1)]
    span = abs(arc.angular_span)
    if span >= math.tau - tolerance.angular:
        return [(0.0, math.tau, 1 if arc.angular_span >= 0.0 else -1)]
    if span <= tolerance.angular:
        return []
    start = _circle_theta(arc.start_xyz, center, radial, basis_1, basis_2)
    end_point = _circle_theta(arc.end_xyz, center, radial, basis_1, basis_2)
    forward = (end_point - start) % math.tau
    backward = (start - end_point) % math.tau
    # angular_span is measured in the carrier-local basis.  Its sign is not
    # invariant under the canonical circle reparameterization, so recover the
    # direction from the endpoint pair in the shared basis.
    if abs(forward - span) <= abs(backward - span):
        direction, begin, end = 1, start, end_point
    else:
        direction, begin, end = -1, end_point, start
    if abs(end - begin) <= tolerance.angular:
        return []
    if begin < end:
        return [(begin, end, direction)]
    return [(begin, math.tau, direction), (0.0, end, direction)]


def _unique_circle_parameters(values: Iterable[float], tolerance: GeometryTolerance) -> list[float]:
    parameters = [0.0, math.tau]
    for value in values:
        value %= math.tau
        if value <= tolerance.angular or math.tau - value <= tolerance.angular:
            value = 0.0
        parameters.append(value)
    parameters.sort()
    unique: list[float] = []
    for value in parameters:
        if not unique or value - unique[-1] > tolerance.angular:
            unique.append(value)
    if unique[-1] < math.tau:
        unique.append(math.tau)
    return unique


def _interval_covers(
    intervals: Sequence[tuple[float, float, int]],
    midpoint: float,
    tolerance: GeometryTolerance,
) -> list[tuple[int, float, float]]:
    return [
        (direction, begin, end)
        for begin, end, direction in intervals
        if begin - tolerance.angular <= midpoint <= end + tolerance.angular
    ]


def _junction_for_point(
    junctions: list[JunctionVertex],
    point: Vec3,
    carrier_ids: Sequence[str],
    carrier_centers: Mapping[str, Vec3],
    circle_id: str,
    tolerance: GeometryTolerance,
) -> str:
    for index, junction in enumerate(junctions):
        if _norm(_sub(junction.xyz, point)) <= max(1e-7, 100.0 * tolerance.distance):
            circles = tuple(sorted(set(junction.source_circle_ids) | {circle_id}))
            carriers = tuple(sorted(set(junction.source_carrier_atom_ids) | set(carrier_ids)))
            junctions[index] = replace(
                junction,
                source_circle_ids=circles,
                source_carrier_atom_ids=carriers,
            )
            return junction.vertex_id
    vertex_id = _token(
        "junction-global",
        tuple(round(value, 8) for value in point),
    )
    direction = _unit(
        _sub(point, carrier_centers.get(min(carrier_ids), point)),
        tolerance.angular,
    ) or (0.0, 0.0, 1.0)
    junctions.append(JunctionVertex(
        vertex_id,
        min(carrier_ids),
        point,
        direction,
        (circle_id,),
        (),
        "PARTIAL",
        tuple(sorted(set(carrier_ids))),
    ))
    return vertex_id


def _patch_retains_direction(
    patch: SphericalPatch,
    direction: Vec3,
    atom_by_id: Mapping[str, AtomBall],
    tolerance: GeometryTolerance,
) -> bool:
    point = _add(patch.carrier_center, _scale(direction, patch.carrier_radius))
    if any(
        _norm(_sub(point, atom_by_id[atom_id].center))
        < atom_by_id[atom_id].radius - tolerance.distance
        for atom_id in patch.self_occluding_atom_ids
        if atom_id in atom_by_id
    ):
        return False
    return any(
        _norm(_sub(point, atom_by_id[atom_id].center))
        <= atom_by_id[atom_id].radius + tolerance.distance
        for atom_id in patch.partner_atom_ids
        if atom_id in atom_by_id
    )


def _exterior_boundary_direction(
    patch: SphericalPatch,
    point: Vec3,
    circle_normal: Vec3,
    theta: float,
    atom_by_id: Mapping[str, AtomBall],
    tolerance: GeometryTolerance,
) -> Optional[int]:
    surface_normal = _unit(_sub(point, patch.carrier_center), tolerance.angular)
    if surface_normal is None:
        return None
    tangent = _canonical_circle_tangent(circle_normal, theta, tolerance)
    if _norm(tangent) <= tolerance.angular:
        return None
    # The positive boundary direction is the one whose left spherical side
    # is retained.  Testing both sides makes the convention independent of
    # carrier ordering and distinguishes outer from inner loops.
    step = max(1e-7, 100.0 * tolerance.distance)
    retained_sides: dict[int, bool] = {}
    for direction, signed_tangent in ((1, tangent), (-1, _scale(tangent, -1.0))):
        left = _unit(_cross(surface_normal, signed_tangent), tolerance.angular)
        if left is None:
            return None
        sample = _unit(_add(surface_normal, _scale(left, step)), tolerance.angular)
        if sample is None:
            return None
        retained_sides[direction] = _patch_retains_direction(
            patch, sample, atom_by_id, tolerance
        )
    if retained_sides[1] == retained_sides[-1]:
        return None
    return 1 if retained_sides[1] else -1


def _refine_shared_circle(
    first: AtomBall,
    second: AtomBall,
    arcs: Sequence[CircularArc],
    arc_patch: Mapping[str, str],
    patch_by_id: Mapping[str, SphericalPatch],
    junctions: list[JunctionVertex],
    tolerance: GeometryTolerance,
    atom_by_id: Optional[Mapping[str, AtomBall]] = None,
) -> tuple[list[RefinedSeamPiece], list[str]]:
    geometry = _canonical_shared_circle(first, second, tolerance)
    pair = tuple(sorted((first.stable_id, second.stable_id)))
    if geometry is None:
        return [], [f"unresolved_shared_circle:{pair[0]}:{pair[1]}"]
    circle_id, center, radial, normal, basis_1, basis_2 = geometry
    mapped: dict[str, list[tuple[float, float, int]]] = {}
    parameters: list[float] = []
    for arc in arcs:
        intervals = _split_common_arc(
            arc, center, radial, basis_1, basis_2, tolerance,
        )
        mapped[arc.arc_id] = intervals
        for begin, end, _ in intervals:
            parameters.extend((begin, end))
    breaks = _unique_circle_parameters(parameters, tolerance)
    pieces: list[RefinedSeamPiece] = []
    for begin, end in zip(breaks, breaks[1:]):
        if end - begin <= tolerance.angular:
            continue
        midpoint = 0.5 * (begin + end)
        first_hits = []
        second_hits = []
        for arc in arcs:
            hits = _interval_covers(mapped[arc.arc_id], midpoint, tolerance)
            if not hits:
                continue
            (first_hits if arc.carrier_atom_id == first.stable_id else second_hits).append((arc, hits[0][0]))
        if not first_hits and not second_hits:
            continue
        classification = (
            "internal" if first_hits and second_hits
            else "exterior_first" if first_hits
            else "exterior_second"
        )
        hits = []
        point = _add(center, _scale(
            _add(_scale(basis_1, math.cos(midpoint)), _scale(basis_2, math.sin(midpoint))),
            radial,
        ))
        for arc, direction in first_hits + second_hits:
            # A closed shared circle is traversed in opposite directions by
            # the two incident spherical patches.
            if (
                classification == "internal"
                and arc.start_vertex_id is None
                and arc.end_vertex_id is None
                and arc.carrier_atom_id == first.stable_id
            ):
                # Interior sides oppose; a retained exterior side traverses
                # its occlusion boundary with the opposite orientation.
                direction = -direction
            elif (
                classification != "internal"
                and arc.arc_id in arc_patch
            ):
                patch = patch_by_id.get(arc_patch[arc.arc_id])
                if patch is not None and atom_by_id is not None:
                    orientation = _exterior_boundary_direction(
                        patch,
                        point,
                        normal,
                        midpoint,
                        atom_by_id,
                        tolerance,
                    )
                    if orientation is None:
                        return pieces, [f"unresolved_boundary_orientation:{circle_id}"]
                    direction = orientation
                elif patch is not None:
                    # Preserve the private diagnostic helper's legacy call
                    # shape; production builds always provide atom geometry.
                    if (
                        classification == "exterior_second"
                        and arc.start_vertex_id is None
                        and arc.end_vertex_id is None
                        and arc.carrier_atom_id == second.stable_id
                    ):
                        direction = -direction
            hits.append((arc, direction))
        source_arc_ids = tuple(sorted(arc.arc_id for arc, _ in hits))
        patch_ids = tuple(sorted({arc_patch[arc.arc_id] for arc, _ in hits if arc.arc_id in arc_patch}))
        source_contacts_set = set()
        for arc, _ in hits:
            source_contacts_set.update(arc.source_contact_ids)
            patch_id = arc_patch.get(arc.arc_id)
            if patch_id in patch_by_id:
                source_contacts_set.update(patch_by_id[patch_id].source_contact_ids)
        source_contacts = tuple(sorted(source_contacts_set))
        orientations = tuple(sorted({
            (arc_patch[arc.arc_id], direction)
            for arc, direction in hits
            if arc.arc_id in arc_patch
        }))
        status = "CERTIFIED"
        if len({patch_id for patch_id, _ in orientations}) != len(orientations):
            status = "PARTIAL"
        start = _add(center, _scale(
            _add(_scale(basis_1, math.cos(begin)), _scale(basis_2, math.sin(begin))),
            radial,
        ))
        end_point = _add(center, _scale(
            _add(_scale(basis_1, math.cos(end)), _scale(basis_2, math.sin(end))),
            radial,
        ))
        centers = {first.stable_id: first.center, second.stable_id: second.center}
        start_id = _junction_for_point(junctions, start, pair, centers, circle_id, tolerance)
        end_id = _junction_for_point(junctions, end_point, pair, centers, circle_id, tolerance)
        seam_id = _token(
            "refined-seam",
            circle_id,
            round(begin, 12),
            round(end, 12),
        )
        pieces.append(RefinedSeamPiece(
            seam_id,
            pair,
            circle_id,
            center,
            radial,
            normal,
            (begin, end),
            (start_id, end_id),
            start,
            end_point,
            patch_ids,
            source_arc_ids,
            source_contacts,
            orientations,
            classification,
            status,
        ))
    return pieces, []


def _classify_local_links(
    junctions: Sequence[JunctionVertex],
    edges: Sequence[tuple[str, str, str, tuple[str, ...]]],
    arcs_by_id: Mapping[str, CircularArc],
) -> tuple[list[JunctionVertex], dict[str, str], int, int]:
    """Classify the combinatorial link of every global junction.

    Each spherical cell contributes one link edge between the two incident
    surface edges.  A cycle is an interior link; a path is a boundary link.
    This deliberately leaves higher-valence links singular instead of
    repairing them by tolerance or by triangulation.
    """
    edge_endpoints = {edge_id: (start, end) for edge_id, start, end, _ in edges}
    edge_patches = {edge_id: patches for edge_id, _, _, patches in edges}
    incident_edges: dict[str, set[str]] = defaultdict(set)
    incident_patches: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for edge_id, start, end, patch_ids in edges:
        for vertex in (start, end):
            incident_edges[vertex].add(edge_id)
            for patch_id in patch_ids:
                incident_patches[vertex][patch_id].add(edge_id)
    arc_contact = {
        arc_id: bool(arc.source_contact_ids)
        for arc_id, arc in arcs_by_id.items()
    }
    result = []
    statuses: dict[str, str] = {}
    manifold = singular = 0
    for junction in junctions:
        vertex = junction.vertex_id
        edge_ids = tuple(sorted(incident_edges.get(vertex, ())))
        if junction.local_link in {"full_sphere", "full_circle_boundary"}:
            link = "interior_cycle" if junction.local_link == "full_sphere" else "boundary_path"
            kind = junction.local_link
            statuses[vertex] = f"{link}:{kind}"
            manifold += 1
            result.append(replace(
                junction,
                incident_edge_ids=(),
                local_link=f"{link}:{kind}",
                status="CERTIFIED",
            ))
            continue
        if (
            len(edge_ids) == 1
            and edge_endpoints[edge_ids[0]] == (vertex, vertex)
            and len(edge_patches[edge_ids[0]]) == 2
        ):
            # One loop edge shared by two patches is a closed seam cycle,
            # not a four-valent branch caused by counting its two ends.
            incident_contact = any(
                arc_contact.get(arc_id, False)
                for arc_id in junction.incident_arc_ids
            )
            if len(junction.source_carrier_atom_ids) >= 3:
                kind = "same_group_triple_junction"
            elif incident_contact:
                kind = "mixed_carrier_partner_junction"
            else:
                kind = "ordinary_seam_subdivision"
            statuses[vertex] = f"interior_cycle:{kind}"
            manifold += 1
            result.append(replace(
                junction,
                incident_edge_ids=edge_ids,
                local_link=f"interior_cycle:{kind}",
                status="CERTIFIED",
            ))
            continue
        if (
            len(edge_ids) == 1
            and edge_endpoints[edge_ids[0]] == (vertex, vertex)
            and len(edge_patches[edge_ids[0]]) == 1
        ):
            # A single closed boundary circle is smooth at its bookkeeping
            # vertex; its local link is a boundary path, not an interior defect.
            incident_contact = any(
                arc_contact.get(arc_id, False)
                for arc_id in junction.incident_arc_ids
            )
            if len(junction.source_carrier_atom_ids) >= 3:
                kind = "same_group_triple_junction"
            elif incident_contact:
                kind = "mixed_carrier_partner_junction"
            else:
                kind = "ordinary_seam_subdivision"
            statuses[vertex] = f"boundary_path:{kind}"
            manifold += 1
            result.append(replace(
                junction,
                incident_edge_ids=edge_ids,
                local_link=f"boundary_path:{kind}",
                status="CERTIFIED",
            ))
            continue
        link_adjacency: dict[str, list[str]] = defaultdict(list)
        invalid = not edge_ids
        for patch_edges in incident_patches.get(vertex, {}).values():
            local = sorted(patch_edges)
            if len(local) == 2:
                left, right = local
                link_adjacency[left].append(right)
                link_adjacency[right].append(left)
            elif len(local) == 1 and edge_endpoints[local[0]] == (vertex, vertex):
                # A full, un-subdivided seam is a legitimate loop 1-cell.
                link_adjacency[local[0]].extend((local[0], local[0]))
            else:
                invalid = True
        degrees = {edge_id: len(link_adjacency.get(edge_id, ())) for edge_id in edge_ids}
        if not invalid and edge_ids:
            seen = set()
            components = 0
            for root in edge_ids:
                if root in seen:
                    continue
                components += 1
                stack = [root]
                while stack:
                    current = stack.pop()
                    if current in seen:
                        continue
                    seen.add(current)
                    stack.extend(link_adjacency.get(current, ()))
            if components == 1 and all(value == 2 for value in degrees.values()):
                link = "interior_cycle"
            elif (
                components == 1
                and sum(value == 1 for value in degrees.values()) == 2
                and all(value in {1, 2} for value in degrees.values())
            ):
                link = "boundary_path"
            else:
                link = "singular_branch"
        else:
            link = "unresolved"
        if link in {"interior_cycle", "boundary_path"}:
            manifold += 1
        else:
            singular += 1
        incident_contact = any(
            arc_contact.get(arc_id, False)
            for arc_id in junction.incident_arc_ids
        )
        if len(junction.source_carrier_atom_ids) >= 3:
            kind = "same_group_triple_junction"
        elif incident_contact:
            kind = "mixed_carrier_partner_junction"
        else:
            kind = "ordinary_seam_subdivision"
        statuses[vertex] = f"{link}:{kind}"
        result.append(replace(
            junction,
            incident_edge_ids=edge_ids,
            local_link=f"{link}:{kind}",
            status="CERTIFIED" if link in {"interior_cycle", "boundary_path"} else "PARTIAL",
        ))
    return result, statuses, manifold, singular


def _ordered_patch_boundary_cycles(
    patch_id: str,
    edges: Sequence[tuple[str, str, str, tuple[str, ...]]],
    refined_by_id: Mapping[str, RefinedSeamPiece],
) -> tuple[tuple[tuple[str, ...], tuple[str, ...]], ...]:
    """Return canonical directed boundary cycles for one spherical region."""
    directed: list[tuple[str, str, str]] = []
    for edge_id, start, end, patch_ids in edges:
        if patch_id not in patch_ids:
            continue
        refined = refined_by_id.get(edge_id)
        if refined is None:
            directed.append((edge_id, start, end))
            continue
        orientation = dict(refined.orientation_on_patches).get(patch_id, 1)
        directed.append((
            edge_id,
            start if orientation >= 0 else end,
            end if orientation >= 0 else start,
        ))
    outgoing: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for edge in directed:
        outgoing[edge[1]].append(edge)
    for values in outgoing.values():
        values.sort(key=lambda edge: edge[0])
    cycles: list[tuple[tuple[str, ...], tuple[str, ...]]] = []
    visited: set[str] = set()
    for edge in sorted(directed, key=lambda value: value[0]):
        if edge[0] in visited:
            continue
        current = edge
        edge_ids: list[str] = []
        vertices: list[str] = [current[1]]
        while current[0] not in visited:
            visited.add(current[0])
            edge_ids.append(current[0])
            vertices.append(current[2])
            candidates = outgoing.get(current[2], ())
            if len(candidates) != 1:
                break
            current = candidates[0]
        if vertices[-1] == vertices[0]:
            vertices.pop()
        cycles.append((tuple(edge_ids), tuple(vertices)))
    cycles.sort(key=lambda cycle: (min(cycle[1]) if cycle[1] else "", len(cycle[0]), cycle[0]))
    return tuple(cycles)


def _make_topological_cuts(
    patches: Sequence[SphericalPatch],
    edges: Sequence[tuple[str, str, str, tuple[str, ...]]],
    refined_by_id: Mapping[str, RefinedSeamPiece],
    junctions: Sequence[JunctionVertex],
) -> tuple[
    tuple[TopologicalCut, ...],
    tuple[tuple[str, int], ...],
    list[str],
]:
    """Make auxiliary edges for multiply-boundary face attachments.

    A cut has no physical area and is not a molecular seam or curvature edge.
    It connects an extra boundary cycle to the deterministic root cycle in
    the attaching word.  The same face traverses that edge twice in opposite
    directions, so the edge is a same-face double incidence, not a zero-face
    edge.  This preserves the physical Euler characteristic while giving the
    global cell model a connected face attachment.
    """
    junction_by_id = {junction.vertex_id: junction for junction in junctions}
    cuts: list[TopologicalCut] = []
    counts: list[tuple[str, int]] = []
    diagnostics: list[str] = []
    for patch in patches:
        cycles = _ordered_patch_boundary_cycles(patch.patch_id, edges, refined_by_id)
        counts.append((patch.patch_id, len(cycles)))
        if len(cycles) <= 1:
            continue
        anchors = [min(cycle[1]) for cycle in cycles if cycle[1]]
        if len(anchors) != len(cycles):
            diagnostics.append(f"unresolved_patch_boundary_cycle:{patch.patch_id}")
            continue
        root = anchors[0]
        for cycle_index, anchor in enumerate(anchors[1:], start=1):
            start = junction_by_id.get(root)
            end = junction_by_id.get(anchor)
            if start is None or end is None:
                diagnostics.append(f"unresolved_patch_cut:{patch.patch_id}:{cycle_index}")
                continue
            cuts.append(TopologicalCut(
                _token("topological-cut", patch.patch_id, 0, cycle_index, root, anchor),
                patch.patch_id,
                patch.carrier_atom_id,
                (0, cycle_index),
                (root, anchor),
                start.xyz,
                end.xyz,
            ))
    return tuple(cuts), tuple(sorted(counts)), diagnostics


def _rank_mod2(rows: Sequence[Sequence[int]], width: int) -> int:
    rows = [sum(1 << index for index in row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((index for index in range(rank, len(rows)) if rows[index] & (1 << column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index] & (1 << column):
                rows[index] ^= rows[rank]
        rank += 1
    return rank


def _cell_edge_incidence_distribution(
    edges: Sequence[tuple[str, str, str, tuple[str, ...]]],
    cycles_by_patch: Mapping[str, Sequence[tuple[tuple[str, ...], tuple[str, ...]]]],
    topological_cuts: Sequence[TopologicalCut],
) -> tuple[tuple[tuple[str, int], ...], tuple[str, ...]]:
    """Audit face-edge attachment, including the two-sided auxiliary cuts.

    A cut is an attachment edge of one face on both sides.  It is therefore
    counted twice for that same face; it is not a zero-face edge.
    """
    occurrences: dict[str, list[str]] = defaultdict(list)
    for patch_id, cycles in cycles_by_patch.items():
        for edge_ids, _vertices in cycles:
            for edge_id in edge_ids:
                occurrences[edge_id].append(patch_id)
    for cut in topological_cuts:
        occurrences[cut.cut_id].extend((cut.patch_id, cut.patch_id))

    edge_ids = {edge_id for edge_id, _start, _end, _patch_ids in edges}
    counts = defaultdict(int)
    issues: list[str] = []
    for edge_id in sorted(edge_ids):
        faces = occurrences.get(edge_id, ())
        if len(faces) == 0:
            category = "zero_face"
            issues.append(f"zero_face:{edge_id}")
        elif len(faces) == 1:
            category = "one_face"
        elif len(faces) == 2 and faces[0] == faces[1]:
            category = "same_face_double"
        elif len(faces) == 2:
            category = "two_face"
        else:
            category = "more_than_two_faces"
            issues.append(f"more_than_two_faces:{edge_id}")
        counts[category] += 1
    return tuple(sorted(counts.items())), tuple(issues)


def _canonical_circle_tangent(normal: Vec3, theta: float, tolerance: GeometryTolerance) -> Vec3:
    basis_1, basis_2 = _circle_basis(normal, tolerance)
    tangent = _add(
        _scale(basis_1, -math.sin(theta)),
        _scale(basis_2, math.cos(theta)),
    )
    return _unit(tangent, tolerance.angular) or (0.0, 0.0, 0.0)


def _curvature_edge_records(surface: MolecularContactSurface) -> tuple[tuple[str, object, tuple[str, ...]], ...]:
    """Return physical edges only; topology-only cuts never enter this list."""
    arc_patch = {
        arc_id: patch.patch_id
        for patch in surface.patches
        for arc_id in patch.boundary_arc_ids
    }
    records: list[tuple[str, object, tuple[str, ...]]] = []
    records.extend(
        ("seam", piece, piece.incident_patch_ids)
        for piece in surface.refined_seams
        if piece.incident_patch_ids
    )
    records.extend(
        ("arc", arc, (arc_patch[arc.arc_id],))
        for arc in surface.arcs
        if arc.source_contact_ids and arc.arc_id in arc_patch
    )
    return tuple(records)


def _curvature_edge_endpoints(kind: str, edge: object) -> tuple[str, str]:
    if kind == "seam":
        piece = edge
        return piece.endpoint_junction_ids
    arc = edge
    return (arc.start_vertex_id, arc.end_vertex_id)


def _curvature_edge_direction(kind: str, edge: object, patch_id: str) -> int:
    if kind == "seam":
        return int(dict(edge.orientation_on_patches).get(patch_id, 0))
    return 1 if edge.angular_span >= 0.0 else -1


def _curvature_edge_directed_start(kind: str, edge: object, patch_id: str) -> str:
    start_id, end_id = _curvature_edge_endpoints(kind, edge)
    direction = _curvature_edge_direction(kind, edge, patch_id)
    return start_id if direction >= 0 else end_id


def _curvature_patch_tangent(
    kind: str,
    edge: object,
    patch: SphericalPatch,
    vertex_id: str,
    tolerance: GeometryTolerance,
) -> Optional[Vec3]:
    """Return the tangent leaving a vertex along a patch boundary."""
    start_id, end_id = _curvature_edge_endpoints(kind, edge)
    if start_id != vertex_id and end_id != vertex_id:
        return None
    if start_id == end_id == vertex_id:
        return None
    if kind == "seam":
        piece = edge
        begin, end = piece.parameter_interval
        direction = _curvature_edge_direction(kind, piece, patch.patch_id)
        directed_start = _curvature_edge_directed_start(kind, piece, patch.patch_id)
        theta = begin if start_id == vertex_id else end
        tangent = _canonical_circle_tangent(piece.circle_normal, theta, tolerance)
        if direction == 0:
            return None
        forward = _scale(tangent, float(direction))
        return forward if directed_start == vertex_id else _scale(forward, -1.0)

    arc = edge
    if arc.start_xyz is None or arc.end_xyz is None:
        return None
    radius = patch.carrier_radius
    direction = 1 if arc.angular_span >= 0.0 else -1
    circle_normal = _unit(arc.circle_normal, tolerance.angular)
    if circle_normal is None or abs(arc.circle_offset) > 1.0 + tolerance.angular:
        return None
    radial = radius * math.sqrt(max(0.0, 1.0 - arc.circle_offset * arc.circle_offset))
    if radial <= tolerance.distance:
        return None
    circle_center = _add(
        patch.carrier_center,
        _scale(circle_normal, radius * arc.circle_offset),
    )
    basis_1, basis_2 = _circle_basis(circle_normal, tolerance)
    point = arc.start_xyz if start_id == vertex_id else arc.end_xyz
    theta = _circle_theta(point, circle_center, radial, basis_1, basis_2)
    tangent = _canonical_circle_tangent(circle_normal, theta, tolerance)
    forward = _scale(tangent, float(direction))
    directed_start = _curvature_edge_directed_start(kind, arc, patch.patch_id)
    return forward if directed_start == vertex_id else _scale(forward, -1.0)


def _curvature_patch_sector_angle(
    patch: SphericalPatch,
    patch_id: str,
    vertex_id: str,
    candidates: Sequence[tuple[str, object]],
    junction: JunctionVertex,
    tolerance: GeometryTolerance,
) -> tuple[Optional[float], Optional[str]]:
    if len(candidates) == 1 and _curvature_edge_endpoints(*candidates[0])[0] == vertex_id == _curvature_edge_endpoints(*candidates[0])[1]:
        return math.pi, None
    if len(candidates) != 2:
        return None, f"curvature_patch_sector:{vertex_id}:{patch_id}:{len(candidates)}"
    tangents = [
        _curvature_patch_tangent(kind, edge, patch, vertex_id, tolerance)
        for kind, edge in candidates
    ]
    if any(tangent is None for tangent in tangents):
        return None, f"curvature_patch_tangent:{vertex_id}:{patch_id}"
    outgoing = []
    incoming = []
    for candidate, tangent in zip(candidates, tangents):
        if _curvature_edge_directed_start(candidate[0], candidate[1], patch_id) == vertex_id:
            outgoing.append(tangent)
        else:
            incoming.append(tangent)
    if len(outgoing) != 1 or len(incoming) != 1:
        return None, f"curvature_patch_direction:{vertex_id}:{patch_id}"
    normal = _unit(_sub(junction.xyz, patch.carrier_center), tolerance.angular)
    if normal is None:
        return None, f"curvature_patch_normal:{vertex_id}:{patch_id}"
    turn = math.atan2(
        _dot(_cross(_scale(incoming[0], -1.0), outgoing[0]), normal),
        _dot(_scale(incoming[0], -1.0), outgoing[0]),
    )
    angle = math.pi - turn
    if angle < -tolerance.angular or angle > math.tau + tolerance.angular:
        return None, f"curvature_patch_angle:{vertex_id}:{patch_id}"
    return _clamp(angle, 0.0, math.tau), None


def _curvature_vertex_terms(
    surface: MolecularContactSurface,
    edges: Sequence[tuple[str, object, tuple[str, ...]]],
    tolerance: GeometryTolerance,
    component_id: Optional[str] = None,
) -> tuple[Optional[float], Optional[float], int, list[str]]:
    """Compute intrinsic vertex defects and exterior boundary corner terms."""
    patch_by_id = {patch.patch_id: patch for patch in surface.patches}
    sectors: dict[str, list[float]] = defaultdict(list)
    incident: dict[tuple[str, str], list[tuple[str, object]]] = defaultdict(list)
    for kind, edge, patch_ids in edges:
        start_id, end_id = _curvature_edge_endpoints(kind, edge)
        for patch_id in patch_ids:
            if patch_id not in patch_by_id:
                continue
            if component_id is not None and patch_by_id[patch_id].component_id != component_id:
                continue
            incident[(patch_id, start_id)].append((kind, edge))
            incident[(patch_id, end_id)].append((kind, edge))

    diagnostics: list[str] = []
    audited = 0
    for junction in surface.junctions:
        vertex_id = junction.vertex_id
        patch_angles: list[float] = []
        patch_ids = sorted({patch_id for (patch_id, vertex) in incident if vertex == vertex_id})
        for patch_id in patch_ids:
            patch = patch_by_id[patch_id]
            candidates = incident[(patch_id, vertex_id)]
            unique = {getattr(edge, "seam_id", getattr(edge, "arc_id", "")): (kind, edge) for kind, edge in candidates}
            candidates = list(unique.values())
            angle, diagnostic = _curvature_patch_sector_angle(
                patch,
                patch_id,
                vertex_id,
                candidates,
                junction,
                tolerance,
            )
            if diagnostic is not None:
                diagnostics.append(diagnostic)
                continue
            patch_angles.append(angle)

        if patch_ids and len(patch_angles) == len(patch_ids):
            audited += 1
            total_angle = sum(patch_angles)
            if junction.local_link.startswith("interior_cycle"):
                sectors["interior"].append(2.0 * math.pi - total_angle)
            elif junction.local_link.startswith("boundary_path"):
                sectors["boundary"].append(math.pi - total_angle)
            else:
                diagnostics.append(f"curvature_junction_link:{vertex_id}:{junction.local_link}")
        elif patch_ids:
            diagnostics.append(f"curvature_junction_incomplete:{vertex_id}")

    if diagnostics:
        return None, None, audited, diagnostics
    return sum(sectors["interior"]), sum(sectors["boundary"]), audited, []


def audit_molecular_contact_surface_gauss_bonnet(
    surface: MolecularContactSurface,
    tolerance: GeometryTolerance | None = None,
) -> dict:
    """Audit Gauss--Bonnet independently on each spherical patch record."""
    tolerance = tolerance or GeometryTolerance()
    patch_by_id = {patch.patch_id: patch for patch in surface.patches}
    junction_by_id = {junction.vertex_id: junction for junction in surface.junctions}
    edges = _curvature_edge_records(surface)
    incident: dict[tuple[str, str], list[tuple[str, object]]] = defaultdict(list)
    for kind, edge, patch_ids in edges:
        start_id, end_id = _curvature_edge_endpoints(kind, edge)
        for patch_id in patch_ids:
            incident[(patch_id, start_id)].append((kind, edge))
            incident[(patch_id, end_id)].append((kind, edge))

    cycle_counts = dict(surface.face_boundary_cycle_counts)
    residuals: list[float] = []
    diagnostics: list[str] = []
    for patch in surface.patches:
        patch_edges = [
            (kind, edge, patch_ids)
            for kind, edge, patch_ids in edges
            if patch.patch_id in patch_ids
        ]
        smooth = patch.area_A2 / (patch.carrier_radius * patch.carrier_radius)
        boundary = 0.0
        vertices: set[str] = set()
        for kind, edge, _patch_ids in patch_edges:
            start_id, end_id = _curvature_edge_endpoints(kind, edge)
            vertices.update((start_id, end_id))
            if kind == "arc":
                boundary += edge.circle_offset * edge.angular_span
                continue
            direction = dict(edge.orientation_on_patches).get(patch.patch_id, 0)
            normal = _unit(edge.circle_normal, tolerance.angular)
            if direction == 0 or normal is None or patch.carrier_radius <= tolerance.distance:
                diagnostics.append(f"patchwise_edge:{patch.patch_id}")
                continue
            offset = _dot(_sub(edge.circle_center, patch.carrier_center), normal) / patch.carrier_radius
            boundary += _clamp(offset) * direction * (edge.parameter_interval[1] - edge.parameter_interval[0])

        corners = 0.0
        for vertex_id in vertices:
            junction = junction_by_id.get(vertex_id)
            if junction is None:
                diagnostics.append(f"patchwise_vertex:{patch.patch_id}:{vertex_id}")
                continue
            candidates = incident.get((patch.patch_id, vertex_id), [])
            unique = {
                getattr(edge, "seam_id", getattr(edge, "arc_id", "")): (kind, edge)
                for kind, edge in candidates
            }
            angle, diagnostic = _curvature_patch_sector_angle(
                patch,
                patch.patch_id,
                vertex_id,
                tuple(unique.values()),
                junction,
                tolerance,
            )
            if diagnostic is not None:
                diagnostics.append(f"patchwise:{diagnostic}")
                continue
            corners += math.pi - angle

        boundaries = cycle_counts.get(patch.patch_id, 1 if vertices else 0)
        chi = 2 - boundaries
        residuals.append(smooth + boundary + corners - 2.0 * math.pi * chi)

    maximum = max((abs(value) for value in residuals), default=0.0)
    return {
        "patch_count": len(surface.patches),
        "max_abs_residual": maximum,
        "sum_residual": sum(residuals),
        "uncertified_patch_count": sum(
            abs(value) > max(1e-8, 100.0 * tolerance.area)
            for value in residuals
        ),
        "diagnostics": tuple(sorted(set(diagnostics))),
        "status": "CERTIFIED" if not diagnostics and maximum <= max(1e-8, 100.0 * tolerance.area) else "PARTIAL",
    }


def compute_molecular_contact_surface_curvature(
    surface: MolecularContactSurface,
    tolerance: GeometryTolerance | None = None,
) -> MolecularContactSurfaceCurvature:
    """Integrate canonical spherical curvature and piecewise GB terms.

    For the oriented patch decomposition,

        ∫K dA + Σ∫_seam k_g ds + Σ∫_boundary k_g ds
            + Σ corner_turn = 2πχ.

    The seam term is intrinsic geodesic curvature on both incident patches;
    it is not the extrinsic dihedral angle used by mean curvature.  Auxiliary
    topology cuts never enter the edge records.
    """
    tolerance = tolerance or GeometryTolerance()
    patch_by_id = {patch.patch_id: patch for patch in surface.patches}
    mean_diagnostics: list[str] = []
    gaussian_diagnostics: list[str] = []

    smooth_h = sum(
        patch.area_A2 / patch.carrier_radius
        for patch in surface.patches
        if patch.carrier_radius > tolerance.distance
    )
    smooth_k = sum(
        patch.area_A2 / (patch.carrier_radius * patch.carrier_radius)
        for patch in surface.patches
        if patch.carrier_radius > tolerance.distance
    )
    if any(patch.carrier_radius <= tolerance.distance for patch in surface.patches):
        mean_diagnostics.append("curvature_invalid_carrier_radius")
        gaussian_diagnostics.append("curvature_invalid_carrier_radius")

    edges = _curvature_edge_records(surface)
    seam_h = 0.0
    physical_seam_count = 0
    for _kind, edge, patch_ids in edges:
        if _kind != "seam" or edge.classification != "internal":
            continue
        physical_seam_count += 1
        if len(patch_ids) != 2 or any(patch_id not in patch_by_id for patch_id in patch_ids):
            mean_diagnostics.append(f"curvature_internal_seam_incidence:{edge.seam_id}")
            gaussian_diagnostics.append(f"curvature_internal_seam_incidence:{edge.seam_id}")
            continue
        ordered = sorted((patch_by_id[patch_ids[0]], patch_by_id[patch_ids[1]]), key=lambda patch: patch.carrier_atom_id)
        first, second = ordered
        distance = _norm(_sub(second.carrier_center, first.carrier_center))
        if distance <= tolerance.distance:
            mean_diagnostics.append(f"curvature_internal_seam_centers:{edge.seam_id}")
            gaussian_diagnostics.append(f"curvature_internal_seam_centers:{edge.seam_id}")
            continue
        begin, end = edge.parameter_interval
        midpoint = 0.5 * (begin + end)
        basis_1, basis_2 = _circle_basis(edge.circle_normal, tolerance)
        point = _add(
            edge.circle_center,
            _scale(
                _add(_scale(basis_1, math.cos(midpoint)), _scale(basis_2, math.sin(midpoint))),
                edge.circle_radius,
            ),
        )
        normal_first = _unit(_sub(point, first.carrier_center), tolerance.angular)
        normal_second = _unit(_sub(point, second.carrier_center), tolerance.angular)
        tangent = _canonical_circle_tangent(edge.circle_normal, midpoint, tolerance)
        if normal_first is None or normal_second is None:
            mean_diagnostics.append(f"curvature_internal_seam_normals:{edge.seam_id}")
            gaussian_diagnostics.append(f"curvature_internal_seam_normals:{edge.seam_id}")
            continue
        beta = math.atan2(
            _dot(tangent, _cross(normal_first, normal_second)),
            _dot(normal_first, normal_second),
        )
        seam_h += 0.5 * beta * edge.circle_radius * (end - begin)

    boundary_geodesic = 0.0
    # Internal seam sides add their intrinsic geodesic curvature.  The
    # outward-normal dihedral beta above is a separate mean-curvature term.
    intrinsic_seam = 0.0
    intrinsic_seam_by_component: dict[str, float] = defaultdict(float)
    boundary_geodesic_by_component: dict[str, float] = defaultdict(float)
    exterior_boundary_count = 0
    for kind, edge, patch_ids in edges:
        if kind == "seam" and edge.classification == "internal":
            if len(patch_ids) != 2 or any(patch_id not in patch_by_id for patch_id in patch_ids):
                gaussian_diagnostics.append(f"curvature_internal_seam_incidence:{edge.seam_id}")
                continue
            begin, end = edge.parameter_interval
            normal = _unit(edge.circle_normal, tolerance.angular)
            if normal is None:
                gaussian_diagnostics.append(f"curvature_internal_seam_circle:{edge.seam_id}")
                continue
            for patch_id in patch_ids:
                patch = patch_by_id[patch_id]
                direction = dict(edge.orientation_on_patches).get(patch_id, 0)
                if direction == 0 or patch.carrier_radius <= tolerance.distance:
                    gaussian_diagnostics.append(f"curvature_internal_seam_orientation:{edge.seam_id}")
                    continue
                offset = _dot(_sub(edge.circle_center, patch.carrier_center), normal) / patch.carrier_radius
                contribution = _clamp(offset) * direction * (end - begin)
                intrinsic_seam += contribution
                intrinsic_seam_by_component[patch.component_id] += contribution
            continue
        if len(patch_ids) != 1 or patch_ids[0] not in patch_by_id:
            gaussian_diagnostics.append(f"curvature_boundary_incidence:{getattr(edge, 'seam_id', getattr(edge, 'arc_id', 'unknown'))}")
            continue
        patch = patch_by_id[patch_ids[0]]
        if kind == "arc":
            contribution = edge.circle_offset * edge.angular_span
            boundary_geodesic += contribution
            boundary_geodesic_by_component[patch.component_id] += contribution
        else:
            begin, end = edge.parameter_interval
            direction = dict(edge.orientation_on_patches).get(patch.patch_id, 0)
            if direction == 0:
                gaussian_diagnostics.append(f"curvature_boundary_orientation:{edge.seam_id}")
                continue
            normal = _unit(edge.circle_normal, tolerance.angular)
            if normal is None or patch.carrier_radius <= tolerance.distance:
                gaussian_diagnostics.append(f"curvature_boundary_circle:{edge.seam_id}")
                continue
            offset = _dot(_sub(edge.circle_center, patch.carrier_center), normal) / patch.carrier_radius
            contribution = _clamp(offset) * direction * (end - begin)
            boundary_geodesic += contribution
            boundary_geodesic_by_component[patch.component_id] += contribution
        exterior_boundary_count += 1

    interior_defect, boundary_corner, audited_junctions, vertex_diagnostics = _curvature_vertex_terms(
        surface,
        edges,
        tolerance,
    )
    gaussian_diagnostics.extend(vertex_diagnostics)
    component_diagnostics: list[str] = []
    component_smooth: dict[str, float] = defaultdict(float)
    component_patches: dict[str, list[SphericalPatch]] = defaultdict(list)
    for patch in surface.patches:
        component_patches[patch.component_id].append(patch)
        if patch.carrier_radius > tolerance.distance:
            component_smooth[patch.component_id] += patch.area_A2 / (patch.carrier_radius * patch.carrier_radius)
    for component_id, patches in sorted(component_patches.items()):
        patch_ids = {patch.patch_id for patch in patches}
        component_edge_ids: set[str] = set()
        component_vertices = {vertex for patch in patches for vertex in patch.junction_ids}
        for kind, edge, edge_patch_ids in edges:
            incident = [patch_id for patch_id in edge_patch_ids if patch_id in patch_ids]
            if not incident:
                continue
            edge_id = getattr(edge, "seam_id", getattr(edge, "arc_id", ""))
            component_edge_ids.add(edge_id)
            component_vertices.update(_curvature_edge_endpoints(kind, edge))
        cycle_counts = [
            dict(surface.face_boundary_cycle_counts).get(
                patch.patch_id, 1 if patch.junction_ids else 0
            )
            for patch in patches
        ]
        requires_cut = any(count > 1 for count in cycle_counts)
        if requires_cut and not surface.topological_cuts:
            # A caller may audit curvature on a representation with its
            # auxiliary cuts stripped.  The physical terms remain valid, but
            # that representation cannot supply a component Euler audit.
            continue
        if requires_cut:
            for cut in surface.topological_cuts:
                if cut.patch_id in patch_ids:
                    component_edge_ids.add(cut.cut_id)
                    component_vertices.update(cut.endpoint_junction_ids)
        component_chi = len(component_vertices) - len(component_edge_ids) + len(patches)
        component_interior, component_corner, _audited, component_vertex_diagnostics = _curvature_vertex_terms(
            surface,
            edges,
            tolerance,
            component_id,
        )
        component_diagnostics.extend(component_vertex_diagnostics)
        if component_interior is None or component_corner is None:
            component_diagnostics.append(f"gauss_bonnet_component_missing:{component_id}")
            continue
        component_completed = (
            component_smooth[component_id]
            + intrinsic_seam_by_component[component_id]
            + boundary_geodesic_by_component[component_id]
            + component_interior
            + component_corner
        )
        component_residual = component_completed - 2.0 * math.pi * component_chi
        if abs(component_residual) > max(1e-7, 100.0 * tolerance.area):
            component_diagnostics.append(
                f"gauss_bonnet_component_residual:{component_id}:{component_residual:.12g}"
            )
    expected = (
        2.0 * math.pi * surface.euler_characteristic
        if surface.euler_characteristic is not None
        else None
    )
    completed = None
    residual = None
    if (
        not gaussian_diagnostics
        and interior_defect is not None
        and boundary_corner is not None
        and expected is not None
    ):
        completed = smooth_k + interior_defect + intrinsic_seam + boundary_geodesic + boundary_corner
        residual = completed - expected
        if abs(residual) > max(1e-7, 100.0 * tolerance.area):
            gaussian_diagnostics.append(f"gauss_bonnet_residual:{residual:.12g}")
    elif not gaussian_diagnostics:
        gaussian_diagnostics.append("gauss_bonnet_missing_expected_topology")

    gaussian_diagnostics.extend(component_diagnostics)

    gaussian_status = "CERTIFIED" if not gaussian_diagnostics else "PARTIAL"
    mean_status = "CERTIFIED" if not mean_diagnostics else "PARTIAL"
    total_h = smooth_h + seam_h if mean_status == "CERTIFIED" else None
    status = "CERTIFIED" if mean_status == gaussian_status == "CERTIFIED" else "PARTIAL"
    diagnostics = tuple(sorted(set(mean_diagnostics + gaussian_diagnostics)))
    return MolecularContactSurfaceCurvature(
        smooth_h,
        seam_h if mean_status == "CERTIFIED" else None,
        total_h,
        smooth_k,
        interior_defect,
        intrinsic_seam,
        interior_defect,
        boundary_geodesic,
        boundary_corner,
        completed,
        expected,
        residual,
        mean_status,
        gaussian_status,
        status,
        "outward molecular normal; canonical spherical patches; physical refined seams; topology-only cuts excluded",
        diagnostics,
        physical_seam_count,
        exterior_boundary_count,
        audited_junctions,
    )


def _cell_complex_invariants(
    junctions: Sequence[JunctionVertex],
    edges: Sequence[tuple[str, str, str, tuple[str, ...]]],
    patches: Sequence[SphericalPatch],
    cycles_by_patch: Mapping[str, Sequence[tuple[tuple[str, ...], tuple[str, ...]]]],
) -> tuple[tuple[int, int, int], int, tuple[str, ...], str]:
    """Return GF(2) Betti numbers and the cellular chain-closure audit."""
    # Topological cuts occur twice with opposite signs in one face attaching
    # word, so their net column in ∂2 is zero; the incidence helper below
    # still records both occurrences explicitly.
    edge_by_id = {edge_id: (start, end) for edge_id, start, end, _ in edges}
    edge_index = {edge_id: index for index, (edge_id, *_rest) in enumerate(edges)}
    vertex_index = {junction.vertex_id: index for index, junction in enumerate(junctions)}
    face_rows: list[list[int]] = []
    chain_defects: list[str] = []
    for patch in patches:
        coefficients: dict[str, int] = defaultdict(int)
        cycles = cycles_by_patch.get(patch.patch_id, ())
        for edge_ids, vertices in cycles:
            if not vertices or len(edge_ids) != len(vertices):
                chain_defects.append(f"face_boundary:{patch.patch_id}")
                continue
            for offset, edge_id in enumerate(edge_ids):
                start = vertices[offset]
                end = vertices[(offset + 1) % len(vertices)]
                if edge_id not in edge_by_id or start not in vertex_index or end not in vertex_index:
                    chain_defects.append(f"face_edge:{patch.patch_id}:{edge_id}")
                    continue
                edge_start, edge_end = edge_by_id[edge_id]
                coefficients[edge_id] += 1 if (start, end) == (edge_start, edge_end) else -1
        row = [edge_index[edge_id] for edge_id, coefficient in coefficients.items() if coefficient % 2]
        face_rows.append(row)
        boundary = defaultdict(int)
        for edge_id, coefficient in coefficients.items():
            if coefficient == 0:
                continue
            start, end = edge_by_id[edge_id]
            boundary[start] -= coefficient
            boundary[end] += coefficient
        if any(value != 0 for value in boundary.values()):
            chain_defects.append(f"boundary_of_boundary:{patch.patch_id}")
    vertex_rows = [
        [
            index
            for index, (_, start, end, _patch_ids) in enumerate(edges)
            # A self-loop has zero cellular boundary, including over GF(2).
            if start != end and junction.vertex_id in {start, end}
        ]
        for junction in junctions
    ]
    rank_1 = _rank_mod2(vertex_rows, len(edges))
    rank_2 = _rank_mod2(face_rows, len(edges))
    beta_0 = len(junctions) - rank_1
    beta_1 = len(edges) - rank_1 - rank_2
    beta_2 = len(patches) - rank_2
    incidence = "CERTIFIED" if not chain_defects else "PARTIAL"
    return (beta_0, beta_1, beta_2), len(set(chain_defects)), tuple(sorted(set(chain_defects))), incidence


def _add_node(nodes: list[_Node], direction: Vec3, circle_ids: Iterable[str], tolerance: GeometryTolerance) -> int:
    direction = _unit(direction, tolerance.angular)
    if direction is None:
        raise ValueError("Spherical arrangement produced a zero direction")
    for index, node in enumerate(nodes):
        if _same_direction(node.direction, direction, tolerance):
            node.source_circle_ids.update(circle_ids)
            return index
    nodes.append(_Node(direction, set(circle_ids), False))
    return len(nodes) - 1


def _constraint_value(constraint: _Constraint, direction: Vec3) -> float:
    return _dot(constraint.normal, direction) - constraint.threshold


def _retained(
    direction: Vec3,
    occluders: Sequence[_Constraint],
    partners: Sequence[_Constraint],
    partner_full_cover: bool,
    tolerance: GeometryTolerance,
) -> bool:
    if any(_constraint_value(constraint, direction) > tolerance.angular for constraint in occluders):
        return False
    return partner_full_cover or any(
        _constraint_value(constraint, direction) >= -tolerance.angular for constraint in partners
    )


def _retained_spherical_path(
    first: Vec3,
    second: Vec3,
    occluders: Sequence[_Constraint],
    partner_constraints: Sequence[_Constraint],
    partner_full_cover: bool,
    tolerance: GeometryTolerance,
) -> bool:
    """Check a minor great-circle witness between two retained samples."""
    cosine = _clamp(_dot(first, second))
    angle = math.acos(cosine)
    if angle <= tolerance.angular:
        return True
    sine = math.sin(angle)
    if abs(sine) <= tolerance.angular:
        return False
    steps = max(16, int(math.ceil(angle / 0.05)))
    for index in range(1, steps):
        fraction = index / steps
        direction = _unit(
            _add(
                _scale(first, math.sin((1.0 - fraction) * angle)),
                _scale(second, math.sin(fraction * angle)),
            ),
            tolerance.angular,
        )
        if direction is None or not _retained(
            direction,
            occluders,
            partner_constraints,
            partner_full_cover,
            tolerance,
        ):
            return False
    return True


def _face_gauge(points: Sequence[Vec3]) -> Vec3:
    candidates = (
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, -1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, -1.0),
    )
    return max(candidates, key=lambda candidate: min(1.0 + _dot(candidate, point) for point in points))


_GL8 = (
    (-0.9602898564975363, 0.1012285362903763),
    (-0.7966664774136267, 0.2223810344533745),
    (-0.5255324099163290, 0.3137066458778873),
    (-0.1834346424956498, 0.3626837833783620),
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
)

_GL16 = (
    (-0.9894009349916499, 0.0271524594117541),
    (-0.9445750230732326, 0.0622535239386479),
    (-0.8656312023878318, 0.0951585116824928),
    (-0.7554044083550030, 0.1246289712555339),
    (-0.6178762444026438, 0.1495959888165767),
    (-0.4580167776572274, 0.1691565193950025),
    (-0.2816035507792589, 0.1826034150449236),
    (-0.0950125098376374, 0.1894506104550685),
    (0.0950125098376374, 0.1894506104550685),
    (0.2816035507792589, 0.1826034150449236),
    (0.4580167776572274, 0.1691565193950025),
    (0.6178762444026438, 0.1495959888165767),
    (0.7554044083550030, 0.1246289712555339),
    (0.8656312023878318, 0.0951585116824928),
    (0.9445750230732326, 0.0622535239386479),
    (0.9894009349916499, 0.0271524594117541),
)


def _arc_integral(circle: _Circle, theta_start: float, delta: float, gauge: Vec3, order: int) -> float:
    """Integrate the spherical-area one-form with adaptive Simpson quadrature."""

    def integrand(fraction: float) -> float:
        theta = theta_start + delta * fraction
        direction = _circle_point(circle, theta)
        derivative = _scale(_circle_tangent(circle, theta), delta)
        denominator = 1.0 + _dot(gauge, direction)
        if denominator <= 1e-11:
            return float("nan")
        return _dot(_cross(gauge, direction), derivative) / denominator

    tolerance = 1e-9 if order == 8 else 1e-12
    a, b = 0.0, 1.0
    c = 0.5
    fa, fb, fc = integrand(a), integrand(b), integrand(c)
    if not all(math.isfinite(value) for value in (fa, fb, fc)):
        return float("nan")
    initial = (b - a) * (fa + 4.0 * fc + fb) / 6.0

    def recurse(left, right, f_left, f_right, f_mid, estimate, remaining, target):
        middle = 0.5 * (left + right)
        left_mid = 0.5 * (left + middle)
        right_mid = 0.5 * (middle + right)
        f_left_mid = integrand(left_mid)
        f_right_mid = integrand(right_mid)
        if not all(math.isfinite(value) for value in (f_left_mid, f_right_mid)):
            return float("nan")
        left_estimate = (middle - left) * (f_left + 4.0 * f_left_mid + f_mid) / 6.0
        right_estimate = (right - middle) * (f_mid + 4.0 * f_right_mid + f_right) / 6.0
        refined = left_estimate + right_estimate
        if remaining <= 0 or abs(refined - estimate) <= 15.0 * target:
            return refined + (refined - estimate) / 15.0
        return recurse(left, middle, f_left, f_mid, f_left_mid, left_estimate, remaining - 1, target / 2.0) + recurse(
            middle, right, f_mid, f_right, f_right_mid, right_estimate, remaining - 1, target / 2.0
        )

    return recurse(a, b, fa, fb, fc, initial, 18, tolerance)


def _face_area(
    face: Sequence[int],
    halfedges: Sequence[_HalfEdge],
    segments: Sequence[_Segment],
    circles: Sequence[_Circle],
    tolerance: GeometryTolerance,
) -> tuple[float, float]:
    if len(face) == 1:
        halfedge = halfedges[face[0]]
        circle = circles[segments[halfedge.segment_id].circle_index]
        if abs(abs(halfedge.delta) - math.tau) <= tolerance.angular:
            # A single complete small circle is analytically a spherical cap
            # (or its complement); avoid a gauge pole in the boundary integral.
            sign = 1.0 if halfedge.delta >= 0.0 else -1.0
            area = 2.0 * math.pi * (1.0 - sign * circle.offset)
            return area, 0.0
    points = []
    for halfedge_id in face:
        halfedge = halfedges[halfedge_id]
        circle = circles[segments[halfedge.segment_id].circle_index]
        points.append(_circle_point(circle, halfedge.theta_start + 0.5 * halfedge.delta))
    gauge = _face_gauge(points)
    integral_8 = 0.0
    integral_16 = 0.0
    for halfedge_id in face:
        halfedge = halfedges[halfedge_id]
        circle = circles[segments[halfedge.segment_id].circle_index]
        integral_8 += _arc_integral(circle, halfedge.theta_start, halfedge.delta, gauge, 8)
        integral_16 += _arc_integral(circle, halfedge.theta_start, halfedge.delta, gauge, 16)
    area_8 = integral_8 % (4.0 * math.pi)
    area_16 = integral_16 % (4.0 * math.pi)
    if abs(area_16 - 4.0 * math.pi) <= tolerance.area:
        area_16 = 0.0
    return area_16, abs(area_16 - area_8)


def _pair_cap_area(carrier: AtomBall, partner: AtomBall, tolerance: GeometryTolerance) -> float:
    displacement = _sub(partner.center, carrier.center)
    distance = _norm(displacement)
    relation = _sphere_relation(carrier.radius, partner.radius, distance, tolerance)
    if relation == "other_contains" or relation == "coincident":
        return 4.0 * math.pi * carrier.radius * carrier.radius
    if relation in {"disjoint", "carrier_contains", "external_tangent", "internal_tangent"}:
        if relation == "internal_tangent" and partner.radius > carrier.radius:
            return 4.0 * math.pi * carrier.radius * carrier.radius
        return 0.0
    height = (partner.radius * partner.radius - (distance - carrier.radius) ** 2) / (2.0 * distance)
    return 2.0 * math.pi * carrier.radius * max(0.0, min(2.0 * carrier.radius, height))


def _build_carrier(
    carrier: AtomBall,
    same_group: Sequence[AtomBall],
    partners: Sequence[tuple[AtomBall, SelectedContact]],
    side: str,
    tolerance: GeometryTolerance,
    source_network_id: str,
    source_interface_id: str,
    definition_version: str,
    contact_selection_id: str,
) -> tuple[list[SphericalPatch], list[CircularArc], list[JunctionVertex], float, float, list[str], int, int, Optional[int], list[str]]:
    diagnostics: list[str] = []
    occluders: list[_Constraint] = []
    partner_constraints: list[_Constraint] = []
    active_occluder_ids: set[str] = set()
    active_partner_ids: set[str] = set()
    partner_full_cover = False
    fully_occluded = False
    naive_area = 0.0

    for other in same_group:
        if other.generator_index == carrier.generator_index:
            continue
        displacement = _sub(other.center, carrier.center)
        distance = _norm(displacement)
        relation = _sphere_relation(carrier.radius, other.radius, distance, tolerance)
        if relation == "other_contains":
            fully_occluded = True
            active_occluder_ids.add(other.stable_id)
            diagnostics.append(f"carrier_fully_occluded:{carrier.stable_id}:{other.stable_id}")
            continue
        if relation == "coincident":
            # Equal coincident same-group balls define a tied carrier
            # boundary.  Keep the geometry, but preserve the unresolved
            # provenance rather than silently treating the tie as generic
            # exposure.
            active_occluder_ids.add(other.stable_id)
            diagnostics.append(f"same_group_coincident:{carrier.stable_id}:{other.stable_id}")
            continue
        if relation in {"partial_overlap", "internal_tangent"}:
            if relation == "internal_tangent":
                diagnostics.append(f"same_group_internal_tangent:{carrier.stable_id}:{other.stable_id}")
                continue
            normal = _unit(displacement, tolerance.angular)
            if normal is None:
                diagnostics.append(f"same_group_coincident:{carrier.stable_id}:{other.stable_id}")
                continue
            threshold = (
                carrier.radius * carrier.radius
                + distance * distance
                - other.radius * other.radius
            ) / (2.0 * carrier.radius * distance)
            occluders.append(_Constraint("occluder", other.generator_index, other.stable_id, normal, threshold))
            active_occluder_ids.add(other.stable_id)

    for partner, contact in partners:
        naive_area += _pair_cap_area(carrier, partner, tolerance)
        displacement = _sub(partner.center, carrier.center)
        distance = _norm(displacement)
        relation = _sphere_relation(carrier.radius, partner.radius, distance, tolerance)
        active_partner_ids.add(partner.stable_id)
        if relation == "other_contains":
            partner_full_cover = True
            continue
        if relation == "coincident":
            partner_full_cover = True
            diagnostics.append(f"partner_coincident:{carrier.stable_id}:{partner.stable_id}")
            continue
        if relation == "partial_overlap":
            normal = _unit(displacement, tolerance.angular)
            if normal is None:
                diagnostics.append(f"partner_coincident:{carrier.stable_id}:{partner.stable_id}")
                continue
            threshold = (
                carrier.radius * carrier.radius
                + distance * distance
                - partner.radius * partner.radius
            ) / (2.0 * carrier.radius * distance)
            partner_constraints.append(
                _Constraint("partner", partner.generator_index, partner.stable_id, normal, threshold, (contact.contact_id,))
            )
            continue
        if relation == "external_tangent" or relation == "internal_tangent":
            diagnostics.append(f"unresolved_contact:{contact.contact_id}:partner_zero_area_contact")
        elif relation == "disjoint" or relation == "carrier_contains":
            diagnostics.append(f"unresolved_contact:{contact.contact_id}:partner_nonoverlap_contact")

    selected_contact_ids = tuple(sorted(contact.contact_id for _, contact in partners))
    if fully_occluded or (not partner_full_cover and not partner_constraints):
        return [], [], [], 0.0, naive_area, list(active_partner_ids | active_occluder_ids), 0, 0, 0, diagnostics

    circles: list[_Circle] = []
    all_constraints = occluders + partner_constraints
    for constraint in all_constraints:
        normal = constraint.normal
        offset = _clamp(constraint.threshold)
        existing = next(
            (
                circle for circle in circles
                if (
                    abs(circle.offset - offset) <= tolerance.angular
                    and _norm(_sub(circle.normal, normal)) <= tolerance.angular
                )
                or (
                    abs(circle.offset + offset) <= tolerance.angular
                    and _norm(_add(circle.normal, normal)) <= tolerance.angular
                )
            ),
            None,
        )
        if existing is not None:
            existing.constraints.append(constraint)
            constraint.circle_id = existing.circle_id
            continue
        circle = _Circle(
            _token("circle", carrier.stable_id, len(circles), constraint.source_atom_id, constraint.kind),
            normal,
            offset,
            [constraint],
        )
        circle.basis_1, circle.basis_2 = _circle_basis(normal, tolerance)
        circle.radial = math.sqrt(max(0.0, 1.0 - offset * offset))
        constraint.circle_id = circle.circle_id
        circles.append(circle)

    if not circles:
        direction = (0.0, 0.0, 1.0)
        if not _retained(direction, occluders, partner_constraints, partner_full_cover, tolerance):
            return [], [], [], 0.0, naive_area, list(active_partner_ids | active_occluder_ids), 0, 0, 0, diagnostics
        patch_id = _token("patch", carrier.stable_id, "full")
        patch = SphericalPatch(
            patch_id, carrier.stable_id, carrier.center, carrier.radius, side,
            tuple(sorted(active_partner_ids)), tuple(sorted(active_occluder_ids)), selected_contact_ids,
            source_network_id, source_interface_id, definition_version, (), (), patch_id,
            4.0 * math.pi * carrier.radius * carrier.radius, 0.0, "full_sphere",
            source_contact_selection_id=contact_selection_id,
        )
        return [patch], [], [], patch.area_A2, naive_area, list(active_partner_ids | active_occluder_ids), 1, 0, 2, diagnostics

    nodes: list[_Node] = []
    circle_pair_intersections: dict[int, list[int]] = defaultdict(list)
    for first_index, first in enumerate(circles):
        for second_index in range(first_index):
            second = circles[second_index]
            for point in _circle_intersections(first, second, tolerance):
                node_index = _add_node(nodes, point, (first.circle_id, second.circle_id), tolerance)
                circle_pair_intersections[first_index].append(node_index)
                circle_pair_intersections[second_index].append(node_index)
    for index, circle in enumerate(circles):
        circle.nodes = sorted(set(circle_pair_intersections.get(index, [])))
        if not circle.nodes:
            synthetic = _circle_point(circle, 0.0)
            nodes.append(_Node(synthetic, {circle.circle_id}, True))
            circle.nodes = [len(nodes) - 1]

    segments: list[_Segment] = []
    for circle_index, circle in enumerate(circles):
        angles = []
        for node_index in circle.nodes:
            direction = nodes[node_index].direction
            angles.append(math.atan2(_dot(direction, circle.basis_2), _dot(direction, circle.basis_1)) % (2.0 * math.pi))
        angles_and_nodes = sorted(zip(angles, circle.nodes), key=lambda item: (item[0], item[1]))
        for offset, (theta, node_index) in enumerate(angles_and_nodes):
            next_theta, next_node = angles_and_nodes[(offset + 1) % len(angles_and_nodes)]
            delta = next_theta - theta
            if offset == len(angles_and_nodes) - 1:
                delta += 2.0 * math.pi
            if delta <= tolerance.angular:
                continue
            segments.append(_Segment(len(segments), circle_index, node_index, next_node, theta, delta))

    halfedges: list[_HalfEdge] = []
    outgoing: dict[int, list[int]] = defaultdict(list)
    for segment in segments:
        circle = circles[segment.circle_index]
        start_tangent = _circle_tangent(circle, segment.theta_start)
        end_theta = segment.theta_start + segment.delta
        end_tangent = _circle_tangent(circle, end_theta)
        forward = _HalfEdge(
            len(halfedges), segment.segment_id, False, segment.start_node, segment.end_node,
            segment.theta_start, segment.delta, start_tangent, end_tangent, len(halfedges) + 1,
        )
        reverse = _HalfEdge(
            len(halfedges) + 1, segment.segment_id, True, segment.end_node, segment.start_node,
            end_theta, -segment.delta, _scale(end_tangent, -1.0), _scale(start_tangent, -1.0), len(halfedges),
        )
        halfedges.extend((forward, reverse))
        outgoing[forward.start_node].append(forward.halfedge_id)
        outgoing[reverse.start_node].append(reverse.halfedge_id)

    def tangent_angle(node_index: int, tangent: Vec3) -> float:
        normal = nodes[node_index].direction
        first, second = _circle_basis(normal, tolerance)
        return math.atan2(_dot(tangent, second), _dot(tangent, first)) % (2.0 * math.pi)

    outgoing_angles = {
        node: sorted((tangent_angle(node, halfedges[edge].tangent_start), edge) for edge in edges)
        for node, edges in outgoing.items()
    }
    for halfedge in halfedges:
        target = tangent_angle(halfedge.end_node, _scale(halfedge.tangent_end, -1.0))
        choices = outgoing_angles[halfedge.end_node]
        angles = [angle for angle, _ in choices]
        position = bisect.bisect_left(angles, target) - 1
        halfedge.next_halfedge = choices[position % len(choices)][1]

    faces: list[list[int]] = []
    for halfedge in halfedges:
        if halfedge.face_id is not None:
            continue
        cycle = []
        current = halfedge.halfedge_id
        while halfedges[current].face_id is None:
            halfedges[current].face_id = -1
            cycle.append(current)
            current = halfedges[current].next_halfedge
            if current == halfedge.halfedge_id:
                break
            if current is None or current in cycle:
                break
        face_id = len(faces)
        for edge_id in cycle:
            halfedges[edge_id].face_id = face_id
        faces.append(cycle)

    retained_faces: set[int] = set()
    face_signatures: dict[int, tuple[bool, ...]] = {}
    face_samples: dict[int, Vec3] = {}
    for face_id, face in enumerate(faces):
        first_edge = halfedges[face[0]]
        circle = circles[segments[first_edge.segment_id].circle_index]
        midpoint = _circle_point(circle, first_edge.theta_start + 0.5 * first_edge.delta)
        tangent = _scale(_circle_tangent(circle, first_edge.theta_start + 0.5 * first_edge.delta), 1.0 if first_edge.delta >= 0 else -1.0)
        left = _unit(_cross(midpoint, tangent), tolerance.angular)
        sample = _unit(_add(midpoint, _scale(left or (0.0, 0.0, 0.0), 1e-7)), tolerance.angular)
        if sample is not None and _retained(sample, occluders, partner_constraints, partner_full_cover, tolerance):
            retained_faces.add(face_id)
        if sample is not None:
            face_samples[face_id] = sample
            face_signatures[face_id] = tuple(
                _constraint_value(constraint, sample) > tolerance.angular
                for constraint in all_constraints
            )

    parent = {face_id: face_id for face_id in retained_faces}

    def find(face_id: int) -> int:
        while parent[face_id] != face_id:
            parent[face_id] = parent[parent[face_id]]
            face_id = parent[face_id]
        return face_id

    def union(first: int, second: int):
        first, second = find(first), find(second)
        if first != second:
            parent[second] = first

    # Disjoint arrangement loops can describe one connected exterior region.
    # Connect retained face witnesses through a small deterministic set of
    # retained spherical witnesses before merging equal-signature cycles.
    # ponytail: the witness set is sufficient, not exhaustive; exotic narrow
    # corridors may remain conservatively over-split until exact cell routing.
    witness_faces = sorted(retained_faces)
    witness_points = [face_samples[face_id] for face_id in witness_faces]
    witness_points.extend(
        direction
        for direction in (
            (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0), (0.0, -1.0, 0.0),
            (0.0, 0.0, 1.0), (0.0, 0.0, -1.0),
        )
        if _retained(direction, occluders, partner_constraints, partner_full_cover, tolerance)
    )
    witness_parent = list(range(len(witness_points)))

    def find_witness(index: int) -> int:
        while witness_parent[index] != index:
            witness_parent[index] = witness_parent[witness_parent[index]]
            index = witness_parent[index]
        return index

    def union_witnesses(first: int, second: int):
        first, second = find_witness(first), find_witness(second)
        if first != second:
            witness_parent[second] = first

    for first_index in range(len(witness_points)):
        for second_index in range(first_index):
            if _retained_spherical_path(
                witness_points[first_index], witness_points[second_index],
                occluders, partner_constraints, partner_full_cover, tolerance,
            ):
                union_witnesses(first_index, second_index)

    by_signature: dict[tuple[bool, ...], list[int]] = defaultdict(list)
    for face_index, face_id in enumerate(witness_faces):
        signature = face_signatures.get(face_id)
        if signature is None:
            continue
        connected_face = next(
            (
                other_id
                for other_index, other_id in enumerate(witness_faces[:face_index])
                if face_signatures.get(other_id) == signature
                and find_witness(face_index) == find_witness(other_index)
            ),
            None,
        )
        if connected_face is not None:
            union(connected_face, face_id)
        by_signature[signature].append(face_id)

    boundary_halfedges: list[int] = []
    for segment in segments:
        first = 2 * segment.segment_id
        second = first + 1
        first_face = halfedges[first].face_id
        second_face = halfedges[second].face_id
        first_retained = first_face in retained_faces
        second_retained = second_face in retained_faces
        if first_retained and second_retained:
            union(first_face, second_face)
        elif first_retained != second_retained:
            boundary_halfedges.append(first if first_retained else second)

    component_faces: dict[int, list[int]] = defaultdict(list)
    for face_id in retained_faces:
        component_faces[find(face_id)].append(face_id)
    component_boundary: dict[int, list[int]] = defaultdict(list)
    for halfedge_id in boundary_halfedges:
        component_boundary[find(halfedges[halfedge_id].face_id)].append(halfedge_id)

    public_node_ids: dict[int, str] = {}
    for node_index, node in enumerate(nodes):
        if node.synthetic:
            continue
        public_node_ids[node_index] = _token("junction", carrier.stable_id, node_index, tuple(sorted(node.source_circle_ids)))

    arcs: list[CircularArc] = []
    arc_ids_by_component: dict[int, list[str]] = defaultdict(list)
    junction_arc_ids: dict[str, list[str]] = defaultdict(list)
    for component_id, halfedge_ids in component_boundary.items():
        for halfedge_id in halfedge_ids:
            halfedge = halfedges[halfedge_id]
            segment = segments[halfedge.segment_id]
            circle = circles[segment.circle_index]
            start_node = halfedge.start_node
            end_node = halfedge.end_node
            start_direction = _circle_point(circle, halfedge.theta_start)
            end_direction = _circle_point(circle, halfedge.theta_start + halfedge.delta)
            start_id = public_node_ids.get(start_node)
            end_id = public_node_ids.get(end_node)
            arc_id = _token("arc", carrier.stable_id, component_id, halfedge_id)
            arc = CircularArc(
                arc_id, carrier.stable_id, circle.circle_id, start_id, end_id,
                _add(carrier.center, _scale(start_direction, carrier.radius)),
                _add(carrier.center, _scale(end_direction, carrier.radius)),
                circle.normal, circle.offset, halfedge.delta,
                circle.source_atom_ids, circle.source_contact_ids,
                1, "CERTIFIED",
            )
            arcs.append(arc)
            arc_ids_by_component[component_id].append(arc_id)
            if start_id:
                junction_arc_ids[start_id].append(arc_id)
            if end_id:
                junction_arc_ids[end_id].append(arc_id)

    junctions = []
    for node_index, vertex_id in public_node_ids.items():
        if vertex_id not in junction_arc_ids:
            continue
        node = nodes[node_index]
        junctions.append(JunctionVertex(
            vertex_id,
            carrier.stable_id,
            _add(carrier.center, _scale(node.direction, carrier.radius)),
            node.direction,
            tuple(sorted(node.source_circle_ids)),
            tuple(sorted(junction_arc_ids[vertex_id])),
            "CERTIFIED",
        ))

    patches: list[SphericalPatch] = []
    total_area = 0.0
    total_error = 0.0
    boundary_loops = 0
    nonmanifold = 0
    for root_face, face_ids in sorted(component_faces.items()):
        component_halfedges = tuple(component_boundary.get(root_face, ()))
        if component_halfedges:
            unit_area, unit_error = _face_area(
                component_halfedges, halfedges, segments, circles, tolerance
            )
        else:
            unit_area, unit_error = 4.0 * math.pi, 0.0
        area = unit_area * carrier.radius * carrier.radius
        error = unit_error * carrier.radius * carrier.radius
        component_arcs = tuple(sorted(arc_ids_by_component.get(root_face, [])))
        component_junctions = tuple(sorted({
            vertex_id
            for arc in arcs if arc.arc_id in component_arcs
            for vertex_id in (arc.start_vertex_id, arc.end_vertex_id)
            if vertex_id is not None
        }))
        outgoing_boundary = defaultdict(int)
        incoming_boundary = defaultdict(int)
        for arc in arcs:
            if arc.arc_id not in component_arcs:
                continue
            if arc.start_vertex_id is not None:
                outgoing_boundary[arc.start_vertex_id] += 1
            if arc.end_vertex_id is not None:
                incoming_boundary[arc.end_vertex_id] += 1
        component_nonmanifold = any(
            outgoing_boundary[key] != 1 or incoming_boundary[key] != 1
            for key in set(outgoing_boundary) | set(incoming_boundary)
        )
        if component_nonmanifold:
            nonmanifold += 1
            diagnostics.append(f"nonmanifold_boundary:{carrier.stable_id}:{root_face}")
        else:
            outgoing_arcs = defaultdict(list)
            for candidate in arcs:
                if candidate.arc_id in component_arcs and candidate.start_vertex_id is not None:
                    outgoing_arcs[candidate.start_vertex_id].append(candidate)
            for candidates in outgoing_arcs.values():
                candidates.sort(key=lambda value: value.arc_id)
            visited = set()
            for arc in arcs:
                if arc.arc_id not in component_arcs or arc.arc_id in visited:
                    continue
                if arc.start_vertex_id is None and arc.end_vertex_id is None:
                    boundary_loops += 1
                    visited.add(arc.arc_id)
                    continue
                current = arc
                while current.arc_id not in visited:
                    visited.add(current.arc_id)
                    next_candidates = outgoing_arcs.get(current.end_vertex_id, [])
                    if not next_candidates:
                        break
                    current = sorted(next_candidates, key=lambda value: value.arc_id)[0]
                boundary_loops += 1
        patch_id = _token("patch", carrier.stable_id, root_face)
        patches.append(SphericalPatch(
            patch_id, carrier.stable_id, carrier.center, carrier.radius, side,
            tuple(sorted(active_partner_ids)), tuple(sorted(active_occluder_ids)), selected_contact_ids,
            source_network_id, source_interface_id, definition_version,
            component_arcs, component_junctions, patch_id, area, error,
            "arrangement_component",
            "CERTIFIED" if error <= tolerance.quadrature and not component_nonmanifold else "PARTIAL",
            contact_selection_id,
        ))
        total_area += area
        total_error += error

    components = len(patches)
    euler = None if nonmanifold else 2 * components - boundary_loops
    public_junctions = junctions
    return (
        patches,
        arcs,
        public_junctions,
        total_area,
        naive_area,
        list(active_partner_ids | active_occluder_ids),
        components,
        boundary_loops,
        euler,
        diagnostics,
    )


def build_molecular_contact_surface(
    atoms: Iterable[AtomBall],
    selected_contacts: Iterable[SelectedContact],
    side: str,
    *,
    source_network_id: str = "",
    source_interface_id: str = "",
    contact_selection_source: Optional[str] = None,
    contact_selection_id: str = "",
    definition_version: str = "selected_partner_union_of_balls_v1",
    tolerance: GeometryTolerance = GeometryTolerance(),
) -> MolecularContactSurface:
    """Build one side of the exact curved selected-contact boundary."""
    if side not in {"A", "B"}:
        raise ValueError("Molecular contact surface side must be A or B")
    atoms = tuple(atoms)
    atom_map = {atom.generator_index: atom for atom in atoms}
    contacts = tuple(selected_contacts)
    if not atom_map:
        raise ValueError("At least one atom ball is required")
    group_atoms = tuple(sorted((atom for atom in atom_map.values() if atom.group == side), key=lambda atom: atom.generator_index))
    carriers: dict[int, list[tuple[AtomBall, SelectedContact]]] = defaultdict(list)
    selected_ids = []
    naive_total = 0.0
    for contact in contacts:
        if contact.status != "SELECTED":
            continue
        a = atom_map.get(contact.atom_a_index)
        b = atom_map.get(contact.atom_b_index)
        if a is None or b is None or a.group != "A" or b.group != "B":
            raise ValueError(f"Selected contact {contact.contact_id} has invalid A/B generators")
        carrier, partner = (a, b) if side == "A" else (b, a)
        carriers[carrier.generator_index].append((partner, contact))
        selected_ids.append(contact.contact_id)
    selected_ids = tuple(sorted(set(selected_ids)))
    patches = []
    arcs = []
    junctions = []
    diagnostics = []
    total_area = 0.0
    total_components = 0
    total_loops = 0
    euler_values = []
    euler_certified = True
    for carrier_index in sorted(carriers):
        carrier = atom_map[carrier_index]
        result = _build_carrier(
            carrier,
            group_atoms,
            sorted(carriers[carrier_index], key=lambda item: (item[0].generator_index, item[1].contact_id)),
            side,
            tolerance,
            source_network_id,
            source_interface_id,
            definition_version,
            contact_selection_id,
        )
        carrier_patches, carrier_arcs, carrier_junctions, area, naive, _, components, loops, euler, carrier_diagnostics = result
        patches.extend(carrier_patches)
        arcs.extend(carrier_arcs)
        junctions.extend(carrier_junctions)
        total_area += area
        naive_total += naive
        total_components += components
        total_loops += loops
        if euler is not None:
            euler_values.append(euler)
        else:
            euler_certified = False
        diagnostics.extend(carrier_diagnostics)

    # Merge physically coincident junction coordinates from neighboring carrier
    # spheres.  Their local unit directions differ, so the first deterministic
    # direction is retained while all carrier provenance is preserved.
    junction_groups: list[list[JunctionVertex]] = []
    for junction in sorted(junctions, key=lambda value: value.vertex_id):
        group = next(
            (
                group for group in junction_groups
                if _norm(_sub(group[0].xyz, junction.xyz)) <= 1e-7
            ),
            None,
        )
        if group is None:
            junction_groups.append([junction])
        else:
            group.append(junction)
    junction_id_map = {}
    merged_junctions = []
    for group in junction_groups:
        first = group[0]
        merged_id = _token("junction-global", tuple(round(value, 8) for value in first.xyz))
        carriers_at_junction = tuple(sorted({item.carrier_atom_id for item in group}))
        source_circles = tuple(sorted({circle for item in group for circle in item.source_circle_ids}))
        incident_arcs = tuple(sorted({arc for item in group for arc in item.incident_arc_ids}))
        for item in group:
            junction_id_map[item.vertex_id] = merged_id
        merged_junctions.append(replace(
            first,
            vertex_id=merged_id,
            source_circle_ids=source_circles,
            incident_arc_ids=incident_arcs,
            source_carrier_atom_ids=carriers_at_junction,
        ))
    junctions = merged_junctions
    arcs = [
        replace(
            arc,
            start_vertex_id=junction_id_map.get(arc.start_vertex_id, arc.start_vertex_id),
            end_vertex_id=junction_id_map.get(arc.end_vertex_id, arc.end_vertex_id),
        )
        for arc in arcs
    ]
    patches = [
        replace(
            patch,
            junction_ids=tuple(sorted(junction_id_map.get(value, value) for value in patch.junction_ids)),
        )
        for patch in patches
    ]
    patch_index = {patch.patch_id: index for index, patch in enumerate(patches)}
    patch_by_id = {patch.patch_id: patch for patch in patches}
    stable_atom_map = {atom.stable_id: atom for atom in atom_map.values()}
    parent = list(range(len(patches)))

    def find_patch(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union_patches(first: int, second: int):
        first, second = find_patch(first), find_patch(second)
        if first != second:
            parent[second] = first

    arc_patch = {
        arc.arc_id: patch.patch_id
        for patch in patches
        for arc_id in patch.boundary_arc_ids
        for arc in arcs
        if arc.arc_id == arc_id
    }
    arc_by_id = {arc.arc_id: arc for arc in arcs}
    self_seam_arcs = [
        arc for arc in arcs
        if not arc.source_contact_ids and len(arc.source_atom_ids) == 1
    ]
    seam_groups: dict[tuple[str, str], list[CircularArc]] = defaultdict(list)
    for arc in self_seam_arcs:
        seam_groups[tuple(sorted((arc.carrier_atom_id, arc.source_atom_ids[0])))].append(arc)
    refined_seams: list[RefinedSeamPiece] = []
    refinement_diagnostics: list[str] = []
    for pair in sorted(seam_groups):
        first = stable_atom_map.get(pair[0])
        second = stable_atom_map.get(pair[1])
        if first is None or second is None:
            refinement_diagnostics.append(f"unresolved_shared_circle:{pair[0]}:{pair[1]}")
            continue
        pieces, piece_diagnostics = _refine_shared_circle(
            first,
            second,
            seam_groups[pair],
            arc_patch,
            patch_by_id,
            junctions,
            tolerance,
            stable_atom_map,
        )
        refined_seams.extend(pieces)
        refinement_diagnostics.extend(piece_diagnostics)
        for piece in pieces:
            if piece.classification != "internal":
                continue
            piece_patches = [patch_index[patch_id] for patch_id in piece.incident_patch_ids]
            for left, right in zip(piece_patches, piece_patches[1:]):
                union_patches(left, right)

    refined_arc_ids_by_junction: dict[str, set[str]] = defaultdict(set)
    for piece in refined_seams:
        for vertex_id in piece.endpoint_junction_ids:
            refined_arc_ids_by_junction[vertex_id].update(piece.source_arc_ids)
    junctions = [
        replace(
            junction,
            incident_arc_ids=tuple(sorted(set(junction.incident_arc_ids) | refined_arc_ids_by_junction.get(junction.vertex_id, set()))),
        )
        for junction in junctions
    ]

    # Canonical global 1-cells are every common-refined seam interval plus the
    # retained partner arcs.  Raw carrier arcs remain in the public geometry
    # for provenance; only this list is used for global topology.
    global_edges: list[tuple[str, str, str, tuple[str, ...]]] = []
    internal_edge_ids: set[str] = set()
    for piece in refined_seams:
        patch_ids = piece.incident_patch_ids
        if not patch_ids:
            continue
        global_edges.append((piece.seam_id, piece.endpoint_junction_ids[0], piece.endpoint_junction_ids[1], patch_ids))
        if piece.classification == "internal":
            internal_edge_ids.add(piece.seam_id)
    for arc in arcs:
        patch_id = arc_patch.get(arc.arc_id)
        if not arc.source_contact_ids or patch_id is None:
            continue
        if arc.start_vertex_id is None or arc.end_vertex_id is None:
            patch = patch_by_id[patch_id]
            basis_1, _ = _circle_basis(arc.circle_normal, tolerance)
            radial = math.sqrt(max(0.0, 1.0 - arc.circle_offset * arc.circle_offset))
            direction = _add(
                _scale(arc.circle_normal, arc.circle_offset),
                _scale(basis_1, radial),
            )
            point = _add(patch.carrier_center, _scale(direction, patch.carrier_radius))
            vertex_id = _token("topological-boundary-vertex", arc.arc_id)
            junctions.append(JunctionVertex(
                vertex_id,
                patch.carrier_atom_id,
                point,
                direction,
                (arc.circle_id,),
                (arc.arc_id,),
                "TOPOLOGY_ONLY",
                (patch.carrier_atom_id,),
                (),
                "full_circle_boundary",
            ))
            arc = replace(
                arc,
                start_vertex_id=vertex_id,
                end_vertex_id=vertex_id,
                start_xyz=point,
                end_xyz=point,
            )
            arcs[arcs.index(arc_by_id[arc.arc_id])] = arc
            patch_by_id[patch_id] = replace(
                patch,
                junction_ids=tuple(sorted(set(patch.junction_ids) | {vertex_id})),
            )
            patches[patch_index[patch_id]] = patch_by_id[patch_id]
        global_edges.append((arc.arc_id, arc.start_vertex_id, arc.end_vertex_id, (patch_id,)))

    arc_by_id = {arc.arc_id: arc for arc in arcs}

    # A complete spherical carrier region has no physical boundary edges.
    # Give its CW model one deterministic auxiliary pole so S² has V=1, E=0,
    # F=1 and χ=2 without inventing any molecular geometry.
    if not global_edges and len(patches) == 1 and patches[0].region_kind == "full_sphere":
        patch = patches[0]
        pole_id = _token("topological-pole", patch.patch_id)
        pole = JunctionVertex(
            pole_id,
            patch.carrier_atom_id,
            _add(patch.carrier_center, _scale((0.0, 0.0, 1.0), patch.carrier_radius)),
            (0.0, 0.0, 1.0),
            (),
            (),
            "TOPOLOGY_ONLY",
            (patch.carrier_atom_id,),
            (),
            "full_sphere",
        )
        junctions.append(pole)
        patches[0] = replace(patch, junction_ids=(pole_id,))

    refined_by_id = {piece.seam_id: piece for piece in refined_seams}
    topological_cuts, face_boundary_cycle_counts, cut_diagnostics = _make_topological_cuts(
        patches,
        global_edges,
        refined_by_id,
        junctions,
    )
    cell_edges = list(global_edges) + [
        (
            cut.cut_id,
            cut.endpoint_junction_ids[0],
            cut.endpoint_junction_ids[1],
            # The incidence audit expands this as two opposite occurrences
            # of the same face; it is not a physical boundary edge.
            (cut.patch_id, cut.patch_id),
        )
        for cut in topological_cuts
    ]
    cycles_by_patch = {
        patch.patch_id: _ordered_patch_boundary_cycles(patch.patch_id, global_edges, refined_by_id)
        for patch in patches
    }
    betti_numbers, chain_defect_count, chain_defects, incidence_status = _cell_complex_invariants(
        junctions,
        cell_edges,
        patches,
        cycles_by_patch,
    )
    edge_incidence_distribution, edge_incidence_issues = _cell_edge_incidence_distribution(
        cell_edges,
        cycles_by_patch,
        topological_cuts,
    )
    if edge_incidence_issues:
        chain_defects = tuple(sorted(set(chain_defects) | set(edge_incidence_issues)))
        incidence_status = "PARTIAL"

    global_roots = {find_patch(index) for index in range(len(patches))}
    global_components = len(global_roots)
    boundary_edges = [edge for edge in global_edges if edge[0] not in internal_edge_ids]
    outgoing: dict[str, list[str]] = defaultdict(list)
    incoming: dict[str, list[str]] = defaultdict(list)
    edge_by_id = {edge_id: (start, end) for edge_id, start, end, _ in boundary_edges}
    edge_root: dict[str, int] = {}
    for edge_id, start, end, patch_ids in boundary_edges:
        for patch_id in patch_ids:
            if patch_id in patch_index:
                edge_root[edge_id] = find_patch(patch_index[patch_id])
                break
        outgoing[start].append(edge_id)
        incoming[end].append(edge_id)
    boundary_degree_defects = {
        vertex
        for vertex in set(outgoing) | set(incoming)
        if len(outgoing[vertex]) != len(incoming[vertex]) or len(outgoing[vertex]) not in {0, 1}
    }
    global_boundary_loops = 0
    boundary_loop_roots: dict[int, int] = defaultdict(int)
    if not boundary_degree_defects:
        visited_edges: set[str] = set()
        for edge_id in sorted(edge_by_id):
            if edge_id in visited_edges:
                continue
            current = edge_id
            closed = False
            while current not in visited_edges:
                visited_edges.add(current)
                _, end = edge_by_id[current]
                candidates = outgoing.get(end, ())
                if len(candidates) != 1:
                    break
                current = candidates[0]
            if current == edge_id:
                closed = True
            if closed:
                global_boundary_loops += 1
                boundary_loop_roots[edge_root.get(edge_id, -1)] += 1

    junctions, link_status, manifold_vertices, singular_vertices = _classify_local_links(
        junctions,
        global_edges,
        arc_by_id,
    )
    topology_diagnostics = list(refinement_diagnostics) + list(cut_diagnostics) + list(chain_defects)
    if boundary_degree_defects:
        topology_diagnostics.append(f"common_refinement_boundary_degree_defects:{len(boundary_degree_defects)}")
    global_nonmanifold = singular_vertices + len(boundary_degree_defects)
    global_euler = None
    orientability = None
    genus_by_component = None
    if (
        not topology_diagnostics
        and global_nonmanifold == 0
        and euler_certified
        and len(junctions) == len({junction.vertex_id for junction in junctions})
    ):
        global_euler = len(junctions) - len(cell_edges) + len(patches)
        # A parity walk over face adjacencies certifies orientability without
        # changing any face orientation in the scientific geometry.
        parity: dict[str, int] = {}
        orientable = True
        for piece in refined_seams:
            if piece.classification != "internal" or len(piece.orientation_on_patches) < 2:
                continue
            left, right = piece.orientation_on_patches[:2]
            required = 0 if left[1] != right[1] else 1
            if left[0] not in parity:
                parity[left[0]] = 0
            expected = parity[left[0]] ^ required
            if right[0] in parity and parity[right[0]] != expected:
                orientable = False
                break
            parity[right[0]] = expected
        orientability = orientable
        if orientable:
            component_genus = []
            for root in sorted(global_roots):
                component_faces = {patch.patch_id for index, patch in enumerate(patches) if find_patch(index) == root}
                component_edges = [edge for edge in cell_edges if any(patch_id in component_faces for patch_id in edge[3])]
                component_vertices = {
                    vertex
                    for patch in patches
                    if find_patch(patch_index[patch.patch_id]) == root
                    for vertex in patch.junction_ids
                }
                component_vertices.update(
                    vertex for edge in component_edges for vertex in edge[1:3]
                )
                component_chi = len(component_vertices) - len(component_edges) + len(component_faces)
                boundary_count = boundary_loop_roots.get(root, 0)
                numerator = 2 - boundary_count - component_chi
                if numerator < 0 or numerator % 2:
                    orientability = False
                    break
                component_genus.append(numerator // 2)
            if orientability:
                genus_by_component = tuple(component_genus)
    if orientability is False:
        topology_diagnostics.append("common_refinement_orientability_conflict")
    diagnostics.extend(topology_diagnostics)
    updated_patches = tuple(
        replace(patch, component_id=_token("component", side, find_patch(index)))
        for index, patch in enumerate(patches)
    )
    unresolved = tuple(sorted(
        diagnostic.split(":", 2)[1]
        for diagnostic in diagnostics
        if diagnostic.startswith("unresolved_contact:")
    ))
    local_nonmanifold_components = sum(
        diagnostic.startswith("nonmanifold_boundary:")
        for diagnostic in diagnostics
    )
    nonmanifold = global_nonmanifold + local_nonmanifold_components
    status = "CERTIFIED" if not diagnostics and not nonmanifold else "PARTIAL"
    euler = (
        global_euler
        if euler_certified and global_nonmanifold == 0 and orientability is True
        else None
    )
    surface = MolecularContactSurface(
        side=side,
        convention="expanded_union_of_balls",
        definition_version=definition_version,
        patches=updated_patches,
        arcs=tuple(arcs),
        junctions=tuple(junctions),
        carrier_atom_ids=tuple(sorted(atom_map[index].stable_id for index in carriers)),
        selected_contact_ids=selected_ids,
        total_area_A2=total_area,
        naive_pairwise_area_A2=naive_total,
        components=global_components,
        boundary_loops=global_boundary_loops,
        euler_characteristic=euler,
        nonmanifold_vertices=nonmanifold,
        unresolved_contact_ids=unresolved,
        area_method="spherical boundary integral, adaptive Simpson convergence check",
        area_error_A2=sum(patch.area_error_A2 for patch in patches),
        status=status,
        diagnostics=tuple(sorted(set(diagnostics))),
        contact_selection_source=contact_selection_source,
        contact_selection_id=contact_selection_id,
        refined_seams=tuple(sorted(refined_seams, key=lambda seam: seam.seam_id)),
        global_vertex_count=len(junctions),
        global_edge_count=len(cell_edges),
        global_face_count=len(patches),
        manifold_vertices=manifold_vertices,
        singular_vertices=singular_vertices,
        junction_link_status=tuple(sorted(link_status.items())),
        orientability=orientability,
        genus_by_component=genus_by_component,
        area_status="CERTIFIED",
        topological_cuts=topological_cuts,
        topological_face_count=len(patches),
        spherical_patch_record_count=len(patches),
        face_boundary_cycle_counts=face_boundary_cycle_counts,
        geometric_edge_count=len(global_edges),
        betti_numbers=betti_numbers,
        incidence_status=incidence_status,
        edge_incidence_distribution=edge_incidence_distribution,
    )
    return replace(surface, curvature=compute_molecular_contact_surface_curvature(surface))


def build_molecular_contact_surface_pair(
    atoms: Iterable[AtomBall],
    selected_contacts: Iterable[SelectedContact],
    **kwargs,
) -> tuple[MolecularContactSurface, MolecularContactSurface]:
    atoms = tuple(atoms)
    contacts = tuple(selected_contacts)
    return (
        build_molecular_contact_surface(atoms, contacts, "A", **kwargs),
        build_molecular_contact_surface(atoms, contacts, "B", **kwargs),
    )
