"""Cazals-style bicolor power-interface curvature on a weighted regular complex.

This analysis is independent of the AW solver and Apollonius complex. The
weighted GUDHI alpha complex supplies regular-simplex incidence and filtration;
power vertices and interface-edge lengths are computed from that incidence.
At alpha zero, alpha selects whole bicolor power facets. Curvature is measured
along each full finite power edge shared by two selected facets, following the
Cazals interior-edge definition. Triangle alpha membership is retained as an
audit field; it does not clip that dual edge or add another selection rule.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
import csv
import math
from pathlib import Path

import numpy as np


DEFAULT_ALPHA = 0.0
DEFAULT_ALPHA_TOLERANCE = 1e-12
DEFAULT_M = 5.0
DEFAULT_POWER_RESIDUAL_TOLERANCE = 1e-8


@dataclass(frozen=True)
class PowerVertex:
    tetrahedron: tuple[int, int, int, int]
    position: tuple[float, float, float]
    power_value: float
    residual: float


@dataclass
class PowerFacetRecord:
    generator_ids: tuple[int, int]
    group_orientation: str
    regular_edge_filtration: float
    alpha_zero_selected: bool
    condition_beta_status: str
    largest_orthogonal_ball_radius: float | None
    smaller_expanded_radius: float
    m_over_r: float | None
    condition_beta_accepted: bool | None
    final_selected: bool
    incident_regular_triangles: tuple[tuple[int, ...], ...]
    status: str
    reason: str = ""


@dataclass
class PowerInterfaceEdgeRecord:
    regular_triangle: tuple[int, int, int]
    group_composition: str
    singleton_generator: int
    singleton_group: str
    bicolor_facets: tuple[tuple[int, int], tuple[int, int]]
    incident_tetrahedra: tuple[tuple[int, ...], ...]
    power_vertex_coordinates: tuple[tuple[float, float, float], ...]
    alpha_zero_facets_selected: tuple[bool, bool]
    condition_beta_facets_accepted: tuple[bool | None, bool | None]
    triangle_alpha_zero_selected: bool
    geometry_status: str
    length: float | None
    full_dual_length: float | None
    clipping_status: str
    beta_radians: float
    beta_degrees: float
    sign: str
    length_beta: float | None
    length_abs_beta: float | None
    status: str
    reason: str = ""


@dataclass
class PowerInterfaceCurvatureResult:
    points: np.ndarray
    expanded_radii: np.ndarray
    weights: np.ndarray
    group_a: frozenset[int]
    group_b: frozenset[int]
    alpha: float
    alpha_tolerance: float
    condition_beta_m: float
    full_simplices: dict[int, set[tuple[int, ...]]]
    alpha_simplices: dict[int, set[tuple[int, ...]]]
    filtration: dict[tuple[int, ...], float]
    power_vertices: dict[tuple[int, ...], PowerVertex]
    unresolved_power_vertices: dict[tuple[int, ...], str]
    facets: list[PowerFacetRecord] = field(default_factory=list)
    edges: list[PowerInterfaceEdgeRecord] = field(default_factory=list)
    fallback_radius_assignments: int = 0

    @property
    def full_counts(self):
        return {d: len(self.full_simplices[d]) for d in range(4)}

    @property
    def alpha_counts(self):
        return {d: len(self.alpha_simplices[d]) for d in range(4)}

    @staticmethod
    def euler_characteristic(counts):
        return counts[0] - counts[1] + counts[2] - counts[3]

    @property
    def full_euler_characteristic(self):
        return self.euler_characteristic(self.full_counts)

    @property
    def alpha_euler_characteristic(self):
        return self.euler_characteristic(self.alpha_counts)

    @property
    def all_bicolor_facets(self):
        return [record for record in self.facets if record.generator_ids]

    @property
    def alpha_zero_facets(self):
        return [record for record in self.facets if record.alpha_zero_selected]

    @property
    def final_facets(self):
        return [record for record in self.facets if record.final_selected]

    @property
    def included_edges(self):
        return [record for record in self.edges if record.status == "included"]

    @property
    def edge_totals(self):
        records = self.included_edges
        length = float(sum(record.length or 0.0 for record in records))
        positive = float(sum(max(record.length_beta or 0.0, 0.0) for record in records))
        negative = float(sum(min(record.length_beta or 0.0, 0.0) for record in records))
        signed = positive + negative
        unsigned = float(sum(record.length_abs_beta or 0.0 for record in records))
        return {
            "length": length,
            "positive": positive,
            "negative": negative,
            "signed": signed,
            "unsigned": unsigned,
            "signed_mean_rad": signed / length if length else float("nan"),
            "unsigned_mean_rad": unsigned / length if length else float("nan"),
            "signed_mean_deg": math.degrees(signed / length) if length else float("nan"),
            "unsigned_mean_deg": math.degrees(unsigned / length) if length else float("nan"),
            "positive_edges": sum(record.beta_radians > 0 for record in records),
            "negative_edges": sum(record.beta_radians < 0 for record in records),
            "zero_edges": sum(record.beta_radians == 0 for record in records),
        }

    @property
    def full_dual_length_total(self):
        return float(sum(record.full_dual_length or 0.0 for record in self.included_edges))

    @property
    def edge_status_counts(self):
        return dict(Counter(record.status for record in self.edges))

    @property
    def edge_geometry_status_counts(self):
        return dict(Counter(record.geometry_status for record in self.edges))

    @property
    def smooth_surface_curvature(self):
        return 0.0

    def summary(self):
        totals = self.edge_totals
        full, alpha = self.full_counts, self.alpha_counts
        return (
            "Power-Diagram Bicolor Interface Curvature\n"
            f"Input atoms: {len(self.points):,}; group A/B: {len(self.group_a):,} / {len(self.group_b):,}\n"
            f"Fallback radius assignments: {self.fallback_radius_assignments:,}\n"
            f"Probe-expanded radii: supplied radius; power weight R^2\n"
            f"Alpha threshold: {self.alpha:g} (tolerance {self.alpha_tolerance:g})\n"
            f"Full regular complex V/E/F/T: {full[0]:,} / {full[1]:,} / {full[2]:,} / {full[3]:,}; "
            f"chi={self.full_euler_characteristic}\n"
            f"Alpha complex V/E/F/T: {alpha[0]:,} / {alpha[1]:,} / {alpha[2]:,} / {alpha[3]:,}; "
            f"chi={self.alpha_euler_characteristic}\n"
            f"Bicolor regular facets: {len(self.all_bicolor_facets):,}\n"
            f"Alpha-zero bicolor facets: {len(self.alpha_zero_facets):,}\n"
            f"Condition-beta M={self.condition_beta_m:g} rejected: "
            f"{sum(f.alpha_zero_selected and f.condition_beta_accepted is False for f in self.facets):,}\n"
            f"Condition-beta unresolved: {sum(f.alpha_zero_selected and f.condition_beta_accepted is None for f in self.facets):,}\n"
            f"Maximum resolved m/r among alpha-zero bicolor facets: "
            f"{max((f.m_over_r for f in self.facets if f.alpha_zero_selected and f.m_over_r is not None), default=float('nan')):.9g}\n"
            f"Final selected facets: {len(self.final_facets):,}\n"
            f"Power-interface edge candidates: {len(self.edges):,}; status counts: {self.edge_status_counts}\n"
            f"Dual geometry status counts: {self.edge_geometry_status_counts}\n"
            f"Included finite full interface edges: {len(self.included_edges):,}\n"
            f"Total edge length: {totals['length']:.9g} A\n"
            f"Full unrestricted length of included dual edges: {self.full_dual_length_total:.9g} A\n"
            f"C_plus = sum(l beta, beta>0): {totals['positive']:.9g} A rad\n"
            f"C_minus = sum(l beta, beta<0): {totals['negative']:.9g} A rad\n"
            f"C_signed = sum(l beta): {totals['signed']:.9g} A rad\n"
            f"C_unsigned = sum(l |beta|): {totals['unsigned']:.9g} A rad\n"
            f"Conventional H=(k1+k2)/2 edge measure, C_signed/2: {0.5 * totals['signed']:.9g} A rad\n"
            f"s_signed: {totals['signed_mean_rad']:.9g} rad / {totals['signed_mean_deg']:.9g} deg\n"
            f"s_unsigned: {totals['unsigned_mean_rad']:.9g} rad / {totals['unsigned_mean_deg']:.9g} deg\n"
            f"Positive/negative/zero edges: {totals['positive_edges']} / "
            f"{totals['negative_edges']} / {totals['zero_edges']}\n"
            f"Smooth facet curvature: {self.smooth_surface_curvature:g} (power facets are planar)\n"
            "Reported C is the raw turning measure sum(l beta); conventional H=(k1+k2)/2 edge measure is C/2.\n"
            "Geometry: full finite dual power-edge length; alpha=0 selects interface facets, not clipped edge segments.\n"
            "Sign: Cazals convention; ABB positive, AAB negative; |beta| is the singleton angle.\n"
            f"Unresolved power vertices: {len(self.unresolved_power_vertices):,}\n"
        )

    def export(self, directory, *, export_filtration=True):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        _write_rows(directory / "power_interface_surfaces.csv", self.facets)
        _write_rows(directory / "power_interface_edges.csv", self.edges)
        (directory / "power_interface_summary.txt").write_text(self.summary(), encoding="utf-8")
        if export_filtration:
            _write_filtration(directory / "power_simplex_filtration.csv", self)
        return directory


def calculate_power_vertex(tetrahedron, points, weights, condition_limit=1e12):
    """Solve the four equal-power equations and verify their common residual."""
    tet = tuple(sorted(int(i) for i in tetrahedron))
    if len(tet) != 4 or len(set(tet)) != 4:
        raise ValueError("A power vertex requires four distinct generators.")
    i = tet[0]
    matrix = np.asarray([2.0 * (points[j] - points[i]) for j in tet[1:]], dtype=float)
    rhs = np.asarray([
        np.dot(points[j], points[j]) - np.dot(points[i], points[i]) - weights[j] + weights[i]
        for j in tet[1:]
    ], dtype=float)
    condition = float(np.linalg.cond(matrix))
    if not np.isfinite(condition) or condition > condition_limit:
        raise np.linalg.LinAlgError(f"Ill-conditioned power vertex {tet}, cond={condition:g}.")
    position = np.linalg.solve(matrix, rhs)
    values = np.asarray([np.dot(position - points[j], position - points[j]) - weights[j]
                         for j in tet], dtype=float)
    common = float(np.mean(values))
    residual = float(np.max(np.abs(values - common)))
    return PowerVertex(tet, tuple(map(float, position)), common, residual)


def cazals_signed_beta(triangle, points, group_a, group_b):
    """Cazals' signed dihedral: ABB positive, AAB negative.

    The paper defines the magnitude as the angle at the singleton generator of
    the mixed-color regular triangle. Canonicalizing IDs makes vertex ordering
    irrelevant; exchanging group labels exchanges the sign.
    """
    tri = tuple(sorted(int(value) for value in triangle))
    a = [i for i in tri if i in group_a]
    b = [i for i in tri if i in group_b]
    if len(tri) != 3 or len(a) + len(b) != 3 or not a or not b:
        raise ValueError("Cazals beta requires a 1:2 bicolor triangle.")
    if len(a) == 1:
        singleton, others, sign, group = a[0], b, 1.0, "A"
    else:
        singleton, others, sign, group = b[0], a, -1.0, "B"
    u = np.asarray(points[others[0]], dtype=float) - points[singleton]
    v = np.asarray(points[others[1]], dtype=float) - points[singleton]
    nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    if min(nu, nv) <= 1e-14:
        raise ValueError("Degenerate regular triangle has coincident generators.")
    angle = float(np.arccos(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0)))
    return singleton, group, sign * angle


def oriented_power_facet_normal(edge, points, group_a, group_b):
    """Unit normal of a power bisector, oriented from group A to group B."""
    i, j = tuple(sorted(map(int, edge)))
    if i in group_a and j in group_b:
        a, b = i, j
    elif i in group_b and j in group_a:
        a, b = j, i
    else:
        raise ValueError("A power-interface facet must join groups A and B.")
    direction = np.asarray(points[b], dtype=float) - np.asarray(points[a], dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm <= 1e-14:
        raise ValueError("Coincident generators do not define a unique power-facet normal.")
    return direction / norm


def power_edge_geometry(triangle, incident_tetrahedra, power_vertices, length_tolerance=1e-10):
    """Return the full dual power-edge classification and finite segment data."""
    tets = tuple(sorted(tuple(sorted(map(int, tet))) for tet in incident_tetrahedra))
    coords = tuple(power_vertices[tet].position for tet in tets if tet in power_vertices)
    if len(tets) == 1:
        return "unbounded_boundary_power_edge", coords, None, "triangle has one incident tetrahedron; dual is a ray"
    if len(tets) != 2:
        status = "no_incident_tetrahedra" if not tets else "nonmanifold_incidence"
        return status, coords, None, f"expected two incident tetrahedra; found {len(tets)}"
    if any(tet not in power_vertices for tet in tets):
        return "unresolved_power_vertex", coords, None, "one or both endpoint power vertices are unresolved"
    length = float(np.linalg.norm(np.asarray(power_vertices[tets[1]].position)
                                  - np.asarray(power_vertices[tets[0]].position)))
    if length <= length_tolerance:
        return "degenerate_zero_length", coords, None, "dual power vertices coincide within length tolerance"
    return "finite", coords, length, ""


def clip_finite_power_edge_to_balls(triangle, endpoint_1, endpoint_2, points, radii,
                                   power_tolerance=1e-8):
    """Measure the exact segment of a finite power edge inside its 3 balls.

    On the power edge dual to ``triangle``, the three power distances agree.
    The restricted alpha-complex edge is therefore the interval where
    q(s)=||x(s)-p_i||^2-R_i^2 <= 0. Since x(s) is linear, q is quadratic and
    its sublevel interval is solved analytically, then intersected with the
    finite full-edge parameter interval [0,1].
    """
    ids = tuple(sorted(map(int, triangle)))
    first, second = np.asarray(endpoint_1, dtype=float), np.asarray(endpoint_2, dtype=float)
    direction = second - first
    norm2 = float(np.dot(direction, direction))
    if len(ids) != 3 or norm2 <= 1e-20:
        return None, "degenerate_power_edge"
    origin = np.asarray(points[ids[0]], dtype=float)
    radius = float(radii[ids[0]])
    projection = -float(np.dot(first - origin, direction)) / norm2
    nearest = first + projection * direction
    minimum_power = float(np.dot(nearest - origin, nearest - origin) - radius * radius)
    scale = max(1.0, radius * radius, float(np.dot(nearest - origin, nearest - origin)))
    if minimum_power > power_tolerance * scale:
        return 0.0, "alpha_triangle_has_no_numerical_ball_intersection"
    half_width = math.sqrt(max(0.0, -minimum_power / norm2))
    lower = max(0.0, projection - half_width)
    upper = min(1.0, projection + half_width)
    fraction = max(0.0, upper - lower)
    q_reference = lambda s: float(np.dot(first + s * direction - origin,
                                         first + s * direction - origin) - radius * radius)
    for atom_id in ids[1:]:
        check_radius = float(radii[atom_id])
        p = np.asarray(points[atom_id], dtype=float)
        q = lambda s: float(np.dot(first + s * direction - p, first + s * direction - p)
                             - check_radius * check_radius)
        error = max(abs(q(lower) - q_reference(lower)),
                    abs(q(upper) - q_reference(upper)))
        if error > 1e-6 * max(1.0, check_radius * check_radius):
            return None, f"power equality residual on clipped boundary ({error:.6g})"
    if fraction <= 1e-14:
        return 0.0, "zero_length_restricted_alpha_edge"
    return float(math.sqrt(norm2) * fraction), "clipped_to_alpha_balls"


def _condition_beta_for_facet(edge, incident_triangles, triangle_to_tetrahedra,
                              power_vertices, radii, threshold):
    """Evaluate the largest empty orthogonal-ball radius on a power facet.

    The power function is convex on the facet plane, so its maximum on a
    bounded polygon occurs at a polygon vertex, i.e. at one of the incident
    power vertices. An unbounded facet has no finite maximum and is rejected.
    ``m/r > M`` is the large-facet rejection; r is the smaller expanded atomic
    ball radius used to construct this weighted diagram.
    """
    edge = tuple(sorted(edge))
    if not incident_triangles:
        return "unresolved_no_facet_boundary", None, float("inf"), False, "no regular triangles in facet star"
    incident_tets = set()
    for tri in incident_triangles:
        count = len(triangle_to_tetrahedra.get(tri, ()))
        if count == 1:
            return "rejected_unbounded", None, float("inf"), False, "facet has a ray on its power-diagram boundary"
        if count == 0:
            return "unresolved_missing_power_vertex", None, float("inf"), False, "facet boundary triangle has no incident tetrahedron"
        if count > 2:
            return "unresolved_nonmanifold_star", None, float("inf"), False, f"boundary triangle has {count} incident tetrahedra"
        incident_tets.update(triangle_to_tetrahedra[tri])
    if any(tet not in power_vertices for tet in incident_tets):
        return "unresolved_power_vertex", None, float("inf"), False, "facet star contains unresolved power vertex"
    max_power = max(power_vertices[tet].power_value for tet in incident_tets)
    radius = math.sqrt(max(0.0, max_power))
    smaller = float(min(radii[edge[0]], radii[edge[1]]))
    ratio = radius / smaller
    accepted = ratio <= threshold
    return ("accepted" if accepted else "rejected_m_over_r", radius, ratio, accepted,
            "" if accepted else f"largest empty orthogonal ball ratio {ratio:.8g} exceeds M={threshold:g}")


def analyze_power_interface_curvature(
        points,
        expanded_radii,
        group_a_ids,
        group_b_ids,
        *,
        alpha=DEFAULT_ALPHA,
        alpha_tolerance=DEFAULT_ALPHA_TOLERANCE,
        condition_beta_m=DEFAULT_M,
        power_residual_tolerance=DEFAULT_POWER_RESIDUAL_TOLERANCE,
        fallback_radius_assignments=0):
    """Construct and measure the power interface from weighted GUDHI output.

    ``expanded_radii`` are the actual ball radii R_i, already including any
    probe radius. The Cazals PDB helper supplies R_i = r_i + probe and weights
    R_i^2. Alpha filtration selects complete regular simplices at ``alpha``;
    only bicolor regular edges use that selection in facet membership.
    """
    points = np.asarray(points, dtype=float)
    radii = np.asarray(expanded_radii, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or radii.shape != (len(points),):
        raise ValueError("points must be (N,3) and expanded_radii must be (N,).")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(radii)) or np.any(radii <= 0):
        raise ValueError("Power input coordinates/radii must be finite and radii positive.")
    group_a = frozenset(map(int, group_a_ids))
    group_b = frozenset(map(int, group_b_ids))
    if not group_a or not group_b or group_a & group_b:
        raise ValueError("Groups must be nonempty and disjoint.")
    if not (group_a | group_b) <= set(range(len(points))):
        raise IndexError("A group generator ID is outside the point array.")
    if condition_beta_m <= 0 or alpha_tolerance < 0:
        raise ValueError("M must be positive and alpha_tolerance nonnegative.")

    # Import lazily so the pure geometric helpers remain usable without GUDHI.
    from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation import (
        cazals_power_interface as power_backend,
    )
    simplex_tree, vertices, edges, triangles, tetrahedra = power_backend.build_regular_complex(points, radii)
    full = {0: {(i,) for i in vertices}, 1: set(edges), 2: set(triangles), 3: set(tetrahedra)}
    filtration = {
        tuple(sorted(map(int, simplex))): float(value)
        for simplex, value in simplex_tree.get_filtration()
    }
    return analyze_power_regular_complex(
        points, radii, group_a, group_b, full, filtration,
        alpha=alpha, alpha_tolerance=alpha_tolerance,
        condition_beta_m=condition_beta_m,
        power_residual_tolerance=power_residual_tolerance,
        fallback_radius_assignments=fallback_radius_assignments,
    )


def analyze_power_regular_complex(
        points, expanded_radii, group_a_ids, group_b_ids, full_simplices, filtration, *,
        alpha=DEFAULT_ALPHA,
        alpha_tolerance=DEFAULT_ALPHA_TOLERANCE,
        condition_beta_m=DEFAULT_M,
        power_residual_tolerance=DEFAULT_POWER_RESIDUAL_TOLERANCE,
        fallback_radius_assignments=0):
    """Analyze audited regular-complex incidence and filtration without rebuilding it.

    This entry point is useful for reproducible audits from a saved GUDHI
    simplex/filtration table. It preserves exactly the supplied incidence and
    does not infer any simplex from coordinates.
    """
    points = np.asarray(points, dtype=float)
    radii = np.asarray(expanded_radii, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or radii.shape != (len(points),):
        raise ValueError("points must be (N,3) and expanded_radii must be (N,).")
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite(radii)) or np.any(radii <= 0):
        raise ValueError("Power input coordinates/radii must be finite and radii positive.")
    group_a, group_b = frozenset(map(int, group_a_ids)), frozenset(map(int, group_b_ids))
    if not group_a or not group_b or group_a & group_b:
        raise ValueError("Groups must be nonempty and disjoint.")
    if not (group_a | group_b) <= set(range(len(points))):
        raise IndexError("A group generator ID is outside the point array.")
    full = {d: {tuple(sorted(map(int, simplex))) for simplex in full_simplices[d]}
            for d in range(4)}
    filtration = {tuple(sorted(map(int, simplex))): float(value)
                 for simplex, value in filtration.items()}
    if condition_beta_m <= 0 or alpha_tolerance < 0:
        raise ValueError("M must be positive and alpha_tolerance nonnegative.")
    edges, triangles, tetrahedra = full[1], full[2], full[3]
    weights = radii ** 2
    alpha_simplices = {
        dimension: {simplex for simplex in full[dimension]
                    if filtration.get(simplex, float("inf")) <= alpha + alpha_tolerance}
        for dimension in range(4)
    }
    full_closure = _closure_issues(full)
    alpha_closure = _closure_issues(alpha_simplices)
    if full_closure or alpha_closure:
        raise RuntimeError(
            "Weighted GUDHI complex is not simplicially closed: "
            f"full={full_closure[:5]}, alpha={alpha_closure[:5]}"
        )

    tri_to_tets = defaultdict(list)
    edge_to_tris = defaultdict(list)
    for tet in tetrahedra:
        for tri in _combinations(tet, 3):
            tri_to_tets[tri].append(tet)
    for tri in triangles:
        for edge in _combinations(tri, 2):
            edge_to_tris[edge].append(tri)

    power_vertices, unresolved = {}, {}
    for tet in sorted(tetrahedra):
        try:
            pv = calculate_power_vertex(tet, points, weights)
            if pv.residual > power_residual_tolerance:
                unresolved[tet] = f"equal-power residual {pv.residual:.8g} exceeds tolerance"
            else:
                power_vertices[tet] = pv
        except (np.linalg.LinAlgError, ValueError) as error:
            unresolved[tet] = str(error)

    bicolor_edges = {edge for edge in edges if _is_bicolor(edge, group_a, group_b)}
    facet_map = {}
    for edge in sorted(bicolor_edges):
        alpha_selected = edge in alpha_simplices[1]
        tris = tuple(sorted(edge_to_tris.get(edge, ())))
        if alpha_selected:
            beta_status, m, ratio, beta_ok, reason = _condition_beta_for_facet(
                edge, tris, tri_to_tets, power_vertices, radii, condition_beta_m
            )
        else:
            beta_status, m, ratio, beta_ok, reason = "not_evaluated_not_alpha0", None, None, None, "facet is not in alpha=0 complex"
        if edge[0] in group_a:
            orientation = f"A:{edge[0]}->B:{edge[1]}"
        else:
            orientation = f"A:{edge[1]}->B:{edge[0]}"
        final = alpha_selected and beta_ok is True
        facet_map[edge] = PowerFacetRecord(
            edge, orientation, filtration[edge], alpha_selected, beta_status,
            m, float(min(radii[list(edge)])), ratio, beta_ok, final, tris,
            "selected" if final else "excluded", reason,
        )

    result = PowerInterfaceCurvatureResult(
        points, radii, weights, group_a, group_b, float(alpha), float(alpha_tolerance),
        float(condition_beta_m), full, alpha_simplices, filtration, power_vertices,
        unresolved, list(facet_map.values()), fallback_radius_assignments=int(fallback_radius_assignments),
    )
    for tri in sorted(triangles):
        bicolor = tuple(sorted(edge for edge in _combinations(tri, 2)
                               if edge in facet_map))
        if len(bicolor) != 2:
            continue
        singleton, singleton_group, beta = cazals_signed_beta(tri, points, group_a, group_b)
        tets = tuple(sorted(tri_to_tets.get(tri, ())))
        geometry_status, coords, length, geom_reason = power_edge_geometry(
            tri, tets, power_vertices
        )
        full_length = length
        clipping_status = "not_clipped_cazals_full_dual_edge"
        alpha_flags = tuple(edge in alpha_simplices[1] for edge in bicolor)
        beta_flags = tuple(facet_map[edge].condition_beta_accepted for edge in bicolor)
        final_facets = all(facet_map[edge].final_selected for edge in bicolor)
        if not all(alpha_flags):
            status, reason = "excluded", "one or both bicolor facets fail the alpha threshold"
        elif not all(flag is True for flag in beta_flags):
            status, reason = "excluded", "one or both bicolor facets fail or cannot resolve condition beta"
        elif geometry_status != "finite":
            status, reason = geometry_status, geom_reason
        elif not final_facets:
            status, reason = "excluded", "incident facet selection is incomplete"
        else:
            status, reason = "included", ""
        result.edges.append(PowerInterfaceEdgeRecord(
            tri, _composition(tri, group_a, group_b), singleton, singleton_group,
            bicolor, tets, coords, alpha_flags, beta_flags,
            tri in alpha_simplices[2], geometry_status, length, full_length,
            clipping_status, beta,
            math.degrees(beta), "positive" if beta > 0 else "negative" if beta < 0 else "zero",
            None if length is None else length * beta,
            None if length is None else length * abs(beta), status, reason,
        ))
    return result


def analyze_prepared_cazals(prepared, **kwargs):
    """Convenience entry point for ``PreparedCazalsPDB`` from cazals_pdb."""
    return analyze_power_interface_curvature(
        prepared.points, prepared.expanded_radii, prepared.group_a, prepared.group_b,
        fallback_radius_assignments=len(prepared.fallback_atoms), **kwargs,
    )


def _combinations(simplex, count):
    from itertools import combinations
    return [tuple(sorted(map(int, values))) for values in combinations(simplex, count)]


def _closure_issues(simplices):
    missing = []
    for dimension in range(1, 4):
        for simplex in simplices[dimension]:
            for face in _combinations(simplex, dimension):
                if face not in simplices[dimension - 1]:
                    missing.append((simplex, face))
    return missing


def _is_bicolor(edge, group_a, group_b):
    i, j = edge
    return (i in group_a and j in group_b) or (i in group_b and j in group_a)


def _composition(simplex, group_a, group_b):
    return "".join(sorted("A" if index in group_a else "B" if index in group_b else "?"
                          for index in simplex))


def _write_rows(path, records):
    rows = [asdict(record) for record in records]
    fields = tuple(rows[0]) if rows else ()
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if fields:
            writer.writeheader()
            for row in rows:
                for key, value in row.items():
                    if isinstance(value, (tuple, list)):
                        row[key] = repr(value)
                writer.writerow(row)


def _write_filtration(path, result):
    with path.open("w", newline="", encoding="utf-8") as handle:
        fields = ("dimension", "simplex", "filtration", "alpha_threshold", "in_alpha")
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for simplex, value in sorted(result.filtration.items(), key=lambda item: (len(item[0]), item[0])):
            writer.writerow({
                "dimension": len(simplex) - 1,
                "simplex": " ".join(map(str, simplex)),
                "filtration": f"{value:.17g}",
                "alpha_threshold": result.alpha,
                "in_alpha": int(value <= result.alpha + result.alpha_tolerance),
            })
