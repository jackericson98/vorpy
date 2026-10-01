"""Bicolor interface curvature on existing additively weighted networks.

Apollonius simplices provide only the incidence filter. All geometry and
curvature come from VorPy's existing AW surface and analytic edge records.
"""

from collections import Counter
from dataclasses import asdict, dataclass, field
import csv
import math
from pathlib import Path

import numpy as np

from vorpy.src.analyze.apollonius import ApolloniusComplex
from vorpy.src.calculations.aw_interface_orientation import aw_interface_edge_orientation
from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge


@dataclass
class InterfaceSurfaceRecord:
    surface_id: int
    generator_ids: tuple[int, int]
    atom_numbers: tuple[int, int]
    group_orientation: str
    area: float | None
    bounded: bool
    feature_complete: bool
    generator_cells_complete: bool
    supported: bool
    included: bool
    dual_edge: tuple[int, int]
    integrated_mean_curvature: float | None
    status: str
    reason: str = ""


@dataclass
class InterfaceEdgeRecord:
    edge_id: int
    generator_ids: tuple[int, int, int]
    atom_numbers: tuple[int, int, int]
    dual_triangle: tuple[int, int, int]
    surface_ids: tuple[int, int]
    geometry_type: str | None
    supported: bool
    feature_complete: bool
    generator_cells_complete: bool
    length: float | None
    signed_contribution: float | None
    unsigned_contribution: float | None
    positive_contribution: float | None
    negative_contribution: float | None
    beta_min: float | None
    beta_max: float | None
    beta_mean: float | None
    beta_statistics_method: str | None
    integration_status: str
    estimated_error: float | None
    quadrature_order: int | None
    status: str
    reason: str = ""


@dataclass
class AWInterfaceCurvatureResult:
    network: object
    complex: ApolloniusComplex
    group_a: frozenset[int]
    group_b: frozenset[int]
    selection_mode: str
    surfaces: list[InterfaceSurfaceRecord] = field(default_factory=list)
    edges: list[InterfaceEdgeRecord] = field(default_factory=list)
    integration_options: dict = field(default_factory=dict)

    @property
    def selected_surfaces(self):
        return [record for record in self.surfaces if record.included]

    @property
    def selected_edges(self):
        return [record for record in self.edges if record.status == "included"]

    @property
    def group_atom_counts(self):
        return {"A": len(self.group_a), "B": len(self.group_b)}

    @property
    def interface_area(self):
        return float(sum(record.area or 0.0 for record in self.selected_surfaces))

    @property
    def surface_curvature(self):
        return float(sum(record.integrated_mean_curvature or 0.0
                         for record in self.selected_surfaces))

    @property
    def edge_curvature(self):
        selected = self.selected_edges
        return {
            "positive": float(sum(record.positive_contribution or 0.0 for record in selected)),
            "negative": float(sum(record.negative_contribution or 0.0 for record in selected)),
            "signed": float(sum(record.signed_contribution or 0.0 for record in selected)),
            "unsigned": float(sum(record.unsigned_contribution or 0.0 for record in selected)),
            "length": float(sum(record.length or 0.0 for record in selected)),
        }

    @property
    def edge_angle_statistics(self):
        measures = self.edge_curvature
        length = measures["length"]
        signed = measures["signed"] / length if length else 0.0
        unsigned = measures["unsigned"] / length if length else 0.0
        positive = negative = mixed = 0
        for edge in self.selected_edges:
            has_pos = (edge.positive_contribution or 0.0) > 1e-10
            has_neg = (edge.negative_contribution or 0.0) < -1e-10
            if has_pos and has_neg:
                mixed += 1
            elif has_pos:
                positive += 1
            elif has_neg:
                negative += 1
        return {
            "signed_radians": signed,
            "signed_degrees": math.degrees(signed),
            "unsigned_radians": unsigned,
            "unsigned_degrees": math.degrees(unsigned),
            "positive_edges": positive,
            "negative_edges": negative,
            "mixed_sign_edges": mixed,
        }

    @property
    def unresolved_reasons(self):
        return dict(Counter(
            record.reason or record.status
            for record in [*self.surfaces, *self.edges]
            if record.status != "included"
        ))

    @property
    def exclusion_diagnostics(self):
        """Count causes independently; a record may have multiple flags."""
        counts = Counter()
        for record in self.surfaces:
            if not record.included:
                if not record.bounded or not record.feature_complete:
                    counts["unbounded_or_incomplete_surface"] += 1
                if not record.generator_cells_complete:
                    counts["incomplete_surface_generator_cells"] += 1
                if not record.supported:
                    counts["unsupported_surface_incidence"] += 1
                if record.integrated_mean_curvature is None:
                    counts["missing_surface_curvature"] += 1
        for record in self.edges:
            if record.status == "included":
                continue
            if not record.supported:
                counts["unsupported_edge_incidence"] += 1
            if not record.feature_complete:
                counts["unbounded_or_incomplete_edge"] += 1
            if not record.generator_cells_complete:
                counts["incomplete_edge_generator_cells"] += 1
            if record.integration_status == "maximum_order_reached":
                counts["quadrature_not_converged"] += 1
            if record.reason.startswith("analytic edge integration failed"):
                counts["analytic_integration_failure"] += 1
        return dict(counts)

    @property
    def combined_curvature(self):
        """Conventional oriented IMC of selected smooth patches and seams.

        The smooth cache stores ``integral H dA`` for H=(k1+k2)/2; edge
        records store raw ``integral beta ds``.  The piecewise-smooth measure
        is consequently their sum with a factor 1/2 on the edge term.
        """
        return float(self.surface_curvature + 0.5*self.edge_curvature["signed"])

    def summary(self):
        edge = self.edge_curvature
        angles = self.edge_angle_statistics
        supported_edges = sum(record.supported for record in self.edges)
        excluded_edges = len(self.edges) - len(self.selected_edges)
        return (
            "Bicolor AW Interface Curvature\n"
            f"Selection mode: {self.selection_mode}\n"
            f"Group A/B atoms: {len(self.group_a):,} / {len(self.group_b):,}\n"
            f"Bicolor AW surfaces selected: {len(self.selected_surfaces):,} / "
            f"{sum(r.bounded and r.feature_complete and r.area is not None for r in self.surfaces):,} geometrically eligible\n"
            f"Interface area: {self.interface_area:.9g} A^2\n"
            f"Bicolor edge candidates: {len(self.edges):,}\n"
            f"Supported edge candidates: {supported_edges:,}\n"
            f"Integrated edges: {len(self.selected_edges):,}; excluded/unresolved: {excluded_edges:,}\n"
            f"Total edge length: {edge['length']:.9g} A\n"
            f"C_edge positive raw integral(beta ds): {edge['positive']:.9g} A rad\n"
            f"C_edge negative raw integral(beta ds): {edge['negative']:.9g} A rad\n"
            f"C_edge signed raw integral(beta ds): {edge['signed']:.9g} A rad\n"
            f"C_edge unsigned raw integral(|beta| ds): {edge['unsigned']:.9g} A rad\n"
            f"s_H,edge signed: {angles['signed_radians']:.9g} rad / {angles['signed_degrees']:.9g} deg\n"
            f"s_H,edge unsigned: {angles['unsigned_radians']:.9g} rad / {angles['unsigned_degrees']:.9g} deg\n"
            f"Positive/negative/mixed edges: {angles['positive_edges']} / "
            f"{angles['negative_edges']} / {angles['mixed_sign_edges']}\n"
            f"C_smooth = integral(H dA): {self.surface_curvature:.9g} A; H=(k1+k2)/2\n"
            "Surface normal: AW clearance gradient, oriented A toward B\n"
            "Edge tangent: induced boundary orientation t=n x m_out; analytic tangent supplies ds only\n"
            "Beta extrema/mean: Gauss-Legendre quadrature-node diagnostics\n"
            f"C_H = C_smooth + 1/2 C_edge,signed: {self.combined_curvature:.9g} A\n"
            f"Quadrature convergence: {sum(r.integration_status == 'converged' for r in self.edges)} converged, "
            f"{sum(r.integration_status == 'maximum_order_reached' for r in self.edges)} reached maximum order\n"
            f"Exclusion diagnostics (independent flags): {self.exclusion_diagnostics}\n"
            f"Unresolved/excluded reasons: {self.unresolved_reasons}\n"
            f"Integration: {self.integration_options}\n"
        )

    def export_csv(self, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        _write_dataclasses(directory / "aw_interface_surfaces.csv", self.surfaces)
        _write_dataclasses(directory / "aw_interface_edges.csv", self.edges)
        (directory / "aw_interface_summary.txt").write_text(self.summary(), encoding="utf-8")


def analyze_aw_interface_curvature(
        network,
        group_a,
        group_b,
        *,
        selection_mode="supported",
        quadrature_order=16,
        max_quadrature_order=256,
        relative_tolerance=1e-7,
        absolute_tolerance=1e-9):
    """Measure a bicolor AW interface using solved surfaces and edges.

    Group arguments are network-local generator IDs. ``raw`` includes finite,
    complete feature geometry while preserving its support flags; ``supported``
    integrates only features in the supported Apollonius subcomplex. ``alpha``
    is intentionally unavailable until the required pair-surface alpha births
    are exact.
    """
    mode = str(selection_mode).lower()
    if mode not in {"raw", "supported", "alpha"}:
        raise ValueError("selection_mode must be 'raw', 'supported', or 'alpha'.")
    if mode == "alpha":
        raise ValueError(
            "Alpha interface selection is unavailable: pair-surface alpha births are unresolved; "
            "no approximation is applied."
        )
    if (getattr(network, "settings", None) or {}).get("net_type", "aw") != "aw":
        raise ValueError("Bicolor interface curvature requires an AW network.")

    group_a = frozenset(int(value) for value in group_a)
    group_b = frozenset(int(value) for value in group_b)
    valid_ids = set(int(index) for index in network.balls.index)
    if not group_a or not group_b:
        raise ValueError("Both interface groups must contain at least one generator.")
    if group_a & group_b:
        raise ValueError("Interface groups A and B must be disjoint.")
    if not (group_a | group_b) <= valid_ids:
        raise ValueError("Interface group contains unknown network generator IDs.")
    if max_quadrature_order < quadrature_order:
        raise ValueError("max_quadrature_order must be >= quadrature_order.")

    complex_ = ApolloniusComplex(network, getattr(network, "group_name", "Network"))
    result = AWInterfaceCurvatureResult(
        network, complex_, group_a, group_b, mode,
        integration_options={
            "initial_order": int(quadrature_order),
            "max_order": int(max_quadrature_order),
            "relative_tolerance": float(relative_tolerance),
            "absolute_tolerance": float(absolute_tolerance),
            "parameterization": "analytic AW edge; ds=|r'(t)|dt",
        },
    )

    surface_curvatures = _oriented_surface_curvatures(network)
    selected_surface_ids = set()
    candidate_surface_ids = set()
    for surface_id, row in network.surfs.iterrows():
        ids = tuple(int(value) for value in row.get("balls", ()))
        if len(ids) != 2:
            continue
        a, b = ids
        if not ((a in group_a and b in group_b) or (a in group_b and b in group_a)):
            continue
        dual_key = tuple(sorted(ids))
        dual_edge = complex_.simplices[1].get(dual_key)
        feature = next((ref for ref in (dual_edge.primal_features if dual_edge else ())
                        if ref.kind == "surf" and ref.feature_id == int(surface_id)), None)
        bounded = bool(feature and feature.bounded)
        feature_complete = bool(feature and feature.complete)
        cells_complete = bool(feature and feature.generator_cells_complete)
        supported = bool(dual_edge and complex_.is_supported(1, dual_key) and feature)
        area = _finite_float(row.get("sa"))
        oriented_cell = a if a in group_a else b
        value = surface_curvatures.get(int(surface_id), {}).get(oriented_cell)
        available = bounded and feature_complete and area is not None and area >= 0.0
        geometry_complete = cells_complete
        candidate = available
        if candidate:
            candidate_surface_ids.add(int(surface_id))
        included = available and geometry_complete and (mode == "raw" or supported)
        reason = "" if included else (
            "surface geometry is unbounded, incomplete, or area is unavailable" if not available
            else "participating AW generator cells are incomplete" if not geometry_complete
            else "surface incidence is not supported"
        )
        if value is None or not np.isfinite(value):
            included = False
            reason = "oriented AW surface mean-curvature value is unavailable"
        if included:
            selected_surface_ids.add(int(surface_id))
            status = "included"
        else:
            status = "excluded"
        record = InterfaceSurfaceRecord(
            int(surface_id), ids, _atom_numbers(network, ids),
            f"A:{oriented_cell}->B:{b if oriented_cell == a else a}",
            area, bounded, feature_complete, cells_complete, supported, included,
            dual_key, None if value is None else float(value), status, reason,
        )
        result.surfaces.append(record)

    for edge_id, row in network.edges.iterrows():
        ids = tuple(int(value) for value in row.get("balls", ()))
        if len(ids) != 3 or not (set(ids) & group_a and set(ids) & group_b):
            continue
        incident_surface_ids = tuple(int(value) for value in row.get("surfs", ()))
        bicolor_ids = tuple(sorted(
            surface_id for surface_id in incident_surface_ids
            if surface_id in candidate_surface_ids
        ))
        if len(bicolor_ids) != 2:
            continue
        tri_key = tuple(sorted(ids))
        dual_triangle = complex_.simplices[2].get(tri_key)
        features = [ref for ref in (dual_triangle.primal_features if dual_triangle else ())
                    if ref.kind == "edge" and ref.feature_id == int(edge_id)]
        feature_complete = any(feature.complete and feature.bounded for feature in features)
        cells_complete = any(feature.generator_cells_complete for feature in features)
        supported = bool(dual_triangle and complex_.is_supported(2, tri_key) and features)
        both_surfaces_included = all(sid in selected_surface_ids for sid in bicolor_ids)
        usable_by_mode = (mode == "raw" or supported) and feature_complete and cells_complete and both_surfaces_included
        if not usable_by_mode:
            reasons = []
            if not bicolor_ids:
                reasons.append("primal edge does not reference two bicolor AW surfaces")
            if mode == "supported" and not supported:
                reasons.append("AW edge incidence is not in the supported Apollonius complex")
            if not feature_complete:
                reasons.append("AW edge geometry is unbounded or incomplete")
            if not cells_complete:
                reasons.append("participating AW generator cells are incomplete")
            if not both_surfaces_included:
                reasons.append("one or both incident bicolor surfaces are excluded")
            reason = "; ".join(reasons)
            result.edges.append(_edge_excluded_record(
                network, edge_id, ids, tri_key, bicolor_ids, supported,
                feature_complete, cells_complete, reason,
            ))
            continue
        try:
            resolved = resolve_aw_network_edge(network, edge_id)
            integral = _integrate_edge_beta(
                network, resolved, bicolor_ids, group_a, group_b,
                order=int(quadrature_order), max_order=int(max_quadrature_order),
                rtol=float(relative_tolerance), atol=float(absolute_tolerance),
            )
            status = "included" if integral["status"] == "converged" else "unresolved"
            reason = "" if status == "included" else "edge quadrature did not converge by maximum order"
            result.edges.append(InterfaceEdgeRecord(
                int(edge_id), ids, _atom_numbers(network, ids), tri_key, bicolor_ids,
                type(resolved.geometry).__name__, supported, feature_complete,
                cells_complete, integral["length"], integral["signed"],
                integral["unsigned"], integral["positive"], integral["negative"],
                integral["beta_min"], integral["beta_max"], integral["beta_mean"],
                "Gauss-Legendre quadrature-node diagnostic",
                integral["status"], integral["error"], integral["order"],
                status, reason,
            ))
        except Exception as error:
            result.edges.append(_edge_excluded_record(
                network, edge_id, ids, tri_key, bicolor_ids, supported,
                feature_complete, cells_complete, f"analytic edge integration failed: {error}",
            ))
    return result


def _oriented_surface_curvatures(network):
    from vorpy.src.calculations.surface_mean_curvature import (
        calculate_aw_network_surface_mean_curvatures,
    )

    if "int_mean_curv_by_ball" in network.surfs.columns:
        values = network.surfs["int_mean_curv_by_ball"].tolist()
    else:
        values = calculate_aw_network_surface_mean_curvatures(network)
    return {int(surface_id): value for surface_id, value
            in zip(network.surfs.index, values, strict=True)}


def _integrate_edge_beta(network, resolved, surface_ids, group_a, group_b,
                         order, max_order, rtol, atol):
    geometry = resolved.geometry
    surfaces = sorted(surface_ids)
    surface_pairs = [tuple(int(value) for value in network.surfs.loc[sid, "balls"])
                     for sid in surfaces]
    edge_row = network.edges.loc[resolved.edge_index]
    generator_ids = tuple(int(value) for value in edge_row["balls"])
    locations = {
        generator_id: np.asarray(network.balls.loc[generator_id, "loc"], dtype=float)
        for generator_id in generator_ids
    }

    def integrate_at(current_order):
        nodes, weights = np.polynomial.legendre.leggauss(current_order)
        lo, hi = float(geometry.t_min), float(geometry.t_max)
        half, center = 0.5 * (hi - lo), 0.5 * (hi + lo)
        signed = absolute = positive = negative = length = 0.0
        beta_values = []
        for node, weight in zip(nodes, weights, strict=True):
            parameter = center + half * float(node)
            point = np.asarray(geometry.point(parameter), dtype=float)
            tangent = np.asarray(geometry.tangent(parameter), dtype=float)
            speed = float(np.linalg.norm(tangent))
            if point.shape != (3,) or not np.all(np.isfinite(point)) or not np.isfinite(speed) or speed <= 1e-12:
                raise ValueError("AW edge point/tangent is non-finite or singular.")
            oriented = aw_interface_edge_orientation(
                point, surface_pairs, generator_ids, group_a, group_b, locations,
            )
            # The analytic parameter direction defines ds but not the sign.
            # It must nevertheless describe the same tangent line as the
            # topology-induced oriented boundary tangent.
            tangent_alignment = abs(float((tangent / speed) @ oriented["tangent"]))
            if not np.isfinite(tangent_alignment) or 1.0 - tangent_alignment > 1e-5:
                raise ValueError(
                    "Analytic edge tangent is not parallel to the topology-induced "
                    f"boundary tangent (absolute dot={tangent_alignment:.9g})."
                )
            beta = oriented["beta"]
            ds_weight = half * float(weight) * speed
            signed += beta * ds_weight
            absolute += abs(beta) * ds_weight
            positive += max(beta, 0.0) * ds_weight
            negative += min(beta, 0.0) * ds_weight
            length += ds_weight
            beta_values.append(beta)
        return {"signed": signed, "unsigned": absolute, "positive": positive,
                "negative": negative, "length": length,
                "beta_min": min(beta_values), "beta_max": max(beta_values),
                "beta_mean": float(np.mean(beta_values))}

    current_order = max(2, int(order))
    previous = integrate_at(current_order)
    error = float("inf")
    while current_order * 2 <= max_order:
        current_order *= 2
        current = integrate_at(current_order)
        keys = ("signed", "unsigned", "positive", "negative", "length")
        error = max(abs(current[key] - previous[key]) for key in keys)
        converged = all(abs(current[key] - previous[key]) <= atol + rtol * abs(current[key])
                        for key in keys)
        previous = current
        if converged:
            return {**current, "status": "converged", "error": error, "order": current_order}
    return {**previous, "status": "maximum_order_reached", "error": error,
            "order": current_order}


def _edge_excluded_record(network, edge_id, ids, tri_key, surface_ids, supported,
                          feature_complete, cells_complete, reason):
    return InterfaceEdgeRecord(
        int(edge_id), ids, _atom_numbers(network, ids), tri_key, tuple(surface_ids),
        None, supported, feature_complete, cells_complete, None, None, None, None,
        None, None, None, None, None, "not_integrated", None, None, "excluded", reason,
    )


def _atom_numbers(network, ids):
    if "num" not in network.balls.columns:
        return tuple(int(value) for value in ids)
    return tuple(int(network.balls.at[index, "num"]) for index in ids)


def _finite_float(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _write_dataclasses(path, rows):
    fields = tuple(asdict(rows[0]).keys()) if rows else ()
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if fields:
            writer.writeheader()
            for row in rows:
                values = asdict(row)
                for key, value in values.items():
                    if isinstance(value, tuple):
                        values[key] = ";".join(map(str, value))
                writer.writerow(values)
