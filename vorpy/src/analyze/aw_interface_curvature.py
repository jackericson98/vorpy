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
from vorpy.src.calculations.edge_geometry import (
    AdditivelyWeightedTrisectorBranch,
    match_aw_trisector_conic_to_samples,
)
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
    surface_ids: tuple[int, ...]
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
    # Conventional H=(k1+k2)/2 seam contribution.  The existing signed and
    # unsigned fields intentionally remain the raw beta ds integrals.
    signed_mean_curvature: float | None = None
    unsigned_mean_curvature: float | None = None


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
    def edge_contribution_complete(self):
        """Whether every retained bicolor seam has a validated C_e."""
        return not self.missing_edge_ids and all(
            record.status == "included" for record in self.edges
        )

    @property
    def missing_edge_ids(self):
        """Surface-referenced edge rows absent from the analyzed network."""
        referenced = set()
        surfaces = getattr(self.network, "surfs", None)
        if surfaces is None:
            return ()
        for record in self.selected_surfaces:
            try:
                values = surfaces.loc[record.surface_id].get("edges", ())
                referenced.update(int(value) for value in values)
            except (KeyError, TypeError, ValueError):
                continue
        represented = {record.edge_id for record in self.edges}
        return tuple(sorted(referenced - represented))

    @property
    def edge_mean_curvature(self):
        """Return the half-scaled seam integral only when edge coverage is complete."""
        if not self.edge_contribution_complete:
            return None
        values = []
        for record in self.edges:
            value = record.signed_mean_curvature
            if value is None and record.signed_contribution is not None:
                value = 0.5 * float(record.signed_contribution)
            if value is None or not np.isfinite(value):
                return None
            values.append(float(value))
        return float(sum(values))

    @property
    def smooth_patch_contribution_complete(self):
        """Whether every retained smooth patch has a finite integral ``H dA``."""
        return all(
            record.status == "included"
            and record.integrated_mean_curvature is not None
            and np.isfinite(record.integrated_mean_curvature)
            for record in self.surfaces
        )

    @property
    def integrated_mean_curvature(self):
        """Return the signed, unit-preserving AW curvature accounting.

        Smooth patches and seams are kept distinct.  A known smooth partial may
        be reported, but an unresolved seam cannot become zero by aggregation.
        Free-boundary and corner terms are not scalar mean-curvature measures
        under the ``H=(k1+k2)/2`` convention, so they remain explicitly
        inapplicable here.  Topology owns certification of the full total.
        """
        smooth_values = [
            float(record.integrated_mean_curvature)
            for record in self.selected_surfaces
            if record.integrated_mean_curvature is not None
            and np.isfinite(record.integrated_mean_curvature)
        ]
        smooth_partial = (
            float(sum(smooth_values))
            if smooth_values or not self.surfaces
            else None
        )
        smooth = smooth_partial if self.smooth_patch_contribution_complete else None

        edge_values = [
            float(record.signed_contribution)
            for record in self.selected_edges
            if record.signed_contribution is not None
            and np.isfinite(record.signed_contribution)
        ]
        raw_seam = (
            float(sum(edge_values))
            if edge_values or (not self.edges and not self.missing_edge_ids)
            else None
        )
        seam = self.edge_mean_curvature
        partial = (
            float(smooth_partial + seam)
            if smooth_partial is not None and seam is not None
            else None
        )
        return {
            # Names shared with the topology-aware open-interface kernel.
            "smooth_H": smooth,
            "internal_edge_H": seam,
            "boundary_H": None,
            "corner_H": None,
            "partial_H": partial,
            "total_H": None,
            "smooth_patch": smooth,
            "smooth_patch_partial": smooth_partial,
            "signed_seam": seam,
            "signed_seam_raw_beta_ds": raw_seam,
            "boundary": None,
            "corner": None,
            "partial": partial,
            "total": None,
            "certified": False,
            "smooth_patch_complete": self.smooth_patch_contribution_complete,
            "signed_seam_complete": self.edge_contribution_complete,
            "boundary_status": "not_applicable_scalar_H",
            "corner_status": "not_applicable_scalar_H",
            "unresolved": self.unresolved_reasons,
            "unsupported": self.exclusion_diagnostics,
            "orientation": "A-to-B",
            "units": "angstrom",
            "convention": "H=(k1+k2)/2; C_e=0.5*integral(beta ds)",
            "provenance": "AW analytic smooth-patch records and validated A-to-B analytic edge records",
        }

    @property
    def total_curvature_certified(self):
        """Total H certification remains owned by the topology accounting pass."""
        return False

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
        reasons = Counter(
            record.reason or record.status
            for record in [*self.surfaces, *self.edges]
            if record.status != "included"
        )
        for edge_id in self.missing_edge_ids:
            reasons[f"missing interface edge record {edge_id}"] += 1
        return dict(reasons)

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
            if len(record.surface_ids) != 2:
                counts["incomplete_bicolor_surface_incidence"] += 1
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
        if self.missing_edge_ids:
            counts["missing_interface_edge_record"] = len(self.missing_edge_ids)
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
            f"Edge contribution complete: {self.edge_contribution_complete}; "
            f"total curvature certified: {self.total_curvature_certified}\n"
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
    complete surface geometry while preserving its support flags; G2b edge
    contributions are support-gated in both modes. ``supported`` additionally
    restricts the selected surfaces to the supported Apollonius subcomplex.
    ``alpha`` is intentionally unavailable until the required pair-surface
    alpha births are exact.
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
        tri_key = tuple(sorted(ids))
        dual_triangle = complex_.simplices[2].get(tri_key)
        features = [ref for ref in (dual_triangle.primal_features if dual_triangle else ())
                    if ref.kind == "edge" and ref.feature_id == int(edge_id)]
        feature_complete = any(feature.complete and feature.bounded for feature in features)
        cells_complete = any(feature.generator_cells_complete for feature in features)
        supported = bool(dual_triangle and complex_.is_supported(2, tri_key) and features)
        if len(bicolor_ids) != 2:
            result.edges.append(_edge_excluded_record(
                network, edge_id, ids, tri_key, bicolor_ids, supported,
                feature_complete, cells_complete,
                "edge does not retain exactly two bicolor incident surfaces "
                f"(found {len(bicolor_ids)})",
            ))
            continue
        both_surfaces_included = all(sid in selected_surface_ids for sid in bicolor_ids)
        # G2b only assigns C_e to a supported, regular analytic seam.  The
        # raw mode still reports candidates, but cannot promote unsupported
        # incidence to a curvature contribution.
        usable_by_mode = supported and feature_complete and cells_complete and both_surfaces_included
        if not usable_by_mode:
            reasons = []
            if not bicolor_ids:
                reasons.append("primal edge does not reference two bicolor AW surfaces")
            if not supported:
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
            converged = integral["status"] == "converged"
            status = "included" if converged else "unresolved"
            reason = "" if status == "included" else "edge quadrature did not converge by maximum order"
            result.edges.append(InterfaceEdgeRecord(
                int(edge_id), ids, _atom_numbers(network, ids), tri_key, bicolor_ids,
                type(resolved.geometry).__name__, supported, feature_complete,
                cells_complete, integral["length"],
                integral["signed"] if converged else None,
                integral["unsigned"] if converged else None,
                integral["positive"] if converged else None,
                integral["negative"] if converged else None,
                integral["beta_min"], integral["beta_max"], integral["beta_mean"],
                "Gauss-Legendre quadrature-node diagnostic",
                integral["status"], integral["error"], integral["order"],
                status, reason,
                0.5 * integral["signed"] if converged else None,
                0.5 * integral["unsigned"] if converged else None,
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


def _validated_edge_geometry(network, resolved):
    """Use a validated nonsingular parameterization when a branch is ill-conditioned."""
    geometry = resolved.geometry
    # A clearance-parameterized branch is exact but becomes poorly conditioned
    # near a conic turning point.  Reparameterize only after independently
    # matching the stored endpoints and samples to a nonsingular conic; if that
    # validation fails, retain the original branch and its fail-closed status.
    branch = getattr(geometry, "source", None)
    if isinstance(branch, AdditivelyWeightedTrisectorBranch):
        edge_row = network.edges.loc[resolved.edge_index]
        generator_ids = tuple(int(value) for value in edge_row["balls"])
        vertex_ids = tuple(int(value) for value in edge_row["verts"])
        try:
            conic_match = match_aw_trisector_conic_to_samples(
                [network.balls.loc[index, "loc"] for index in generator_ids],
                [network.balls.loc[index, "rad"] for index in generator_ids],
                vertex_ids,
                [network.verts.loc[index, "loc"] for index in vertex_ids],
                edge_row["points"],
                tolerance=1e-6,
            )
            if conic_match.status == "matched_nonsingular":
                geometry = conic_match.geometry
        except (TypeError, ValueError, np.linalg.LinAlgError):
            pass
    return geometry


def _integrate_edge_beta(network, resolved, surface_ids, group_a, group_b,
                         order, max_order, rtol, atol):
    geometry = _validated_edge_geometry(network, resolved)
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
        values = (signed, absolute, positive, negative, length, *beta_values)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("AW interface edge curvature integral is non-finite.")
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
        edge_id=int(edge_id),
        generator_ids=ids,
        atom_numbers=_atom_numbers(network, ids),
        dual_triangle=tri_key,
        surface_ids=tuple(surface_ids),
        geometry_type=None,
        supported=supported,
        feature_complete=feature_complete,
        generator_cells_complete=cells_complete,
        length=None,
        signed_contribution=None,
        unsigned_contribution=None,
        positive_contribution=None,
        negative_contribution=None,
        beta_min=None,
        beta_max=None,
        beta_mean=None,
        beta_statistics_method=None,
        integration_status="unresolved",
        estimated_error=None,
        quadrature_order=None,
        status="excluded",
        reason=reason,
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
