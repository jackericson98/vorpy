"""Experimental AW restricted-cell nerve over a solved VorPy network."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from math import isfinite

import numpy as np

from vorpy.src.analyze.apollonius import ApolloniusComplex

from .births import edge_birth, generator_birth, surface_birth, vertex_birth


@dataclass
class AWAlphaRecord:
    dimension: int
    generator_tuple: tuple
    geometric_birth: float | None
    filtration_birth: float | None = None
    birth_source: str = "unresolved"
    closure_adjustment: float | None = None
    supported: bool = False
    residual: float | None = None
    notes: str = ""
    diagnostics: dict = field(default_factory=dict)

    @property
    def simplex_id(self):
        return f"{self.dimension}:" + ",".join(map(str, self.generator_tuple))


class AWAlphaSubcomplex:
    """Birth-filtered subcomplex, retaining blocked-simplex diagnostics."""

    def __init__(
        self, source, alpha, records, blocked, tolerance, truncate_negative=False
    ):
        self.source = source
        self.alpha = float(alpha)
        self.records = records
        self.blocked = list(blocked)
        self.tolerance = float(tolerance)
        self.truncate_negative = bool(truncate_negative)

    @property
    def counts(self):
        return {
            dimension: sum(row.dimension == dimension for row in self.records)
            for dimension in range(4)
        }

    @property
    def euler_characteristic(self):
        counts = self.counts
        return counts[0] - counts[1] + counts[2] - counts[3]

    def validate(self):
        keys = {row.generator_tuple for row in self.records}
        errors = []
        for row in self.records:
            for face in combinations(row.generator_tuple, row.dimension):
                if tuple(face) not in keys:
                    errors.append(
                        f"{row.simplex_id} missing face {tuple(face)} at alpha={self.alpha:g}"
                    )
        return errors


class AWAlphaFiltration:
    """Birth table over the existing AW incidence; no tessellation is built."""

    def __init__(self, network, name=None, tolerance=1e-6):
        self.network = network
        self.name = str(name or getattr(network, "group_name", None) or "Network")
        self.tolerance = float(tolerance)
        if not isfinite(self.tolerance) or self.tolerance <= 0:
            raise ValueError("AW alpha tolerance must be finite and positive.")
        if (getattr(network, "settings", None) or {}).get("net_type", "aw") != "aw":
            raise ValueError(
                "Experimental AW alpha filtration requires a solved AW network."
            )
        self.incidence = ApolloniusComplex(network, name=self.name)
        self.records = {dimension: {} for dimension in range(4)}
        self.issues = []
        self._built = False
        self._max_dimension_built = -1
        # Births share boundary features heavily.  Keep these caches scoped
        # to this filtration so no state escapes the solved network/query.
        self._edge_birth_cache = {}
        self._vertex_birth_cache = {}
        self._geometry_cache = {}

    def _record(
        self,
        dimension,
        ids,
        geometric_birth,
        *,
        supported,
        source,
        residual=None,
        notes="",
        diagnostics=None,
    ):
        self.records[dimension][tuple(ids)] = AWAlphaRecord(
            dimension=dimension,
            generator_tuple=tuple(ids),
            geometric_birth=None if geometric_birth is None else float(geometric_birth),
            birth_source=source,
            supported=bool(supported),
            residual=None if residual is None else float(residual),
            notes=str(notes),
            diagnostics=diagnostics or {},
        )

    def calculate_births(self, max_dimension=3):
        max_dimension = int(max_dimension)
        if max_dimension < 0 or max_dimension > 3:
            raise ValueError("max_dimension must be between 0 and 3.")
        if self._built or max_dimension <= self._max_dimension_built:
            return self
        if max_dimension >= 0 and not self.records[0]:
            for key, simplex in self.incidence.simplices[0].items():
                radius = float(self.network.balls.loc[key[0], "rad"])
                self._record(
                    0,
                    key,
                    generator_birth(radius),
                    supported=True,
                    source="generator_radius_convention",
                    residual=0.0,
                    notes="Formal union-of-balls birth alpha=-r; values remain negative when r>0.",
                )

        if max_dimension >= 1 and not self.records[1]:
            simplices = list(self.incidence.simplices[1].items())
            total = len(simplices)
            progress_step = max(1, total // 20)
            for index, (key, simplex) in enumerate(simplices, start=1):
                if index == 1 or index == total or index % progress_step == 0:
                    system = getattr(self.network, 'sys', None)
                    updater = getattr(system, 'update_progress', None)
                    if updater is not None:
                        updater(
                            process=f'Interface analysis | Alpha selection {index}/{total}',
                            progress=20.0 + 15.0 * index / max(total, 1),
                            network=self.network,
                        )
                if not self.incidence.is_supported(1, key):
                    self._record(
                        1,
                        key,
                        None,
                        supported=False,
                        source="unresolved",
                        notes="No complete bounded AW surface patch in the solved network.",
                    )
                    continue
                try:
                    value, diagnostic = surface_birth(
                        self.network, simplex, self.tolerance,
                        edge_birth_cache=self._edge_birth_cache,
                        vertex_birth_cache=self._vertex_birth_cache,
                        geometry_cache=self._geometry_cache,
                    )
                except (
                    TypeError,
                    ValueError,
                    KeyError,
                    IndexError,
                    np.linalg.LinAlgError,
                ) as error:
                    value, diagnostic = (
                        None,
                        {"birth_source": "unresolved", "reason": str(error)},
                    )
                self._record(
                    1,
                    key,
                    value,
                    supported=value is not None,
                    source=diagnostic.get("birth_source", "unresolved"),
                    residual=diagnostic.get("residual_A"),
                    notes=diagnostic.get("reason", diagnostic.get("notes", "")),
                    diagnostics=diagnostic,
                )

        if max_dimension >= 2 and not self.records[2]:
            for key, simplex in self.incidence.simplices[2].items():
                if not self.incidence.is_supported(2, key):
                    self._record(
                        2,
                        key,
                        None,
                        supported=False,
                        source="unresolved",
                        notes="AW edge is not a complete bounded incidence-closed primal feature.",
                    )
                    continue
                minima, diagnostics = [], []
                for feature in simplex.primal_features:
                    try:
                        value, diagnostic = edge_birth(
                            self.network, feature.feature_id, self.tolerance
                        )
                        minima.append(value)
                        diagnostics.append({"feature": feature.key, **diagnostic})
                    except (
                        TypeError,
                        ValueError,
                        KeyError,
                        IndexError,
                        np.linalg.LinAlgError,
                    ) as error:
                        diagnostics.append(
                            {
                                "feature": feature.key,
                                "birth_source": "unresolved",
                                "reason": str(error),
                            }
                        )
                if minima and len(minima) == len(simplex.primal_features):
                    value = min(minima)
                    chosen = min(
                        diagnostics,
                        key=lambda item: min(
                            item.get("candidate_clearances_A", (float("inf"),))
                        ),
                    )
                    self._record(
                        2,
                        key,
                        value,
                        supported=True,
                        source=chosen.get("birth_source", "other"),
                        residual=chosen.get("residual_A"),
                        notes="Minimum over all actual bounded AW edge components.",
                        diagnostics={"components": diagnostics},
                    )
                else:
                    self._record(
                        2,
                        key,
                        None,
                        supported=False,
                        source="unresolved",
                        notes="At least one actual AW edge component has no certified minimum.",
                        diagnostics={"components": diagnostics},
                    )

        if max_dimension >= 3 and not self.records[3]:
            for key, simplex in self.incidence.simplices[3].items():
                if not self.incidence.is_supported(3, key):
                    self._record(
                        3,
                        key,
                        None,
                        supported=False,
                        source="unresolved",
                        notes="AW vertex lacks complete bounded incidence-closed primal support.",
                    )
                    continue
                minima, diagnostics = [], []
                for feature in simplex.primal_features:
                    try:
                        value, diagnostic = vertex_birth(
                            self.network, feature.feature_id, key, self.tolerance
                        )
                        minima.append(value)
                        diagnostics.append({"feature": feature.key, **diagnostic})
                    except (TypeError, ValueError, KeyError, IndexError) as error:
                        diagnostics.append(
                            {
                                "feature": feature.key,
                                "birth_source": "unresolved",
                                "reason": str(error),
                            }
                        )
                if minima and len(minima) == len(simplex.primal_features):
                    index = int(np.argmin(minima))
                    self._record(
                        3,
                        key,
                        minima[index],
                        supported=True,
                        source="apollonius_vertex",
                        residual=diagnostics[index]["residual_A"],
                        notes="Common AW clearance verified independently at all four generators.",
                        diagnostics={"components": diagnostics},
                    )
                else:
                    self._record(
                        3,
                        key,
                        None,
                        supported=False,
                        source="unresolved",
                        notes="At least one corresponding AW vertex failed four-distance verification.",
                        diagnostics={"components": diagnostics},
                    )
        self._close_filtration()
        self._max_dimension_built = max_dimension
        self._built = max_dimension == 3
        return self

    def _close_filtration(self):
        """Apply the nerve's face-maximum closure, preserving raw births."""
        self.issues.clear()
        for dimension in range(4):
            for record in self.records[dimension].values():
                if record.geometric_birth is None or not record.supported:
                    record.filtration_birth = None
                    record.closure_adjustment = None
                    continue
                face_records = (
                    [
                        self.records[dimension - 1].get(tuple(face))
                        for face in combinations(record.generator_tuple, dimension)
                    ]
                    if dimension
                    else []
                )
                if any(
                    face is None
                    or face.geometric_birth is None
                    or face.filtration_birth is None
                    or not face.supported
                    for face in face_records
                ):
                    record.filtration_birth = None
                    record.closure_adjustment = None
                    record.notes = (
                        record.notes
                        + " Filtration birth unresolved because a proper face birth is unresolved."
                    ).strip()
                    self.issues.append(
                        {
                            "kind": "unresolved_face",
                            "simplex": record.simplex_id,
                            "message": "One or more proper-face births are unresolved.",
                        }
                    )
                    continue
                face_max = max(
                    (face.filtration_birth for face in face_records),
                    default=record.geometric_birth,
                )
                record.filtration_birth = max(record.geometric_birth, face_max)
                record.closure_adjustment = (
                    record.filtration_birth - record.geometric_birth
                )
                if record.closure_adjustment > self.tolerance:
                    self.issues.append(
                        {
                            "kind": "closure_adjustment",
                            "simplex": record.simplex_id,
                            "adjustment_A": record.closure_adjustment,
                            "message": "Geometric birth violates face monotonicity; closure was recorded, not hidden.",
                        }
                    )

    def validate_filtration(self, tolerance=None):
        self.calculate_births()
        tolerance = self.tolerance if tolerance is None else float(tolerance)
        issues = []
        for dimension in range(1, 4):
            for record in self.records[dimension].values():
                if record.geometric_birth is None:
                    continue
                for face_tuple in combinations(record.generator_tuple, dimension):
                    face = self.records[dimension - 1].get(tuple(face_tuple))
                    if face is None:
                        issues.append(
                            {
                                "kind": "missing_face",
                                "simplex": record.simplex_id,
                                "face": tuple(face_tuple),
                            }
                        )
                    elif face.geometric_birth is None:
                        issues.append(
                            {
                                "kind": "unresolved_face",
                                "simplex": record.simplex_id,
                                "face": face.simplex_id,
                            }
                        )
                    elif face.geometric_birth > record.geometric_birth + tolerance:
                        issues.append(
                            {
                                "kind": "monotonicity_violation",
                                "simplex": record.simplex_id,
                                "face": face.simplex_id,
                                "face_birth_A": face.geometric_birth,
                                "coface_birth_A": record.geometric_birth,
                            }
                        )
        return issues

    def alpha_complex(self, alpha, *, truncate_negative=False, max_dimension=3):
        self.calculate_births(max_dimension=max_dimension)
        alpha = float(alpha)
        if not isfinite(alpha):
            raise ValueError("AW alpha must be finite.")
        if truncate_negative and alpha < 0.0:
            raise ValueError(
                "A filtration truncated to alpha>=0 cannot be queried below zero."
            )
        query_alpha = alpha
        included, blocked = [], []
        for dimension in range(4):
            for record in self.records[dimension].values():
                if (
                    record.filtration_birth is None
                    or not record.supported
                    or record.filtration_birth > query_alpha + self.tolerance
                ):
                    continue
                face_keys = (
                    [
                        tuple(face)
                        for face in combinations(record.generator_tuple, dimension)
                    ]
                    if dimension
                    else []
                )
                if any(
                    not any(row.generator_tuple == face for row in included)
                    for face in face_keys
                ):
                    blocked.append(
                        {
                            "simplex": record.simplex_id,
                            "reason": "face missing, unresolved, or above alpha",
                        }
                    )
                    continue
                included.append(record)
        return AWAlphaSubcomplex(
            self, alpha, included, blocked, self.tolerance, truncate_negative
        )

    def bicolor_surface_records(self, alpha, group_a, group_b, *, max_dimension=3):
        """Map accepted alpha edges to their existing network AW surface rows."""
        subcomplex = self.alpha_complex(alpha, max_dimension=max_dimension)
        group_a, group_b = set(map(int, group_a)), set(map(int, group_b))
        accepted = {
            row.generator_tuple: row
            for row in subcomplex.records
            if row.dimension == 1
            and (
                (
                    row.generator_tuple[0] in group_a
                    and row.generator_tuple[1] in group_b
                )
                or (
                    row.generator_tuple[1] in group_a
                    and row.generator_tuple[0] in group_b
                )
            )
        }
        rows, edge_sets = [], {}
        for pair, birth in accepted.items():
            simplex = self.incidence.simplices[1][pair]
            seen_surfaces = set()
            for feature in simplex.primal_features:
                if feature.feature_id in seen_surfaces:
                    continue
                seen_surfaces.add(feature.feature_id)
                surf = self.network.surfs.loc[feature.feature_id]
                edge_ids = tuple(sorted(int(value) for value in surf.get("edges", ())))
                edge_sets[feature.feature_id] = set(edge_ids)
                atom_rows = [self.network.balls.loc[index] for index in pair]
                atom_i = _network_atom_metadata(self.network, pair[0], atom_rows[0])
                atom_j = _network_atom_metadata(self.network, pair[1], atom_rows[1])
                atom_i["group"] = "A" if pair[0] in group_a else "B"
                atom_j["group"] = "A" if pair[1] in group_a else "B"
                rows.append(
                    {
                        "generator_tuple": pair,
                        "surface_id": feature.feature_id,
                        **{f"atom_i_{key}": value for key, value in atom_i.items()},
                        **{f"atom_j_{key}": value for key, value in atom_j.items()},
                        "geometric_birth_A": birth.geometric_birth,
                        "filtration_birth_A": birth.filtration_birth,
                        "birth_source": birth.birth_source,
                        "surface_area_A2": _optional_float(surf.get("sa")),
                        "surface_integrated_mean_curvature_A": _optional_float(
                            surf.get("int_mean_curv")
                        ),
                        "surface_mean_curvature": _optional_float(
                            surf.get("mean_curv")
                        ),
                        "feature_complete": bool(feature.complete and feature.bounded),
                        "supported": birth.supported,
                        "confidence_status": birth.birth_source,
                        "primal_feature": feature.key,
                    }
                )
        components = _components(edge_sets)
        return rows, {
            "pair_count": len(accepted),
            "mapped_surface_count": len(rows),
            "interface_atom_count": len({i for pair in accepted for i in pair}),
            "total_surface_area_A2": sum(row["surface_area_A2"] or 0.0 for row in rows),
            "connected_component_count": components,
            "unsupported_birth_count": sum(
                record.birth_source == "unresolved"
                for record in self.records[1].values()
            ),
        }

    def all_records(self):
        self.calculate_births()
        return [
            record
            for dimension in range(4)
            for record in self.records[dimension].values()
        ]

    def summary(self):
        self.calculate_births()
        resolved = sum(
            record.filtration_birth is not None for record in self.all_records()
        )
        unresolved = sum(
            record.filtration_birth is None for record in self.all_records()
        )
        violations = [
            issue
            for issue in self.validate_filtration()
            if issue["kind"] == "monotonicity_violation"
        ]
        adjustments = [
            record
            for record in self.all_records()
            if record.closure_adjustment is not None
            and record.closure_adjustment > self.tolerance
        ]
        return (
            "Experimental additively weighted alpha filtration\n"
            f"System: {self.name}\n"
            "Distance: d_i(x)=||x-p_i||-r_i; alpha is an additive radial expansion in A.\n"
            "alpha=0 uses the supplied network radii as-is; if probe-expanded radii were supplied, they remain probe-expanded.\n"
            "No Power-alpha formula, GUDHI alpha complex, or Cazals M=5 condition is used.\n"
            f"Resolved filtration births: {resolved}; unresolved: {unresolved}\n"
            f"Filtration-order violations: {len(violations)}; closure adjustments > tolerance: {len(adjustments)}\n"
            "Pair surfaces use the analytic pair-contact minimum only when feasible in the actual AW patch; otherwise the recorded boundary is required.\n"
        )


def _optional_float(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _atom_metadata(row, atom_index):
    chain = _first_nonempty(row, ("chain_name", "chain"))
    residue_number = _first_nonempty(row, ("auth_seq_id", "res_seq", "residue_number"))
    insertion_code = _first_nonempty(row, ("pdb_ins_code", "insertion_code", "icode"))
    residue_name = _first_nonempty(row, ("auth_comp_id", "res_name", "residue_name"))
    atom_name = _first_nonempty(row, ("auth_atom_id", "name", "atom_name"))
    stable_id = _first_nonempty(row, ("stable_id",)) or (
        f"{chain}|{residue_number}{insertion_code}|{residue_name}|{atom_name}"
        if any((chain, residue_number, residue_name, atom_name))
        else ""
    )
    return {
        "atom_index": int(row.get("num", atom_index)),
        "generator_id": int(atom_index),
        "stable_id": stable_id,
        "metadata_status": "available" if stable_id else "unavailable",
        "chain": chain,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
        "residue_name": residue_name,
        "atom_name": atom_name,
        "radius_A": _optional_float(row.get("rad")),
        "group": str(row.get("group", "")),
    }


def _first_nonempty(row, keys):
    for key in keys:
        value = row.get(key, "")
        if value is not None and str(value).strip() not in {"", "nan", "None"}:
            return str(value).strip()
    return ""


def _network_atom_metadata(network, generator_id, geometry_row):
    """Use the parent molecular atom table when the network geometry table is compact."""
    source = geometry_row
    system = getattr(network, "sys", None)
    system_balls = getattr(system, "balls", None)
    if system_balls is not None and hasattr(system_balls, "loc"):
        source_id = geometry_row.get("num", generator_id)
        try:
            candidate = system_balls.loc[int(source_id)]
            if np.allclose(
                np.asarray(candidate["loc"], dtype=float),
                np.asarray(geometry_row["loc"], dtype=float),
                rtol=0.0,
                atol=1e-8,
            ):
                source = candidate
        except (KeyError, TypeError, ValueError, IndexError):
            pass
    metadata = _atom_metadata(source, generator_id)
    metadata["atom_index"] = int(
        source.get("num", geometry_row.get("num", generator_id))
    )
    metadata["generator_id"] = int(generator_id)
    metadata["radius_A"] = _optional_float(geometry_row.get("rad"))
    return metadata


def _components(edge_sets):
    if not edge_sets:
        return 0
    parents = {key: key for key in edge_sets}

    def find(key):
        while parents[key] != key:
            parents[key] = parents[parents[key]]
            key = parents[key]
        return key

    owners = {}
    for surface_id, edge_ids in edge_sets.items():
        for edge_id in edge_ids:
            if edge_id in owners:
                parents[find(surface_id)] = find(owners[edge_id])
            else:
                owners[edge_id] = surface_id
    return len({find(key) for key in parents})
