"""Common query API for scheme-specific alpha filtrations.

Primitive and Power births are read from GUDHI simplex trees. AW births are
delegated to VorPy's experimental restricted-cell filtration. Native values
and units are preserved; a shared API does not imply shared alpha semantics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite

import numpy as np


@dataclass(frozen=True)
class AlphaSimplex:
    dimension: int
    generator_tuple: tuple[int, ...]
    geometric_birth: float | None
    filtration_birth: float | None
    birth_source: str
    supported: bool
    residual: float | None = None
    notes: str = ""
    diagnostics: dict = field(default_factory=dict)

    @property
    def simplex_id(self):
        return f"{self.dimension}:" + ",".join(map(str, self.generator_tuple))


class AlphaSubcomplex:
    """The resolved simplices present at one scheme-native alpha value."""

    def __init__(self, scheme, alpha, records, blocked=(), metadata=None):
        self.scheme = scheme
        self.alpha = float(alpha)
        self.records = list(records)
        self.blocked = list(blocked)
        self.metadata = dict(metadata or {})

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

    def bicolor_edges(self, group_a, group_b):
        """Return selected 1-simplices crossing the supplied generator groups."""
        a, b = set(map(int, group_a)), set(map(int, group_b))
        if a & b:
            raise ValueError("Bicolor groups must be disjoint.")
        return [
            row for row in self.records
            if row.dimension == 1
            and ((row.generator_tuple[0] in a and row.generator_tuple[1] in b)
                 or (row.generator_tuple[1] in a and row.generator_tuple[0] in b))
        ]


class AlphaFiltration:
    """Scheme-specific births behind a small common inspection API."""

    def __init__(self, network, scheme, records, metadata, *, tolerance=1e-12,
                 blocked=(), aw_source=None):
        self.network = network
        self.scheme = str(scheme)
        self.records = dict(records)
        self.metadata = dict(metadata)
        self.tolerance = float(tolerance)
        self.blocked = list(blocked)
        self._aw_source = aw_source

    def simplex_birth(self, generator_tuple):
        key = tuple(sorted(map(int, generator_tuple)))
        dimension = len(key) - 1
        if dimension < 0 or dimension > 3:
            return None
        return self.records.get(dimension, {}).get(key)

    def simplices_at(self, alpha, *, truncate_negative=False):
        alpha = float(alpha)
        if not isfinite(alpha):
            raise ValueError("alpha must be finite")
        if truncate_negative and alpha < 0:
            raise ValueError("A filtration truncated to alpha >= 0 cannot be queried below zero")
        rows, blocked = [], {record.simplex_id: record for record in self.blocked}
        for dimension in range(4):
            for record in self.records.get(dimension, {}).values():
                birth = record.filtration_birth
                if not record.supported or birth is None:
                    blocked[record.simplex_id] = record
                elif birth <= alpha + self.tolerance:
                    rows.append(record)
        return AlphaSubcomplex(
            self.scheme, alpha, rows, blocked.values(),
            metadata={**self.metadata, "alpha_query_native": alpha,
                      "negative_alpha_truncated": bool(truncate_negative)},
        )

    def interface_at(self, alpha, group_a, group_b, *, truncate_negative=False):
        """Return alpha-selected bicolor dual edges (surface mapping is Phase 5)."""
        subcomplex = self.simplices_at(alpha, truncate_negative=truncate_negative)
        return subcomplex.bicolor_edges(group_a, group_b)


def _records_from_simplex_tree(simplex_tree, id_by_index, source):
    records = {dimension: {} for dimension in range(4)}
    for simplex, birth in simplex_tree.get_simplices():
        if len(simplex) > 4:
            continue
        ids = tuple(sorted(id_by_index[int(index)] for index in simplex))
        dimension = len(ids) - 1
        value = float(birth)
        records[dimension][ids] = AlphaSimplex(
            dimension, ids, value, value, source, True, residual=0.0,
        )
    return records


def _generator_inputs(network):
    balls = network.balls
    ids, points, radii = [], [], []
    for index, row in balls.iterrows():
        ids.append(int(index))
        points.append(np.asarray(row["loc"], dtype=float))
        radii.append(float(row["rad"]))
    if len(set(ids)) != len(ids):
        raise ValueError("Network generator IDs must be unique.")
    if not points or any(point.shape != (3,) for point in points):
        raise ValueError("Network balls must provide 3D loc coordinates.")
    return ids, np.asarray(points, dtype=float), np.asarray(radii, dtype=float)


def _gudhi_module():
    try:
        import gudhi
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "GUDHI is required for Primitive and Power alpha filtrations; "
            "install a GUDHI version supporting AlphaComplex."
        ) from error
    if not hasattr(gudhi, "AlphaComplex"):
        raise RuntimeError("Installed GUDHI does not expose AlphaComplex.")
    return gudhi


def _build_primitive(network, tolerance):
    ids, points, _radii = _generator_inputs(network)
    gudhi = _gudhi_module()
    tree = gudhi.AlphaComplex(points=points.tolist(), precision="safe").create_simplex_tree()
    records = _records_from_simplex_tree(tree, ids, "gudhi_ordinary_alpha")
    return AlphaFiltration(
        network, "prm", records,
        {
            "filtration_kind": "ordinary alpha complex",
            "native_alpha": "squared radius",
            "native_alpha_units": "A^2",
            "physical_offset_A": None,
            "effective_radius_rule": "sqrt(alpha) for alpha >= 0",
            "simplex_tree_dimension": int(tree.dimension()),
            "generator_ids": ids,
        }, tolerance=tolerance,
    )


def _build_power(network, tolerance):
    ids, points, radii = _generator_inputs(network)
    # Reuse the frozen weighted GUDHI constructor: weights remain exactly R_i^2.
    from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_power_interface import (
        build_regular_complex,
    )

    tree, *_ = build_regular_complex(points, radii)
    records = _records_from_simplex_tree(tree, ids, "gudhi_weighted_alpha")
    weights = radii ** 2
    return AlphaFiltration(
        network, "pow", records,
        {
            "filtration_kind": "weighted Power alpha complex",
            "native_alpha": "power distance threshold",
            "native_alpha_units": "A^2",
            "power_sublevel": "||x-p_i||^2 - w_i <= alpha",
            "weight_rule": "w_i = exact Network.balls.rad_i squared",
            "weights_A2": weights.tolist(),
            "physical_offset_A": None,
            "effective_radius_rule": "sqrt(max(w_i + alpha, 0)) per generator",
            "simplex_tree_dimension": int(tree.dimension()),
            "generator_ids": ids,
        }, tolerance=tolerance,
    )


def _build_aw(network, tolerance, max_dimension):
    from vorpy.src.geometry.aw_alpha import AWAlphaFiltration

    source = AWAlphaFiltration(network, tolerance=max(float(tolerance), 1e-6))
    source.calculate_births(max_dimension=max_dimension)
    records = {dimension: {} for dimension in range(4)}
    blocked = []
    for dimension, dimension_records in source.records.items():
        for ids, row in dimension_records.items():
            record = AlphaSimplex(
                dimension, tuple(map(int, ids)), row.geometric_birth,
                row.filtration_birth, row.birth_source, row.supported,
                row.residual, row.notes, dict(row.diagnostics),
            )
            records[dimension][record.generator_tuple] = record
            if not row.supported or row.filtration_birth is None:
                blocked.append(record)
    return AlphaFiltration(
        network, "aw", records,
        {
            "filtration_kind": "experimental AW restricted-cell nerve",
            "experimental": True,
            "distance": "d_i(x)=||x-p_i||-r_i",
            "native_alpha": "additive radial offset",
            "native_alpha_units": "A",
            "physical_offset_A": "native alpha",
            "effective_radius_rule": "r_i + alpha",
            "higher_dimension_limitations": "unsupported geometric births remain unresolved",
            "birth_dimensions_calculated": max_dimension,
            "closure_issues": list(source.issues),
            "generator_ids": [int(index) for index in network.balls.index],
        }, tolerance=source.tolerance, blocked=blocked, aw_source=source,
    )


def build_alpha_filtration(network, *, tolerance=1e-12, max_dimension=3):
    """Build the alpha adapter matching an already solved network scheme."""
    max_dimension = int(max_dimension)
    if max_dimension < 0 or max_dimension > 3:
        raise ValueError("max_dimension must be between 0 and 3")
    scheme = str((getattr(network, "settings", None) or {}).get("net_type", "aw"))
    if scheme == "prm":
        result = _build_primitive(network, tolerance)
    elif scheme == "pow":
        result = _build_power(network, tolerance)
    elif scheme == "aw":
        return _build_aw(network, tolerance, max_dimension)
    else:
        raise ValueError(f"Unsupported network scheme {scheme!r}; expected prm, pow, or aw.")
    if max_dimension < 3:
        result.records = {dimension: rows for dimension, rows in result.records.items()
                          if dimension <= max_dimension}
    result.metadata["birth_dimensions_calculated"] = max_dimension
    return result
