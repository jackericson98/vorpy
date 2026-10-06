"""Common data model for scheme-specific abstract dual incidence complexes.

The records describe generator tuples and references to the solved Voronoi
network. Their optional center-connected geometry is only a visualization of
the abstract dual; physical geometry remains in the referenced network.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


@dataclass(frozen=True)
class DualSimplex:
    dimension: int
    generator_ids: tuple[int, ...]
    primal_features: tuple
    bounded: bool
    complete: bool
    supported: bool
    status: str
    generator_coordinates: tuple[tuple[float, float, float], ...]
    generator_radii: tuple[float, ...]
    generator_weights: tuple[float, ...] | None = None
    metadata: dict = field(default_factory=dict)

    @property
    def generator_tuple(self):
        return self.generator_ids

    @property
    def simplex_id(self):
        return f"{self.dimension}:" + ",".join(map(str, self.generator_ids))

    @property
    def primal_feature_ids(self):
        return tuple(feature.key for feature in self.primal_features)

    @property
    def visualization_geometry(self):
        """Generator centers; not the physical realization for curved AW."""
        return self.generator_coordinates


@dataclass(frozen=True)
class DualIncidenceAudit:
    scheme: str
    generator_count: int
    simplex_counts: dict[int, int]
    feature_counts: dict[str, int]
    feature_status_counts: dict[str, dict[str, int]]
    unmapped_feature_count: int
    reverse_mapping_errors: tuple[str, ...]

    @property
    def valid(self):
        return self.unmapped_feature_count == 0 and not self.reverse_mapping_errors


class DualComplex:
    """A read-only scheme-specific dual view over one solved network."""

    def __init__(self, network, scheme, simplices, feature_to_simplices,
                 metadata=None, unresolved_features=()):
        self.network = network
        self.scheme = str(scheme)
        self.generators = tuple(
            {
                "generator_id": int(index),
                "coordinates": tuple(float(x) for x in row["loc"]),
                "radius": float(row["rad"]),
                "weight": float(row["rad"]) ** 2 if self.scheme == "pow" else None,
            }
            for index, row in network.balls.iterrows()
        )
        self.simplices = {dimension: {} for dimension in range(4)}
        for simplex in simplices:
            self.simplices[simplex.dimension][simplex.generator_ids] = simplex
        self.feature_to_simplices = dict(feature_to_simplices)
        self.unresolved_features = tuple(unresolved_features)
        self.metadata = dict(metadata or {})
        self.metadata.update({
            "scheme": self.scheme,
            "generator_id_scope": "network-local Network.balls index",
            "physical_geometry_source": "solved VorPy network",
            "visualization_geometry": "straight generator-center realization of abstract dual",
        })

    def simplex(self, dimension, generator_ids):
        return self.simplices[int(dimension)].get(tuple(sorted(map(int, generator_ids))))

    @property
    def simplex_counts(self):
        return {dimension: len(rows) for dimension, rows in self.simplices.items()}

    def audit(self):
        feature_counts = Counter()
        feature_status = {}
        seen = set()
        errors = []
        feature_tables = ((0, "balls", "cell"), (1, "surfs", "surf"),
                          (2, "edges", "edge"), (3, "verts", "vert"))
        for dimension, table, kind in feature_tables:
            records = getattr(self.network, table, None)
            if records is None:
                continue
            for feature_id in records.index:
                key = (kind, int(feature_id))
                mapped = self.feature_to_simplices.get(key, ())
                if mapped:
                    seen.add(key)
                else:
                    errors.append(f"{kind}:{feature_id} has no dual simplex mapping")
                feature_counts[kind] += 1
                row = records.loc[feature_id]
                status = "unresolved"
                ids = (int(feature_id),) if dimension == 0 else tuple(
                    int(value) for value in row.get("balls", ())
                )
                simplex = self.simplex(dimension, ids) if len(ids) == dimension + 1 else None
                if simplex is not None:
                    ref = next((f for f in simplex.primal_features if f.kind == kind and f.feature_id == int(feature_id)), None)
                    if ref is not None:
                        if not ref.bounded:
                            status = "unbounded"
                        elif not ref.complete or not ref.generator_cells_complete:
                            status = "incomplete"
                        elif not simplex.supported:
                            status = "unresolved_incidence"
                        else:
                            status = "bounded_complete_supported"
                feature_status.setdefault(kind, Counter())[status] += 1
                for simplex_id in mapped:
                    simplex = next((s for s in self.simplices[dimension].values()
                                    if s.simplex_id == simplex_id), None)
                    if simplex is None or not any(f.kind == kind and f.feature_id == int(feature_id)
                                                  for f in simplex.primal_features):
                        errors.append(f"Reverse mapping mismatch for {kind}:{feature_id} -> {simplex_id}")
        expected = sum(feature_counts.values())
        unmapped = expected - len(seen)
        for feature in self.unresolved_features:
            errors.append(
                f"{feature.kind}:{feature.feature_id} unresolved: {feature.reason}"
            )
        return DualIncidenceAudit(
            self.scheme, len(self.generators), self.simplex_counts, dict(feature_counts),
            {kind: dict(counts) for kind, counts in feature_status.items()}, unmapped,
            tuple(errors),
        )
