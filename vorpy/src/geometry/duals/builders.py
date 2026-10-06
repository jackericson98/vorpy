"""Scheme dispatch for duals derived from existing solved network incidence."""

from __future__ import annotations

from vorpy.src.analyze.apollonius import ApolloniusComplex

from .model import DualComplex, DualSimplex


class _NetworkIncidenceDualBuilder:
    scheme = None

    def __init__(self, network):
        self.network = network
        actual = str((getattr(network, "settings", None) or {}).get("net_type", "aw"))
        if actual != self.scheme:
            raise ValueError(f"{type(self).__name__} requires net_type={self.scheme!r}; got {actual!r}.")

    def build(self):
        # This extractor only reads the network. Its historical name is kept
        # for compatibility with experimental AW-alpha clients.
        incidence = ApolloniusComplex(self.network)
        rows, feature_to_simplices = [], {}
        for dimension, table in incidence.simplices.items():
            for ids, source in table.items():
                features = tuple(source.primal_features)
                supported = incidence.is_supported(dimension, ids)
                radii = tuple(float(value) for value in source.generator_radii)
                rows.append(DualSimplex(
                    dimension=dimension,
                    generator_ids=tuple(ids),
                    primal_features=features,
                    bounded=bool(source.bounded),
                    complete=bool(source.complete),
                    supported=bool(supported),
                    status="supported" if supported else "incomplete_or_unbounded",
                    generator_coordinates=tuple(tuple(float(x) for x in point)
                                                for point in source.generator_coordinates),
                    generator_radii=radii,
                    generator_weights=tuple(radius * radius for radius in radii) if self.scheme == "pow" else None,
                    metadata={
                        "geometry_kind": {"prm": "ordinary_voronoi", "pow": "power_laguerre", "aw": "additively_weighted_apollonius"}[self.scheme],
                        "dual_kind": {"prm": "delaunay", "pow": "regular_weighted_delaunay", "aw": "apollonius_abstract_incidence"}[self.scheme],
                        "physical_realization": "primal network feature; dual center geometry is abstract",
                    },
                ))
                simplex_id = f"{dimension}:" + ",".join(map(str, ids))
                for feature in features:
                    feature_to_simplices.setdefault((feature.kind, feature.feature_id), []).append(simplex_id)
        return DualComplex(
            self.network, self.scheme, rows, feature_to_simplices,
            metadata={"construction": "read-only extraction from solved VorPy incidence",
                      "input_radii": "exact Network.balls.rad values"},
            unresolved_features=incidence.unsupported_features,
        )


class PrimitiveDualBuilder(_NetworkIncidenceDualBuilder):
    """Ordinary Delaunay incidence from a primitive Voronoi network."""
    scheme = "prm"


class PowerDualBuilder(_NetworkIncidenceDualBuilder):
    """Regular-triangulation incidence from a Power/Laguerre network."""
    scheme = "pow"


class AWDualBuilder(_NetworkIncidenceDualBuilder):
    """Abstract Apollonius incidence from an existing AW network."""
    scheme = "aw"


def build_dual(network):
    """Build the read-only dual adapter matching a network's solved scheme."""
    scheme = str((getattr(network, "settings", None) or {}).get("net_type", "aw"))
    builders = {"prm": PrimitiveDualBuilder, "pow": PowerDualBuilder, "aw": AWDualBuilder}
    try:
        builder = builders[scheme]
    except KeyError as error:
        raise ValueError(f"Unsupported network scheme {scheme!r}; expected prm, pow, or aw.") from error
    return builder(network).build()
