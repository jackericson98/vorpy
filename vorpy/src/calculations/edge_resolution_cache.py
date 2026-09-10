"""Shared analytic AW edge-geometry cache for network curvature calculations."""

from time import perf_counter

from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge


_CACHE_ATTR = "_aw_edge_geometry_cache"
_META_ATTR = "_aw_edge_geometry_cache_meta"


def clear_aw_edge_geometry_cache(net):
    """Remove any cached analytic AW edge geometries from ``net``."""
    for name in (_CACHE_ATTR, _META_ATTR):
        if hasattr(net, name):
            delattr(net, name)


def get_aw_edge_geometry_cache(net, tolerance=1e-7):
    """Return one resolved analytic geometry object per physical network edge.

    The cache is stored directly on the Network instance. The first caller
    resolves all edges at a strict tolerance; later curvature stages reuse the
    same objects. A stricter later request automatically rebuilds the cache.

    Returns
    -------
    tuple(dict, bool, float)
        ``(cache, built_now, build_seconds)``.
    """
    requested = float(tolerance)
    if requested <= 0:
        raise ValueError("Edge-geometry cache tolerance must be positive.")

    # Use at least 1e-7 precision for the shared production cache so Edge-M,
    # Edge-G, and Vertex-G can normally share the same resolved geometries.
    cache_tolerance = min(requested, 1e-7)
    edges = getattr(net, "edges", None)

    if edges is None:
        raise ValueError("The network has no edge table.")

    meta = getattr(net, _META_ATTR, None)
    cache = getattr(net, _CACHE_ATTR, None)

    valid = (
        isinstance(cache, dict)
        and isinstance(meta, dict)
        and meta.get("edges_id") == id(edges)
        and meta.get("edge_count") == len(edges)
        and float(meta.get("tolerance", float("inf"))) <= cache_tolerance
        and len(cache) == len(edges)
    )

    if valid:
        return cache, False, 0.0

    start = perf_counter()
    cache = {
        int(edge_index): resolve_aw_network_edge(
            net,
            edge_index,
            tolerance=cache_tolerance,
        )
        for edge_index in edges.index
    }
    elapsed = perf_counter() - start

    setattr(net, _CACHE_ATTR, cache)
    setattr(net, _META_ATTR, {
        "edges_id": id(edges),
        "edge_count": len(edges),
        "tolerance": cache_tolerance,
    })

    return cache, True, elapsed
