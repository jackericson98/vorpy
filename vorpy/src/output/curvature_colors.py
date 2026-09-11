"""Shared curvature coloring for VorPy geometry exports.

Three visualization interpretations are supported for integrated curvature:

``boundary``
    Signed curvature of the boundary of a selected union of cells. This is
    used by shell_* exports and is the appropriate view for morphometric
    molecular-boundary interpretation.

``magnitude``
    Orientation-free component magnitude. This is used by generic surfs,
    edges, and verts exports, where a physical component has no unique cell
    orientation. Values occupy only the neutral-to-positive half of a
    diverging colormap, so arbitrary blue/red sign flips are avoided.

``cell``
    Curvature from the perspective of an individual AW cell. This is used by
    atom-cell exports.
"""

import numpy as np
import matplotlib as mpl


INTEGRATED_CURVATURE_SCHEMES = {"int_mean_curv", "int_gauss_curv"}

SCHEME_ALIASES = {
    "int_mean": "int_mean_curv",
    "integrated_mean": "int_mean_curv",
    "integrated mean curvature": "int_mean_curv",
    "int_gauss": "int_gauss_curv",
    "integrated_gauss": "int_gauss_curv",
    "integrated gaussian curvature": "int_gauss_curv",
}

SUPPORTED_COMPONENTS = {
    "int_mean_curv": ("surface", "edge"),
    "int_gauss_curv": ("surface", "edge", "vertex"),
}


def canonical_curvature_scheme(scheme):
    if scheme is None:
        return None
    name = str(scheme).strip().lower()
    name = SCHEME_ALIASES.get(name, name)
    return name if name in INTEGRATED_CURVATURE_SCHEMES else None


def canonical_color_mode(mode):
    name = "boundary" if mode is None else str(mode).strip().lower()
    aliases = {
        "shell": "boundary", "signed": "boundary", "boundary": "boundary",
        "mag": "magnitude", "abs": "magnitude", "absolute": "magnitude",
        "nondirectional": "magnitude", "non_directional": "magnitude",
        "magnitude": "magnitude",
        "atom": "cell", "cell": "cell", "cell_relative": "cell",
    }
    return aliases.get(name, "boundary")


def _cmap(name):
    name = "coolwarm" if name is None else str(name).strip()
    try:
        return mpl.colormaps.get_cmap(name)
    except Exception:
        try:
            return mpl.cm.get_cmap(name)
        except Exception:
            return mpl.colormaps.get_cmap("coolwarm")


def _target_set(target_cells):
    if target_cells is None:
        return None
    return {int(value) for value in target_cells}


def _row_balls(row):
    try:
        return tuple(int(value) for value in row["balls"])
    except (KeyError, TypeError, ValueError):
        return ()


def _mapping_value(mapping, cell_index):
    if not isinstance(mapping, dict):
        return None
    for key in (cell_index, str(cell_index)):
        if key in mapping:
            try:
                value = float(mapping[key])
            except (TypeError, ValueError):
                return None
            return value if np.isfinite(value) else None
    return None


def _mapping_values(mapping, target_cells=None):
    if not isinstance(mapping, dict):
        return []
    targets = _target_set(target_cells)
    values = []
    for key, value in mapping.items():
        try:
            cell = int(key)
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if targets is not None and cell not in targets:
            continue
        if np.isfinite(numeric):
            values.append(numeric)
    return values


def _sum_mapping(mapping, target_cells=None):
    values = _mapping_values(mapping, target_cells)
    return float(sum(values)) if values else 0.0


def _magnitude_mapping(mapping, target_cells=None):
    """Orientation-free representative magnitude for one physical component."""
    values = _mapping_values(mapping, target_cells)
    if not values:
        return 0.0
    return float(np.mean(np.abs(values)))


def _boundary_surface_mean_value(row, target_cells):
    """Signed smooth-surface M of the boundary of a cell union."""
    targets = _target_set(target_cells)
    balls = _row_balls(row)
    if targets is None or len(balls) != 2:
        return None

    inside = [ball for ball in balls if ball in targets]
    if len(inside) != 1:
        return 0.0  # internal or unrelated pairwise surface

    value = _mapping_value(row.get("int_mean_curv_by_ball", None), inside[0])
    if value is not None:
        return value

    try:
        value = float(row["int_mean_curv"])
    except (KeyError, TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _boundary_edge_mean_value(row, target_cells):
    """Signed singular edge-M for the boundary of a cell union.

    A regular AW edge has three generating cells. One selected cell produces
    a convex shell ridge. Two selected cells produce the complementary,
    re-entrant ridge and therefore carry the negative turning magnitude of the
    excluded cell. Zero/three selected generators are not union-boundary edges.
    """
    targets = _target_set(target_cells)
    balls = _row_balls(row)
    if targets is None or len(balls) != 3:
        return None

    inside = [ball for ball in balls if ball in targets]
    outside = [ball for ball in balls if ball not in targets]
    mapping = row.get("int_mean_curv_by_ball", None)

    if len(inside) == 1:
        value = _mapping_value(mapping, inside[0])
        return 0.0 if value is None else value
    if len(inside) == 2:
        value = _mapping_value(mapping, outside[0])
        return 0.0 if value is None else -value
    return 0.0


def component_value(row, component, scheme, target_cells=None, mode="boundary"):
    """Return one scalar used to color a geometric curvature component."""
    scheme = canonical_curvature_scheme(scheme)
    mode = canonical_color_mode(mode)
    if scheme is None:
        return None

    mapping_key = f"{scheme}_by_ball"
    mapping = row.get(mapping_key, None)

    # Non-directional physical geometry: magnitude only, never arbitrary sign.
    if mode == "magnitude":
        if isinstance(mapping, dict):
            return _magnitude_mapping(mapping, target_cells)
        if component == "surface" and scheme in row.index:
            try:
                value = float(row[scheme])
            except (TypeError, ValueError):
                return None
            return abs(value) if np.isfinite(value) else None
        return None

    # Individual-cell geometry: use that cell's oriented value directly.
    if mode == "cell":
        targets = _target_set(target_cells)
        if isinstance(mapping, dict) and targets:
            values = [
                value for cell in targets
                if (value := _mapping_value(mapping, cell)) is not None
            ]
            if values:
                return float(sum(values))
        if component == "surface" and scheme in row.index:
            try:
                value = float(row[scheme])
            except (TypeError, ValueError):
                return None
            return value if np.isfinite(value) else None
        return None

    # Signed molecular/group boundary. Mean curvature needs union-boundary
    # orientation rather than a sum of positive cell turning magnitudes.
    if scheme == "int_mean_curv":
        if component == "surface":
            return _boundary_surface_mean_value(row, target_cells)
        if component == "edge":
            return _boundary_edge_mean_value(row, target_cells)
        return None

    # Gaussian boundary coloring currently retains the validated cell-relative
    # data. Exact union-boundary G bookkeeping remains a separate theory step.
    if component not in SUPPORTED_COMPONENTS[scheme]:
        return None
    if isinstance(mapping, dict):
        return _sum_mapping(mapping, target_cells)
    if component == "surface" and scheme in row.index:
        try:
            value = float(row[scheme])
        except (TypeError, ValueError):
            return None
        return value if np.isfinite(value) else None
    return None


def component_values(table, component, scheme, target_cells=None, indices=None,
                     mode="boundary"):
    if table is None or len(table) == 0:
        return np.empty(0, dtype=float)
    positions = range(len(table)) if indices is None else list(indices)
    values = []
    for position in positions:
        try:
            row = table.iloc[int(position)]
        except (IndexError, TypeError, ValueError):
            continue
        value = component_value(row, component, scheme, target_cells, mode=mode)
        if value is not None and np.isfinite(value):
            values.append(float(value))
    return np.asarray(values, dtype=float)


def _cell_mapping_values(table, component, scheme, target_cells=None):
    if table is None or len(table) == 0:
        return []
    key = f"{scheme}_by_ball"
    values = []
    for _, row in table.iterrows():
        if key in row.index:
            values.extend(_mapping_values(row[key], target_cells))
        elif component == "surface" and scheme in row.index:
            try:
                value = float(row[scheme])
            except (TypeError, ValueError):
                continue
            if np.isfinite(value):
                values.append(value)
    return values


def curvature_color_limit(net, scheme, target_cells=None, mode="boundary"):
    """Return a common absolute scale for one visualization interpretation."""
    scheme = canonical_curvature_scheme(scheme)
    mode = canonical_color_mode(mode)
    if scheme is None:
        return None

    arrays = []
    support = SUPPORTED_COMPONENTS[scheme]

    # For atom-cell output the global scale should contain individual cell
    # contributions, not sums across all requested atoms.
    if mode == "cell":
        if "surface" in support:
            arrays.extend(_cell_mapping_values(net.surfs, "surface", scheme, target_cells))
        if "edge" in support:
            arrays.extend(_cell_mapping_values(net.edges, "edge", scheme, target_cells))
        if "vertex" in support:
            arrays.extend(_cell_mapping_values(net.verts, "vertex", scheme, target_cells))
        values = np.asarray(arrays, dtype=float)
    else:
        collected = []
        if "surface" in support:
            collected.append(component_values(net.surfs, "surface", scheme, target_cells, mode=mode))
        if "edge" in support:
            collected.append(component_values(net.edges, "edge", scheme, target_cells, mode=mode))
        if "vertex" in support:
            collected.append(component_values(net.verts, "vertex", scheme, target_cells, mode=mode))
        collected = [array for array in collected if array.size]
        values = np.concatenate(collected) if collected else np.empty(0, dtype=float)

    values = values[np.isfinite(values)]
    if values.size == 0:
        return 1.0
    limit = float(np.max(np.abs(values)))
    return limit if np.isfinite(limit) and limit > 0.0 else 1.0


def _signed_color_coordinate(value, limit, factor="log", log_gain=1000.0):
    value = float(value)
    limit = abs(float(limit))
    if not np.isfinite(value):
        value = 0.0
    if not np.isfinite(limit) or limit <= 0.0:
        limit = 1.0

    ratio = np.clip(value / limit, -1.0, 1.0)
    factor = "log" if factor is None else str(factor).strip().lower()
    if factor in {"lin", "linear", "none"}:
        return float(ratio)
    if factor in {"log", "signed_log", "symlog"}:
        gain = float(log_gain)
        if not np.isfinite(gain) or gain <= 0.0:
            gain = 1000.0
        transformed = np.log1p(gain * abs(ratio)) / np.log1p(gain)
        return float(np.sign(ratio) * transformed)
    return float(ratio)


def curvature_color(value, color_map="coolwarm", limit=1.0, factor="log",
                    log_gain=1000.0, mode="boundary"):
    """Map curvature to RGB using signed or magnitude-only normalization."""
    mode = canonical_color_mode(mode)
    try:
        if mode == "magnitude":
            coordinate = abs(_signed_color_coordinate(
                abs(float(value)), limit, factor=factor, log_gain=log_gain,
            ))
            # Use only the neutral -> positive half of the diverging map.
            normalized = np.clip(0.5 + 0.5 * coordinate, 0.5, 1.0)
        else:
            coordinate = _signed_color_coordinate(
                value, limit, factor=factor, log_gain=log_gain,
            )
            normalized = np.clip(0.5 + 0.5 * coordinate, 0.0, 1.0)
    except (TypeError, ValueError):
        normalized = 0.5
    return np.asarray(_cmap(color_map)(normalized)[:3], dtype=float)


def component_color(row, component, scheme, target_cells=None,
                    color_map="coolwarm", limit=1.0, factor="log",
                    log_gain=1000.0, mode="boundary"):
    value = component_value(row, component, scheme, target_cells, mode=mode)
    if value is None:
        return None
    return curvature_color(
        value, color_map=color_map, limit=limit, factor=factor,
        log_gain=log_gain, mode=mode,
    )


def mean_vertex_display_value(net, vertex_row, target_cells=None, mode="boundary"):
    """Display-only local mean-curvature proxy for a geometric vertex.

    Integrated mean curvature has no independent vertex term. For visualizing
    local convex/re-entrant character, use the average incident edge-M value.
    This value MUST NOT be included in M totals or free-energy sums.
    """
    try:
        edge_indices = [int(value) for value in vertex_row["edges"]]
    except (KeyError, TypeError, ValueError):
        return None

    values = []
    for edge_index in edge_indices:
        try:
            edge = net.edges.iloc[edge_index]
        except (IndexError, TypeError, ValueError):
            continue
        value = component_value(
            edge, "edge", "int_mean_curv", target_cells, mode=mode,
        )
        if value is not None and np.isfinite(value):
            # Ignore edges that are not part of the boundary in boundary mode.
            if mode == "boundary" and abs(value) < 1e-15:
                continue
            values.append(float(value))

    return float(np.mean(values)) if values else 0.0


def mean_vertex_display_color(net, vertex_row, target_cells=None,
                              color_map="coolwarm", limit=1.0, factor="log",
                              log_gain=1000.0, mode="boundary"):
    value = mean_vertex_display_value(net, vertex_row, target_cells, mode=mode)
    if value is None:
        return None
    return curvature_color(
        value, color_map=color_map, limit=limit, factor=factor,
        log_gain=log_gain, mode=mode,
    )
