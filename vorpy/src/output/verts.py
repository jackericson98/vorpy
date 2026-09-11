"""Prepare and export drawable Voronoi vertex geometry."""

import numpy as np

from vorpy.src.output.colors import color_dict
from vorpy.src.output.curvature_colors import (
    canonical_curvature_scheme,
    component_color,
    curvature_color_limit,
    mean_vertex_display_color,
)
from vorpy.src.output.draw import DEFAULT_VERTEX_RADIUS, draw_joint
from vorpy.src.output.mesh import combine_mesh_parts, write_mesh


def _resolve_color(color):
    if color is None:
        return color_dict.get("red", [1.0, 0.0, 0.0])
    if isinstance(color, str):
        return color_dict.get(color, color_dict.get("red", [1.0, 0.0, 0.0]))
    return color


def prepare_verts(net, verts, color=None, vert_rad=DEFAULT_VERTEX_RADIUS,
                  subdivisions=0, color_scheme=None, color_map=None,
                  color_limit=None, target_cells=None, color_mode="boundary"):
    """Prepare selected vertices with optional integrated-G coloring.

    Vertices carry Gaussian angular-defect curvature but no mean-curvature
    term in the current piecewise-smooth decomposition. Therefore
    ``int_mean_curv`` intentionally leaves the normal fixed vertex color.
    """
    if verts is None or len(verts) == 0:
        return None

    fixed_color = np.asarray(_resolve_color(color), dtype=float)

    settings = getattr(net, "settings", None) or {}
    scheme = canonical_curvature_scheme(
        settings.get("surf_scheme") if color_scheme is None else color_scheme
    )
    cmap = settings.get("surf_col", "coolwarm") if color_map is None else color_map

    if scheme is not None and color_limit is None:
        color_limit = curvature_color_limit(net, scheme, target_cells, mode=color_mode)

    point_parts, triangle_parts, color_parts, index_parts = [], [], [], []

    for index in list(verts):
        vertex = net.verts.iloc[index]
        location = np.asarray(vertex["loc"], dtype=float)
        points, triangles = draw_joint(
            location,
            radius=vert_rad,
            subdivisions=subdivisions,
        )

        vertex_color = fixed_color
        if scheme == "int_mean_curv":
            # Display-only local crease proxy. Mean curvature has no separate
            # vertex term, so this color is never used in curvature totals.
            mapped = mean_vertex_display_color(
                net, vertex, target_cells, color_map=cmap, limit=color_limit,
                factor=settings.get("scheme_factor", "log"), mode=color_mode,
            )
            if mapped is not None:
                vertex_color = mapped
        elif scheme == "int_gauss_curv":
            mapped = component_color(
                vertex, "vertex", scheme, target_cells,
                color_map=cmap, limit=color_limit,
                factor=settings.get("scheme_factor", "log"), mode=color_mode,
            )
            if mapped is not None:
                vertex_color = mapped

        point_parts.append(points)
        triangle_parts.append(triangles)
        color_parts.append([vertex_color] * len(triangles))
        index_parts.append(np.full(len(triangles), index, dtype=np.int64))

    return combine_mesh_parts(
        point_parts,
        triangle_parts,
        color_parts,
        {"vertex_index": index_parts},
    )


def write_verts(net, verts, file_name, atom_type=None, directory=None, color=None,
                vert_rad=DEFAULT_VERTEX_RADIUS, file_type="off", chunk_size=10000,
                subdivisions=0, color_scheme=None, color_map=None,
                color_limit=None, target_cells=None, color_mode="boundary"):
    """Prepare selected vertices once and write OFF, PLY, or VTP."""
    mesh = prepare_verts(
        net, verts, color, vert_rad, subdivisions,
        color_scheme=color_scheme, color_map=color_map,
        color_limit=color_limit, target_cells=target_cells, color_mode=color_mode,
    )
    if mesh is None:
        return None
    return write_mesh(mesh, file_name, file_type, directory, chunk_size)


def write_off_verts(net, verts, file_name, atom_type=None, directory=None, color=None,
                    vert_rad=DEFAULT_VERTEX_RADIUS, file_type="off", chunk_size=10000,
                    subdivisions=0, color_scheme=None, color_map=None,
                    color_limit=None, target_cells=None, color_mode="boundary"):
    """Backward-compatible vertex-export entry point."""
    return write_verts(
        net, verts, file_name, atom_type=atom_type, directory=directory,
        color=color, vert_rad=vert_rad, file_type=file_type,
        chunk_size=chunk_size, subdivisions=subdivisions,
        color_scheme=color_scheme, color_map=color_map,
        color_limit=color_limit, target_cells=target_cells, color_mode=color_mode,
    )


def write_off_verts1(verts, file_name, atom_type=None, directory=None, color=None,
                     vert_rad=DEFAULT_VERTEX_RADIUS, file_type="off", chunk_size=10000,
                     subdivisions=0):
    """Legacy DataFrame vertex writer retained for compatibility."""
    if verts is None or len(verts) == 0:
        return None

    color = _resolve_color(color)
    point_parts, triangle_parts, color_parts = [], [], []

    for _, vertex in verts.iterrows():
        points, triangles = draw_joint(
            vertex["loc"],
            radius=vert_rad,
            subdivisions=subdivisions,
        )
        point_parts.append(points)
        triangle_parts.append(triangles)
        color_parts.append([color] * len(triangles))

    mesh = combine_mesh_parts(point_parts, triangle_parts, color_parts)
    return write_mesh(mesh, file_name, file_type, directory, chunk_size)
