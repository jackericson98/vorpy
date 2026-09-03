"""Prepare and export drawable Voronoi vertex geometry."""

import numpy as np

from vorpy.src.output.colors import color_dict
from vorpy.src.output.draw import DEFAULT_VERTEX_RADIUS, draw_joint
from vorpy.src.output.mesh import combine_mesh_parts, write_mesh


def _resolve_color(color):
    if color is None:
        return color_dict.get('red', [1.0, 0.0, 0.0])
    if isinstance(color, str):
        return color_dict.get(color, color_dict.get('red', [1.0, 0.0, 0.0]))
    return color


def prepare_verts(net, verts, color=None, vert_rad=DEFAULT_VERTEX_RADIUS,
                  subdivisions=0):
    """Prepare selected Voronoi vertices as icosahedral markers.

    The markers use the same geometry generator as edge junctions. Standalone
    vertices default to a slightly larger radius so they remain visible within
    the connected edge network.
    """
    if verts is None or len(verts) == 0:
        return None

    color = _resolve_color(color)
    point_parts = []
    triangle_parts = []
    color_parts = []
    index_parts = []

    for index in list(verts):
        location = np.asarray(net.verts['loc'][index], dtype=float)
        points, triangles = draw_joint(
            location,
            radius=vert_rad,
            subdivisions=subdivisions,
        )
        point_parts.append(points)
        triangle_parts.append(triangles)
        color_parts.append([color] * len(triangles))
        index_parts.append(np.full(len(triangles), index, dtype=np.int64))

    return combine_mesh_parts(
        point_parts,
        triangle_parts,
        color_parts,
        {'vertex_index': index_parts},
    )


def write_verts(net, verts, file_name, atom_type=None, directory=None, color=None,
                vert_rad=DEFAULT_VERTEX_RADIUS, file_type='off', chunk_size=10000,
                subdivisions=0):
    """Prepare selected vertices once and write OFF, PLY, or VTP."""
    mesh = prepare_verts(net, verts, color, vert_rad, subdivisions)
    if mesh is None:
        return None
    return write_mesh(mesh, file_name, file_type, directory, chunk_size)


def write_off_verts(net, verts, file_name, atom_type=None, directory=None, color=None,
                    vert_rad=DEFAULT_VERTEX_RADIUS, file_type='off', chunk_size=10000,
                    subdivisions=0):
    """Backward-compatible vertex-export entry point; OFF remains the default."""
    return write_verts(
        net,
        verts,
        file_name,
        atom_type=atom_type,
        directory=directory,
        color=color,
        vert_rad=vert_rad,
        file_type=file_type,
        chunk_size=chunk_size,
        subdivisions=subdivisions,
    )


def write_off_verts1(verts, file_name, atom_type=None, directory=None, color=None,
                     vert_rad=DEFAULT_VERTEX_RADIUS, file_type='off', chunk_size=10000,
                     subdivisions=0):
    """Legacy DataFrame vertex writer retained for compatibility."""
    if verts is None or len(verts) == 0:
        return None

    color = _resolve_color(color)
    point_parts = []
    triangle_parts = []
    color_parts = []

    for _, vertex in verts.iterrows():
        points, triangles = draw_joint(
            vertex['loc'],
            radius=vert_rad,
            subdivisions=subdivisions,
        )
        point_parts.append(points)
        triangle_parts.append(triangles)
        color_parts.append([color] * len(triangles))

    mesh = combine_mesh_parts(point_parts, triangle_parts, color_parts)
    return write_mesh(mesh, file_name, file_type, directory, chunk_size)
