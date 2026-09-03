"""Prepare and export drawable Voronoi edge geometry."""

import os

import numpy as np

from vorpy.src.output.colors import color_dict
from vorpy.src.output.draw import (
    DEFAULT_EDGE_RADIUS,
    DEFAULT_JOINT_RADIUS_FACTOR,
    draw_edge,
    draw_joint,
)
from vorpy.src.output.mesh import combine_mesh_parts, write_mesh


def _resolve_color(color):
    if color is None:
        return color_dict.get('gray', [0.5, 0.5, 0.5])
    if isinstance(color, str):
        return color_dict.get(color, color_dict.get('gray', [0.5, 0.5, 0.5]))
    return color


def _cache_columns(net):
    if 'draw_points' not in net.edges:
        net.edges['draw_points'] = [[] for _ in range(len(net.edges))]
    if 'draw_tris' not in net.edges:
        net.edges['draw_tris'] = [[] for _ in range(len(net.edges))]


def prepare_edges(net, edges, color=None, radius=DEFAULT_EDGE_RADIUS, add_joints=True,
                  joint_radius=None, joint_subdivisions=0, endpoint_tolerance=1e-6):
    """Prepare selected edge tubes and deduplicated endpoint junctions.

    Junctions slightly overlap every incident tube. This produces visually
    continuous intersections even though each tube has an independently
    transported triangular cross-section.
    """
    if edges is None or len(edges) == 0:
        return None
    if endpoint_tolerance <= 0:
        raise ValueError('endpoint_tolerance must be positive')

    _cache_columns(net)
    edge_indices = list(edges)
    color = _resolve_color(color)
    joint_radius = DEFAULT_JOINT_RADIUS_FACTOR * radius if joint_radius is None else joint_radius

    point_parts = []
    triangle_parts = []
    color_parts = []
    index_parts = []
    kind_parts = []
    endpoints = {}

    for index in edge_indices:
        edge = net.edges.iloc[index]
        draw_points = edge['draw_points']
        draw_tris = edge['draw_tris']

        if draw_points is None or draw_tris is None or len(draw_points) == 0 or len(draw_tris) == 0:
            draw_points, draw_tris = draw_edge(edge, radius=radius)
            net.edges.at[net.edges.index[index], 'draw_points'] = draw_points
            net.edges.at[net.edges.index[index], 'draw_tris'] = draw_tris

        point_parts.append(draw_points)
        triangle_parts.append(draw_tris)
        color_parts.append([color] * len(draw_tris))
        index_parts.append(np.full(len(draw_tris), index, dtype=np.int64))
        kind_parts.append(np.zeros(len(draw_tris), dtype=np.int64))

        if add_joints and len(edge['points']) > 0:
            for endpoint in (edge['points'][0], edge['points'][-1]):
                endpoint = np.asarray(endpoint, dtype=float)
                key = tuple(np.rint(endpoint / endpoint_tolerance).astype(np.int64))
                endpoints.setdefault(key, endpoint)

    for endpoint in endpoints.values():
        joint_points, joint_tris = draw_joint(
            endpoint,
            radius=joint_radius,
            subdivisions=joint_subdivisions,
        )
        point_parts.append(joint_points)
        triangle_parts.append(joint_tris)
        color_parts.append([color] * len(joint_tris))
        index_parts.append(np.full(len(joint_tris), -1, dtype=np.int64))
        kind_parts.append(np.ones(len(joint_tris), dtype=np.int64))

    return combine_mesh_parts(
        point_parts,
        triangle_parts,
        color_parts,
        {'edge_index': index_parts, 'geometry_kind': kind_parts},
    )


def write_edges(net, edges, file_name, color=None, directory=None, profile=True,
                file_type='off', chunk_size=10000, radius=DEFAULT_EDGE_RADIUS,
                add_joints=True, joint_radius=None, joint_subdivisions=0):
    """Prepare selected edges once and write OFF, PLY, or VTP."""
    mesh = prepare_edges(
        net,
        edges,
        color=color,
        radius=radius,
        add_joints=add_joints,
        joint_radius=joint_radius,
        joint_subdivisions=joint_subdivisions,
    )
    if mesh is None:
        return None
    return write_mesh(mesh, file_name, file_type, directory, chunk_size)


def write_edges1(edges, file_name, color=None, directory=None, radius=DEFAULT_EDGE_RADIUS,
                 file_type='off', chunk_size=10000):
    """Legacy DataFrame edge writer retained for compatibility."""
    if edges is None or len(edges) == 0:
        return None
    color = _resolve_color(color)
    point_parts, triangle_parts, color_parts = [], [], []
    for _, edge in edges.iterrows():
        draw_points, draw_tris = draw_edge(edge, radius=radius)
        point_parts.append(draw_points)
        triangle_parts.append(draw_tris)
        color_parts.append([color] * len(draw_tris))
    mesh = combine_mesh_parts(point_parts, triangle_parts, color_parts)
    return write_mesh(mesh, file_name, file_type, directory, chunk_size)
