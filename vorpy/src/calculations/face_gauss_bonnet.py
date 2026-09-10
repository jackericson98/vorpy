import numpy as np

from vorpy.src.calculations.edge_geometry import (
    network_surface_boundary_cycle,
    oriented_aw_face_edge_geodesic_curvature,
)
from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge


def aw_cell_gauss_bonnet(
        net,
        cell_index,
        quadrature_order=32,
        tolerance=1e-7):
    """Evaluate intrinsic Gaussian-curvature terms for one complete AW cell."""
    from vorpy.src.calculations.edge_geometry import (
        oriented_aw_face_edge_geodesic_curvature,
    )
    from vorpy.src.calculations.vertex_geometry import (
        aw_vertex_face_angle,
    )
    from vorpy.src.network.edge_geometry_diagnostics import (
        resolve_aw_network_edge,
    )

    cell_index = int(cell_index)

    # Locate the ball robustly; ball number need not equal DataFrame position.
    matches = net.balls.index[net.balls["num"].astype(int) == cell_index]

    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one ball numbered {cell_index}; found {len(matches)}."
        )

    ball = net.balls.loc[matches[0]]
    surface_indices = [int(v) for v in ball["surfs"]]
    edge_indices = [int(v) for v in ball["edges"]]
    vertex_indices = [int(v) for v in ball["verts"]]

    # --------------------------------------------------------------
    # Smooth Gaussian curvature
    # --------------------------------------------------------------
    surface_sum = sum(
        float(net.surfs.loc[surface_index, "int_gauss_curv"])
        for surface_index in surface_indices
    )

    # Resolve each cell edge once.
    resolved_edges = {
        edge_index: resolve_aw_network_edge(
            net,
            edge_index,
            tolerance=tolerance,
        )
        for edge_index in edge_indices
    }

    # --------------------------------------------------------------
    # Edge singular contribution
    #
    # Every cell edge bounds exactly two cell faces. Each face contributes
    # its own signed geodesic-curvature integral.
    # --------------------------------------------------------------
    edge_sum = 0.0

    for edge_index in edge_indices:
        resolved = resolved_edges[edge_index]

        incident_surfaces = [
            int(surface_index)
            for surface_index in net.edges.loc[edge_index, "surfs"]
            if cell_index in {
                int(v)
                for v in net.surfs.loc[int(surface_index), "balls"]
            }
        ]

        if len(incident_surfaces) != 2:
            raise ValueError(
                f"Cell {cell_index}, edge {edge_index} has "
                f"{len(incident_surfaces)} incident cell surfaces; expected 2."
            )

        for surface_index in incident_surfaces:
            surface_balls = tuple(
                int(v)
                for v in net.surfs.loc[surface_index, "balls"]
            )

            other_index = (
                surface_balls[1]
                if surface_balls[0] == cell_index
                else surface_balls[0]
            )

            edge_sum += oriented_aw_face_edge_geodesic_curvature(
                edge_geometry=resolved.geometry,
                generator_indices=resolved.ball_indices,
                generator_locations=resolved.locations,
                cell_index=cell_index,
                face_other_index=other_index,
                order=quadrature_order,
            )

    # --------------------------------------------------------------
    # Vertex angular defects
    # --------------------------------------------------------------
    vertex_sum = 0.0

    for vertex_index in vertex_indices:
        incident_surfaces = [
            int(surface_index)
            for surface_index in net.verts.loc[vertex_index, "surfs"]
            if cell_index in {
                int(v)
                for v in net.surfs.loc[int(surface_index), "balls"]
            }
        ]

        if len(incident_surfaces) < 3:
            raise ValueError(
                f"Cell {cell_index}, vertex {vertex_index} has only "
                f"{len(incident_surfaces)} incident surfaces."
            )

        interior_angles = [
            aw_vertex_face_angle(
                net=net,
                vertex_index=vertex_index,
                cell_index=cell_index,
                surface_index=surface_index,
                resolved_edges=resolved_edges,
                tolerance=tolerance,
            )
            for surface_index in incident_surfaces
        ]

        vertex_sum += (
            2.0 * np.pi
            - float(np.sum(interior_angles))
        )

    total = surface_sum + edge_sum + vertex_sum
    target = 4.0 * np.pi

    return {
        "surface": float(surface_sum),
        "edge": float(edge_sum),
        "vertex": float(vertex_sum),
        "total": float(total),
        "target": float(target),
        "error": float(total - target),
        "surfaces": len(surface_indices),
        "edges": len(edge_indices),
        "vertices": len(vertex_indices),
    }