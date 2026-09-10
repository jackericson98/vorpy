"""Network-level vertex contributions to integrated Gaussian curvature."""

from time import perf_counter as now

import numpy as np

from vorpy.src.calculations.vertex_geometry import (
    aw_cell_vertex_angular_defect,
)
from vorpy.src.calculations.edge_resolution_cache import get_aw_edge_geometry_cache


def calculate_aw_network_vertex_gaussian_curvatures(
        net,
        tolerance=1e-7):
    """Calculate angular-defect Gaussian curvature for target AW cells.

    Each network vertex stores a dictionary mapping participating target
    cells to their angular defect at that vertex.

    Only cells in ``net.group`` are evaluated. Surrounding/support balls may
    remain in the network for construction but are not treated as closed cells.

    Analytic edges incident to each vertex are resolved once and reused for
    every target cell evaluated at that vertex.
    """
    settings = getattr(net, "settings", None) or {}

    if settings.get("net_type", "aw") != "aw":
        raise ValueError("AW vertex Gaussian curvature requires an AW network.")
    if getattr(net, "verts", None) is None:
        raise ValueError("The network has no vertex table.")
    if getattr(net, "edges", None) is None:
        raise ValueError("The network has no edge table.")
    if getattr(net, "surfs", None) is None:
        raise ValueError("The network has no surface table.")
    if net.group is None:
        raise ValueError("AW vertex Gaussian curvature requires a target network group.")

    total_start = now()
    timing = {"edge_resolution": 0.0, "vertex_geometry": 0.0, "storage": 0.0}
    target_cells = {int(value) for value in net.group}
    contributions = []
    resolved_edges, cache_built, cache_time = get_aw_edge_geometry_cache(net, tolerance=tolerance)
    timing["edge_resolution"] += cache_time
    resolved_edge_calls = 0
    defect_values = 0

    for vertex_index, vertex in net.verts.iterrows():
        vertex_balls = {int(value) for value in vertex["balls"]}
        participating_cells = sorted(vertex_balls.intersection(target_cells))

        if not participating_cells:
            t = now()
            contributions.append({})
            timing["storage"] += now() - t
            continue

        vertex_edges = tuple(int(value) for value in vertex["edges"])

        if len(vertex_edges) < 3:
            raise ValueError(
                f"Target vertex {vertex_index} has only {len(vertex_edges)} incident edges."
            )

        # Reuse the network-wide analytic cache; only select this vertex's edges.
        vertex_resolved_edges = {
            edge_index: resolved_edges[edge_index]
            for edge_index in vertex_edges
        }
        resolved_edge_calls += len(vertex_edges)

        vertex_values = {}

        for cell_index in participating_cells:
            t = now()
            defect = aw_cell_vertex_angular_defect(
                net=net,
                vertex_index=vertex_index,
                cell_index=cell_index,
                resolved_edges=vertex_resolved_edges,
                tolerance=tolerance,
            )
            timing["vertex_geometry"] += now() - t

            defect = float(defect)

            if not np.isfinite(defect):
                raise ValueError(
                    "Non-finite AW angular defect for "
                    f"vertex {vertex_index}, ball {cell_index}: {defect}"
                )

            t = now()
            vertex_values[cell_index] = defect
            timing["storage"] += now() - t
            defect_values += 1

        if len(vertex_values) != len(participating_cells):
            raise ValueError(
                f"Vertex {vertex_index} did not produce one Gaussian-curvature "
                "defect for every participating target cell."
            )

        t = now()
        contributions.append(vertex_values)
        timing["storage"] += now() - t

    if len(contributions) != len(net.verts):
        raise ValueError(
            "Vertex Gaussian-curvature calculation did not return one result for every network vertex."
        )

    if settings.get("verbose", False):
        total = now() - total_start
        measured = sum(timing.values())
        overhead = max(total - measured, 0.0)

        print("\n" + "=" * 70)
        print("VERTEX GAUSSIAN CURVATURE TIMING")
        print("=" * 70)
        for label, value in (
            ("Shared edge cache build" if cache_built else "Shared edge cache reuse", timing["edge_resolution"]),
            ("Vertex angular defects", timing["vertex_geometry"]),
            ("Storage", timing["storage"]),
            ("Other / loop overhead", overhead),
        ):
            pct = 100.0 * value / total if total > 0 else 0.0
            print(f"{label:<30} {value:10.4f} s  {pct:6.2f} %")

        print("-" * 70)
        print(f"{'TOTAL':<30} {total:10.4f} s  100.00 %")
        print()
        print(f"Vertices:            {len(net.verts):,}")
        print(f"Resolved edge calls: {resolved_edge_calls:,}")
        print(f"Cell defects:        {defect_values:,}")
        print("=" * 70)

    return contributions
