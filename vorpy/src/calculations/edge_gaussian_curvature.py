"""Network-level edge contribution to integrated Gaussian curvature."""

from time import perf_counter as now

import numpy as np

from vorpy.src.calculations.edge_geometry import (
    oriented_aw_face_edge_geodesic_curvature,
)
from vorpy.src.network.edge_geometry_diagnostics import (
    resolve_aw_network_edge,
)


def calculate_aw_network_edge_gaussian_curvatures(
        net,
        quadrature_order=32,
        tolerance=1e-6):
    """Calculate geodesic-curvature contribution for each AW cell edge.

    For a given cell, an edge belongs to two of that cell's boundary faces.
    The cell-relative edge contribution is therefore the sum of the signed
    geodesic-curvature integrals from those two face-edge incidences.

    Returns
    -------
    list[dict]
        One dictionary per edge:

            {cell_index: integrated_edge_G, ...}
    """
    settings = getattr(net, "settings", None) or {}

    if net.group is None:
        raise ValueError(
            "AW edge Gaussian curvature requires a target network group."
        )

    target_cells = {int(value) for value in net.group}

    if settings.get("net_type", "aw") != "aw":
        raise ValueError("AW edge Gaussian curvature requires an AW network.")

    total_start = now()
    timing = {"resolution": 0.0, "surface_lookup": 0.0, "integration": 0.0, "storage": 0.0}
    contributions = []
    integration_calls = 0
    cell_values = 0

    for edge_index, edge in net.edges.iterrows():
        t = now()
        resolved = resolve_aw_network_edge(net, edge_index, tolerance=tolerance)
        timing["resolution"] += now() - t

        edge_balls = tuple(int(v) for v in resolved.ball_indices)

        if len(edge_balls) != 3:
            raise ValueError(f"Edge {edge_index} does not have exactly three generators.")

        edge_surfaces = [int(v) for v in edge["surfs"]]
        edge_values = {}

        for cell_index in edge_balls:
            if cell_index not in target_cells:
                continue

            t = now()
            incident_surfaces = [
                surface_index
                for surface_index in edge_surfaces
                if cell_index in {int(v) for v in net.surfs.loc[surface_index, "balls"]}
            ]
            timing["surface_lookup"] += now() - t

            # Support/incomplete cells may not retain both adjacent surfaces.
            if len(incident_surfaces) != 2:
                continue

            value = 0.0

            for surface_index in incident_surfaces:
                surface_balls = tuple(int(v) for v in net.surfs.loc[surface_index, "balls"])
                other_index = surface_balls[1] if surface_balls[0] == cell_index else surface_balls[0]

                t = now()
                value += oriented_aw_face_edge_geodesic_curvature(
                    edge_geometry=resolved.geometry,
                    generator_indices=resolved.ball_indices,
                    generator_locations=resolved.locations,
                    cell_index=cell_index,
                    face_other_index=other_index,
                    order=quadrature_order,
                )
                timing["integration"] += now() - t
                integration_calls += 1

            if not np.isfinite(value):
                raise ValueError(
                    f"Non-finite edge Gaussian curvature for edge {edge_index}, "
                    f"ball {cell_index}: {value}"
                )

            t = now()
            edge_values[cell_index] = float(value)
            timing["storage"] += now() - t
            cell_values += 1

        contributions.append(edge_values)

    if len(contributions) != len(net.edges):
        raise ValueError(
            "Edge Gaussian-curvature calculation did not return one result for every network edge."
        )

    if settings.get("verbose", False):
        total = now() - total_start
        measured = sum(timing.values())
        overhead = max(total - measured, 0.0)

        print("\n" + "=" * 70)
        print("EDGE GAUSSIAN CURVATURE TIMING")
        print("=" * 70)
        for label, value in (
            ("Analytic edge resolution", timing["resolution"]),
            ("Surface / topology lookup", timing["surface_lookup"]),
            ("Geodesic integration", timing["integration"]),
            ("Storage", timing["storage"]),
            ("Other / loop overhead", overhead),
        ):
            pct = 100.0 * value / total if total > 0 else 0.0
            print(f"{label:<30} {value:10.4f} s  {pct:6.2f} %")

        print("-" * 70)
        print(f"{'TOTAL':<30} {total:10.4f} s  100.00 %")
        print()
        print(f"Edges:              {len(net.edges):,}")
        print(f"Quadrature order:   {quadrature_order}")
        print(f"Face integrations:  {integration_calls:,}")
        print(f"Cell values:        {cell_values:,}")
        print("=" * 70)

    return contributions
