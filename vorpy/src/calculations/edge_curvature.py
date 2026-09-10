"""Fused network-level AW edge curvature calculation."""

from time import perf_counter

import numpy as np

from vorpy.src.calculations.edge_geometry import aw_edge_curvature_measures
from vorpy.src.calculations.edge_resolution_cache import get_aw_edge_geometry_cache


_CACHE_ATTR = "_aw_fused_edge_curvature_cache"


def _cache_key(net, quadrature_order, tolerance):
    group = () if net.group is None else tuple(sorted(int(value) for value in net.group))
    return (id(net.edges), len(net.edges), int(quadrature_order), float(tolerance), group)


def calculate_aw_network_edge_curvatures(net, quadrature_order=32, tolerance=1e-6):
    """Calculate edge M and edge G together in one traversal per physical edge.

    Returns
    -------
    tuple[list[dict], list[dict]]
        ``(mean_contributions, gaussian_contributions)`` with one dictionary
        per network edge.
    """
    settings = getattr(net, "settings", None) or {}

    if settings.get("net_type", "aw") != "aw":
        raise ValueError("Fused AW edge curvature requires an AW network.")
    if getattr(net, "edges", None) is None:
        raise ValueError("The network has no edge table.")
    if getattr(net, "surfs", None) is None:
        raise ValueError("The network has no surface table.")
    if net.group is None:
        raise ValueError("Fused AW edge curvature requires a target network group.")

    key = _cache_key(net, quadrature_order, tolerance)
    cached = getattr(net, _CACHE_ATTR, None)
    if isinstance(cached, dict) and cached.get("key") == key:
        return cached["mean"], cached["gaussian"]

    verbose = settings.get("verbose", False)
    total_start = perf_counter()
    target_cells = {int(value) for value in net.group}
    n_edges = len(net.edges)
    quadrature = np.polynomial.legendre.leggauss(int(quadrature_order))

    resolved_edges, cache_built, cache_time = get_aw_edge_geometry_cache(
        net, tolerance=tolerance,
    )

    topology_time = 0.0
    integration_time = 0.0
    kernel_timing = {
        "orientation": 0.0,
        "point_tangent": 0.0,
        "second_derivative": 0.0,
        "radial": 0.0,
        "normals": 0.0,
        "mean": 0.0,
        "gaussian": 0.0,
        "other": 0.0,
    }
    reduction_time = 0.0
    mean_contributions = []
    gaussian_contributions = []
    face_values = 0
    gauss_cell_values = 0

    for count, (edge_index, edge) in enumerate(net.edges.iterrows(), start=1):
        resolved = resolved_edges[int(edge_index)]
        edge_balls = tuple(int(value) for value in resolved.ball_indices)

        if len(edge_balls) != 3:
            raise ValueError(f"Edge {edge_index} does not have exactly three generators.")

        # Identify only the G face incidences needed by requested cells.
        t = perf_counter()
        edge_surfaces = [int(value) for value in edge["surfs"]]
        cell_faces = {}

        for cell_index in edge_balls:
            if cell_index not in target_cells:
                continue

            incident_surfaces = [
                surface_index for surface_index in edge_surfaces
                if cell_index in {
                    int(value) for value in net.surfs.loc[surface_index, "balls"]
                }
            ]
            if len(incident_surfaces) != 2:
                continue

            others = []
            for surface_index in incident_surfaces:
                surface_balls = tuple(
                    int(value) for value in net.surfs.loc[surface_index, "balls"]
                )
                other = surface_balls[1] if surface_balls[0] == cell_index else surface_balls[0]
                others.append(int(other))
            cell_faces[cell_index] = tuple(others)

        face_pairs = tuple(
            (cell_index, other_index)
            for cell_index, others in cell_faces.items()
            for other_index in others
        )
        topology_time += perf_counter() - t

        t = perf_counter()
        measures = aw_edge_curvature_measures(
            edge_geometry=resolved.geometry,
            generator_indices=resolved.ball_indices,
            generator_locations=resolved.locations,
            face_pairs=face_pairs,
            order=quadrature_order,
            quadrature=quadrature,
            tol=tolerance,
            timing=kernel_timing if verbose else None,
        )
        integration_time += perf_counter() - t

        mean_values = measures["mean"]
        face_integrals = measures["gaussian"]

        if len(mean_values) != 3:
            raise ValueError(
                f"Edge {edge_index} produced {len(mean_values)} mean-curvature "
                "cell contributions instead of 3."
            )

        t = perf_counter()
        gauss_values = {}
        for cell_index, others in cell_faces.items():
            value = sum(face_integrals[(cell_index, other)] for other in others)
            if not np.isfinite(value):
                raise ValueError(
                    f"Non-finite edge Gaussian curvature for edge {edge_index}, "
                    f"ball {cell_index}: {value}"
                )
            gauss_values[cell_index] = float(value)
            gauss_cell_values += 1

        for cell_index, value in mean_values.items():
            if not np.isfinite(value):
                raise ValueError(
                    f"Non-finite edge mean curvature for edge {edge_index}, "
                    f"ball {cell_index}: {value}"
                )

        mean_contributions.append(mean_values)
        gaussian_contributions.append(gauss_values)
        face_values += len(face_integrals)
        reduction_time += perf_counter() - t

        if count == n_edges or count % 100 == 0:
            net.update_progress(
                f"Calculating edge curvature: {count:,} / {n_edges:,}",
                100.0 * count / max(n_edges, 1),
            )

    total_time = perf_counter() - total_start
    other_time = max(
        total_time - cache_time - topology_time - integration_time - reduction_time,
        0.0,
    )

    result = {
        "key": key,
        "mean": mean_contributions,
        "gaussian": gaussian_contributions,
    }
    setattr(net, _CACHE_ATTR, result)

    if verbose:
        print("\n" + "=" * 70)
        print("FUSED EDGE CURVATURE TIMING")
        print("=" * 70)
        cache_label = "Shared edge cache build" if cache_built else "Shared edge cache reuse"
        for label, value in (
            (cache_label, cache_time),
            ("Surface / topology lookup", topology_time),
            ("Fused M + G integration", integration_time),
            ("Storage / G reduction", reduction_time),
            ("Other / loop overhead", other_time),
        ):
            pct = 100.0 * value / total_time if total_time > 0 else 0.0
            print(f"{label:<30} {value:10.4f} s  {pct:6.2f} %")
        print("-" * 70)
        print(f"{'TOTAL':<30} {total_time:10.4f} s  100.00 %")
        print()
        print(f"Edges:                  {n_edges:,}")
        print(f"Quadrature order:       {quadrature_order}")
        print(f"Edge quadrature points: {n_edges * quadrature_order:,}")
        print(f"Mean cell values:       {3 * n_edges:,}")
        print(f"Gaussian face values:   {face_values:,}")
        print(f"Gaussian cell values:   {gauss_cell_values:,}")
        print("=" * 70)

        kernel_total = sum(kernel_timing.values())
        print("\n" + "=" * 70)
        print("FUSED EDGE KERNEL TIMING")
        print("=" * 70)
        for label, value in (
            ("Boundary orientation", kernel_timing["orientation"]),
            ("Point + tangent", kernel_timing["point_tangent"]),
            ("Second derivative", kernel_timing["second_derivative"]),
            ("Radial vectors / normalize", kernel_timing["radial"]),
            ("Pairwise normal construction", kernel_timing["normals"]),
            ("Mean-curvature accumulation", kernel_timing["mean"]),
            ("Gaussian accumulation", kernel_timing["gaussian"]),
            ("Kernel loop / other", kernel_timing["other"]),
        ):
            pct = 100.0 * value / kernel_total if kernel_total > 0 else 0.0
            print(f"{label:<30} {value:10.4f} s  {pct:6.2f} %")
        print("-" * 70)
        print(f"{'KERNEL TOTAL':<30} {kernel_total:10.4f} s  100.00 %")
        print("=" * 70)

    return mean_contributions, gaussian_contributions
