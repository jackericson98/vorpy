import time
from itertools import combinations
import numpy as np
from numpy import pi, sqrt
from time import perf_counter as now
from vorpy.src.calculations import calc_sphericity
from vorpy.src.calculations import calc_isoperimetric_quotient
from vorpy.src.calculations import calc_contacts_cached
from vorpy.src.calculations import calc_cell_point_properties_cached
from vorpy.src.calculations import calc_cell_mass_properties_cached


def _row_ids(row, field):
    """Read an incidence list without changing the producer's tables."""
    value = row.get(field, ())
    if value is None:
        return ()
    try:
        return tuple(int(item) for item in value)
    except (TypeError, ValueError, OverflowError):
        return ()


def _weighted_cell_boundedness(net, ball_id, *, tolerance=1e-9):
    """Check AW-cell boundedness from the exact weighted half-space recession cone.

    For generator ``i``, the additive-weighted Voronoi cell inequalities have
    recession constraints ``2(p_j-p_i) dot d <= 0``.  A non-zero bounded-box
    solution is therefore a proof that the cell is unbounded.  This diagnostic
    is deliberately AW-only; it does not infer boundedness for Power or
    Primitive networks.
    """
    settings = getattr(net, "settings", None) or {}
    if str(settings.get("net_type", "aw")).lower() != "aw":
        return None, "unsupported_scheme"

    try:
        from scipy.optimize import linprog
    except ImportError:
        return None, "boundedness_solver_unavailable"

    balls = net.balls
    try:
        generator_ids = tuple(int(value) for value in balls["num"])
        location_by_id = {
            generator_id: np.asarray(location, dtype=float)
            for generator_id, location in zip(generator_ids, balls["loc"])
        }
        center = location_by_id[int(ball_id)]
        locations = np.asarray(list(location_by_id.values()), dtype=float)
    except (KeyError, TypeError, ValueError):
        return None, "invalid_generator_geometry"
    if center.shape != (3,) or locations.ndim != 2 or locations.shape[1] != 3:
        return None, "invalid_generator_geometry"
    if not np.isfinite(center).all() or not np.isfinite(locations).all():
        return None, "invalid_generator_geometry"

    rows = [
        2.0 * (location - center)
        for generator_id, location in location_by_id.items()
        if generator_id != int(ball_id)
    ]
    if not rows:
        return False, "unbounded_recession_cone"
    matrix = np.asarray(rows, dtype=float)
    zero = np.zeros(len(matrix), dtype=float)
    bounds = [(-1.0, 1.0)] * 3
    for axis in range(3):
        direction = np.zeros(3, dtype=float)
        direction[axis] = 1.0
        maximum = linprog(
            -direction,
            A_ub=matrix,
            b_ub=zero,
            bounds=bounds,
            method="highs",
        )
        minimum = linprog(
            direction,
            A_ub=matrix,
            b_ub=zero,
            bounds=bounds,
            method="highs",
        )
        if not maximum.success or not minimum.success:
            return None, "boundedness_solver_failed"
        if -float(maximum.fun) > tolerance or float(minimum.fun) < -tolerance:
            return False, "unbounded_recession_cone"
    return True, "bounded_recession_cone"


def diagnose_geometry_completeness(net, *, tolerance=1e-9):
    """Return deterministic, read-only per-cell closure evidence.

    ``missing_faces`` contains generator-pair relations expected by an
    incident three-generator edge.  ``missing_edges`` contains generator-
    triple relations expected by an incident four-generator vertex.
    ``missing_vertices`` contains invalid edge endpoint references.  The
    ``closure_relations`` mapping reports one-sided edge/surface incidence and
    vertex-link counts.  AW boundedness is classified independently using the
    weighted recession cone, so an open cell is not mistaken for a producer
    omission.  The function never builds or mutates geometry.
    """
    if getattr(net, "balls", None) is None:
        raise ValueError("Geometry completeness diagnostics require Network.balls.")

    balls = net.balls
    edges = getattr(net, "edges", None)
    verts = getattr(net, "verts", None)
    surfs = getattr(net, "surfs", None)
    if edges is None or verts is None or surfs is None:
        raise ValueError("Geometry completeness diagnostics require edges, vertices, and surfaces.")

    ball_ids = tuple(int(value) for value in balls["num"])
    ball_rows = {
        int(row["num"]): row
        for _, row in balls.iterrows()
    }
    boundary_positions = set(
        int(value) for value in (getattr(net, "boundary_indices", ()) or ())
    )
    boundary_ids = {
        ball_ids[position]
        for position in boundary_positions
        if 0 <= position < len(ball_ids)
    }
    group_ids = tuple(int(value) for value in (getattr(net, "group", ()) or ()))
    surface_keys = {
        tuple(sorted(_row_ids(row, "balls")))
        for _, row in surfs.iterrows()
        if len(_row_ids(row, "balls")) == 2
    }
    edge_keys = {
        tuple(sorted(_row_ids(row, "balls")))
        for _, row in edges.iterrows()
        if len(_row_ids(row, "balls")) == 3
    }
    valid_edge_ids = {int(value) for value in edges.index}
    valid_vertex_ids = {int(value) for value in verts.index}
    valid_surface_ids = {int(value) for value in surfs.index}
    surface_by_edge = {int(edge_id): [] for edge_id in edges.index}
    for surface_id, row in surfs.iterrows():
        for edge_id in _row_ids(row, "edges"):
            if edge_id in surface_by_edge:
                surface_by_edge[edge_id].append(int(surface_id))

    unresolved_vertices = {}
    for entry in getattr(net, "_aw_unresolved_vertex_curvature", ()) or ():
        try:
            vertex_id, cell_id, reason = entry
            unresolved_vertices.setdefault(int(cell_id), []).append(
                (int(vertex_id), str(reason))
            )
        except (TypeError, ValueError):
            continue

    cells = []
    for cell_id in group_ids:
        if cell_id not in ball_rows:
            continue
        row = ball_rows[cell_id]
        cell_surface_ids = _row_ids(row, "surfs")
        cell_edge_ids = _row_ids(row, "edges")
        cell_vertex_ids = _row_ids(row, "verts")
        missing_faces = set()
        missing_edges = set()
        missing_vertices = set()
        open_edges = []
        bad_vertices = []

        for edge_id in cell_edge_ids:
            if edge_id not in valid_edge_ids:
                missing_edges.add((cell_id, edge_id))
                continue
            edge_row = edges.loc[edge_id]
            edge_balls = _row_ids(edge_row, "balls")
            for other_id in edge_balls:
                if other_id == cell_id:
                    continue
                pair = tuple(sorted((cell_id, other_id)))
                if pair not in surface_keys:
                    missing_faces.add(pair)
            incident_surfaces = tuple(
                surface_id
                for surface_id in surface_by_edge.get(edge_id, ())
                if cell_id in _row_ids(surfs.loc[surface_id], "balls")
            )
            if len(incident_surfaces) != 2:
                open_edges.append((edge_id, len(incident_surfaces)))
            endpoints = _row_ids(edge_row, "verts")
            if len(endpoints) != 2:
                missing_vertices.add((edge_id, -1))
            else:
                missing_vertices.update(
                    (edge_id, vertex_id)
                    for vertex_id in endpoints
                    if vertex_id not in valid_vertex_ids
                )

        for vertex_id in cell_vertex_ids:
            if vertex_id not in valid_vertex_ids:
                missing_vertices.add((-1, vertex_id))
                continue
            vertex_row = verts.loc[vertex_id]
            vertex_balls = _row_ids(vertex_row, "balls")
            for other_pair in combinations(
                (value for value in vertex_balls if value != cell_id), 2
            ):
                triple = tuple(sorted((cell_id,) + tuple(other_pair)))
                if triple not in edge_keys:
                    missing_edges.add(triple)
            vertex_edge_ids = _row_ids(vertex_row, "edges")
            cell_vertex_edges = tuple(
                edge_id
                for edge_id in vertex_edge_ids
                if edge_id in valid_edge_ids
                and cell_id in _row_ids(edges.loc[edge_id], "balls")
            )
            vertex_surface_ids = tuple(
                surface_id
                for surface_id in _row_ids(vertex_row, "surfs")
                if surface_id in valid_surface_ids
                and cell_id in _row_ids(surfs.loc[surface_id], "balls")
            )
            if len(cell_vertex_edges) != 3 or len(vertex_surface_ids) < 3:
                bad_vertices.append(
                    (vertex_id, len(cell_vertex_edges), len(vertex_surface_ids))
                )

        if len(cell_vertex_ids) < 3 or len(cell_edge_ids) < 4 or len(cell_surface_ids) < 3:
            minimum_feature_count = True
        else:
            minimum_feature_count = False

        bounded, bounded_reason = _weighted_cell_boundedness(
            net, cell_id, tolerance=tolerance
        )
        analytic_issues = tuple(sorted(unresolved_vertices.get(cell_id, ())))
        topology_issues = bool(
            minimum_feature_count
            or missing_faces
            or missing_edges
            or missing_vertices
            or open_edges
            or bad_vertices
        )
        if bool(row.get("complete", False)):
            classification = "COMPLETE"
        elif bounded is False:
            classification = "UNBOUNDED"
        elif topology_issues:
            classification = "INCOMPLETE_CONSTRUCTION"
        elif analytic_issues or bounded is None:
            classification = "UNSUPPORTED_ANALYTIC_GEOMETRY"
        else:
            classification = "OVERLY_CONSERVATIVE"

        stable_id = row.get("stable_id")
        if stable_id is None or str(stable_id) == "":
            stable_id = f"system:{cell_id}"
        system_id = row.get("system_num")
        try:
            system_id = None if system_id is None else int(system_id)
        except (TypeError, ValueError, OverflowError):
            system_id = None
        cells.append({
            "local_id": cell_id,
            "system_id": system_id,
            "stable_id": str(stable_id),
            "boundary_generator": cell_id in boundary_ids,
            "complete": bool(row.get("complete", False)),
            "classification": classification,
            "bounded": bounded,
            "bounded_reason": bounded_reason,
            "feature_counts": {
                "surfaces": len(cell_surface_ids),
                "edges": len(cell_edge_ids),
                "vertices": len(cell_vertex_ids),
            },
            "missing_faces": tuple(sorted(missing_faces)),
            "missing_edges": tuple(sorted(missing_edges)),
            "missing_vertices": tuple(sorted(missing_vertices)),
            "closure_relations": {
                "open_edges": tuple(sorted(open_edges)),
                "bad_vertices": tuple(sorted(bad_vertices)),
            },
            "unsupported_analytic_geometry": analytic_issues,
        })

    cells.sort(key=lambda item: item["local_id"])
    classifications = {}
    for cell in cells:
        key = cell["classification"]
        classifications[key] = classifications.get(key, 0) + 1
    return {
        "schema_version": 1,
        "scheme": str((getattr(net, "settings", None) or {}).get("net_type", "aw")),
        "counts": {
            "cells": len(cells),
            "complete": sum(cell["complete"] for cell in cells),
            "classifications": dict(sorted(classifications.items())),
        },
        "cells": tuple(cells),
    }


def _rmsd(values, average):
    """Return the RMS deviation of values from average."""
    if not values:
        return 0.0
    return sqrt(sum((value - average) ** 2 for value in values) / len(values))


def _print_analysis_timing(timer, total):
    """Print a detailed analysis timing breakdown."""
    print('\n' + '=' * 70)
    print('ANALYSIS TIMING')
    print('=' * 70)
    labels = [
        ('setup', 'Setup / cache construction'),
        ('surface_gather', 'Surface gathering'),
        ('completeness', 'Cell completeness'),
        ('basic', 'Surface area + volume'),
        ('curvs', 'Curvature'),
        ('geometric', 'Geometric metrics'),
        ('neighbors_1', 'First neighbors'),
        ('neighbors_2', 'Second neighbors'),
        ('spikes', 'Spikes + bounding box'),
        ('contacts', 'Contacts'),
        ('com', 'COM + moment of inertia'),
        ('moi', 'Moment of inertia only'),
        ('b_box', 'Bounding box only'),
        ('surface_assign', 'Surface assignment'),
        ('ball_assign', 'Ball assignment'),
    ]
    for key, label in labels:
        elapsed = timer.get(key, 0.0)
        pct = 100.0 * elapsed / total if total > 0 else 0.0
        print(f'{label:<28} {elapsed:10.4f} s  {pct:6.2f} %')
    measured = sum((timer.get(key, 0.0) for key, _ in labels))
    other = max(total - measured, 0.0)
    pct = 100.0 * other / total if total > 0 else 0.0
    print(f"{'Other / loop overhead':<28} {other:10.4f} s  {pct:6.2f} %")
    print('-' * 70)
    print(f"{'TOTAL':<28} {total:10.4f} s  100.00 %")

def analyze_geometry_only(net):
    """Populate only cell completeness needed by dual/interface selection."""
    ball_index = net.balls.index.to_numpy()
    ball_verts = net.balls['verts'].to_numpy()
    ball_edges = net.balls['edges'].to_numpy()
    ball_surfs = net.balls['surfs'].to_numpy()
    edge_ball_sets = [frozenset(balls) for balls in net.edges['balls'].to_numpy()]
    vert_edges = net.verts['edges'].to_numpy()
    vert_balls = net.verts['balls'].to_numpy()
    num_to_pos = {int(num): pos for pos, num in enumerate(ball_index)}
    complete = [False] * len(net.balls)

    for ball_num in net.group:
        ball_num = int(ball_num)
        ball_pos = num_to_pos.get(ball_num)
        if ball_pos is None:
            continue
        verts = ball_verts[ball_pos]
        edges = ball_edges[ball_pos]
        surfs = ball_surfs[ball_pos]
        if len(verts) < 3 or len(edges) < 4 or len(surfs) < 3:
            continue
        valid = True
        for vertex in verts:
            vertex = int(vertex)
            owning = sum(
                ball_num in edge_ball_sets[int(edge)]
                for edge in vert_edges[vertex]
            )
            if owning != 3:
                vertex_balls = set(vert_balls[vertex])
                if sum(
                    edge_ball_sets[int(edge)].issubset(vertex_balls)
                    for edge in edges
                ) < 3:
                    valid = False
                    break
        complete[ball_pos] = valid

    net.balls['complete'] = complete


def analyze(
    net,
    complicated=True,
    spikes=None,
    contacts=None,
    second_neighbors=None,
    com=None,
    moi=None,
    bounding_box=None,
):
    """Analyze cells in ``net.group`` and assemble geometric cell properties.

    For AW networks:
        M = M_surface + M_edge
        G = G_surface + G_edge + G_vertex

    Canonical ``int_mean_curv`` and ``int_gauss_curv`` store the complete
    piecewise-smooth cell measures. Component fields are retained for
    validation, interpretation, and regression testing.
    """
    analysis_start = now()
    net.update_progress("Analysis | Setup", 0.0)

    spikes = complicated if spikes is None else spikes
    contacts = complicated if contacts is None else contacts
    second_neighbors = complicated if second_neighbors is None else second_neighbors
    com = complicated if com is None else com
    moi = complicated if moi is None else moi
    bounding_box = complicated if bounding_box is None else bounding_box

    timer = {
        'setup': 0.0, 'surface_gather': 0.0, 'completeness': 0.0,
        'basic': 0.0, 'curvs': 0.0, 'geometric': 0.0,
        'neighbors_1': 0.0, 'neighbors_2': 0.0, 'spikes': 0.0,
        'contacts': 0.0, 'com': 0.0, 'moi': 0.0, 'b_box': 0.0,
        'surface_assign': 0.0, 'ball_assign': 0.0,
    }

    setup_start = now()

    # Cache DataFrame columns once; the hot loop should not repeatedly slice DataFrames.
    n_balls = len(net.balls)
    n_surfs = len(net.surfs)

    ball_nums = net.balls['num'].to_numpy()
    ball_surfs_all = net.balls['surfs'].to_numpy()
    ball_verts_all = net.balls['verts'].to_numpy()
    ball_edges_all = net.balls['edges'].to_numpy()
    ball_locs = net.balls['loc'].to_numpy()
    ball_locs_matrix = np.asarray([np.asarray(loc, dtype=float) for loc in ball_locs])
    ball_rads = net.balls['rad'].to_numpy(dtype=float)
    ball_index = net.balls.index.to_numpy()

    surf_sa = net.surfs['sa'].to_numpy()
    surf_vols = net.surfs['vols'].to_numpy()
    surf_balls = net.surfs['balls'].to_numpy()
    surf_mean_curv = net.surfs['mean_curv'].to_numpy()
    surf_gauss_curv = net.surfs['gauss_curv'].to_numpy()
    surf_int_mean = net.surfs['int_mean_curv'].to_numpy()
    surf_int_mean_sq = net.surfs['int_mean_curv_sq'].to_numpy()
    surf_int_gauss = net.surfs['int_gauss_curv'].to_numpy()

    settings = getattr(net, 'settings', None) or {}
    is_aw = settings.get('net_type', 'aw') == 'aw'

    # AW uses cell-relative surface M, edge M, and vertex G caches. Missing
    # production caches are errors rather than silently dropping curvature.
    # AW complete curvature uses surface/edge M and edge/vertex G caches.
    if is_aw:
        required = (
            (net.surfs, 'int_mean_curv_by_ball', 'surface mean-curvature'),
            (net.edges, 'int_mean_curv_by_ball', 'edge mean-curvature'),
            (net.edges, 'int_gauss_curv_by_ball', 'edge Gaussian-curvature'),
            (net.verts, 'int_gauss_curv_by_ball', 'vertex Gaussian-curvature'),
        )

        for table, column, description in required:
            if column not in table.columns:
                raise ValueError(
                    f"AW analysis requires {description} cache '{column}'."
                )

        surf_int_mean_by_ball = net.surfs['int_mean_curv_by_ball'].to_numpy()
        edge_mean_curv_by_ball = net.edges['int_mean_curv_by_ball'].to_numpy()
        edge_gauss_curv_by_ball = net.edges['int_gauss_curv_by_ball'].to_numpy()
        vert_gauss_curv_by_ball = net.verts['int_gauss_curv_by_ball'].to_numpy()

    else:
        surf_int_mean_by_ball = None
        edge_mean_curv_by_ball = None
        edge_gauss_curv_by_ball = None
        vert_gauss_curv_by_ball = None

    # Normalize mesh arrays once for compiled geometry helpers.
    raw_surf_points = net.surfs['points'].to_numpy()
    raw_surf_tris = net.surfs['tris'].to_numpy()
    surf_points = np.empty(n_surfs, dtype=object)
    surf_tris = np.empty(n_surfs, dtype=object)

    for surf_id in range(n_surfs):
        points = np.asarray(raw_surf_points[surf_id], dtype=np.float64)
        if points.size == 0:
            points = np.empty((0, 3), dtype=np.float64)
        elif points.ndim != 2:
            points = points.reshape((-1, 3))
        surf_points[surf_id] = np.ascontiguousarray(points)

        tris = np.asarray(raw_surf_tris[surf_id], dtype=np.int64)
        if tris.size == 0:
            tris = np.empty((0, 3), dtype=np.int64)
        elif tris.ndim != 2:
            tris = tris.reshape((-1, 3))
        surf_tris[surf_id] = np.ascontiguousarray(tris)

    edge_balls = net.edges['balls'].to_numpy()
    vert_edges = net.verts['edges'].to_numpy()
    vert_balls = net.verts['balls'].to_numpy()

    num_to_pos = {int(num): pos for pos, num in enumerate(ball_nums)}
    group_nums = [int(num) for num in net.group if int(num) in num_to_pos]
    group_set = set(group_nums)
    n_group = len(group_nums)
    edge_ball_sets = [frozenset(balls) for balls in edge_balls]

    # Build first-neighbor adjacency once from solved pairwise surfaces.
    ball_neighbors = [[] for _ in range(n_balls)]
    neighbor_sets = [set() for _ in range(n_balls)]

    for surf_id, balls in enumerate(surf_balls):
        if balls is None or len(balls) < 2:
            continue

        for i, ball_num in enumerate(balls):
            ball_num = int(ball_num)
            ball_pos = num_to_pos.get(ball_num)
            if ball_pos is None:
                continue

            for j, neighbor_num in enumerate(balls):
                if i == j:
                    continue
                neighbor_num = int(neighbor_num)
                if neighbor_num not in num_to_pos:
                    continue
                ball_neighbors[ball_pos].append((neighbor_num, surf_id))
                neighbor_sets[ball_pos].add(neighbor_num)

    # Full-table outputs. Balls outside net.group retain zero/default values.
    b_vols = [0.0] * n_balls
    b_sas = [0.0] * n_balls
    b_cell = [0] * n_balls

    b_max_mean_curvs = [0.0] * n_balls
    b_avg_mean_surf_curvs = [0.0] * n_balls
    b_max_gauss_curvs = [0.0] * n_balls
    b_avg_gauss_surf_curvs = [0.0] * n_balls

    # Complete integrated curvature measures plus explicit decompositions.
    b_int_mean_curvs = [0.0] * n_balls
    b_int_mean_curv_surfaces = [0.0] * n_balls
    b_int_mean_curv_edges = [0.0] * n_balls
    b_int_mean_curv_totals = [0.0] * n_balls
    b_int_mean_curv_sqs = [0.0] * n_balls

    b_int_gauss_curvs = [0.0] * n_balls
    b_int_gauss_curv_surfaces = [0.0] * n_balls
    b_int_gauss_curv_edges = [0.0] * n_balls
    b_int_gauss_curv_vertices = [0.0] * n_balls
    b_int_gauss_curv_totals = [0.0] * n_balls

    b_sphrctys = [0.0] * n_balls
    b_isopmqs = [0.0] * n_balls

    num_nbors = [0] * n_balls
    near_nbors = [0] * n_balls
    near_nbor_dists = [0.0] * n_balls
    nbor_lyr_rmsds = [0] * n_balls
    nbor_dst_avgs = [0] * n_balls
    b_inner = [0] * n_balls

    b_min_spikes = [0.0] * n_balls
    b_max_spikes = [0.0] * n_balls

    contact_areas = [0.0] * n_balls
    non_olap_vols = [0.0] * n_balls
    olap_vols = [0.0] * n_balls
    num_olaps = [0] * n_balls

    coms = [0.0] * n_balls
    mois = [0.0] * n_balls
    b_boxs = [0.0] * n_balls

    surface_contact_area = np.zeros(n_surfs, dtype=float)
    surface_overlap = np.zeros(n_surfs, dtype=float)

    timer['setup'] += now() - setup_start
    count = 0
    last_update = now()

    # ==================================================================
    # HOT LOOP: only requested cells in net.group
    # ==================================================================
    for ball_num in group_nums:
        ball_pos = num_to_pos[ball_num]
        surf_ids = ball_surfs_all[ball_pos]

        if surf_ids is None or len(surf_ids) == 0:
            count += 1
            continue

        # --------------------------------------------------------------
        # Surface quantities
        # --------------------------------------------------------------
        t = now()

        sa = 0.0
        volume = 0.0
        max_mean_curv = -float('inf')
        max_gauss_curv = -float('inf')
        int_mean_curv_surface = 0.0
        int_mean_curv_sq = 0.0
        int_gauss_curv_surface = 0.0

        for surf_id in surf_ids:
            surf_id = int(surf_id)
            sa += surf_sa[surf_id]

            vols_i = surf_vols[surf_id]
            try:
                volume += vols_i[ball_num]
            except (IndexError, KeyError, TypeError):
                volume += vols_i[ball_pos]

            max_mean_curv = max(max_mean_curv, surf_mean_curv[surf_id])
            max_gauss_curv = max(max_gauss_curv, surf_gauss_curv[surf_id])

            # H changes sign with orientation, so AW uses the cell-relative cache.
            if is_aw:
                surface_values = surf_int_mean_by_ball[surf_id]
                if not isinstance(surface_values, dict):
                    raise ValueError(
                        f"Invalid AW surface mean-curvature cache for surface {surf_id}: expected dict."
                    )
                if ball_num not in surface_values:
                    raise ValueError(
                        f"Missing AW surface mean-curvature contribution for ball {ball_num}, surface {surf_id}."
                    )
                surface_value = float(surface_values[ball_num])
                if not np.isfinite(surface_value):
                    raise ValueError(
                        f"Non-finite AW surface mean-curvature contribution for ball {ball_num}, "
                        f"surface {surf_id}: {surface_value}"
                    )
                int_mean_curv_surface += surface_value
            else:
                int_mean_curv_surface += surf_int_mean[surf_id]

            # H^2 and K are invariant under normal reversal.
            int_mean_curv_sq += surf_int_mean_sq[surf_id]
            int_gauss_curv_surface += surf_int_gauss[surf_id]

        timer['surface_gather'] += now() - t

        if sa == 0:
            count += 1
            continue

        # --------------------------------------------------------------
        # Cell completeness
        # --------------------------------------------------------------
        t = now()

        ball_verts = ball_verts_all[ball_pos]
        ball_edges = ball_edges_all[ball_pos]
        complete = True
        ball_index_value = ball_index[ball_pos]
        unresolved_vertices = getattr(net, "_aw_unresolved_vertex_curvature", set())
        if any(int(cell) == int(ball_num) for _, cell, _ in unresolved_vertices):
            complete = False

        for vert in ball_verts:
            vert = int(vert)
            owning_edge_count = sum(
                ball_index_value in edge_ball_sets[int(edge)]
                for edge in vert_edges[vert]
            )

            if owning_edge_count != 3:
                vert_ball_set = set(vert_balls[vert])
                new_count = sum(
                    edge_ball_sets[int(edge)].issubset(vert_ball_set)
                    for edge in ball_edges
                )
                if new_count < 3:
                    complete = False
                    break

        if len(ball_verts) < 3 or len(ball_edges) < 4 or len(surf_ids) < 3:
            complete = False

        b_cell[ball_pos] = complete
        timer['completeness'] += now() - t

        # --------------------------------------------------------------
        # Basic geometry
        # --------------------------------------------------------------
        t = now()
        b_sas[ball_pos] = sa
        b_vols[ball_pos] = volume
        timer['basic'] += now() - t

        # --------------------------------------------------------------
        # Complete piecewise-smooth curvature measures
        # --------------------------------------------------------------
        t = now()

        int_mean_curv_edge = 0.0
        int_gauss_curv_edge = 0.0
        int_gauss_curv_vertex = 0.0

        if is_aw:
            # Edge contributions: M turning and intrinsic G boundary curvature.
            for edge_id in ball_edges:
                edge_id = int(edge_id)

                mean_values = edge_mean_curv_by_ball[edge_id]
                gauss_values = edge_gauss_curv_by_ball[edge_id]

                if not isinstance(mean_values, dict):
                    raise ValueError(
                        f"Invalid AW edge mean-curvature cache for edge {edge_id}."
                    )
                if not isinstance(gauss_values, dict):
                    raise ValueError(
                        f"Invalid AW edge Gaussian-curvature cache for edge {edge_id}."
                    )
                if ball_num not in mean_values:
                    raise ValueError(
                        f"Missing AW edge mean-curvature contribution for "
                        f"ball {ball_num}, edge {edge_id}."
                    )
                if ball_num not in gauss_values:
                    if complete:
                        raise ValueError(
                            f"Missing AW edge Gaussian-curvature contribution for "
                            f"complete ball {ball_num}, edge {edge_id}."
                        )
                    continue

                mean_value = float(mean_values[ball_num])
                gauss_value = float(gauss_values[ball_num])

                if not np.isfinite(mean_value):
                    raise ValueError(
                        f"Non-finite AW edge mean curvature for "
                        f"ball {ball_num}, edge {edge_id}: {mean_value}"
                    )
                if not np.isfinite(gauss_value):
                    raise ValueError(
                        f"Non-finite AW edge Gaussian curvature for "
                        f"ball {ball_num}, edge {edge_id}: {gauss_value}"
                    )

                int_mean_curv_edge += mean_value
                int_gauss_curv_edge += gauss_value

            # Intrinsic Gaussian angular defects at vertices.
            for vert_id in ball_verts:
                vert_id = int(vert_id)
                vertex_values = vert_gauss_curv_by_ball[vert_id]

                if not isinstance(vertex_values, dict):
                    raise ValueError(
                        f"Invalid AW vertex Gaussian-curvature cache for vertex {vert_id}."
                    )
                if ball_num not in vertex_values:
                    if not complete:
                        continue
                    raise ValueError(
                        f"Missing AW vertex Gaussian-curvature contribution for "
                        f"ball {ball_num}, vertex {vert_id}."
                    )

                value = float(vertex_values[ball_num])

                if not np.isfinite(value):
                    raise ValueError(
                        f"Non-finite AW vertex Gaussian curvature for "
                        f"ball {ball_num}, vertex {vert_id}: {value}"
                    )

                int_gauss_curv_vertex += value

        int_mean_curv_total = int_mean_curv_surface + int_mean_curv_edge
        int_gauss_curv_total = (
                int_gauss_curv_surface
                + int_gauss_curv_edge
                + int_gauss_curv_vertex
        )

        b_max_mean_curvs[ball_pos] = max_mean_curv
        b_max_gauss_curvs[ball_pos] = max_gauss_curv

        # Canonical fields store the complete cell measures.
        b_int_mean_curvs[ball_pos] = int_mean_curv_total
        b_int_gauss_curvs[ball_pos] = int_gauss_curv_total

        b_int_mean_curv_surfaces[ball_pos] = int_mean_curv_surface
        b_int_mean_curv_edges[ball_pos] = int_mean_curv_edge
        b_int_mean_curv_totals[ball_pos] = int_mean_curv_total
        b_int_mean_curv_sqs[ball_pos] = int_mean_curv_sq

        b_int_gauss_curv_surfaces[ball_pos] = int_gauss_curv_surface
        b_int_gauss_curv_edges[ball_pos] = int_gauss_curv_edge
        b_int_gauss_curv_vertices[ball_pos] = int_gauss_curv_vertex
        b_int_gauss_curv_totals[ball_pos] = int_gauss_curv_total

        # These remain smooth-surface averages; singular terms are not area densities.
        b_avg_mean_surf_curvs[ball_pos] = int_mean_curv_surface / sa
        b_avg_gauss_surf_curvs[ball_pos] = int_gauss_curv_surface / sa

        timer['curvs'] += now() - t

        # --------------------------------------------------------------
        # Shape metrics
        # --------------------------------------------------------------
        t = now()
        b_sphrctys[ball_pos] = calc_sphericity(volume=volume, surface_area=sa)
        b_isopmqs[ball_pos] = calc_isoperimetric_quotient(volume=volume, surface_area=sa)
        timer['geometric'] += now() - t

        # --------------------------------------------------------------
        # First neighbors
        # --------------------------------------------------------------
        t = now()

        ball_loc = ball_locs[ball_pos]
        ball_rad = ball_rads[ball_pos]
        adjacency = ball_neighbors[ball_pos]
        neighbors_nums = []
        neighbor_dists = []

        for neighbor_num, surf_id in adjacency:
            neighbor_pos = num_to_pos[neighbor_num]
            delta = ball_locs_matrix[neighbor_pos] - ball_locs_matrix[ball_pos]
            center_dist = sqrt(float(np.dot(delta, delta)))
            neighbor_dist = center_dist - ball_rad - ball_rads[neighbor_pos]

            neighbors_nums.append(neighbor_num)
            neighbor_dists.append(neighbor_dist)

            overlap_dist = max(-neighbor_dist, 0.0)
            if overlap_dist > surface_overlap[surf_id]:
                surface_overlap[surf_id] = overlap_dist

        if neighbor_dists:
            b_inner[ball_pos] = group_set.issuperset(neighbors_nums)
            num_nbors[ball_pos] = len(neighbors_nums)

            min_index = min(range(len(neighbor_dists)), key=neighbor_dists.__getitem__)
            near_nbor_dists[ball_pos] = neighbor_dists[min_index]
            near_nbors[ball_pos] = neighbors_nums[min_index]

            nbor_dist_avg = sum(neighbor_dists) / len(neighbor_dists)
            nbor_dst_avgs[ball_pos] = [nbor_dist_avg]
            nbor_lyr_rmsds[ball_pos] = [_rmsd(neighbor_dists, nbor_dist_avg)]
        else:
            b_inner[ball_pos] = True
            nbor_dst_avgs[ball_pos] = [0.0]
            nbor_lyr_rmsds[ball_pos] = [0.0]

        timer['neighbors_1'] += now() - t

        # --------------------------------------------------------------
        # Point properties
        # --------------------------------------------------------------
        if spikes or bounding_box:
            t = now()
            min_spike, max_spike, box = calc_cell_point_properties_cached(
                ball_loc, surf_ids, surf_points
            )

            if spikes:
                b_min_spikes[ball_pos] = min_spike
                b_max_spikes[ball_pos] = max_spike
            if bounding_box:
                b_boxs[ball_pos] = box

            timer['spikes'] += now() - t

        # --------------------------------------------------------------
        # Contacts / overlap volume
        # --------------------------------------------------------------
        if contacts:
            t = now()
            contact_area, vdw_vol = calc_contacts_cached(
                ball_loc, ball_rad, surf_ids, surf_points, surf_tris
            )

            num_olaps[ball_pos] = sum(dist < 0 for dist in neighbor_dists)
            contact_areas[ball_pos] = sum(contact_area.values())

            for surf_id in surf_ids:
                surf_id = int(surf_id)
                surface_contact_area[surf_id] = contact_area[surf_id]

            non_olap_vols[ball_pos] = vdw_vol
            olap_vols[ball_pos] = (4.0 / 3.0) * pi * ball_rad ** 3 - vdw_vol
            timer['contacts'] += now() - t

        # --------------------------------------------------------------
        # Second neighbors
        # --------------------------------------------------------------
        if second_neighbors:
            t = now()
            first_set = neighbor_sets[ball_pos]
            previous = set(first_set)
            previous.add(ball_num)

            second_set = set()
            for first_num in previous:
                first_pos = num_to_pos.get(first_num)
                if first_pos is not None:
                    second_set.update(neighbor_sets[first_pos])

            second_set.difference_update(previous)

            if second_set:
                ball_loc_array = ball_locs_matrix[ball_pos]
                layer2_dists = []

                for ball_2 in second_set:
                    pos_2 = num_to_pos[ball_2]
                    delta = ball_locs_matrix[pos_2] - ball_loc_array
                    layer2_dists.append(sqrt(float(np.dot(delta, delta))))

                lyr2_dist_avg = sum(layer2_dists) / len(layer2_dists)
                nbor_dst_avgs[ball_pos].append(lyr2_dist_avg)
                nbor_lyr_rmsds[ball_pos].append(_rmsd(layer2_dists, lyr2_dist_avg))
            else:
                nbor_dst_avgs[ball_pos].append(0.0)
                nbor_lyr_rmsds[ball_pos].append(0.0)

            timer['neighbors_2'] += now() - t

        # --------------------------------------------------------------
        # Center of mass / moment of inertia
        # --------------------------------------------------------------
        if com and moi:
            t = now()
            com_val, moi_val = calc_cell_mass_properties_cached(
                ball_loc, surf_ids, surf_points, surf_tris, volume
            )
            coms[ball_pos] = com_val
            mois[ball_pos] = moi_val
            timer['com'] += now() - t

        elif com:
            t = now()
            com_val, _ = calc_cell_mass_properties_cached(
                ball_loc, surf_ids, surf_points, surf_tris, volume
            )
            coms[ball_pos] = com_val
            timer['com'] += now() - t

        elif moi:
            t = now()
            _, moi_val = calc_cell_mass_properties_cached(
                ball_loc, surf_ids, surf_points, surf_tris, volume
            )
            mois[ball_pos] = moi_val
            timer['moi'] += now() - t

        count += 1
        current_time = now()

        if current_time - last_update >= 0.25 or count == n_group:
            net.update_progress("Analysis", 100.0 * count / max(n_group, 1))
            last_update = current_time

    # ------------------------------------------------------------------
    # Bulk assignments
    # ------------------------------------------------------------------
    t = now()

    net.balls = net.balls.assign(
        vol=b_vols,
        sa=b_sas,
        complete=b_cell,
        max_mean_curv=b_max_mean_curvs,
        max_gauss_curv=b_max_gauss_curvs,
        avg_mean_surf_curv=b_avg_mean_surf_curvs,
        avg_gauss_surf_curv=b_avg_gauss_surf_curvs,

        int_mean_curv=b_int_mean_curvs,
        int_mean_curv_surface=b_int_mean_curv_surfaces,
        int_mean_curv_edge=b_int_mean_curv_edges,
        int_mean_curv_total=b_int_mean_curv_totals,
        int_mean_curv_sq=b_int_mean_curv_sqs,

        int_gauss_curv=b_int_gauss_curvs,
        int_gauss_curv_surface=b_int_gauss_curv_surfaces,
        int_gauss_curv_edge=b_int_gauss_curv_edges,
        int_gauss_curv_vertex=b_int_gauss_curv_vertices,
        int_gauss_curv_total=b_int_gauss_curv_totals,

        sphericity=b_sphrctys,
        isometric_quotient=b_isopmqs,
        ball_inside=b_inner,
        number_of_neighbors=num_nbors,
        nearest_neighbor=near_nbors,
        nearest_neighbor_distance=near_nbor_dists,
        neighbor_distance_average=nbor_dst_avgs,
        neighbor_distance_rmsd=nbor_lyr_rmsds,
        number_of_olaps=num_olaps,
        min_spike=b_min_spikes,
        max_spike=b_max_spikes,
        contact_area=contact_areas,
        olap_vol=olap_vols,
        vdw_vol=non_olap_vols,
        com=coms,
        moi=mois,
        bounding_box=b_boxs,
    )

    timer['ball_assign'] += now() - t

    t = now()
    net.surfs = net.surfs.assign(contact_area=surface_contact_area, overlap=surface_overlap)
    timer['surface_assign'] += now() - t

    analysis_total = now() - analysis_start
    net.metrics['anal'] = (
        now() - net.metrics['start'] - net.metrics['surf']
        - net.metrics['con'] - net.metrics['vert']
    )

    net.analysis_timing = timer.copy()
    net.analysis_timing['total'] = analysis_total

    if net.settings.get('verbose', False):
        _print_analysis_timing(timer, analysis_total)
