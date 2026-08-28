import time
import numpy as np
from shapely import Polygon, Point
from vorpy.src.calculations import calc_surf_func
from vorpy.src.calculations.curvature import calc_surf_tri_curvs_both
from vorpy.src.network.perimeter import build_perimeter
from vorpy.src.network.fill import calc_surf_point
from vorpy.src.network.fill import calc_surf_point_from_plane
from vorpy.src.calculations import calc_com
from vorpy.src.calculations import project_to_plane
from vorpy.src.calculations import calc_dist
from vorpy.src.calculations import unproject_to_3d
from vorpy.src.network.triangulate import triangulate_2D_Surface
from vorpy.src.network.triangulate import is_within


def _record_timing(timing, key, elapsed):
    """Add elapsed time to a timing dictionary when profiling is enabled."""
    if timing is not None:
        timing[key] = timing.get(key, 0.0) + elapsed


def get_com(locs, rads, perimeter, surf_loc, surf_norm, func, flat,
            net_type='aw', timing=None):
    """
    Calculate a valid representative center point for a surface.

    The historical fallback order and Shapely ``contains`` semantics are
    preserved. The projected perimeter polygon is built once per surface and
    reused for all containment queries.
    """
    # Normalize perimeter once for all COM calculations.
    t0 = time.perf_counter()
    perimeter_array = np.asarray(perimeter, dtype=float)
    if timing is not None:
        timing['com_setup'] = timing.get('com_setup', 0.0) + time.perf_counter() - t0

    # Historical fast path for flat network types.
    if net_type in {'del', 'pow'}:
        t0 = time.perf_counter()
        result = calc_com(points=perimeter_array)
        if timing is not None:
            timing['com_centroid'] = timing.get('com_centroid', 0.0) + time.perf_counter() - t0
            timing['com_flat_returns'] = timing.get('com_flat_returns', 0) + 1
        return result, False

    # --------------------------------------------------------------
    # Build the same 2D Shapely perimeter representation used by
    # triangulate.is_within(), but only once for this surface.
    # --------------------------------------------------------------
    t0 = time.perf_counter()
    try:
        if perimeter_array.ndim == 2 and perimeter_array.shape[1] == 3:
            perimeter_2d = project_to_plane(
                perimeter_array,
                surf_loc,
                surf_norm
            )
        else:
            perimeter_2d = perimeter_array

        com_polygon = Polygon(perimeter_2d)
    except TypeError:
        com_polygon = None

    if timing is not None:
        timing['com_contains_setup'] = (
            timing.get('com_contains_setup', 0.0)
            + time.perf_counter() - t0
        )

    def _contains(point):
        t_query = time.perf_counter()

        if com_polygon is None:
            result = False
        else:
            try:
                # calc_surf_point() uses np.roots(), which may return values
                # with a numerically-zero complex component. Handle those
                # explicitly instead of relying on a ComplexWarning-producing
                # cast to float.
                point_array = np.asarray(point)

                if np.iscomplexobj(point_array):
                    real_point = np.real_if_close(point_array, tol=1000)

                    # If a genuinely non-negligible imaginary component remains,
                    # the point is not a valid real surface point.
                    if np.iscomplexobj(real_point):
                        result = False
                        if timing is not None:
                            timing['com_complex_rejects'] = (
                                timing.get('com_complex_rejects', 0) + 1
                            )
                        return result

                    point_array = real_point

                point_array = np.asarray(point_array, dtype=float)

                if point_array.shape[0] == 3:
                    point_2d = project_to_plane(
                        [point_array],
                        surf_loc,
                        surf_norm
                    )[0]
                else:
                    point_2d = point_array

                result = com_polygon.contains(Point(point_2d))
            except (TypeError, ValueError):
                result = False

        if timing is not None:
            timing['com_contains_query'] = (
                timing.get('com_contains_query', 0.0)
                + time.perf_counter() - t_query
            )

        return result

    # 1) Prefer supplied surface location if inside perimeter.
    if _contains(surf_loc):
        if timing is not None:
            timing['com_surf_loc_returns'] = timing.get('com_surf_loc_returns', 0) + 1
        return surf_loc, False

    # 2) Project true perimeter centroid onto curved surface.
    t0 = time.perf_counter()
    true_com = calc_com(points=perimeter_array)
    if timing is not None:
        timing['com_centroid'] = timing.get('com_centroid', 0.0) + time.perf_counter() - t0

    t0 = time.perf_counter()
    my_com = calc_surf_point(locs, point=true_com, func=func)
    if timing is not None:
        timing['com_surface_projection'] = (
            timing.get('com_surface_projection', 0.0)
            + time.perf_counter() - t0
        )

    if my_com is not None and _contains(my_com):
        if timing is not None:
            timing['com_true_returns'] = timing.get('com_true_returns', 0) + 1
        return my_com, False

    # 3) Historical every-fifth-perimeter-point centroid.
    t0 = time.perf_counter()
    sampled_com = calc_com(points=perimeter_array[::5])
    if timing is not None:
        timing['com_centroid'] = timing.get('com_centroid', 0.0) + time.perf_counter() - t0

    t0 = time.perf_counter()
    my_com = calc_surf_point(locs, point=sampled_com, func=func)
    if timing is not None:
        timing['com_surface_projection'] = (
            timing.get('com_surface_projection', 0.0)
            + time.perf_counter() - t0
        )

    if my_com is not None and _contains(my_com):
        if timing is not None:
            timing['com_sample_returns'] = timing.get('com_sample_returns', 0) + 1
        return my_com, False

    # 4) Historical hard fallback: closest perimeter point to true centroid.
    t0 = time.perf_counter()
    deltas = perimeter_array - true_com
    dist_sq = np.einsum('ij,ij->i', deltas, deltas)
    my_point = perimeter_array[int(np.argmin(dist_sq))]
    if timing is not None:
        timing['com_nearest_fallback'] = (
            timing.get('com_nearest_fallback', 0.0)
            + time.perf_counter() - t0
        )
        timing['com_hard_returns'] = timing.get('com_hard_returns', 0) + 1

    return my_point, True


def project_to_hyperboloid(twoD_points, small_ball_loc, surf_func, plane_normal,
                           plane_location, timing=None):
    """Project triangulated plane points back onto the curved surface."""
    # Normalize values that are constant for the entire surface once rather
    # than reallocating them for every projected point.
    plane_normal = np.asarray(plane_normal, dtype=float)
    surf_func = np.asarray(surf_func, dtype=float)
    small_ball_loc = np.asarray(small_ball_loc, dtype=float)

    plane_points = np.asarray(
        unproject_to_3d(twoD_points, plane_location, plane_normal),
        dtype=float
    )

    new_points = []
    for point in plane_points:
        new_point = calc_surf_point_from_plane(
            point,
            plane_normal,
            surf_func,
            small_ball_loc,
            timing=timing,
        )
        if new_point is not None:
            new_points.append(new_point)
    return new_points


def build_surf(locs, rads, epnts, res, net_type, sfunc=None, perimeter=None,
               surf_loc=None, surf_norm=None, timing=None):
    """
    Build one surface and calculate its mesh and curvature descriptors.

    The optional ``timing`` dictionary is populated in-place with cumulative
    wall-clock timings for the major stages. The normal return structure is
    unchanged, so existing callers that do not pass ``timing`` continue to work.
    """
    total_start = time.perf_counter()

    # Surface function
    stage_start = time.perf_counter()
    if sfunc is None:
        sfunc = calc_surf_func(np.array(locs[0]), rads[0], np.array(locs[1]), rads[1])
    _record_timing(timing, 'surf_func', time.perf_counter() - stage_start)

    # Flatness check
    stage_start = time.perf_counter()
    flat = net_type in {'prm', 'pow'} or rads[0] == rads[1]
    _record_timing(timing, 'flat_check', time.perf_counter() - stage_start)

    # Perimeter construction
    stage_start = time.perf_counter()
    if perimeter is None:
        perimeter, surf_loc, surf_norm = build_perimeter(locs, rads, epnts=epnts, net_type=net_type)
    _record_timing(timing, 'perimeter', time.perf_counter() - stage_start)

    if surf_loc is None or surf_norm is None:
        _record_timing(timing, 'total', time.perf_counter() - total_start)
        return

    # Surface center
    stage_start = time.perf_counter()
    surf_com, filter_hard = get_com(
        locs, rads, perimeter=perimeter, surf_loc=surf_loc, surf_norm=surf_norm,
        flat=flat, func=sfunc, net_type=net_type, timing=timing
    )
    _record_timing(timing, 'get_com', time.perf_counter() - stage_start)

    # Project perimeter to plane
    stage_start = time.perf_counter()
    flat_points = project_to_plane(
        np.array(perimeter), plane_normal=surf_norm, plane_point=surf_loc
    )
    _record_timing(timing, 'project_perimeter', time.perf_counter() - stage_start)

    # Project COM / surface location to plane
    stage_start = time.perf_counter()
    flat_com, flat_loc = project_to_plane(
        np.array([surf_com, surf_loc]), plane_normal=surf_norm, plane_point=surf_loc
    )
    _record_timing(timing, 'project_com', time.perf_counter() - stage_start)

    # Triangulation
    stage_start = time.perf_counter()
    my_2d_points, surf_tris = triangulate_2D_Surface(flat_points, res=res, center=flat_loc, timing=timing)
    _record_timing(timing, 'triangulate', time.perf_counter() - stage_start)

    if not flat:
        # Project mesh onto hyperboloid
        stage_start = time.perf_counter()
        spoints = project_to_hyperboloid(my_2d_points, locs[0], sfunc, surf_norm, surf_loc, timing=timing)
        _record_timing(timing, 'project_hyperboloid', time.perf_counter() - stage_start)

        # Mean + Gaussian curvature in one shared triangle pass.
        stage_start = time.perf_counter()
        (
            mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv,
            gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv,
            int_mean_curv, int_mean_curv_sq, int_gauss_curv,
            curvature_area,
        ) = calc_surf_tri_curvs_both(sfunc, spoints, surf_tris)
        _record_timing(timing, 'combined_curvature', time.perf_counter() - stage_start)
    else:
        # Flat surfaces need only unprojection; curvature is zero.
        stage_start = time.perf_counter()
        spoints = unproject_to_3d(my_2d_points, surf_loc, surf_norm)
        _record_timing(timing, 'unproject_flat', time.perf_counter() - stage_start)

        stage_start = time.perf_counter()
        n_tris = len(surf_tris)
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv = [0 for _ in range(n_tris)], 0, 0
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv = [0 for _ in range(n_tris)], 0, 0
        int_mean_curv, int_mean_curv_sq, int_gauss_curv = 0.0, 0.0, 0.0
        curvature_area = None
        _record_timing(timing, 'flat_curvature_init', time.perf_counter() - stage_start)

    _record_timing(timing, 'total', time.perf_counter() - total_start)

    return (
        spoints, surf_tris,
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv,
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv,
        int_mean_curv, int_mean_curv_sq, int_gauss_curv,
        curvature_area,
        sfunc, surf_com, flat, surf_loc
    )
