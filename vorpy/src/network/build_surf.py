import time
import numpy as np
from vorpy.src.calculations import calc_surf_func
from vorpy.src.calculations import calc_surf_tri_curvs
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


def get_com(locs, rads, perimeter, surf_loc, surf_norm, func, flat, net_type='aw'):
    """Calculate a valid representative center point for a surface."""
    if net_type in {'del', 'pow'}:
        return calc_com(points=np.array(perimeter)), False

    if is_within(perimeter, surf_loc, surf_loc, surf_norm):
        return surf_loc, False

    true_com = calc_com(points=np.array(perimeter))
    my_com = calc_surf_point(locs, point=true_com, func=func)
    if my_com is not None:
        if is_within(perimeter, my_com, surf_loc, surf_norm):
            return my_com, False

    my_com = calc_surf_point(locs, point=calc_com(points=np.array(perimeter[::5])), func=func)
    if my_com is not None:
        if is_within(perimeter, my_com, surf_loc, surf_norm):
            return my_com, False

    min_dist, my_point = np.inf, None
    for point in perimeter:
        dist = calc_dist(point, true_com)
        if dist < min_dist:
            my_point, min_dist = point, dist
    return my_point, True


def project_to_hyperboloid(twoD_points, small_ball_loc, surf_func, plane_normal, plane_location):
    """Project triangulated plane points back onto the curved surface."""
    plane_points = unproject_to_3d(twoD_points, plane_location, plane_normal)

    new_points = []
    for point in plane_points:
        new_point = calc_surf_point_from_plane(point, plane_normal, surf_func, small_ball_loc)
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
        flat=flat, func=sfunc, net_type=net_type
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
    my_2d_points, surf_tris = triangulate_2D_Surface(flat_points, res=res, center=flat_loc)
    _record_timing(timing, 'triangulate', time.perf_counter() - stage_start)

    if not flat:
        # Project mesh onto hyperboloid
        stage_start = time.perf_counter()
        spoints = project_to_hyperboloid(my_2d_points, locs[0], sfunc, surf_norm, surf_loc)
        _record_timing(timing, 'project_hyperboloid', time.perf_counter() - stage_start)

        # Mean curvature
        stage_start = time.perf_counter()
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv = calc_surf_tri_curvs(
            sfunc, spoints, surf_tris, curvature_type='mean'
        )
        _record_timing(timing, 'mean_curvature', time.perf_counter() - stage_start)

        # Gaussian curvature
        stage_start = time.perf_counter()
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv = calc_surf_tri_curvs(
            sfunc, spoints, surf_tris, curvature_type='gauss'
        )
        _record_timing(timing, 'gauss_curvature', time.perf_counter() - stage_start)
    else:
        # Flat surfaces need only unprojection; curvature is zero.
        stage_start = time.perf_counter()
        spoints = unproject_to_3d(my_2d_points, surf_loc, surf_norm)
        _record_timing(timing, 'unproject_flat', time.perf_counter() - stage_start)

        stage_start = time.perf_counter()
        n_tris = len(surf_tris)
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv = [0 for _ in range(n_tris)], 0, 0
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv = [0 for _ in range(n_tris)], 0, 0
        _record_timing(timing, 'flat_curvature_init', time.perf_counter() - stage_start)

    _record_timing(timing, 'total', time.perf_counter() - total_start)

    return (
        spoints, surf_tris,
        mean_tri_curvs, mean_surf_curv, avg_mean_surf_curv,
        gauss_tri_curvs, gauss_surf_curv, avg_gauss_surf_curv,
        sfunc, surf_com, flat, surf_loc
    )