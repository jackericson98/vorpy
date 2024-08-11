from System.sys_funcs.calcs.calcs import calc_dist_numba
from numba import jit


@jit(nopython=True)
def verify_aw(loc, rad, test_locs, test_rads):
    # Go through the balls in the overlap test ball list
    for i, b_loc in enumerate(test_locs):
        # Get the ball's location and radius
        b_rad = test_rads[i]
        if calc_dist_numba(b_loc, loc) - b_rad < rad:
            return False
    return True


@jit(nopython=True)
def verify_del(loc, rad, test_locs):
    # Go through the balls in the overlap test ball list
    for i, b_loc in enumerate(test_locs):
        if calc_dist_numba(b_loc, loc) < rad:
            return False
    return True


@jit(nopython=True)
def verify_pow(loc, rad, test_locs, test_rads):
    # Go through the balls in the overlap test ball list
    for i, b_loc in enumerate(test_locs):
        # Get the ball's location and radius
        b_rad = test_rads[i]
        if calc_dist_numba(b_loc, loc) ** 2 - b_rad ** 2 < rad:
            return False
    return True


def verify_site(loc, rad, test_locs, test_rads, net_type='aw'):
    """
    Compares a vertex to the balls around to see if they overlap, balls pre-gathered
    """
    # Verification for a voronoi network
    if net_type == 'aw':
        return verify_aw(loc, rad, test_locs, test_rads)
    # Verify Delaunay
    elif net_type == 'del':
        return verify_del(loc, rad, test_locs)
    # Verify power network
    elif net_type == 'pow':
        return verify_pow(loc, rad, test_locs, test_rads)
