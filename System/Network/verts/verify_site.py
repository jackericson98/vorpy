from System.sys_funcs.calcs.calcs import calc_dist
from numba import jit


@jit(nopython=True)
def verify_vor(loc, rad, test_locs, test_rads):
    # Go through the atoms in the overlap test atom list
    for i, aloc in enumerate(test_locs):
        # Get the atom's location and radius
        arad = test_rads[i]
        if calc_dist(aloc, loc) - arad < rad:
            return False
    return True


@jit(nopython=True)
def verify_del(loc, rad, test_locs):
    # Go through the atoms in the overlap test atom list
    for i, aloc in enumerate(test_locs):
        if calc_dist(aloc, loc) < rad:
            return False
    return True


@jit(nopython=True)
def verify_pow(loc, rad, test_locs, test_rads):
    # Go through the atoms in the overlap test atom list
    for i, aloc in enumerate(test_locs):
        # Get the atom's location and radius
        arad = test_rads[i]
        if calc_dist(aloc, loc) ** 2 - arad ** 2 < rad:
            return False
    return True


def verify_site(loc, rad, test_locs, test_rads, net_type='vor'):
    """
    Compares a vertex to the atoms around to see if they overlap, atoms pre-gathered
    """
    # Verification for a voronoi network
    if net_type == 'vor':
        return verify_vor(loc, rad, test_locs, test_rads)
    # Verify Delaunay
    elif net_type == 'del':
        return verify_del(loc, rad, test_locs)
    # Verify power network
    elif net_type == 'pow':
        return verify_pow(loc, rad, test_locs, test_rads)

