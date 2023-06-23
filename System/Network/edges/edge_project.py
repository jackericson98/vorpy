import numpy as np
from System.sys_funcs.calcs.calcs import calc_circ, calc_angle, calc_surf_func
from numba import jit


# Find projection values. Calculates the correct end and projection points for the edge
@jit(nopython=True)
def calc_edge_proj_pt(pv0, pv1, loc):
    # Get the projection point
    # Find the point in between the two vertex points
    r01 = pv1 - pv0  # Vector between vertices
    r_mag = np.linalg.norm(r01)  # Magnitude of the vector between the two vertex points
    rn01 = r01 / r_mag  # Normal to the vector between the vertices
    pc01 = pv0 + 0.5 * rn01 * r_mag  # Center point

    # Determine if the theoretical center of the edge is inside the vertices or not
    dr = 1
    if np.sqrt(sum(np.square(loc - pv0))) < r_mag or np.sqrt(sum(np.square(loc - pv1))) < r_mag:
        dr = -1

    # Find the vector normal to the projection plane
    p_norm = dr * np.cross(loc - pc01, pv1 - pc01)
    # Find the vector perpendicular to the plane's normal (i.e. in the plane) and the vector between vertices
    r_pcr = - np.cross(p_norm, rn01)
    rn_pcr = r_pcr / np.linalg.norm(r_pcr)
    # Calculate the reference point
    return pc01 + 0.5 * r_mag * rn_pcr


# Project method. Projects a point onto the surface using a reference point
@jit(nopython=True)
def edge_project(rn, pa, func, ep_1, ep_2=None):
    # Get the function values
    f = func
    # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
    a = f[0] * rn[0] ** 2 + f[1] * rn[1] ** 2 + f[2] * rn[2] ** 2 + f[3] * rn[0] * rn[1] + f[4] * rn[
        1] * rn[2] + f[5] * rn[2] * rn[0]
    b = 2 * f[0] * rn[0] * pa[0] + 2 * f[1] * rn[1] * pa[1] + 2 * f[2] * rn[2] * pa[2] + f[3] \
        * (rn[0] * pa[1] + rn[1] * pa[0]) + f[4] * (rn[1] * pa[2] + rn[2] * pa[1]) + f[5] \
        * (rn[2] * pa[0] + rn[0] * pa[2]) + f[6] * rn[0] + f[7] * rn[1] + f[8] * rn[2]
    c = f[0] * pa[0] ** 2 + f[1] * pa[1] ** 2 + f[2] * pa[2] ** 2 + f[3] * pa[0] * pa[1] + f[4] * pa[1] * pa[
        2] + f[5] * pa[2] * pa[0] + f[6] * pa[0] + f[7] * pa[1] + f[8] * pa[2] + f[9]
    # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
    # and add that point to our surface list of points
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        # Calculate the roots
        roots = np.roots(np.array([a, b, c]))
        # If one root exists return it
        if len(roots) == 1:
            return pa + roots[0] * rn
        else:
            p1 = pa + min(roots) * rn
            p2 = pa + max(roots) * rn
        # If the point we are calculating is the first in the edge choose the one closest to the vertex
        if ep_2 is None:
            point = p1
            if np.sqrt(sum(np.square(p2 - ep_1))) <= np.sqrt(sum(np.square(p1 - ep_1))):
                point = p2
        # If we have 2 points to choose from, choose the one that makes the angle closer to 180
        else:
            point = p1
            if calc_angle(ep_1, ep_2, p2) >= calc_angle(ep_1, ep_2, p1):
                point = p2
        # Return the point we choose
        return point