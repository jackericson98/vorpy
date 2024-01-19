from numba import jit
import numpy as np


@jit(nopython=True)
def calc_circ_coefs(l0, l1, l2, r0, r1, r2):
    # Move the other atoms to the location of the first
    x2, y2, z2 = l1[0] - l0[0], l1[1] - l0[1], l1[2] - l0[2]
    x3, y3, z3 = l2[0] - l0[0], l2[1] - l0[1], l2[2] - l0[2]
    # Calculate coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (r0 - r1), r0 ** 2 - r1 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (r0 - r2), r0 ** 2 - r2 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3 = y2 * z3 - z2 * y3, z2 * x3 - x2 * z3, x2 * y3 - y2 * x3
    abcs = [[a1, a1, a3], [b1, b2, b3], [c1, c2, c3]]
    # More coefficients
    F = a3 * b2 * c1 - a2 * b3 * c1 - a3 * b1 * c2 + a1 * b3 * c2 + a2 * b1 * c3 - a1 * b2 * c3
    Fx0 = b3 * c2 * f1 - b2 * c3 * f1 - b3 * c1 * f2 + b1 * c3 * f2
    Fx1 = b3 * c2 * d1 - b2 * c3 * d1 - b3 * c1 * d2 + b1 * c3 * d2
    Fy0 = - a3 * c2 * f1 + a2 * c3 * f1 + a3 * c1 * f2 - a1 * c3 * f2
    Fy1 = - a3 * c2 * d1 + a2 * c3 * d1 + a3 * c1 * d2 - a1 * c3 * d2
    Fz0 = a3 * b2 * f1 - a2 * b3 * f1 - a3 * b1 * f2 + a1 * b3 * f2
    Fz1 = a3 * b2 * d1 - a2 * b3 * d1 - a3 * b1 * d2 + a1 * b3 * d2
    Fs = F, Fx0, Fx1, Fy0, Fy1, Fz0, Fz1

    return Fs, abcs


@jit(nopython=True)
def calc_circ_abcs(Fs, r0):
    F, Fx0, Fx1, Fy0, Fy1, Fz0, Fz1 = Fs
    # Find the radius of the tangential circle using the quadratic formula
    a = (Fx1 ** 2 + Fy1 ** 2 + Fz1 ** 2) / F ** 2 - 1
    b = 2 * (Fx0 * Fx1 + Fy0 * Fy1 + Fz0 * Fz1) / F ** 2 - 2 * r0
    c = (Fx0 ** 2 + Fy0 ** 2 + Fz0 ** 2) / F ** 2 - r0 ** 2
    return a, b, c


def calc_circ(l0, l1, l2, r0, r1, r2):
    """
    Takes in 3 atoms, calculates the center and radius of inscribed circle
    :param : Locations and radii for the circle
    :return: Center and radius of the inscribed circle
    """
    # Make sure the locations are arrays
    l0, l1, l2 = np.array(l0), np.array(l1), np.array(l2)

    Fs, abcs = calc_circ_coefs(l0, l1, l2, r0, r1, r2)
    # Catch for F=0 (i.e. no circle exists)
    if Fs[0] == 0:
        return
    a, b, c = calc_circ_abcs(Fs, r0)
    # Calculate the discriminant.
    disc = b ** 2 - 4 * a * c
    # If the discriminant is negative then the tangential circle does not exist.
    if round(disc, 10) > 0:
        # Grab the two roots
        rs = [_ for _ in np.roots(np.array([a, b, c])) if np.isreal(_)]
        # If there is only one root return it
        if len(rs) == 1:
            r = rs[0]
        # If there are 2 roots choose between them
        else:
            # If the smaller of the two roots is negative return the other root
            if min(rs) < 0:
                r = max(rs)
            # If they're both positive, return the smaller of the two
            elif rs[0] > 0 and rs[1] > 0:
                r = min(rs)
            # If they're both negative return
            else:
                return
        F, Fx0, Fx1, Fy0, Fy1, Fz0, Fz1 = Fs
        # Calculate the vertex based off of our coefficient values and the sphere's radius
        x = Fx0 / F + r * Fx1 / F + l0[0]
        y = Fy0 / F + r * Fy1 / F + l0[1]
        z = Fz0 / F + r * Fz1 / F + l0[2]
        return np.array([x, y, z]), r