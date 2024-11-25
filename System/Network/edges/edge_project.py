import numpy as np
from System.sys_funcs.calcs.calcs import calc_dist
from numba import jit


# @jit(nopython=True)
def edge_project(rn, pa, func, elocs, erads, direction=None):
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
    # Given a positive discriminant, find the root closer to the sphere, corresponding to the 181L surface
    # and add that point to our surface list of points
    if round(b ** 2 - 4 * a * c, 10) >= 0:
        # Calculate the roots
        roots = np.roots(np.array([a, b, c]))
        # If one root exists return it
        if len(roots) == 1:
            return pa + roots[0] * rn
        else:
            p1 = pa + roots[0] * rn
            p2 = pa + roots[1] * rn
        # Check the distance from the edge balls
        d11 = round(calc_dist(p1, elocs[0]) - erads[0], 5)
        d12 = round(calc_dist(p1, elocs[1]) - erads[1], 5)
        d13 = round(calc_dist(p1, elocs[2]) - erads[2], 5)
        d21 = round(calc_dist(p2, elocs[0]) - erads[0], 5)
        d22 = round(calc_dist(p2, elocs[1]) - erads[1], 5)
        d23 = round(calc_dist(p2, elocs[2]) - erads[2], 5)

        # First find out if either of the roots is negative for later flagging
        if roots[0] < 0:
            point = p2
        elif roots[1] < 0:
            point = p1
        elif roots[0] < roots[1]:
            point = p1
        else:
            point = p2

        # First test for the distance between the points and the surfaces of the edge balls
        if d11 == d12 == d13 and d21 == d22 == d23:
            pass
        elif d11 == d12 == d13:
            if all([point[i] == p1[i] for i in range(3)]):
                print("1 INCOMPATIBLE SOLUTION:   ", p1, p2, point)
                print(rn, pa, func, elocs, erads, direction)
            return p1
        elif d21 == d22 == d23:
            if not all([point[i] == p2[i] for i in range(3)]):
                print("2 INCOMPATIBLE SOLUTION:   ", p1, p2, point)
                print(rn, pa, func, elocs, erads, direction)
            return p2
        # Next test to see if either of the points is heading the wrong direction
        if direction is not None:
            # Calculate the dot product of the vectors between the middle point and the direction
            p1_dot = np.dot(p1 - pa, direction)
            p2_dot = np.dot(p2 - pa, direction)
            if p1_dot > 0 > p2_dot:
                if not all([point[i] == p1[i] for i in range(3)]):
                    print("3 INCOMPATIBLE SOLUTION:   ", p1, p2, point)
                    print(rn, pa, func, elocs, erads, direction)
                return p1
            elif p2_dot > 0 > p1_dot:
                if not all([point[i] == p2[i] for i in range(3)]):
                    print("4 INCOMPATIBLE SOLUTION:   ", p1, p2, point)
                    print(rn, pa, func, elocs, erads, direction)
                return p2
        # Last go with the smaller of the two distances to the balls
        if d11 < d22:
            if not all([point[i] == p1[i] for i in range(3)]):
                print("5 INCOMPATIBLE SOLUTION:   ", p1, p2, point)
                print(rn, pa, func, elocs, erads, direction)
            return p1
        if not all([point[i] == p2[i] for i in range(3)]):
            print("6 INCOMPATIBLE SOLUTION:   ", p1, p2, point)
            print(rn, pa, func, elocs, erads, direction)
        return p2
