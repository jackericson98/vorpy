from numba import jit
import numpy as np


@jit(nopython=True, cache=True)
def edge_project(rn, pa, f):
    """
    Projects an edge onto a surface defined by a polynomial equation. Calculates the intersection points of a line,
    defined by a point and a direction, with a surface defined by a polynomial equation.

    Args:
        rn (numpy.ndarray): The direction vector of the line.
        pa (numpy.ndarray): A point on the line.
        f (numpy.ndarray): Coefficients of the polynomial defining the surface.

    Returns:
        numpy.ndarray: The point on the surface closest to the starting point 'pa' along the direction 'rn' that
                       intersects the surface defined by 'func'.

    Notes:
        - The polynomial is assumed to be of the form:
          f(x, y, z) = f[0]*x^2 + f[1]*y^2 + f[2]*z^2 + f[3]*xy + f[4]*yz + f[5]*zx +
                       f[6]*x + f[7]*y + f[8]*z + f[9]
        - The function calculates the coefficients of the quadratic equation in t, obtained by substituting
          the parametric equations of the line into the polynomial equation.
        - It then solves this quadratic equation to find the values of t that correspond to the intersection points.
    """

    # Calculate coefficients of the quadratic equation at**2 + bt + c = 0 for the line
    a = (f[0] * rn[0] ** 2 + f[1] * rn[1] ** 2 + f[2] * rn[2] ** 2 + f[3] * rn[0] * rn[1] + f[4] * rn[1] * rn[2] +
         f[5] * rn[2] * rn[0])
    b = (2 * f[0] * rn[0] * pa[0] + 2 * f[1] * rn[1] * pa[1] + 2 * f[2] * rn[2] * pa[2] +
         f[3] * (rn[0] * pa[1] + rn[1] * pa[0]) + f[4] * (rn[1] * pa[2] + rn[2] * pa[1]) +
         f[5] * (rn[2] * pa[0] + rn[0] * pa[2]) + f[6] * rn[0] + f[7] * rn[1] + f[8] * rn[2])
    c = (f[0] * pa[0] ** 2 + f[1] * pa[1] ** 2 + f[2] * pa[2] ** 2 + f[3] * pa[0] * pa[1] + f[4] * pa[1] * pa[2] +
         f[5] * pa[2] * pa[0] + f[6] * pa[0] + f[7] * pa[1] + f[8] * pa[2] + f[9])

    # Calculate the discriminant of the quadratic equation
    discriminant = b ** 2 - 4 * a * c

    # Proceed only if the discriminant is non-negative, indicating real roots
    if round(discriminant, 10) >= 0:
        # Solve the quadratic equation to find intersection points
        roots = np.roots(np.array([a, b, c]))

        # Determine the valid root closest to the given point
        if len(roots) == 1:
            return pa + roots[0] * rn
        else:
            p1 = pa + roots[0] * rn
            p2 = pa + roots[1] * rn

            # Evaluate which point is valid based on the direction and proximity
            if roots[0] < 0:
                return p2
            elif roots[1] < 0:
                return p1
            elif roots[0] < roots[1]:
                return p1
            else:
                return p2

