from System.Network.Edges.edge_calcs import *


# Calculate edge point function. Takes in a surface and a point and returns the intersection point of the vector
# from the center of the smallest of the surfaces 2 atoms through the point into the surface
def calc_edge_point(edge, s0):
    # Grab the function's coefficients
    f = s0.func
    # Grab the vertex points
    pv0, pv1 = np.array(edge.verts[0].loc), np.array(edge.verts[1].loc)
    # Find the point in between the two vertex points
    r01 = pv1 - pv0
    r01_mag = np.linalg.norm(r01)
    rn01 = r01 / r01_mag
    # Get the center point of the vertices
    cp = pv0 + 0.5 * rn01 * r01_mag

    # Finding the a, b, c, values that satisfy at**2 + bt + c = 0
    a = f[0] * rn01[0] ** 2 + f[1] * rn01[1] ** 2 + f[2] * rn01[2] ** 2 + f[3] * rn01[0] * rn01[1] + f[4] * rn01[1] * rn01[2] + f[5] \
        * rn01[2] * rn01[0]
    b = 2 * f[0] * rn01[0] * cp[0] + 2 * f[1] * rn01[1] * cp[1] + 2 * f[2] * rn01[2] * cp[2] + f[3] \
        * (rn01[0] * cp[1] + rn01[1] * cp[0]) + f[4] * (rn01[1] * cp[2] + rn01[2] * cp[1]) + f[5] \
        * (rn01[2] * cp[0] + rn01[0] * cp[2]) + f[6] * rn01[0] + f[7] * rn01[1] + f[8] * rn01[2]
    c = f[0] * cp[0] ** 2 + f[1] * cp[1] ** 2 + f[2] * cp[2] ** 2 + f[3] * cp[0] * cp[1] + f[4] * cp[1] * cp[2] + \
        f[5] * cp[2] * cp[0] + f[6] * cp[0] + f[7] * cp[1] + f[8] * cp[2] + f[9]

    # Given a positive discriminant, find the root closer to the sphere, corresponding to the correct surface
    # and add that point to our surface list of points
    if round(b ** 2 - 4 * a * c, 4) >= 0:
        # If the projection point on a0's surface is outside a1's surface take the smallest of the roots
        roots = np.roots([a, b, c])
        mag = min(abs(roots))
        return cp + mag * rn01

