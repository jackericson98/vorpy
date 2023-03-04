import numpy as np


def draw_line(points, radius):
    # Initiate the draw attributes
    draw_points, draw_tris = [], []
    # Go through the points
    for i in range(len(points)):
        # If we are at the end of the points list, use the previous point for calibration
        if i == len(points) - 1:
            p0, p1 = np.array(points[i]), np.array(points[i - 1])
            r = p0 - p1
        else:
            p0, p1 = np.array(points[i]), np.array(points[i + 1])
            r = p1 - p0
        # Find the vector and its normal between the two points
        rn = r / np.linalg.norm(r)
        # In the case that the vector between the points is in the z direction only, move it
        if rn[0] == 0 and rn[1] == 0:
            r = r + np.array([0.001, 0.001, 0])
            rn = r / np.linalg.norm(r)
        # Take the cross product with the +z direction and normalize it
        v0_0x = np.cross(rn, [0, 0, 1])
        v0_0n = v0_0x / np.linalg.norm(v0_0x)
        # Calculate the location of the first point
        p0_0 = v0_0n * radius + p0
        # Take the cross product of the edge vector and the vector to the first point and normalize it
        v0_1x = np.cross(rn, v0_0n)
        v0_1nx = v0_1x / np.linalg.norm(v0_1x)
        # Find the vectors for the other two points (30/60/90 triangle)
        v0_1 = - 0.5 * radius * v0_0n + 0.5 * np.sqrt(3) * radius * v0_1nx
        v0_2 = - 0.5 * radius * v0_0n - 0.5 * np.sqrt(3) * radius * v0_1nx
        # Get the points and add them to the list of draw points
        p0_1, p0_2 = v0_1 + p0, v0_2 + p0
        draw_points += [p0_0, p0_1, p0_2]
    # Go through the points
    for i in range(len(points) - 1):
        # List the points
        p0_0, p0_1, p0_2, p1_0, p1_1, p1_2 = range(3 * i, 3 * (i + 2))
        # Create the triangles
        draw_tris += [[p0_0, p0_1, p1_0], [p1_0, p1_1, p0_1],
                           [p0_1, p0_2, p1_1], [p1_1, p1_2, p0_2],
                           [p0_2, p0_0, p1_2], [p1_2, p1_0, p0_0]]
    # Return the points and triangles
    return draw_points, draw_tris


def draw_grid(net):
    # Set up the grid of points
    grid_points = [[[]]]


# Draw Edge Function. Takes in an edge and updates its attributes draw_points, draw_tris
def draw_edge(edge, radius=0.01):
    # Make sure the edge is built already
    if edge.points is None:
        edge.build()
    # Calculate the lines
    edge.draw_points, edge.draw_tris = draw_line(edge.points, radius)


# Draw Vertices Function. Takes in a vertex and updates the loc_points, loc_tris, sphere_points, sphere_points attribute
def draw_vert(vert, sphere=False, radius=0.05, resolution=0.1):
    loc = np.array(vert.loc)
    # Draw the point
    xp, xn = loc + np.array([radius, 0, 0]), loc - np.array([radius, 0, 0])
    yp, yn = loc + np.array([0, radius, 0]), loc - np.array([0, radius, 0])
    zp, zn = loc + np.array([0, 0, radius]), loc - np.array([0, 0, radius])
    # Connect the points
    vert.loc_points = [xp, xn, yp, yn, zp, zn]
    vert.loc_tris = [[0, 2, 4], [0, 2, 5], [0, 3, 4], [0, 3, 5], [1, 2, 4], [1, 2, 5], [1, 3, 4], [1, 3, 5]]
    # If th user wants a sphere
    if sphere:
        # Set the resolution of the spheres
        circ = 2 * np.pi * vert.rad
        f = max(int(circ / resolution), 3)
        # Find u, v values that span phi and theta
        u, v = np.mgrid[0:2 * np.pi:f * resolution * 2j, 0:np.pi:f * resolution * 1j]
        # Plot each sphere
        # Get x, y, z data for the wireframe
        x = vert.rad * np.cos(u) * np.sin(v) + loc[0]
        y = vert.rad * np.sin(u) * np.sin(v) + loc[1]
        z = vert.rad * np.cos(v) + loc[2]

