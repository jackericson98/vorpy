import numpy as np
from scipy.spatial import Delaunay as dl


# Calculate tetrahedron volume function. Calculated the volume of a tetrahedron defined by its vertices
def calc_tetra_vol(p0, p1, p2, p3):
    # Choose a base point (p0) and find the vectors between it and other points
    r01 = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    r02 = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    r03 = p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]
    # Formula for tetrahedron volume: 1/6 * r03 dot (r01 cross r02)
    vol = (1/6)*abs(np.dot(r03, np.cross(r01, r02)))
    return vol


# Calculate triangle are function. Takes in 3 points in 3 space and returns the area of the triangle created by them
def calc_tri(points):
    # Get the two triangles vectors
    AB = np.array(points[0]) - np.array(points[1])
    AC = np.array(points[0]) - np.array(points[2])
    # Return half the cross product between the two vectors
    return 0.5 * np.linalg.norm((np.cross(AB, AC)))


# Calculate surface area function. Takes in a
def calc_sa(points):
    # Calculate the Delaunay simplexes
    tri = dl(points)
    sa = 0
    # Go through each simplex
    for simplex in tri.simplices:
        spoints = []
        # Grab the points of the simplex
        for ndx in simplex:
           spoints.append(points[ndx])
        # Calculate the area of the triangle made by the simplex and add it to the total surface area
        sa += calc_tri(spoints)
    # Return the total surface area
    return sa


# Calculate interface function. Takes in a set of surfaces and calculates the total surface area
def calc_interface(surfs):
    sa = 0
    for surf in surfs:
        sa += calc_sa(surf.points + surf.edge_points)
    return sa


# Calculate cell volume function. Grabs the points in a cell and calculates the volume made by the tetrahedrons
def calc_vol(atom):
    points = []
    # Grab the points for the surfaces in the cell
    for surf in atom.surfs:
        points += surf.points
    # Grab the points for the edges in the cell
    for edge in atom.edges:
        points += edge.points
    # Grab the vertices locations for the cell
    for vert in atom.verts:
        points += [vert.loc]


# Analyze system function. Finds the surfaces and volumes of the system
def analyze(sys):
    # Go through each surface in the system and find the simplices and the surface area
    for surf in sys.net.surf:
        # Get the surfaces simplices
        surf.simps = surf.find_simps()
        # Get the surface area of the surface
        surf.sa = calc_sa(surf)

    # Go through each atom in the system
    for atom in sys.atoms:
        atom.cell_vol = calc_vol(atom)
