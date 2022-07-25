from System.sys_calcs import *
from System.Network.vertex import Vertex


# Calculate edge function. Takes in an edge and the network and returns the closest vertex along that edge
def calc_edge(edge, net):
    # Grab the vertex of the edge for later use
    v0 = edge.verts[0]
    # Find the center of the edge atoms
    circ = calc_circ(edge.atoms)
    # If the circle made between the atom edges is Nonetype return
    if not circ:
        return
    else:
        circ = circ[0][0]
    # Find the distance between the old vertex and the center of the bottleneck
    d1 = calc_dist(v0.loc, circ)
    # Make a list of the closest neighbors
    neighbors = sortbyDist(edge.atoms, net)
    # Find the other atom from the old vertex
    an = None
    # Find the atom from the v0 that is not in the edge's atom list
    for atom in v0.atoms:
        # Set our changing atom to the outlier atom found above and change the radius by 5%
        if atom not in edge.atoms:
            an = atom
    # Set the closest distance to infinity and the vertex to None
    c_dist = np.inf
    myVert = None
    # Go through each atom in the neighbor list
    for n in neighbors[:15]:
        # Check to see if the neighbor is the old vertex's other atom. If so continue
        if n == an:
            continue
        # Calculate the vertex of the edge atoms with the neighbor atom
        vn = Vertex(edge.atoms + [n])
        # Make sure that the vertex exists and does not overlap with the old vertex
        if not vn or check_vert(set(vn.atoms), net.verts):
            continue
        # Check to see if it is a real vertex or not
        i = 0
        overlap = False
        while i < len(neighbors):
            if {neighbors[i]}.issubset(edge.atoms + [n]):
                continue
            if calc_dist(vn.loc, neighbors[i].loc) - (vn.rad + neighbors[i].rad) < 0:
                overlap = True

        if not overlap:
            # Calculate the distance between the new vertex and the center
            d2 = calc_dist(vn.loc, circ)
            # Calculate the distance between the new vertex and the old vertex
            d3 = calc_dist(vn.loc, v0.loc)
            # Mevdevev's edge site finding checks. Find the shortest relative distance from v0 to vn
            if d1 <= d3 or d2 <= d3:
                r_len = d1 + d2
            else:
                r_len = d3
            if r_len < c_dist:
                myVert = vn
                c_dist = r_len
    return myVert


# Calculate edge function. Chases an edge toward the next vertex
def calc_edge1(edge, net, dt=None):
    # Find the closest neighbors
    neighbors = sortbyDist(edge.atoms, net, length=50)
    # Get the old vertex
    v0 = edge.verts[0]
    # Estimate a working dt
    if dt is None:
        dt = edge.atoms[0].rad / 20
    an = None
    # Find the atom from the v0 that is not in the edge's atom list
    for atom in v0.atoms:
        # Set our changing atom to the outlier atom found above and change the radius by 5%
        if atom not in edge.atoms:
            an = atom
    # Calculate the bottleneck of the edge
    bn = calc_circ(edge.atoms)[0][1]
    # Adjust the size of the atom to fit through the bottleneck
    if bn < 1.05*an.rad:
        an.rad = 0.95*bn
    # Find the vertex between the edge atoms and the adjusted atom
    vn = Vertex(edge.atoms + [an])
    # If we get a None vertex the shrink went too far. Keep increasing radius until a vertex is found.
    while vn is None:
        an.rad = an.rad * 1.01
        vn = Vertex(edge.atoms + [an])
    # If the radius is larger than the bottleneck, continue and hope that the other side of the edge will be able to
    if vn.rad > bn:
        return
    # Find the initial direction by getting the vector between the new vertex formed after the atom got smaller
    dr = np.array([v0.loc[0] - vn.loc[0], v0.loc[1] - vn.loc[1], v0.loc[2] - vn.loc[2]])
    elen = 0
    vfound = False
    vert = None
    # Keep adding points to the edge until the next vertex is found or the edge left the network
    while not vfound:
        # Normalize the direction vector
        dr_mag = np.sqrt(dr.dot(dr))
        dr = dr / dr_mag
        # Add up the length of the edge
        elen += dr_mag
        if elen > net.rad:
            edge.verts.append(None)
            return None
        # Record vns location before changing it
        vn_1 = vn
        # Move the atom along the direction of the edge by dt increments
        an.loc = an.loc[0] + dt*dr[0], an.loc[1] + dt*dr[1], an.loc[2] + dt*dr[2]
        # Calculate the new vertex
        vn = Vertex(edge.atoms + [an])
        # Add the vertex location to the edges points
        edge.points.append(vn.loc)
        # Find the new move direction by finding the direction from vn-1 to vn
        dr = np.array([vn.loc[0] - vn_1.loc[0], vn.loc[1] - vn_1.loc[1], vn.loc[2] - vn_1.loc[2]])
        # Check to see if we have passed a vertex
        for vert in neighbors:
            # Calculate the vectors between the vertex and the new and old edge points
            d1 = np.array([vn_1.loc[0] - vert.loc[0], vn_1.loc[0] - vert.loc[0], vn_1.loc[0] - vert.loc[0]])
            d2 = np.array([vn.loc[0] - vert.loc[0], vn.loc[0] - vert.loc[0], vn.loc[0] - vert.loc[0]])
            # Check to see if the vertex is in between the new and old edge points
            if np.sqrt(d1.dot(d1)) <= dr_mag and np.sqrt(d2.dot(d2)) <= dr_mag:
                # If so, we have found our vert and exit
                vfound = True
    # Add the vertex to the edge and the network
    edge.verts.append(vert)
    net.verts.append(vert)
    return vert