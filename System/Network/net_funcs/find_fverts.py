from System.system import System
from System.Network.network import Network, Vertex
from System.sys_funcs.calcs import calc_angle, np, ndx_search
from System.Network.net_objs.surface import Surface
from System.Network.net_objs.edge import Edge


def fget_circ_rad(dist, r0, r1):
    # First get the distance between the locations
    rad = (1/(2*dist)) * np.sqrt(4 * dist ** 2 * r0 ** 2 - (dist ** 2 - r1 ** 2 + r0 ** 2))
    return rad


def ffind_near_atoms(net, a0, inc=0):
    # Get the closest atoms
    near_atoms = net.get_atoms(cells=a0.box, reach=inc)
    dists = []
    a0_array = np.array(a0.loc)
    # Get the atoms distance from a0
    for atom in near_atoms:
        dists.append(np.sqrt(sum(np.square(a0_array - np.array(atom.loc)))) - (a0.rad + atom.rad))
    # Sort the atoms
    sorted_atoms = [[x, _] for _, x in sorted(zip(dists, near_atoms), key=lambda pair: pair[0]) if x.num != a0.num]
    return sorted_atoms


# Expand atom function. Expands into the surrounding atoms, creating intersecting spheres -> circles
def expand_atom(net, a0, atoms, dists, max_dist=None):
    # Set up the centers and rn lists
    surfs = []
    # Go through the atoms in the list of atoms
    for i in range(len(atoms)):
        # Set up the atom
        an = atoms[i]
        # Check that the atom is not the one we are checking against
        if a0.num == an.num:
            continue
        # If the center is out of range return
        if max_dist is not None and 0.5 * dists[i] > max_dist:
            return surfs
        # Get the vector between the atoms, the normal of that vector and the center
        r = np.array(a0.loc) - np.array(an.loc)
        rn = r / np.linalg.norm(r)
        center = np.array(an.loc) + (a0.rad + 0.5 * dists[i]) * rn
        # Add the centers and rns
        surfs.append(Surface(net=net, atoms=[a0, an], center=center, rn=rn))
    # Return a tuple
    return surfs


# Expand circle function. Keeps expanding the given circle into surrounding circles providing overlapping lines
def expand_circle(net, s0, surfs, max_dist=None):
    """
    Expands the intersecting circle between two balloons until it meets another circle around the same atom
    :param net:
    :param s0:
    :param surfs:
    :param max_dist:
    :return:
    """
    # Create the line vectors and centers lists
    l_vctrs, l_cntrs, l_dists, l_atoms = [], [], [], []
    # Go through the circles finding where they intersect
    for i in range(len(surfs)):
        # Check to see if the c_vctr is the same as c0
        if surfs[i].ndx == s0.ndx:
            continue
        # Get the edge atom
        ndx = [_ for _ in surfs[i].ndx if _ not in s0.ndx] + s0.ndx
        ndx.sort()
        # Get the plane coefficients
        a0, b0, c0 = s0.rn
        a1, b1, c1 = surfs[i].rn
        # Get the offset
        d0, d1 = np.dot(s0.rn, s0.center), np.dot(surfs[i].rn, surfs[i].center)
        # Get the parameterized line equations
        denominator = a1 * b0 - a0 * b1
        xt = [(-b1 * c0 + b0 * c1) / denominator, (-b1 * d0 + b0 * d1) / denominator]
        yt = [(a1 * c0 - a0 * c1) / denominator, (a1 * d0 - a0 * d1) / denominator]
        zt = [1, 0]

        # Get the normal to the line
        r = np.array([_[0] for _ in [xt, yt, zt]])
        l_vctr = r / np.linalg.norm(r)
        # Find an arbitrary point on the line (t = 0)
        pa = np.array([_[1] for _ in [xt, yt, zt]])
        # Get the vector between this and the center of the circle
        pac = np.array(c0[1]) - pa
        # Get the center of the line by dotting the atom's location onto it
        l_cntr = pa + l_vctr * np.dot(pac, l_vctr)
        # Check to see if the line center is too far away from the atom's location
        l_dist = np.sqrt(sum(np.square(np.array(c0[0]) - np.array(l_cntr))))

        if max_dist is not None and l_dist > max_dist:
            continue
        # Load up the variables
        l_dists.append(l_dist)
        l_vctrs.append(l_vctr)
        l_cntrs.append(l_cntr)
        l_atoms.append(ndx)
    # Sort the vectors and the centers by their distance from the atom
    ls = [[d, v, c, a] for d, v, c, a in sorted(zip(l_dists, l_vctrs, l_cntrs, l_atoms), key=lambda quad: quad[0])]
    edges = [Edge(atoms=[net.atoms[_] for _ in ls[i][3]], rn=ls[i][1], loc=ls[i][2], dist=ls[i][0]) for i in range(len(ls))]
    return edges


# Expand line function. Keeps expanding the given line until both sides intersect another line
def expand_line(l0, lines, atoms, max_dist=np.inf):
    """

    :param l0: iterable - baseline with a normalized vector along the line from the point closest to the base atom
    :param lines: list - connecting check lines
    :return: two vertices or 1 vertex and None
    """
    # Set up the vertex locations and radii lists
    neg_vert, pos_vert, neg_vert_dist, pos_vert_dist = None, None, max_dist, -max_dist
    # Go through the lines in the list and
    for i in range(len(lines)):
        # Get the vector between the line's center points from l0
        rx0, ry0, rz0, x0, y0, z0 = l0[0] + l0[1]
        rx1, ry1, rz1, x1, y1, z1 = l0[0] + l0[1]
        # Calculate the intersecting point between the two lines
        t1 = (ry0 * x0 + ry0 * x1 - rx0 * y0 - rx0 * y1)
        t0 = (x0 + x1) / rx0 + rx1 * t1 / rx0
        # Get our x, y, z, variables
        x, y, z = rx0 * t0 + x0, ry0 * t0 + y0, rz0 * t0 + z0
        # Calculate the distance between this and l0
        my_dist = np.sqrt(sum(np.square(np.array([x, y, z]) - np.array(l0[1]))))
        # Add it to the correct bin
        if 0 <= t0 and my_dist < pos_vert_dist:
            # Check to see if the vertex is closer to the center of the edge or not
            pos_vert = [x, y, z]
            pos_vert_dist = my_dist
        elif 0 > t0 and my_dist < neg_vert_dist:
            neg_vert = [x, y, z]
            neg_vert_dist = my_dist
    # We should have the closest two vertices to this edge
    return pos_vert, neg_vert


def ffind_verts(net, a0=None, max_dist=10):
    # Get the starting atom
    if a0 is None:
        # Get a random atom
        a0 = net.atoms[np.random.randint(0, len(net.atoms) - 1)]
    # Get the maximum increment to search against for each atom
    max_inc = int(net.max_vert / max(net.sub_box_size)) + 1
    # Set up the atoms stack
    atoms_stack = [a0]
    # Keep searching until the atoms stack is empty
    while atoms_stack:
        # Pull the most recent atom
        my_atom = atoms_stack.pop()
        cell_complete = False
        # Get the closest atoms based on the current inc
        search_atoms = ffind_near_atoms(net, my_atom, inc=max_inc)
        # Set up the atoms and distances lists
        dists = [_[1] for _ in search_atoms]
        atoms = [_[0] for _ in search_atoms]
        # Get the surfaces
        my_surfs = expand_atom(net=net, a0=my_atom, dists=dists, atoms=atoms, max_dist=max_dist)
        inc = 0
        # Keep looking for close atoms until the atom is complete or the distance is less than the max allowed
        while not cell_complete or len(my_surfs) > 0:

            # Get one of the
            my_surf = my_surfs.pop(0)
            ndx = sorted([_.num for _ in atoms])
            # Search to see if the atoms have been found before
            surf_ndx = ndx_search(net.surf_ndxs, ndx)
            # Skip if found before
            if len(net.surf_ndxs) < surf_ndx and net.surf_ndxs[surf_ndx] == ndx:
                continue
            # Get the edges of the current surface
            lines = expand_circle(net=net, s0=my_surf, surfs=my_surfs, max_dist=max_dist)
            # Start the lines while loop
            loop_complete = False
            found_edges_ndxs = []
            while not loop_complete or len(lines) == 0:
                # Get the edge from the lines
                l_vctr, l_cntr, l_atoms = lines.pop(0)
                edge_ndx = sorted([_.ndx for _ in l_atoms])
                # Search for previously found edge
                line_ndx = ndx_search(net.edge_ndxs, edge_ndx)
                if len(net.edge_ndxs) > line_ndx and net.edge_ndxs[line_ndx] == edge_ndx:
                    found_edges_ndxs.append(edge_ndx)
                    continue
                # Get the vertices for the edge
                vert_locs, vert_atoms = expand_line(a0, [l_vctr, l_cntr], my_atom.edges)



            # Test to see if the cell is complete


if __name__ == '__main__':
    mySys = System(file='./Data/test_data/cube.pdb')
    mySys.net = Network(atoms=mySys.atoms, sys=mySys)
    ffind_verts(mySys.net)

