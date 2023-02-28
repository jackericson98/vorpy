import numpy as np
from System.system import System
from System.Network.network import Network
from System.Network.net_objs.surface import Surface
from System.Network.net_objs.edge import Edge
from System.Network.net_objs.vertex import Vertex
from System.sys_funcs.calcs import ndx_search


def ffind_near_atoms(net, a0, max_dist=20):
    # Get the closest atoms to a0
    max_atom_dist = 0
    dists = []
    a0_array = np.array(a0.loc)
    max_inc = int(max_dist/min(net.sub_box_size)) + 1
    # Get the atoms
    my_atoms = net.get_atoms([a0.box], reach=max_inc)
    # Get the atoms distance from a0
    for atom in my_atoms[len(dists):]:
        # Calculate the distance between the atoms
        my_dist = np.sqrt(sum(np.square(a0_array - np.array(atom.loc)))) - (a0.rad + atom.rad)
        # Replace the maximum distance if needed
        if my_dist > max_atom_dist:
            max_atom_dist = my_dist
            # Check that the atom's distance isn't larger than the specified maximum distance
            if my_dist < max_dist:
                dists.append(my_dist)
        else:
            dists.append(my_dist)
    # Sort the atoms
    sorted_atoms = [[x, _] for _, x in sorted(zip(dists, my_atoms), key=lambda pair: pair[0]) if x.num != a0.num]
    # Make sure the return doesn't trigger an index error
    return sorted_atoms


def ffind_near_surfs(net, s0):
    # First get the first atom's surfaces
    my_surfs = [_ for _ in s0.atoms[0].surfs if _.ndx != s0.ndx]
    # Add the surfaces from a1
    my_surfs += [_ for _ in s0.atoms[1].surfs if _ not in my_surfs and _.ndx != s0.ndx]
    my_edges, my_dists = [], []
    # Find the intersection of the closest surfaces with s0 and calculate center of the edge
    for i in range(len(my_surfs)):
        # Get the x, y, z variables for the check surfaces normal
        a0, b0, c0 = s0.norm
        a1, b1, c1 = my_surfs[i].norm
        # Get the offset
        d0, d1 = np.dot(s0.norm, s0.loc), np.dot(my_surfs[i].norm, my_surfs[i].loc)
        # Get the parameterized line equations
        denominator = a1 * b0 - a0 * b1
        if denominator == 0:
            continue
        xt = [(-b1 * c0 + b0 * c1) / denominator, (-b1 * d0 + b0 * d1) / denominator]
        yt = [(a1 * c0 - a0 * c1) / denominator, (a1 * d0 - a0 * d1) / denominator]
        zt = [1, 0]
        # Get the normal to the line
        r = np.array([_[0] for _ in [xt, yt, zt]])
        l_vctr = r / np.linalg.norm(r)
        # Find an arbitrary point on the line (t = 0)
        pa = np.array([_[1] for _ in [xt, yt, zt]])
        # Get the vector between this and the center of the circle
        pac = np.array(s0.loc) - pa
        # Get the center of the line by dotting the atom's location onto it
        l_cntr = pa + l_vctr * np.dot(pac, l_vctr)
        # Get the atoms by finding the outlier in the my_surf.atoms list and adding the check_surf.atoms
        edge_atoms = [_ for _ in s0.atoms if _.num not in my_surfs[i].ndx] + my_surfs[i].atoms
        edge_ndx = [_.num for _ in edge_atoms]
        edge_ndx.sort()
        edge_dist = np.sqrt(sum(np.square(np.array(s0.loc) - np.array(l_cntr))))
        my_edges.append(Edge(net=net, atoms=edge_atoms, ndx=edge_ndx, center=l_cntr, normal=l_vctr, dist=edge_dist))
        my_dists.append(edge_dist)
    # Return the list of surfaces sorted by maximum dot product with s0's normal
    return [[x, _] for _, x in sorted(zip(my_dists, my_surfs), key=lambda pair: pair[0])]


def calc_flat_vert(atoms):
    # Get the plane equations
    coeffs = []
    # Go through the atoms to make the planes
    for an in atoms[1:]:
        # Get the point between the atoms
        r = np.array(an.loc) - np.array(atoms[0].loc)
        norm = np.linalg.norm(r)
        rn = r / norm
        atom_surf_dist = norm - an.rad - atoms[0].rad
        center = (0.5 * atom_surf_dist + atoms[0].rad) * rn + np.array(atoms[0].loc)
        coeffs.append(rn.tolist() + [np.dot(rn, center)])

    a, b, c, d = coeffs[0]
    e, f, g, h = coeffs[1]
    i, j, k, m = coeffs[2]

    disc = c*f*i - b*g*i - c*e*j + a*g*j + b*e*k - a*f*k
    x_numerator = d*g*j - c*h*j - d*f*k + b*h*k + c*f*m - b*g*m
    y_numerator = - d*g*i + c*h*i + d*e*k - a*h*k - c*e*m + a*g*m
    z_numerator = d*f*i - b*h*i - d*e*j + a*h*j + b*e*m - a*f*m
    x, y, z = x_numerator / disc, y_numerator / disc, z_numerator / disc
    return [x, y, z]


def ffind_all_verts(surf, net):
    # Make sure the initial vertices are found
    if len(surf.verts) <= 3:
        ffind_init_edges(net, surf)
    j = 3
    surf_vert_dists = [_.dist for _ in surf.edges]
    # Keep looking for edges until the maximum vertex distance is less than the closest outside edge
    while True:
        check_edge, check_edge_dist = surf.edge_bank[j]
        if check_edge_dist < max(surf_vert_dists):
            # Create a variable to track the minimum set of edges
            min_edges = []
            min_edge_dists = []
            # We need to add the edge to the surface where it is still sorted
            for k in range(len(surf.edges)):
                check_dist = np.sqrt(
                    sum(np.square(np.array(surf.edge_bank[j][0].loc) - np.array(surf.edges[k].loc))))
                if len(min_edges) < 2 or check_dist < min_edge_dists[1]:
                    min_edges[1], min_edge_dists[1] = check_edge, check_dist
                # Quick sort
                if min_edge_dists[0] > min_edge_dists[1]:
                    min_edges[0], min_edges[1] = min_edges[1], min_edges[0]
                    min_edge_dists[0], min_edge_dists[1] = min_edge_dists[1], min_edge_dists[0]
            nsrt_ndx = min([surf.edges.index(_) for _ in min_edges])
            surf.edges.insert(nsrt_ndx, check_edge)
        else:
            break
        j += 1


def ffind_init_edges(net, my_surf):
    # Create the sorted list of surfaces
    my_surf.edge_bank = ffind_near_surfs(net=net, s0=my_surf)
    # Add the closest 3 edges to the surface's list of surfaces
    my_surf.edges = [my_surf.edge_bank[i][0] for i in range(3)]
    # Get the initial vertices for each of the edge combinations
    surf_verts = []
    for j in range(3):
        # Get the edges for this vertex
        e1, e2 = my_surf.edges[j], my_surf.edges[(j + 1) % 3]
        # Pull the atoms from the edges and get their indices
        my_vert_atoms = [_ for _ in e1.atoms if _.num not in e2.ndx] + e2.atoms
        my_vert_ndx = [_.num for _ in my_vert_atoms]
        my_vert_ndx.sort()
        # Create the vertex object and add it to the list of surface vertices
        surf_verts.append(Vertex(net=net, atoms=my_vert_atoms, ndx=my_vert_ndx, edges=[e1, e2],
                          location=calc_flat_vert(my_vert_atoms)))
    # Calculate the distances from the surface's center
    surf_vert_dists = [np.sqrt(sum(np.square(np.array(_.loc) - np.array(my_surf.loc)))) for _ in surf_verts]




def ffind_init_surfs(net, my_atom, max_dist=20):
    # Get the sorted closest max_surfs number of atoms, atoms
    my_atom.surf_bank = ffind_near_atoms(net=net, a0=my_atom, max_dist=max_dist)
    # Add the atoms to the net surfs list
    for j in range(4):
        # Get the check atom variable
        check_atom, check_dist = my_atom.surf_bank[j]
        # Get the index to search for
        surf_ndx = [my_atom.num, check_atom.num]
        surf_ndx.sort()
        surf_net_ndx = ndx_search(net.surf_ndxs, surf_ndx)
        if surf_net_ndx < len(net.surfs) and net.surf_ndxs[surf_net_ndx] == surf_ndx:
            continue
        # Create a new surface object
        my_surf = Surface(atoms=[my_atom, check_atom], distance=check_dist)
        my_surf.ndx = surf_ndx
        net.surfs.insert(surf_net_ndx, my_surf)
        net.surf_ndxs.insert(surf_net_ndx, surf_ndx)
        # Find where to insert the surface
        k = 0
        while k < len(my_atom.surfs) and check_dist > my_atom.surfs[k].dist:
            k += 1
        # Insert the surface's other atom's index and distance
        my_atom.surfs.insert(k, my_surf)
        # Find where to insert my_atom into the close atom's list
        m = 0
        while m < len(check_atom.surfs) and check_dist > check_atom.surfs[m].dist:
            m += 1
        check_atom.surfs.insert(m, my_surf)


if __name__ == '__main__':
    my_sys = System(file='C:/Users/jacke/PycharmProjects/vorpy/Data/test_data/Na5.pdb', root_dir='C:/Users/jacke/PycharmProjects/vorpy')
    my_sys.net = Network(sys=my_sys, atoms=my_sys.atoms)
    my_sys.net.sort_atoms()
    for atom in my_sys.net.atoms:
        ffind_init_surfs(my_sys.net, atom)
    for surf in my_sys.net.surfs:
        ffind_init_edges(my_sys.net, surf)
