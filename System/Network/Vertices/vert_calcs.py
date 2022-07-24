from System.sys_calcs import *
from System.system import *
"""Calculator functions"""


# Calculate direction function. Takes in a vertex and an edge and returns True if it is facing the center
def calc_dir(edge):
    # Grab the previous vertex
    vn_1 = edge.verts[0]
    # Find ak and copy it
    ak = None
    for atom in vn_1.atoms:
        if not {atom}.issubset(edge.atoms):
            ak = atom
    akp = ak
    # Find the direction toward the center of the edge
    r0 = [edge.loc[0] - vn_1.loc[0], edge.loc[1] - vn_1.loc[1], edge.loc[2] - vn_1.loc[2]]
    r0_mag = np.sqrt(r0[0]**2 + r0[1]**2 + r0[2]**2)
    r0_hat = [r0[0]/r0_mag, r0[1]/r0_mag, r0[2]/r0_mag]
    # Move the copy toward the center of the edge.
    akp.loc = [akp.loc[0] + r0_hat[0]*0.1, akp.loc[1] + r0_hat[1]*0.1, akp.loc[2] + r0_hat[2]*0.1]
    # Calculate the new vertex made by akp
    vkp = calc_vert(edge.atoms + [akp])
    while not vkp:
        akp.loc = [akp.loc[0] - r0_hat[0]*0.01, akp.loc[1] - r0_hat[1]*0.1, akp.loc[2] - r0_hat[2]*0.1]
        vkp = calc_vert(edge.atoms + [akp])
    # If the new inscribed sphere overlaps with ak, flip the direction of tang_hat
    if calc_dist(ak.loc, vkp.loc) - (ak.rad + vkp.rad) < 0:
        return False
    return True


# Calculate relative length function. Takes in 3 points and returns a float value for the relative distance
def calc_rel_dist(v0, v1, edge):
    # Grab the center
    c = np.array(edge.loc)
    # Find the distances between the 3 points
    r0, r1, r2 = np.linalg.norm(c - np.array(v0.loc)), np.linalg.norm(c - np.array(v1.loc)), \
                 np.linalg.norm(np.array(v0.loc) - np.array(v1.loc))
    # Cases 1 and 2: r0 > r1 > r2 and r0 > r2 > r1
    if r0 >= r1 and r0 > r2 and edge.dir:
        rel_dist = r2
    # Cases 3 and 4: r1 > r0 > r2 and r1 > r2 > r0
    elif r1 > r0 and r1 > r2 and not edge.dir:
        rel_dist = r2
    # Cases 5 and 6: r2 > r0 > r1 and r2 > r1 > r0
    elif r2 > r0 and r2 > r1 and edge.dir:
        rel_dist = r0 + r1
    # All other cases should not give a distance
    else:
        rel_dist = np.inf
    # Return the relative distance
    return rel_dist


# Move sphere function. Takes in a location, an Atom object and a direction and updates the Atom's location
def move(loc, atom, to_home=False):
    # Change whether we are adding or subtracting the location to the sphere's location.
    d = 1
    if not to_home:
        d = -1
    # Update the atom's location
    atom.loc[0] = atom.loc[0] + d * loc[0]
    atom.loc[1] = atom.loc[1] + d * loc[1]
    atom.loc[2] = atom.loc[2] + d * loc[2]


# Calculate vertex function. Takes in 4 atoms, calculates the center and radius of the inscribed sphere and returns them
def calc_vert(atoms):
    # The real location and radius of the base sphere
    l1, R1 = atoms[0].loc, atoms[0].rad
    # Set the radii and x, y, z values for the 3 spheres
    R2, R3, R4 = atoms[1].rad, atoms[2].rad, atoms[3].rad
    x2, y2, z2 = atoms[1].loc[0] - l1[0], atoms[1].loc[1] - l1[1], atoms[1].loc[2] - l1[2]
    x3, y3, z3 = atoms[2].loc[0] - l1[0], atoms[2].loc[1] - l1[1], atoms[2].loc[2] - l1[2]
    x4, y4, z4 = atoms[3].loc[0] - l1[0], atoms[3].loc[1] - l1[1], atoms[3].loc[2] - l1[2]

    # Calculate our System of linear equations coefficients
    a1, b1, c1, d1, f1 = 2 * x2, 2 * y2, 2 * z2, 2 * (R2 - R1), R1 ** 2 - R2 ** 2 + x2 ** 2 + y2 ** 2 + z2 ** 2
    a2, b2, c2, d2, f2 = 2 * x3, 2 * y3, 2 * z3, 2 * (R3 - R1), R1 ** 2 - R3 ** 2 + x3 ** 2 + y3 ** 2 + z3 ** 2
    a3, b3, c3, d3, f3 = 2 * x4, 2 * y4, 2 * z4, 2 * (R4 - R1), R1 ** 2 - R4 ** 2 + x4 ** 2 + y4 ** 2 + z4 ** 2

    A, B, C, d, f = [a1, a2, a3], [b1, b2, b3], [c1, c2, c3], [d1, d2, d3], [f1, f2, f2]

    # Calculate the ranks of the matrices
    ABC_rank = np.linalg.matrix_rank([A, B, C])
    m_rank = np.linalg.matrix_rank([A, B, C, d])
    f_rank = np.linalg.matrix_rank([A, B, C, d, f])

    # Calculate the F values
    F, F10, F11, F20, F21, F30, F31 = a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1, \
                                      b1*c2*f3 - b1*c3*f2 - b2*c1*f3 + b2*c3*f1 + b3*c1*f2 - b3*c2*f1, \
                                      -b1*c2*d3 + b1*c3*d2 + b2*c1*d3 - b2*c3*d1 - b3*c1*d2 + b3*c2*d1, \
                                      -a1*c2*f3 + a1*c3*f2 + a2*c1*f3 - a2*c3*f1 - a3*c1*f2 + a3*c2*f1, \
                                      a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1, \
                                      a1*b2*f3 - a1*b3*f2 - a2*b1*f3 + a2*b3*f1 + a3*b1*f2 - a3*b2*f1, \
                                      -a1*b2*d3 + a1*b3*d2 + a2*b1*d3 - a2*b3*d1 - a3*b1*d2 + a3*b2*d1
    # Catch for F = 0.
    if F == 0:
        return
    # Instantiate our root arrays
    xs, ys, zs, Rs = [], [], [], []
    verts = []
    # Case 1:
    if ABC_rank == 3 and m_rank == 3 and f_rank == 3:
        # Calculate the radius polynomial coefficients
        a = ((F11 ** 2 + F21 ** 2 + F31 ** 2) / F ** 2) - 1
        b = (2 * (F10 * F11 + F20 * F21 + F30 * F31) / F ** 2) - 2 * R1
        c = ((F10 ** 2 + F20 ** 2 + F30 ** 2) / F ** 2) - R1 ** 2
        # If the discriminant is positive, find the real positive roots of the quadratic
        if -4*a*c + b**2 > 0:
            Rs = [R for R in np.roots([a, b, c]) if np.isreal(R) and R > 0]
        # Instantiate the verts array
        verts = []
        # Go through each radius and calculate the vertex
        for R in Rs:
            x, y, z = F10/F + R*F11/F, F20/F + R*F21/F, F30/F + R*F31/F
            # Move the vertex back to the actual location of the atoms
            verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))

    # Case 2:
    elif ABC_rank == 2 and m_rank == 3 and f_rank == 3:
        # Case 2 subcases:
        # Case 2.1
        if np.linalg.matrix_rank([A, B, d]) == 3:
            # Calculate the z value polynomial coefficients
            a = F**2 + F11**2 + F21**2 - F31**2
            b = 2*(F10*F11 + F20*F21 - F30*F31 - F*F31*R1)
            c = F10**2 + F20**2 - (F30 + F*R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                zs = [z for z in np.roots([a, b, c]) if np.isreal(z) and z > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for z in zs:
                x, y, R = F10 / F + z * F11 / F, F20 / F + z * F21 / F, F30 / F + z * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))

        # Case 2.2
        elif np.linalg.matrix_rank([A, d, C]) == 3:
            # Calculate the z value polynomial coefficients
            a = F ** 2 + F11 ** 2 - F21 ** 2 + F31 ** 2
            b = 2 * (F10 * F11 - F20 * F21 + F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F30 ** 2 - (F20 + F * R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                ys = [y for y in np.roots([a, b, c]) if np.isreal(y) and y > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for y in ys:
                x, R, z = F10 / F + y * F11 / F, F20 / F + y * F21 / F, F30 / F + y * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))

        # Case 2.3
        elif np.linalg.matrix_rank([d, B, C]):
            # Calculate the z value polynomial coefficients
            a = F ** 2 + F11 ** 2 + F21 ** 2 - F31 ** 2
            b = 2 * (F10 * F11 + F20 * F21 - F30 * F31 - F * F31 * R1)
            c = F10 ** 2 + F20 ** 2 - (F30 + F * R1)
            # If the discriminant is positive, find the real positive roots of the quadratic
            if -4 * a * c + b ** 2 > 0:
                xs = [x for x in np.roots([a, b, c]) if np.isreal(x) and x > 0]
            # Instantiate the verts array
            verts = []
            # Go through each radius and calculate the vertex
            for x in xs:
                R, y, z = F10 / F + x * F11 / F, F20 / F + x * F21 / F, F30 / F + x * F31 / F
                # Move the vertex back to the actual location of the atoms
                verts.append(Vertex([x + l1[0], y + l1[1], z + l1[2]], R, atoms=atoms))
    # If no verts are found return None
    if not verts:
        return
    else:
        if len(verts) == 2 and verts[0].rad > verts[1].rad:
            return verts[1]
        return verts[0]
