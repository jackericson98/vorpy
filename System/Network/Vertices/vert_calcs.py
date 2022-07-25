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

