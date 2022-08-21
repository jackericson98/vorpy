import numpy as np

from System.system import *

# Create 2 sets of 10 subcases for comparison
radii = [[0.4, 1], [0.35, 1], [0.3, 1], [0.275, 1], [0.25, 1], [0.225, 1], [0.2, 1], [0.175, 1], [0.15, 1], [0.1, 1]]
locs = [[0, 0.05], [0, 0.1], [0, 0.2], [0, 0.25], [0, 0.3], [0, 0.35], [0, 0.5], [0, 1], [0, 2], [0, 5]]

# Create Atoms and Surfaces for the different subcases
atoms, surfs, syss = [], [], []
for j in range(len(locs)):
    for i in range(len(radii)):
        myAtoms = [Atom([locs[i][0], 0, 0], radii[j][0]), Atom([locs[i][1], 0, 0], radii[j][1])]
        sys = System(myAtoms)
        syss.append(sys)
        atoms.append(myAtoms)
        surfs.append(Surface(myAtoms, sys.net))

ang_incs = np.linspace(0, np.pi, 1000)
for surf in surfs:
    root_diff = np.inf
    crit_ang = 0
    for ang in ang_incs:
        point = [surf.atoms[0].rad*np.cos(ang), surf.atoms[0].rad*np.sin(ang), 0]
        roots = surf.calc_surf_point1(point)
        my_diff = np.inf
        if roots is not None and len(roots) > 1:
            my_diff = abs(abs(roots[0]) - abs(roots[1]))
        if my_diff < root_diff:
            crit_ang = ang
            root_diff = my_diff
    print(surfs.index(surf), crit_ang, crit_ang/np.pi, surf.atoms[0].loc, surf.atoms[0].rad, surf.atoms[1].loc, surf.atoms[1].rad, "\n")