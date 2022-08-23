import numpy as np
import matplotlib.pyplot as plt

from System.system import *

# Create 2 sets of 10 subcases for comparison
radii = [[0.35, 1], [0.325, 1], [0.3, 1], [0.275, 1], [0.25, 1], [0.225, 1], [0.2, 1], [0.175, 1], [0.15, 1], [0.45, 1]]
locs = [[0, x] for x in np.linspace(0.5, 1, 99)]
print(locs)
# Create Atoms and Surfaces for the different subcases
atoms, surfs, syss = [], [], []
for j in range(len(locs)):
    for i in range(len(radii)):
        myAtoms = [Atom([locs[j][0], 0, 0], radii[i][0]), Atom([locs[j][1], 0, 0], radii[i][1])]
        sys = System(myAtoms)
        syss.append(sys)
        atoms.append(myAtoms)
        surfs.append(Surface(myAtoms, sys.net))

xs, ys, yts = [], [], []
ang_incs = np.linspace(0, np.pi, 1000)
for surf in surfs:
    root_diff = np.inf
    crit_ang = 0
    for ang in ang_incs:
        point = [surf.atoms[0].rad*np.cos(ang), surf.atoms[0].rad*np.sin(ang), 0]
        roots = surf.calc_surf_point(point, roots=True)
        my_diff = np.inf
        if roots is not None and len(roots) > 1:
            my_diff = abs(abs(roots[0]) - abs(roots[1]))
        if my_diff < root_diff:
            crit_ang = ang
            root_diff = my_diff
    xs.append(surf.atoms[1].loc[0])
    ys.append(crit_ang)
    yts.append(surf.atoms[1].loc[0])

for i in range(len(ys)):
    r0, r1 = surfs[i].atoms[0].rad, surfs[i].atoms[1].rad
    yts[i] = np.tan(xs[i]*20 + np.pi/2 + r0/r1) + (r0 ** 0.5/r1 + 1)
    if abs(ys[i]) > 20:
        ys[i] = 0

# Seperate the plots
plots = []
for i in range(10):
    plots.append([xs[i*100:(i+1)*100], ys[i*100:(i+1)*100]])

for i in range(10):
    plt.plot(plots[i][0], plots[i][1])
plt.legend(radii)

plt.show()
