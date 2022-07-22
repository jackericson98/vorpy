import os
from Meshes.build_meshes import make_mesh
from Network.build_network import *
from System.system import Surface, Atom
import matplotlib.pyplot as plt
from Presentation.Visualize.visualize import *
os.chdir("../..")

cases = []
# Case 0: No overlap - distance(a0, a1) < a0.rad + a1.rad
cases.append([Atom([-2, 0, 0], 1.5), Atom([2, 0, 0], 0.5)])
# Case 1: Minimal overlap - a0.rad, a1.rad < distance(a0, a1) < a0.rad + a1.rad
cases.append([Atom([-1, 0, 0], 1.5), Atom([1, 0, 0], 1)])
# Case 2: More than half of one atom - a0.rad < distance(a0, a1) < a1.rad
cases.append([Atom([-1, 0, 0], 2.25), Atom([1, 0, 0], 1)])
# Case 3: Overlap is more than half of both atoms - distance(a0, a1) < a0.rad, a1.rad
cases.append([Atom([-1, 0, 0], 4), Atom([.56, 0, 0], 2.5)])
# Case 4: Full encapsulation of one atom - distance(a0, a1) + a0.rad < a1.rad
cases.append([Atom([-.25, 0, 0], 1.5), Atom([.25, 0, 0], 0.5)])

surfs = []
for case_atoms in cases:
    mySurf = Surface(case_atoms)
    surfs.append(mySurf)
    make_mesh(mySurf, 0.1)


fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No Overlap", "Case 1: Minimal overlap", "Case 2:"]
for i in range(len(cases)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['r', 'b'], alpha=.1)
    plot_surfs([surfs[i]], simps=True, fig=fig, ax=axn, dfo=5)

plt.show()
