import os
from System.system import Surface
from Presentation.Visualize.visualize import *
os.chdir("../..")

cases = []
# Case 0: No overlap (equal radii) - distance(a0, a1) < a0.rad + a1.rad
cases.append([Atom([-2, 0, 0], 1.4999999999999), Atom([2, 0, 0], 1.5)])
# Case 1: No overlap - distance(a0, a1) < a0.rad + a1.rad
cases.append([Atom([-2, 0, 0], 1.5), Atom([2, 0, 0], 0.5)])
# Case 2: Minimal overlap - a0.rad, a1.rad < distance(a0, a1) < a0.rad + a1.rad
cases.append([Atom([-0.95, 0, 0], 1.5), Atom([0.95, 0, 0], 0.5)])
# Case 3: More than half of one atom - a0.rad < distance(a0, a1) < a1.rad
cases.append([Atom([-0.74, 0, 0], 1.5), Atom([0.74, 0, 0], 0.5)])
# Case 4: Overlap is more than half of both atoms - distance(a0, a1) < a0.rad, a1.rad
cases.append([Atom([-0.65, 0, 0], 1.5), Atom([.65, 0, 0], 0.5)])
# Case 5: Full encapsulation of one atom - distance(a0, a1) + a0.rad < a1.rad
cases.append([Atom([0, 0, 0], 1.5), Atom([0.1, 0, 0], 0.25)])

surfs = []
for case_atoms in cases:
    mySurf = Surface(case_atoms)
    mySurf.build()
    surfs.append(mySurf)



fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No overlap (equal radii)",
          "Case 1: No overlap",
          "Case 2: Minimal overlap",
          "Case 3: More than half of one atom",
          "Case 4: Overlap is more than half of both atoms",
          "Case 5: Full encapsulation of one atom"]
for i in range(len(cases)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['r', 'b'], alpha=.1)
    plot_surfs([surfs[i]], simps=True, fig=fig, ax=axn, dfo=5)

plt.show()
