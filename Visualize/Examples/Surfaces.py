import os
from System.system import System
from Visualize.visualize import *
os.chdir("../..")

# Create atoms for vertices and edges
vert_atoms = [Atom([1, 5, 5], 1), Atom([1, -5, 5], 1), Atom([1, 5, -5], 1), Atom([1, -5, -5], 1)]

cases = []
# Case 0: No overlap (equal radii) - distance(a0, a1) < a0.rad + a1.rad
cases.append(System([Atom([-2, 0, 0], 1.4999999999999), Atom([2, 0, 0], 1.5)] + vert_atoms))
# Case 1: No overlap - distance(a0, a1) < a0.rad + a1.rad
cases.append(System([Atom([-2, 0, 0], 1.5), Atom([2, 0, 0], 0.5)] + vert_atoms))
# Case 2: Minimal overlap - a0.rad, a1.rad < distance(a0, a1) < a0.rad + a1.rad
cases.append(System([Atom([-0.95, 0, 0], 1.5), Atom([0.95, 0, 0], 0.5)] + vert_atoms))
# Case 3: More than half of one atom - a0.rad < distance(a0, a1) < a1.rad
cases.append(System([Atom([-0.75, 0, 0], 1.5), Atom([0.75, 0, 0], 0.5)] + vert_atoms))
# Case 4: Overlap is more than half of both atoms - distance(a0, a1) < a0.rad, a1.rad
cases.append(System([Atom([-0.65, 0, 0], 1.5), Atom([.65, 0, 0], 0.5)] + vert_atoms))
# Case 5: Full encapsulation of one atom - distance(a0, a1) + a0.rad < a1.rad
# cases.append(System([Atom([0, 0, 0], 1.5), Atom([0.1, 0, 0], 0.25)] + vert_atoms))


for sys in cases:
    sys.build_network(0.05)


fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No overlap (equal radii)",
          "Case 1: No overlap",
          "Case 2: Minimal overlap",
          "Case 3: More than half of one atom",
          "Case 4: Overlap is more than half of both atoms",
          "Case 5: Full encapsulation of one atom"]
for i in range(len(cases)):
    sys = cases[i]
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    plot_atoms(sys.atoms[:2], fig=fig, ax=axn, colors=['r', 'b'], alpha=.1)
    plot_surfs(sys.net.surfs, simps=True, fig=fig, ax=axn, dfo=5)

plt.show()