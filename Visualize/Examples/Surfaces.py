import os
from System.system import System, Atom
from Visualize.visualize import *
os.chdir("../..")

# cases = []
# # Case 0: No overlap (equal radii) - distance(a0, a1) < a0.rad + a1.rad
# cases.append(System([Atom([-1, 0, 0], 1.5), Atom([3, 0, 0], 1.5)] + vert_atoms))
# # Case 1: No overlap - distance(a0, a1) < a0.rad + a1.rad
# cases.append(System([Atom([-2, 0, 0], 1.5), Atom([2, 0, 0], 0.5)] + vert_atoms))
# # Case 2: Minimal overlap - a0.rad, a1.rad < distance(a0, a1) < a0.rad + a1.rad
# cases.append(System([Atom([-0.95, 0, 0], 1.5), Atom([0.95, 0, 0], 0.5)] + vert_atoms))
# # Case 3: More than half of one atom - a0.rad < distance(a0, a1) < a1.rad
# cases.append(System([Atom([-1, 0, 0], 2), Atom([1.5, 0, 0], 0.5)] + vert_atoms))
# # Case 4: Overlap is more than half of both atoms - distance(a0, a1) < a0.rad, a1.rad
# cases.append(System([Atom([-0.65, 0, 0], 1.5), Atom([.65, 0, 0], 0.5)] + vert_atoms))
# # Case 5: Full encapsulation of one atom - distance(a0, a1) + a0.rad < a1.rad
# cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([-0.24, 0, 0], 0.25)] + vert_atoms))
# Create atoms for vertices and edges
vert_atoms = [Atom([2, 5, 5], 2.5), Atom([2, -5, 5], 2.5), Atom([2, 5, -5], 2.5), Atom([2, -5, -5], 2.5)]


cases = []
cases.append(System([Atom([-1, 0, 0], .5), Atom([3, 0, 0], .5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([1, 0, 0], 0.5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([0.6, 0, 0], 0.5)] + vert_atoms))
# cases.append(System([Atom([-15, 0, 0], 15), Atom([0.5, 0, 0], 1), Atom([10, 5, 5], 2.5), Atom([10, -5, 5], 2.5), Atom([10, 5, -5], 2.5), Atom([10, -5, -5], 2.5)]))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([0.4, 0, 0], 0.5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([0.25, 0, 0], 0.5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([0.1, 0, 0], 0.5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([-.1, 0, 0], 0.5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([-0.25, 0, 0], 0.5)] + vert_atoms))
cases.append(System([Atom([-1.5, 0, 0], 1.5), Atom([-0.4, 0, 0], 0.5)] + vert_atoms))


for i in range(len(cases)):
    cases[i].net.build_surfs(0.05)


fig = plt.figure(figsize=(20, 40))
titles = ["Equal Radii",
          "No overlap",
          "Very close",
          "Minimal overlap",
          "Minimal Overlap",
          "Minimal Overlap",
          "Case 4: Overlap is more than half of both atoms",
          "Case 5: Full encapsulation of one atom", "", ""]
for i in range(len(cases)):

    sys = cases[i]
    axn = fig.add_subplot(int("33" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    try:
        plot_verts(sys.net.surfs[0].verts, fig=fig, ax=axn)
        plot_atoms(sys.net.surfs[0].atoms, fig=fig, ax=axn, colors=['r', 'b'], alpha=.5, res=3)
        plot_edges(sys.net.surfs[0].edges, fig=fig, ax=axn)
        plot_surfs(sys.net.surfs[:1], simps=True, fig=fig, ax=axn, dfo=1.5)
    except IndexError:
        continue
plt.show()
