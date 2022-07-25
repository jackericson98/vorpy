from Visualize.visualize import *
from System.Network.vertex import Vertex
"""Example for vertices"""

# With vertices, we have __ different cases
cases = []
# Case 0: No overlap
cases.append([Atom([1, 0, 0], 0.25), Atom([0, 1, 0], 0.25), Atom([-1, 0, 0], 0.25), Atom([0, -1, 0], 0.25)])
# Case 1: One overlapping ball
cases.append([Atom([0.25, 1, 0], 0.25), Atom([0, 1, 0], 0.25), Atom([-1, 0, 0], 0.25), Atom([0, -1, 0], 0.25)])
# Case 2: Two overlapping balls
cases.append([Atom([0.3, 1, 0], 0.25), Atom([0, 1, 0], 0.25), Atom([-0.3, 1, 0], 0.25), Atom([0, -1, 0], 0.25)])
# Case 3: Three overlapping balls
cases.append([Atom([0.25, .75, 0], 0.25), Atom([0, 1, 0], 0.25), Atom([-.25, .75, 0], 0.25), Atom([.4, 0.5, 0], 0.25)])
# Case 4: Four overlapping balls
cases.append([Atom([1, 0, 0], 1), Atom([0, 1, 0], 1), Atom([-1, 0, 0], 1), Atom([0, -1, 0], 1)])


verts = []
for i in range(len(cases)):
    myVert = Vertex(atoms=cases[i])
    verts.append(myVert)


fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No Overlap", "Case 1: One overlapping ball", "Case 2: Two overlapping balls",
          "Case 3: Three overlapping balls", "Case 4: Four overlapping balls"]
for i in range(len(cases)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['r', 'b'], alpha=.1, Show=True, dfo=5)
    # plot_verts([verts[i]], fig=fig, ax=axn, dfo=5, Show=True)