import os
from Visualize.visualize import *
from System.Network.edge import Edge
from System.Network.vertex import Vertex
from System.system import System, Atom
os.chdir("../..")


cases = []
# Case 0a: No overlap equal sizes
cases.append([Atom([0.5, -0.75, 0], 0.5), Atom([0.5, 0.75, 0], 0.5), Atom([-0.5, 0, 0], 0.5)])
# Case 0b: No overlap unequal sizes
cases.append([Atom([0.5, -0.75, 0], 0.75), Atom([0.5, 0.75, 0], 0.5), Atom([-0.5, 0, 0], 0.25)])
# Case 1: One overlapping set
cases.append([Atom([0.5, -0.5, 0], 0.75), Atom([0.5, 0.5, 0], 0.5), Atom([-0.5, 0, 0], 0.25)])
# Case 2: Two overlapping sets
# cases.append([Atom([0.5, -1, 0], 1), Atom([0.5, 1, 0], 0.75), Atom([0, 0, 0], 0.5)])
# Case 3a: All overlapping not going through the atoms
cases.append([Atom([0.5, -1, 0], 0.95), Atom([0.5, 0.5, 0], 0.65), Atom([-0.5, 0, 0], 0.5)])
# Case 3b: All overlapping going through the atoms
cases.append([Atom([0.5, -0.5, 0], 0.75), Atom([0.5, 0.5, 0], 0.5), Atom([0, 0.1, 0], 0.25)])

edges, verts = [], []
# Create 2 dummy atoms for the vertices
va0, va1 = Atom([0, 0, 5], 0.5), Atom([0, 0, -5], 0.5)
syss = []
# Create vertices, Edge object and calculate the edges points for each edge
for case_atoms in cases:
    syss.append(System(case_atoms + [va0, va1]))
    # Calculate the two vertices
    v0 = Vertex([va0] + case_atoms)
    v0.calc_vert()
    v1 = Vertex([va1] + case_atoms)
    v1.calc_vert()
    verts.append([v0, v1])
    # Create the edge
    myEdge = Edge(case_atoms, verts=[v0, v1])
    # Get the edges points
    myEdge.build(min_dist=.05)
    # Add the edge to the list of
    edges.append(myEdge)


fig = plt.figure(figsize=(10, 5))
fig.tight_layout()
titles = ["Case 0a: No overlap, all equal", "Case 0b: No Overlap, different radii", "Case 1: One overlapping",
          "Case 2: Two overlapping", "Case 3a: All overlapping, hole", "Case 3b: All overlapping, no hole"]
for i in range(len(cases)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    plot_verts(verts[i], fig=fig, ax=axn, grid=True)
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['b', 'b', 'b'], alpha=.1, dfo=3, grid=True, res=6)
    plot_edges([edges[i]], fig=fig, ax=axn, dfo=3, grid=False, alpha=1)
fig.suptitle('Edges', fontsize=20)
plt.show()
