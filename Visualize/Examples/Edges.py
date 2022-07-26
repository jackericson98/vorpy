import os
from Visualize.visualize import *
from System.Network.edge import Edge
from System.Network.vertex import Vertex
from System.Network.surface import Surface
os.chdir("../..")


cases = []
cases.append([Atom([0.5, -0.5, 0], 0.25), Atom([0.5, 0.5, 0], 0.25), Atom([-0.5, 0, 0], 0.5)])
# Case 0: No overlap, equal atom sizes
# Case 0bi: No overlap, one bigger, two equal smaller
# Case 0bii: No overlap, one smaller (equivalent to two larger)
# Case 0c: No overlap, all different sizes

# Case 1: One overlapping set
cases.append([Atom([0.5, -np.sqrt(3)*np.pi/16, 0], 0.25), Atom([0.5, np.sqrt(3)*np.pi/16, 0], 0.5), Atom([-0.5, 0, 0], 0.5)])
# Case 1a: One overlapping set, all equal sizes
# Case 1bi: One overlapping set, two equal smaller/larger atoms overlapping
# Case 1bii: One overlapping set, two equal smaller/larger atoms not overlapping
# Case 1c: One overlapping set, two non equal smaller/larger atoms overlapping

# Case 2: Two overlapping sets
cases.append([Atom([0.25, -.3, 0], .25), Atom([0.25, .3, 0], .25), Atom([-.1, 0, 0], .25)])

# Case 2a: two overlapping sets, all atoms equal size
# Case 2bi: two overlapping sets, two equal smaller/larger atoms overlapping
# Case 2bii: two overlapping sets, two equal smaller/larger atoms not overlapping
# Case 2ci: two overlapping sets, all different sizes [s, m l]. (s-m-l, l-s-m, m-l-s)

# Case 3: All overlapping
cases.append([Atom([0.25, -np.sqrt(3)*np.pi/16, 0], .5), Atom([0.25, np.sqrt(3)*np.pi/16, 0], .5), Atom([-.25, 0, 0], .5)])
# Case 3a: All overlapping, equal sizes
# Case 3b: All overlapping, 2 equal smaller/larger atoms
# Case 3c: All Overlapping, all different sizes


edges, verts = [], []
# Create 2 dummy atoms for the vertices
va0, va1 = Atom([0, 0, 10], 0.5), Atom([0, 0, -10], 0.5)
# Create vertices, Edge object and calculate the edges points for each edge
for case_atoms in cases:
    # Calculate the two vertices
    v0 = Vertex([va0] + case_atoms)
    v1 = Vertex([va1] + case_atoms)
    verts.append([v0, v1])
    # Create the edge
    myEdge = Edge(case_atoms, [v0, v1])
    # Get the edges points
    myEdge.calc_points(min_dist=.1)
    # Add the edge to the list of
    edges.append(myEdge)


fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No Overlap", "Case 1: One overlap", "Case 2: Two overlaps", "Case 3: all overlapping"]
for i in range(len(cases)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    axn.set_title(titles[i])
    plot_verts(verts[i], fig=fig, ax=axn, grid=True)
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['r', 'b'], alpha=.1, dfo=2, grid=True)
    plot_edges([edges[i]], fig=fig, ax=axn, dfo=5, grid=True)

plt.show()
