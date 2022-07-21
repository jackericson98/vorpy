import os
from Meshes.build_meshes import edge_trace, edge_trace1
from Network.build_network import *
from System.system import Surface, Atom, Edge
import matplotlib.pyplot as plt
from Presentation.Visualize.visualize import *
os.chdir("../..")


cases = [[], [], [], []]
# Case 0a: No overlap, equal atom sizes
cases[0] = [Atom([0.5, -np.sqrt(3)*np.pi/2, 0], 0.5), Atom([0.5, np.sqrt(3)*np.pi/2, 0], 0.5), Atom([1, 0, 0], 0.5)]
# Case 0bi: No overlap, one bigger, two equal smaller
# Case 0bii: No overlap, one smaller (equivalent to two larger)
# Case 0c: No overlap, all different sizes

# Case 1a: One overlapping set, all equal sizes
# Case 1bi: One overlapping set, two equal smaller/larger atoms overlapping
# Case 1bii: One overlapping set, two equal smaller/larger atoms not overlapping
# Case 1c: One overlapping set, two non equal smaller/larger atoms overlapping

# Case 2a: two overlapping sets, all atoms equal size
# Case 2bi: two overlapping sets, two equal smaller/larger atoms overlapping
# Case 2bii: two overlapping sets, two equal smaller/larger atoms not overlapping
# Case 2ci: two overlapping sets, all different sizes [s, m l]. (s-m-l, l-s-m, m-l-s)

# Case 3a: All overlapping, equal sizes
# Case 3b: All overlapping, 2 equal smaller/larger atoms
# Case 3c: All Overlapping, all different sizes


edges = []
# Create 2 dummy atoms for the vertices
va0, va1 = Atom([0, 0, 10], 0.5), Atom([0, 0, -10], 0.5)
for case_atoms in cases:
    # Calculate the two vertices
    v0 = calc_vert([va0] + case_atoms)
    v1 = calc_vert([va1] + case_atoms)
    # Create the edge
    myEdge = Edge(case_atoms, [v0, v1])
    # Get the edges points

    edges.append(myEdge)


fig = plt.figure(figsize=(20, 40))
titles = ["Case 0: No Overlap", "Case 1: Minimal overlap", "Case 2:"]
for i in range(len(cases)):
    axn = fig.add_subplot(int("23" + str(i + 1)), projection="3d", xlim=10)
    plot_atoms(cases[i], fig=fig, ax=axn, colors=['r', 'b'], alpha=.1)
    plot_surfs([surfs[i]], fig=fig, ax=axn, dfo=5)

plt.show()
