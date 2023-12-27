from numpy import array as ar
import matplotlib.pyplot as plt
from Visualize.mpl_visualize import plot_atoms, plot_verts, plot_edges, plot_surfs
from System.Network.verts.calc_vert import calc_vert
from System.Network.edges.build_edge import build_edge


"""
Edge plotting code. Choose an edge type below.

    1. Flat Edge: All balls equal
    2. Semi Curved Edge 1: No overlap, 2 balls the same size, relatively close
    3. Semi Curved Edge 2: No overlap, all balls different sizes, relatively close
    4. Semi Curved Edge 3: Some overlap, 2 balls the same size, very close
    5. Semi curved Edge 4. Some overlap, all balls different sizes, very close
    ...
"""

# Choose here
edge_choice = 5

# Flat edge
if edge_choice == 1:
    rads = 1.0, 1.0, 1.0
# Semi Curved edge 1
elif edge_choice == 2:
    rads = 1.0, 0.75, 0.75
# Semi Curved Edge 2
elif edge_choice == 3:
    rads = 1.0, 0.75, 0.5
# Semi Curved Edge 3
elif edge_choice == 4:
    rads = 3.0, 2.5, 2.5
# Semi Curved Edge 4
elif edge_choice == 5:
    rads = 3.5, 3.0, 2.5


# Set the radii and distances for the surrounding atoms
r, d = 0.5, 15.0
my_vert_atoms = [(ar([0.0, d, 0.0]), r), (ar([0.0, -d, 0.0]), r)]
# Create the edge Atoms
edge_atoms = [(ar([2.5, 0.0, 0.0]), rads[0]), (ar([-2.5, 0.0, 0.0]), rads[1]), (ar([0.0, 0.0, 2.5]), rads[2])]

# Calculate the vertices
vert_atoms = [edge_atoms + [_] for _ in my_vert_atoms]
my_verts = [calc_vert(ar([_[0] for _ in my_atoms]), ar([_[1] for _ in my_atoms])) for my_atoms in vert_atoms]

# Calculate the Edge
my_edge = build_edge(alocs=ar([_[0] for _ in edge_atoms]), arads=ar([_[1] for _ in edge_atoms]),
                     vlocs=ar([_[0] for _ in my_verts]), res=0.5)


# Plot everything
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot_edges([my_edge[0]], fig=fig, ax=ax, colors=['k'], thickness=1)
plot_verts([_[0] for _ in my_verts], [_[1] for _ in my_verts], fig=fig, ax=ax)
plot_atoms(alocs=[_[0] for _ in edge_atoms], arads=[_[1] for _ in edge_atoms], alpha=0.3, fig=fig, ax=ax, Show=True)


