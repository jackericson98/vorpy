from visualize import *
from build_network import calc_vertex, build_network
import matplotlib.patches as mpatches

# Create the figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# Create atom objects from sets of points
atoms = [Atom([0, 0, 0], 1.5)]
dist = 5
rad = .5

atoms += [Atom([dist, 0, 0], rad+1), Atom([-dist, 0, 0], rad), Atom([0, dist, 0], rad), Atom([0, -dist, 0], rad),
          Atom([0, 0, dist], rad), Atom([0, 0, -dist], rad)]
# Plot the atom objects
plot_atoms(atoms, fig=fig, ax=ax)
# Create a system of from the atoms
mySys = System()
mySys.atoms = atoms
mySys.net.atoms = atoms

# Calculate the vertices
vert_nums = [[1, 3, 5], [2, 3, 5], [2, 4, 5], [1, 4, 5], [1, 3, 6], [2, 3, 6], [2, 4, 6], [1, 4, 6]]
verts = []
for i in range(8):
    vn = calc_vertex([atoms[0], atoms[vert_nums[i][0]], atoms[vert_nums[i][1]], atoms[vert_nums[i][2]]])
    verts.append(vn)
    print(vn.rad)


# Build the network of vertices
build_network(mySys)

# Plot the vertices
red_patch = mpatches.Patch(color='red', label='Vertices')
blue_patch = mpatches.Patch(color='blue', label='Atoms')
black_patch = mpatches.Patch(color='black', label='Interstitial Spheres')
ax.legend(handles=[blue_patch, red_patch, black_patch], loc='upper right')
plot_verts(verts, fig=fig, ax=ax, plot_spheres=True, Show=True)


