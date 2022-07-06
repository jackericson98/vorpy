from visualize import *
from build_network import calc_vertex, build_network
import matplotlib.patches as mpatches
from objects import System, Surface

# Create the figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# Create atom objects from sets of points
atoms = [[[0, 0, 0], 2.5]]
dist = 5
rad = .5

atoms += [[[dist, 0, 0], rad + 1], [[-dist - 4, 0, 0], rad], [[0, dist -1, 0], rad],[[0, -dist - 15, 0], rad], [[0, 0, dist], rad],
          [[0, 0, -dist], rad]]
# Create a system of from the atoms
mySys = System(atoms)

# Plot the atom objects
plot_atoms(mySys.atoms, fig=fig, ax=ax)

# Calculate the vertices
vert_nums = [[1, 3, 5], [2, 3, 5], [2, 4, 5], [1, 4, 5], [1, 3, 6], [2, 3, 6], [2, 4, 6], [1, 4, 6]]
verts = []
for i in range(8):
    vn = calc_vertex([mySys.atoms[0], mySys.atoms[vert_nums[i][0]], mySys.atoms[vert_nums[i][1]], mySys.atoms[vert_nums[i][2]]])
    if vn:
        verts.append(vn)


# Build the network of vertices
build_network(mySys)
for vert in mySys.net.verts:
    print(vert.loc)
# Plot the vertices
red_patch = mpatches.Patch(color='red', label='Vertices')
blue_patch = mpatches.Patch(color='blue', label='Atoms')
black_patch = mpatches.Patch(color='black', label='Interstitial Spheres')
ax.legend(handles=[blue_patch, red_patch, black_patch], loc='upper right')
plot_verts(mySys.net.verts, fig=fig, ax=ax, dfo=10, Show=True, colors=['r' for i in range(4)], grid=True)
f = calc_surf(mySys.atoms[:2])
mySurf = Surface(f, [mySys.atoms[0], mySys.atoms[1]])

# surf = make_mesh(mySurf)
