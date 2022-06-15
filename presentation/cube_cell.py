from visualize import *
from build_network import calc_vertex

# Create atom objects from sets of points
atoms = [Atom([0, 0, 0], 13), Atom([5, 0, 0], 1), Atom([-5, 0, 0], 1), Atom([0, 5, 0], 1), Atom([0, -5, 0], 1),
         Atom([0, 0, 5], 1), Atom([0, 0, -5], 1)]

# Create a system of from the atoms
mySys = System()
mySys.atoms = atoms
mySys.net.atoms = atoms

# Calculate the vertices
vert_nums = [[1, 3, 5], [2, 3, 5], [2, 4, 5], [1, 4, 5], [1, 3, 6], [2, 3, 6], [2, 4, 6], [1, 4, 6]]
verts = []
for i in range(8):
    verts.append(calc_vertex([atoms[0], atoms[vert_nums[i][0]], atoms[vert_nums[i][1]], atoms[vert_nums[i][2]]]))

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
# Plot the vertices
plot_verts(verts, fig=fig, ax=ax)
# Plot the atom objects
plot_atoms(mySys.atoms, colors=['r', 'b'], fig=fig, ax=ax, Show=True)
