from System.system import *
from Visualize.visualize import *

# Create the atoms for the mesh
atoms = [Atom([-3, 0, 0], 3), Atom([2, 0, 0], 1)]
mySurf = Surface(atoms)
# Make the mesh
mySurf.build(min_dist=0.5)


# Plot the surface and the simplices
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# ax.scatter()
plot_atoms(mySurf.atoms, fig=fig, ax=ax, grid=True, alpha=0.1)
plot_simps(mySurf, Show=True, fig=fig, ax=ax, dfo=5)

