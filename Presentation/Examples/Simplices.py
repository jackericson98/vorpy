from System.system import *
from Meshes.build_meshes import *
from Presentation.Visualize.visualize import *

# Create the atoms for the mesh
atoms = [Atom([-5, 0, 0], 3), Atom([5, 0, 0], 1)]
mySurf = Surface(atoms)
# Make the mesh
make_mesh(surf=mySurf, min_dist=.1, circ_mesh=True, radius=10)
# Find the simplices
mySurf.simps = find_simps(mySurf)

# Plot the surface and the simplices
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# ax.scatter()
plot_surfs([mySurf], simps=True, fig=fig, ax=ax, Show=True)
