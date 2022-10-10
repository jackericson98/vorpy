import os
from System.system import System
from Visualize.visualize import plot_edges, plot_verts
import matplotlib.pyplot as plt
os.chdir("../..")


# Files
m_file = os.getcwd() + "./Data/test_data/Na_W_cluster5.pdb"
v_file = os.getcwd() + "./Data/test_data/Na5_verts1.txt"

# Get the System
sys = System(m_file, box_size=1.1, min_dist=0.05)
# sys.add_verts(v_file)
sys.net.find_verts()
sys.export_verts()
sys.net.build_surfs()

surfs = sys.net.surfs
for i in range(len(surfs)):
    # Calculate and print the running percentage for mesh calculations
    surfs[i].build_surfs()
    # Calculate and print the running percentage for mesh calculations
    percentage = int((i + 1) / len(surfs) * 100)
    print("\rBuilding Surfaces: ",
          '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')

# Plot the System
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# plot_atoms(sys.atoms[:1], fig=fig, ax=ax, alpha=.1, colors=['w' for i in range(len(sys.atoms))])
plot_verts(sys.net.verts, fig=fig, ax=ax, colors=['r' for i in range(len(sys.net.verts))])
plot_edges(sys.net.edges, fig=fig, ax=ax)
# plot_surfs(sys.atoms[0].surfs, fig=fig, ax=ax, alpha=1, simps=True)
plt.show()
sys.analyze()
sys.export(export_all=True)
