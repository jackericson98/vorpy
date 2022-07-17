from Presentation.Visualize.visualize import plot_atoms
from Network.build_network import find_v0, find_edges
from System.objects import System

mySys = System()
plot_atoms(mySys.atoms)
v0 = find_v0(mySys.net)

v1 = find_edges(v0, mySys.net)
print(v1)
