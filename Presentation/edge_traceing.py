from visualize import plot_atoms
from build_network import find_v0, find_edges, calc_edge1
from objects import System

mySys = System()
plot_atoms(mySys.atoms)
v0 = find_v0(mySys.net)

v1 = find_edges(v0, mySys.net)
print(v1)
