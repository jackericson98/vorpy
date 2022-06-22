from load_system import random_system
from visualize import plot_atoms
from build_network import find_v0, find_edges, calc_edge1

mySys = random_system(anums=20, )
plot_atoms(mySys.atoms)
v0 = find_v0(mySys.net)

v1 = find_edges(v0, mySys.net, nedges=1)
print(v1)
