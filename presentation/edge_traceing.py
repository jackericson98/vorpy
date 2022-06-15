from load_system import random_system
from visualize import plot_atoms
from build_network import find_v0, find_edges

mySys = random_system(anums=20)
# plot_atoms(mySys.atoms)
v1 = find_v0(mySys.net)
find_edges(v1, mySys.net)
