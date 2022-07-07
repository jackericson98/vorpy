import os
from objects import System
from build_network import build_network
os.chdir("..")

# Get the system
sys = System("./Data/test_data/Na_W_cluster5.pdb")
build_network(sys)
# print(len(sys.net.verts))
for vert in sys.net.verts:
    ndxs = []
    for atom in vert.atoms:
        # print(atom.loc)
        ndxs.append(sys.atoms.index(atom))
    ndxs.sort()
    # print(ndxs, vert.loc)
