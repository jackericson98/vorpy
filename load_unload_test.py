from System.system import System
import os

print(list("[]"))
sys1 = System(os.getcwd() + "/Data/test_data/Na5.pdb")

sys1.load_net(file=os.getcwd() + "/Data/User_data/Na5176/Na5_net.csv")

sys2 = System(os.getcwd() + "./Data/test_data/Na5.pdb")

sys2.build_network()

