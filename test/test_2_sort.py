from System.system import System
from System.Network.network import Network


"""
____________________________________________Atom Sorting________________________________________________________________
Tests the atom sorting schemes for a loaded network
"""


def test_net_calc_box():
    my_sys = System('../Data/test_data/cambrin.pdb')
    my_net = Network(sys=my_sys, atoms=my_sys.atoms)
    my_sys.net = my_net
    my_sys.net.sort_atoms()
    assert my_net.box == [[4.678000000000004, 8.424, -16.369249999999994], [71.356, 68.952, 52.170249999999996]]
    assert my_net.atoms_box == [[15.791, 18.512, -4.946], [60.243, 58.864, 40.747]]


def test_net_sort_atoms():
    my_sys = System('../Data/test_data/cambrin.pdb')
    my_net = Network(sys=my_sys, atoms=my_sys.atoms)
    my_sys.net = my_net
    my_sys.net.sort_atoms()
    assert my_net.atoms[45].box == [32, 31, 29]
    assert my_net.sub_box_size == [1.1113, 1.0088, 1.1423249999999998]
    assert len(my_net.sub_boxes) == 60
