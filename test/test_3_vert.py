from System.system import System
from System.Network.network import Network
from System.Network.net_funcs.find_verts import *

"""
___________________________________________________Voronoi Calculations_________________________________________________
Tests: calc_vert, verify_site, find_verts
"""


def test_calc_vor_vert():
    my_sys = System("../Data/test_data/Na5.pdb")
    loc, rad, loc2, rad2 = calc_vert(locs=[_.loc for _ in [my_sys.atoms[_] for _ in [0, 1, 2, 23]]],
                                     rads=[_.rad for _ in [my_sys.atoms[_] for _ in [0, 1, 2, 23]]])
    assert [round(_, 3) for _ in loc] == [25.682, 31.695, 8.389]
    assert round(rad, 3) == 0.635


def test_verify_vor_site():
    my_sys = System("../Data/test_data/Na5.pdb")
    my_net = Network(my_sys, my_sys.atoms)
    my_net.sort_atoms()
    loc, rad, loc2, rad2 = calc_vert(locs=[_.loc for _ in [my_sys.atoms[_] for _ in [0, 1, 2, 23]]],
                                     rads=[_.rad for _ in [my_sys.atoms[_] for _ in [0, 1, 2, 23]]])
    assert verify_site(loc, rad, [0, 1, 2, 23], my_net)
