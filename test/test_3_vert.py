from System.system import System
from System.Network.network import Network
from System.Group.group import Group
from System.Network.net_funcs.find_verts import *

"""
___________________________________________________Voronoi Calculations_________________________________________________
Tests: calc_vert, verify_site, find_verts
"""


def test_calc_vor_vert():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_atoms = [my_sys.atoms[_] for _ in [15, 255, 256, 332]]
    loc, rad, loc2, rad2 = calc_vert(locs=[_.loc for _ in my_atoms], rads=[_.rad for _ in my_atoms])
    assert [round(_, 3) for _ in loc] == [19.939, 18.897, 12.831]
    assert round(rad, 3) == 0.116


def test_verify_vor_site():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_atoms = [my_sys.atoms[_] for _ in [15, 255, 256, 332]]
    loc, rad, loc2, rad2 = calc_vert(locs=[_.loc for _ in my_atoms], rads=[_.rad for _ in my_atoms])
    assert verify_site(loc, rad, [15, 255, 256, 332], my_sys.net)


def test_find_vor_verts():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_group = Group(sys=my_sys, atoms=[my_sys.atoms[15]])
    my_sys.net.build(output=False, build_surfs=False, print_actions=False, my_group=my_group)
    assert len(my_sys.net.verts) == 32


"""
________________________________________________________Delaunay Calculations___________________________________________
"""


def test_calc_del_vert():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb", )
    my_sys.net.type = 'del'
    loc, rad = calc_flat_vert(locs=[_.loc for _ in [my_sys.atoms[_] for _ in [0, 1, 3, 4]]],
                                          rads=[_.rad for _ in [my_sys.atoms[_] for _ in [0, 1, 3, 4]]], power=False)
    assert [round(_, 3) for _ in loc] == [21.353, 21.368, 9.513]
    assert round(rad, 3) == 1.449


def test_verify_del_site():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_sys.net.type = 'del'
    loc, rad = calc_flat_vert(locs=[_.loc for _ in [my_sys.atoms[_] for _ in [0, 1, 3, 4]]],
                                     rads=[_.rad for _ in [my_sys.atoms[_] for _ in [0, 1, 3, 4]]], power=False)
    assert verify_site(loc, rad, [0, 1, 3, 4], my_sys.net)


def test_find_del_verts():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_net = Network(my_sys, my_sys.atoms, net_type='del')
    my_group = Group(sys=my_sys, atoms=[my_sys.atoms[0]])
    my_net.build(output=False, build_surfs=False, print_actions=False, my_group=my_group)
    assert len(my_net.verts) == 32


"""
__________________________________________________________Power Calculations____________________________________________
"""


def test_calc_pow_vert():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_sys.net.type = 'pow'
    my_atoms = [my_sys.atoms[i] for i in [0, 1, 3, 4]]
    loc, rad = calc_flat_vert(locs=[_.loc for _ in my_atoms], rads=[_.rad for _ in my_atoms], power=True)
    assert [round(_, 3) for _ in loc] == [21.233, 21.301, 9.761]
    assert round(rad, 3) == -0.372


def test_verify_pow_site():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_sys.net.type = 'pow'
    my_atoms = [my_sys.atoms[i] for i in [0, 1, 3, 4]]
    loc, rad = calc_flat_vert(locs=[_.loc for _ in my_atoms], rads=[_.rad for _ in my_atoms], power=True)
    assert verify_site(loc, rad, [0, 1, 3, 4], my_sys.net)


def test_find_pow_verts():
    my_sys = System("../Data/test_data/EDTA_Mg.pdb")
    my_sys.net.type = 'pow'
    my_group = Group(sys=my_sys, atoms=[my_sys.atoms[0]])
    my_sys.net.build(output=False, build_surfs=False, print_actions=False, my_group=my_group)
    assert len(my_sys.net.verts) == 20

