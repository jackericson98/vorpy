import pytest

from System.Group.group import Group
from System.Network.net_funcs.find_verts import verify_site
from test.test_1_load import pdb_sys_vor, pdb_sys_pow, pdb_sys_del


"""
___________________________________________________Voronoi Calculations_________________________________________________
Tests: calc_vert, verify_site, find_verts
"""


@pytest.fixture(scope='session')
def vor_net(pdb_sys_vor):
    return pdb_sys_vor.net


# Test the calculate voronoi vertex functionality
def test_calc_vor_vert(vor_net):
    assert round(vor_net.verts[0].rad, 3) == 0.186
    assert [round(_, 3) for _ in vor_net.verts[0].loc] == [20.114, 20.532, 10.307]


# Test the doublet functionality
def test_doublet(vor_net):
    assert len([0 for _ in vor_net.verts if _.doublet is not None]) == 2


# Test the ability to verify a vertex site
def test_verify_vor_site(vor_net):
    assert verify_site([18.695, 19.126, 12.54], 0.135, [15, 255, 256, 332], vor_net, vor_net.type)


def test_find_vor_verts(vor_net):
    assert len(vor_net.verts) == 32


"""
________________________________________________________Delaunay Calculations___________________________________________
"""


@pytest.fixture(scope='session')
def del_net(pdb_sys_del):
    return pdb_sys_del.net


def test_calc_del_vert(del_net):
    assert [round(_, 3) for _ in del_net.verts[0].loc] == [21.353, 21.368, 9.513]
    assert round(del_net.verts[0].rad, 3) == 1.449


def test_verify_del_site(del_net):
    assert verify_site([21.353, 21.368, 9.513], 1.449, [0, 1, 3, 4], del_net, del_net.type)


def test_find_del_verts(del_net):
    assert len(del_net.verts) == 32


"""
__________________________________________________________Power Calculations____________________________________________
"""


@pytest.fixture(scope='session')
def pow_net(pdb_sys_pow):
    return pdb_sys_pow.net


def test_calc_pow_vert(pow_net):
    assert [round(_, 3) for _ in pow_net.verts[0].loc] == [21.233, 21.301, 9.761] and round(pow_net.verts[0].rad, 3) == -0.372


def test_verify_pow_site(pow_net):
    assert verify_site([21.233, 21.301, 9.761], -0.372, [0, 1, 3, 4], pow_net, pow_net.type)


def test_find_pow_verts(pow_net):
    assert len(pow_net.verts) == 20

