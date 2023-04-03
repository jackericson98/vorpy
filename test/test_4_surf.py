from test.test_2_vert import vor_net, pow_net, del_net
from test.test_1_load import pdb_sys_vor, pdb_sys_pow, pdb_sys_del


"""
______________________________________________Surface Values____________________________________________________________
"""


def test_vor_surf_num_points(vor_net):
    assert len(vor_net.surfs[0].points) == 48


def test_vor_surf_sa(vor_net):
    assert round(vor_net.surfs[0].sa, 3) == 0.684


def test_vor_surf_curve(vor_net):
    assert round(vor_net.surfs[0].curv, 3) == 21.367


def test_del_surf_num_points(del_net):
    assert len(del_net.surfs[0].points) == 45


def test_del_surf_sa(del_net):
    assert round(del_net.surfs[0].sa, 3) == 0.667


def test_pow_surf_num_points(pow_net):
    assert len(pow_net.surfs[0].points) == 67


def test_pow_surf_sa(pow_net):
    assert round(pow_net.surfs[0].sa, 3) == 0.875


"""
______________________________________________Network Connections_______________________________________________________
"""


def test_vor_num_surfs(vor_net):
    assert len(vor_net.surfs) == 18
