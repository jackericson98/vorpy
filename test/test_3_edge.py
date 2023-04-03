from test.test_2_vert import vor_net, pow_net, del_net
from test.test_1_load import pdb_sys_vor, pdb_sys_pow, pdb_sys_del


"""
This holds edge building tests and network connections
"""

"""
________________________________________Edge Building___________________________________________________________________
Test for issues in 
"""


def test_vor_edge_length(vor_net):
    assert len(vor_net.edges[0].points) == 7


def test_vor_edge_loc(vor_net):
    assert [round(_, 3) for _ in vor_net.edges[0].loc] == [20.007, 20.473, 10.206]


def test_pow_edge_length(pow_net):
    assert len(pow_net.edges[0].points) == 6


def test_del_edge_length(del_net):
    assert len(del_net.edges[0].points) == 5


"""
_______________________________________Network Connections______________________________________________________________
"""


def test_num_vor_edges(vor_net):
    assert len(vor_net.edges) == 49
