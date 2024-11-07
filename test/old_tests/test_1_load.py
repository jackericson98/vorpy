import pytest

from System.system import System


"""
__________________________________________________Fixtures______________________________________________________________
"""


@pytest.fixture(scope='session')
def pdb_sys_vor():
    my_sys = System('../../Data/test_data/EDTA_Mg.pdb')
    my_sys.create_group(atoms=[my_sys.atoms[15]])
    my_sys.net.build(my_group=my_sys.groups[0], surf_res=0.5)
    return my_sys


@pytest.fixture(scope='session')
def pdb_sys_pow():
    my_sys = System('../../Data/test_data/EDTA_Mg.pdb')
    my_sys.create_group(atoms=[my_sys.atoms[0]])
    my_sys.net.build(my_group=my_sys.groups[0], net_type='pow')
    return my_sys


@pytest.fixture(scope='session')
def pdb_sys_del():
    my_sys = System('../../Data/test_data/EDTA_Mg.pdb')
    my_sys.create_group(atoms=[my_sys.atoms[0]])
    my_sys.net.build(my_group=my_sys.groups[0], net_type='del')
    return my_sys


@pytest.fixture(scope='module')
def gro_sys():
    my_sys = System('../Data/test_data/EDTA_Mg.gro')
    return my_sys


@pytest.fixture(scope='module')
def mol_sys():
    my_sys = System('../Data/test_data/EDTA_Mg.mol')
    return my_sys


@pytest.fixture(scope='module')
def cif_sys():
    my_sys = System('../Data/test_data/EDTA_Mg.cif')
    return my_sys


"""
____________________________________________________Atom Files__________________________________________________________
Test the functionality of read_pdb, read_gro, read_mol, read_cif
"""


# Load Pdb
def test_read_pdb(pdb_sys_vor):
    assert len(pdb_sys_vor.atoms) == 612
    assert pdb_sys_vor.atoms[30].loc == [20.39, 18.3, 8.23]

#
# # Load GRO
# def test_read_gro(gro_sys):
#     assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]
#
#
# # Load MOL
# def test_read_mol(mol_sys):
#     assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]
#
#
# # Load CIF
# def test_read_cif(cif_sys):
#     assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]


"""
___________________________________________________Network File_________________________________________________________
Load a network file for a large atom file and test the number of network elements, the values of random network elements
"""


# def test_load_net(pdb_sys):
#     pdb_sys.load_net(file="../Data/test_data/cambrin_net.csv")
#     assert len(pdb_sys.net.verts) == 4679

"""
____________________________________________________Index File__________________________________________________________
Load a GROMACS index file and test that the groupings are 181L
"""


# def test_read_ndx(pdb_sys):
#     pdb_sys.load_ndx(file="../Data/test_data/cambrin_ndx.ndx")
#     assert len(pdb_sys.groups[0].atoms) == 45


"""
_____________________________________________________Voronota Files_____________________________________________________
"""

# def test_read_vta_data(pdb_sys):
#     pdb_sys.load_verts(verts_file="../Data/test_data/cambrin_verts.txt",
#                        balls_file="../Data/test_data/cambrin_balls.txt")
#     assert pdb_sys.verts is not None


"""
____________________________________________Atom Sorting________________________________________________________________
Tests the atom sorting schemes for a loaded network
"""


def test_net_calc_box(pdb_sys_vor):
    assert pdb_sys_vor.net.box == [[4.754, 2.664, -5.988], [36.891, 38.341, 25.568]]
    assert pdb_sys_vor.net.ball_box == [[10.11, 8.61, -0.729], [31.535, 32.395, 20.309]]


def test_net_sort_atoms(pdb_sys_vor):
    assert pdb_sys_vor.net.atoms[45].box == [14, 8, 14]
    assert pdb_sys_vor.net.sub_box_size == [1.285, 1.427, 1.262]
    assert len(pdb_sys_vor.net.sub_boxes) == 25
