from System.sys_funcs.input import *
from System.system import System


"""
____________________________________________________Atom Files__________________________________________________________
Test the functionality of read_pdb, read_gro, read_mol, read_cif
"""


# Load Pdb
def test_read_pdb():
    sys = System('../Data/test_data/Na5.pdb')
    read_pdb(sys)
    assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]

#
# # Load GRO
# def test_read_gro():
#     sys = System('../Data/test_data/Na5.gro')
#     read_gro(sys)
#     assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]
#
#
# # Load MOL
# def test_read_mol():
#     sys = System('../Data/test_data/Na5.mol')
#     read_mol(sys)
#     assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]
#
#
# # Load CIF
# def test_read_cif():
#     sys = System('../Data/test_data/Na5.cif')
#     read_cif(sys)
#     assert len(sys.atoms) == 61 and sys.atoms[30].loc == [28.854, 34.009, 7.627]
#
#
# # Load System
# def test_read_full_sys():
#     sys = System('../Data/test_data/test_sys.pdb')
#     assert len(sys.atoms) == 20

"""
___________________________________________________Network File_________________________________________________________
Load a network file for a large atom file and test the number of network elements, the values of random network elements
"""


# def test_load_net():
#     my_sys = System("../Data/test_data/cambrin.pdb", network_file="../Data/test_data/cambrin_net.csv")
#     assert len(my_sys.net.verts) == 0

"""
____________________________________________________Index File__________________________________________________________
Load a GROMACS index file and test that the groupings are correct
"""


# def test_read_ndx():
#     my_sys = System("../Data/test_data/cambrin.pdb", index_file="../Data/test_data/cambrin_ndx.ndx")
#     assert len(my_sys.groups[0].atoms) == 45


"""
_____________________________________________________Voronota Files_____________________________________________________
"""

# def test_read_vta_data():
#     my_sys = System("../Data/test_data/cambrin.pdb", verts_file="../Data/test_data/cambrin_verts.txt",
#                     balls_file="../Data/test_data/cambrin_balls.txt")
#     assert my_sys.verts is not None
