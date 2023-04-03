from test.test_1_load import pdb_sys_vor


def test_group_atoms(pdb_sys_vor):
    assert len(pdb_sys_vor.groups[0].atoms) == 1


def test_group_vertices(pdb_sys_vor):
    assert len(pdb_sys_vor.groups[0].verts) == 32


def test_group_edges(pdb_sys_vor):
    assert len(pdb_sys_vor.groups[0].edges) == 49


def test_group_surfs(pdb_sys_vor):
    assert len(pdb_sys_vor.groups[0].surfs) == 17
