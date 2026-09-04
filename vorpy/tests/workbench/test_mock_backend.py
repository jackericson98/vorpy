from vorpy.workbench.services.mock_backend import MockBackend


def test_mock_backend_returns_displayable_result():
    updates = []
    result = MockBackend().solve(None, lambda label, value: updates.append((label, value)), lambda: False)

    assert result.atoms
    assert result.bonds
    assert result.layers
    assert updates[-1][1] == 100


def test_pdb_loader_reads_atoms_metadata_and_conect(tmp_path):
    from vorpy.workbench.services.structure_loader import load_pdb

    pdb = tmp_path / "tiny.pdb"
    pdb.write_text(
        "ATOM      1  N   GLY A   7       0.000   0.000   0.000  1.00  0.00           N  \n"
        "ATOM      2  CA  GLY A   7       1.450   0.000   0.000  1.00  0.00           C  \n"
        "CONECT    1    2\nEND\n",
        encoding="utf-8",
    )
    result = load_pdb(pdb)
    assert len(result.atoms) == 2
    assert len(result.bonds) == 1
    assert result.atoms[0].residue_name == "GLY"
    assert result.atoms[0].residue_sequence == "7"
    assert result.atoms[0].chain == "A"
