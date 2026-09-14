import numpy as np

from vorpy.workbench.domain import Atom
from vorpy.workbench.services.structure_loader import COVALENT_RADII, _infer_bonds


def test_inferred_bonds_match_pairwise_distances():
    rng = np.random.default_rng(42)
    elements = list(COVALENT_RADII) + ["UNKNOWN"]
    atoms = [Atom(i, i + 1, "X", elements[i % len(elements)], tuple(pos), "ILE", "1", "A")
             for i, pos in enumerate(rng.uniform(-5, 5, (200, 3)))]
    expected = set()
    for i, atom in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            other = atoms[j]
            distance = np.linalg.norm(np.asarray(atom.position) - other.position)
            cutoff = COVALENT_RADII.get(atom.element, 0.76) + COVALENT_RADII.get(other.element, 0.76) + 0.45
            if 0.4 <= distance <= cutoff:
                expected.add((i, j))
    assert _infer_bonds(atoms) == expected
    assert _infer_bonds([]) == set()
    assert _infer_bonds(atoms[:1]) == set()


def test_progressive_loading_preserves_file_order_and_preview_independence(tmp_path):
    from dataclasses import replace
    from vorpy.workbench.services.structure_loader import load_pdb
    from vorpy.workbench.chemistry import is_solvent
    path = tmp_path / 'interleaved.pdb'
    def record(serial, residue, element, x):
        return (f'ATOM  {serial:5d} {element:^4s} {residue:>3s} A   1    '
                f'{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2s}\n')
    path.write_text(record(1, 'HOH', 'O', 10.) + record(2, 'ILE', 'C', 0.)
                    + record(3, 'HOH', 'H', 11.) + record(4, 'ILE', 'N', 1.)
                    + 'CONECT    2    4\nEND\n')
    previews = []
    staged = load_pdb(path, preview=previews.append)
    ordinary = load_pdb(path)
    assert staged.atoms == ordinary.atoms
    assert staged.bonds == ordinary.bonds
    assert len(previews) == 1
    primary = previews[0]
    assert [atom.index for atom in primary.atoms] == [0, 1]
    assert [atom.serial for atom in primary.atoms] == [2, 4]
    assert primary.atoms == [replace(atom, index=i) for i, atom in enumerate(
        atom for atom in staged.atoms if not is_solvent(atom))]
    # The callback's list must remain unchanged as the full list is assembled.
    assert len(primary.atoms) == 2
    assert primary.bonds == [type(staged.bonds[0])(0, 1)]


def test_solvent_only_and_invalid_solvent_do_not_corrupt_indices(tmp_path):
    from vorpy.workbench.services.structure_loader import load_pdb
    path = tmp_path / 'water.pdb'
    water = 'ATOM      1  O   HOH A   1      10.000   0.000   0.000  1.00  0.00           O\n'
    path.write_text(water)
    previews = []
    result = load_pdb(path, preview=previews.append)
    assert previews == []
    assert len(result.atoms) == 1
    invalid = water[:30] + 'INVALID ' + water[38:]
    primary = water.replace('HOH', 'ILE').replace('10.000', '20.000')
    path.write_text(invalid + primary + water)
    result = load_pdb(path, preview=previews.append)
    assert result.atoms == load_pdb(path).atoms
    assert [atom.index for atom in result.atoms] == [0, 1]
    assert len(previews[-1].atoms) == 1


def test_ion_contacts_are_not_inferred_but_polyatomic_bonds_remain():
    atoms = [
        Atom(0, 1, "O", "O", (0., 0., 0.), "MOL", "1", "A"),
        Atom(1, 2, "MG", "MG", (2., 0., 0.), "MG", "2", "A"),
        Atom(2, 3, "S", "S", (10., 0., 0.), "SO4", "3", "A"),
        Atom(3, 4, "O", "O", (11.4, 0., 0.), "SO4", "3", "A"),
        Atom(4, 5, "O", "O", (10., 1.4, 0.), "HOH", "4", "A"),
    ]
    assert _infer_bonds(atoms) == {(2, 3)}


def test_bundled_ions_have_no_inferred_ligand_bonds():
    from pathlib import Path
    from vorpy.workbench.services.structure_loader import load_pdb
    from vorpy.workbench.chemistry import ION_RESIDUES
    for name in ("Na5", "EDTA"):
        result = load_pdb(Path(__file__).parents[2] / "data" / f"{name}.pdb")
        ions = {a.index for a in result.atoms if a.residue_name.upper() in ION_RESIDUES}
        assert ions
        assert result.bonds
        assert all(b.atom_a not in ions and b.atom_b not in ions for b in result.bonds)
