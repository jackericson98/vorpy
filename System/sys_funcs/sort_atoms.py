from System.sys_objs.molcule import Molecule
from System.sys_objs.residue import Residue


def sort_atoms(sys):
    # Set up the chain names list
    sys.mols, sys.mol_names, sys.atom_names, = [], [], []
    # Go through each of the atoms in the system adding the atoms to their respective chains
    for atom in sys.atoms:

        # Set the atom's name
        sys.atom_names.append(atom.element + str(sys.atoms.index(atom)))
        # Add the solution
        if atom.mol_class.lower() == 'sol':
            if sys.sol is None:
                sys.sol = Molecule(atoms=[atom])
                atom.mol = sys.sol
            else:
                sys.sol.atoms.append(atom)
                atom.mol = sys.sol
        else:
            # If no chain is specified, set the chain to 'None'
            mol_name = atom.mol_class + atom.chain
            # If the atom's chain does not exist add it to the list of chains
            if mol_name not in sys.mol_names:
                my_mol = Molecule(atoms=[atom], name=mol_name)
                sys.mols.append(my_mol)
                sys.mol_names.append(mol_name)
                atom.mol = my_mol
            else:
                my_mol = sys.mols[sys.mol_names.index(mol_name)]
                my_mol.atoms.append(atom)
                atom.mol = my_mol
    # Add the solution to the molecules list
    if sys.sol is not None:
        sys.mols.append(Molecule(atoms=sys.sol))
        sys.mol_names.append("SOL")
    # Set up the residues names list
    sys.residues, sys.res_names = [], []
    # Set up the residues
    for atom in sys.atoms:
        res_name = atom.mol_class + atom.res_seq
        # If the residue name does not exist, add it
        if res_name not in sys.res_names:
            my_res = Residue(atoms=[atom], sequence=atom.res_seq, seg_id=atom.seg_id, mol=atom.mol, name=res_name)
            sys.residues.append(my_res)
            sys.res_names.append(res_name)
            atom.res = my_res
            atom.mol.resids.append(my_res)
        else:
            my_res = sys.residues[sys.res_names.index(res_name)]
            my_res.atoms.append(atom)
            atom.res = my_res
