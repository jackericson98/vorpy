from System.sys_objs.residue import Residue
from System.sys_objs.chain import Chain, Sol


def sort_atoms(sys):
    # Go through the atoms in the system
    resids, chains, sys.residues, sys.chains = {}, {}, [], []
    # Go through the atoms in the system
    for atom in sys.atoms:
        # Create the chain and residue dictionaries
        res_name, chn_name = atom.residue + atom.res_seq, atom.chain
        # If the chain has been made before
        if chn_name in chains:
            # Get the chain from the dictionary and add the atom
            my_chn = chains[chn_name]
            my_chn.add_atom(atom)
            atom.chn = my_chn
        # Create the chain
        else:
            # If the chain is the sol chain
            if atom.residue.lower() == 'sol':
                my_chn = Sol(atoms=[atom], residues=[], name=atom.chain)
                sys.sol = my_chn
            # If the chain is not sol create a regular chain object
            else:
                my_chn = Chain(atoms=[atom], residues=[], name=atom.chain)
                sys.chains.append(my_chn)
            # Set the chain in the dictionary and give the atom it's chain
            chains[chn_name] = my_chn
            atom.chn = my_chn
        # Assign the atoms and create the residues
        if res_name in resids:
            my_res = resids[res_name]
            my_res.atoms.append(atom)
            atom.res = my_res
        else:
            my_res = Residue(atoms=[atom], name=atom.residue, sequence=atom.res_seq, chain=atom.chn)
            atom.chn.residues.append(my_res)
            resids[res_name] = my_res
            atom.res = my_res
            if atom.residue.lower() != 'sol':
                sys.residues.append(my_res)
            else:
                sys.sol.residues.append(my_res)
