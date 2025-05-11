import numpy as np
from vorpy.src.calculations import calc_dist
from vorpy.src.objects import Residue


def fix_sol(sys, residue):
    """
    Reorganizes water molecules within a residue by correctly grouping oxygen and hydrogen atoms.

    This function processes a residue containing water molecules and ensures proper organization
    of atoms into complete water molecules (H2O). It handles various edge cases including:
    - Missing hydrogen atoms
    - Extra hydrogen atoms
    - Incorrect initial grouping
    - Atoms from different water molecules being mixed together

    Parameters:
    -----------
    residue : Residue
        The residue object containing water molecule atoms to be reorganized.
        Expected to contain oxygen and hydrogen atoms that need to be properly grouped.

    Returns:
    --------
    list
        A list of Residue objects, each representing a complete water molecule with:
        - One oxygen atom
        - Two hydrogen atoms
        - Proper residue assignments
        - Correct chain and sequence information

    Notes:
    ------
    - Uses a distance-based approach to assign hydrogens to nearest oxygen
    - Maximum O-H bond length is assumed to be 1.5 Å
    - Incomplete water molecules are attempted to be fixed by reassigning hydrogens
    """

    # Initialize containers for oxygen and hydrogen atoms
    oxy_res = []
    hydrogens = []

    # Separate oxygen and hydrogen atoms into different lists
    for a in residue.atoms:
        # Get the atom
        atom = sys.balls.iloc[a]
        # If the atom is an oxygen atom
        if atom['element'].lower() == 'o':
            # Create a new residue for each oxygen atom
            oxy_res.append(Residue(sys=residue.sys, atoms=[a], name=atom['residue'],
                                   sequence=atom['res_seq'], chain=atom['chn']))
        # If the atom is a hydrogen atom
        elif atom['element'].lower() == 'h':
            # Add the atom to the hydrogens list
            hydrogens.append(atom['num'])

    # Assign hydrogens to the nearest oxygen atom to form water molecules
    for h in hydrogens:
        # Initialize the closest residue and the minimum distance
        closest_res, min_dist = None, np.inf
        # Go through the oxygen residues
        for res in oxy_res:
            # Calculate the distance between the oxygen and the hydrogen
            dist = calc_dist(sys.balls['loc'][res.atoms[0]], sys.balls['loc'][h])
            # If the distance is less than the minimum distance
            if dist < min_dist:
                # Update the minimum distance and the closest residue
                min_dist = dist
                closest_res = res
        # If the closest residue is not None and the minimum distance is less than 1.5 Å
        if closest_res and min_dist < 1.5:
            # Add the hydrogen to the closest residue
            closest_res.atoms.append(h)
            # Remove the hydrogen from the hydrogens list
            hydrogens.remove(h)

    # Check the integrity of newly formed residues
    good_resids = []
    incomplete_resids = []
    # Go through the oxygen residues
    for res in oxy_res:
        # If the residue has 3 atoms
        if len(res.atoms) == 3:
            # Add the residue to the good residues
            good_resids.append(res)
            # Go through the atoms in the residue
            for a in res.atoms:
                # Set the residue for the atom
                sys.balls.loc[a, 'res'] = res
        else:
            incomplete_resids.append(res)

    # Attempt to correct incomplete residues
    for res in incomplete_resids:
        # If the residue has less than 3 atoms
        if len(res.atoms) < 3:
            # This block tries to find hydrogens that can be moved to this residue
            for h in hydrogens:
                dist = calc_dist(sys.balls['loc'][res.atoms[0]], sys.balls['loc'][h])
                if dist < 1.5:  # Assumed maximum bond length for O-H
                    res.atoms.append(h)
                    hydrogens.remove(h)
                if len(res.atoms) == 3:
                    break
            # print([(sys.balls['name'][_], sys.balls['res_seq'][_], sys.balls['loc'][_][0]) for _ in res.atoms])
        # Add the residue to the good residues
        good_resids.append(res)

    # Las sort the hydrogens
    if len(hydrogens) == 1:
        h = hydrogens[0]
        good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=sys.balls['name'][h],
                                   sequence=sys.balls['res_seq'][h], chain=sys.balls['chn'][h]))
    elif len(hydrogens) == 2:
        h1, h2 = sys.balls.iloc[hydrogens[0]], sys.balls.iloc[hydrogens[1]]
        if calc_dist(h1['loc'], h2['loc']) < 2 and h1['name'] != h2['name']:
            good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=h1['residue'],
                                       sequence=h1['res_seq'], chain=h1['chn']))
        else:
            for h in hydrogens:
                good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=sys.balls['name'][h],
                                           sequence=sys.balls['res_seq'][h], chain=sys.balls['chn'][h]))
    else:
        for h in hydrogens:
            good_resids.append(Residue(sys=residue.sys, atoms=hydrogens, name=sys.balls['name'][h],
                                       sequence=sys.balls['res_seq'][h], chain=sys.balls['chn'][h]))
    # print([(sys.balls['name'][_], sys.balls['res_seq'][_], sys.balls['loc'][_][0]) for _ in hydrogens])

    return good_resids