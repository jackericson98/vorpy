from vorpy.src.chemistry import special_radii


def set_pymol_atoms(sys):

    """
    Generates a PyMOL script to configure atomic radii for visualization.

    Explicit-radius sphere systems use each PDB B-factor as the sphere radius,
    covering every loaded sphere in one PyMOL command. Molecular PDBs retain
    their element and residue-specific radius rules.

    Args:
        sys: System object containing atomic information including:
            - element_radii: Dictionary mapping element symbols to their radii
            - residues: List of residue objects with atom information
            - balls: DataFrame containing atom names and radii
            - type: System type ('foam' or 'coarse' for special handling)

    Returns:
        None: Creates a PyMOL script file in the system's output directory
    """
    # Explicit-radius systems store each sphere's radius in the PDB B-factor.
    # Apply it to the whole loaded object so every sphere is covered, including
    # repeated atom/residue names and serial-number rollovers.
    if sys.type in {'foam', 'coarse', 'balls'}:
        with open('set_atoms.pml', 'w') as file:
            file.write(f"alter {sys.name}, vdw=b\n\nrebuild\n")
        return
    # Check to see if the atoms in the system are all accounted for
    for i, res in enumerate(sys.residues):
        if res.name not in special_radii:
            special_radii[res.name] = {sys.balls['name'][j]: round(sys.balls['rad'][j], 2) for j in res.atoms}
    # Create the file
    with open('set_atoms.pml', 'w') as file:
        # Write the change radii script for the system's set atomic radii
        for radius in sys.element_radii:
            if radius != '':
                file.write("alter {} and e. {}, vdw={}\n".format(sys.name, radius, sys.element_radii[radius]))
        # Change the radii for special atoms
        for res in special_radii:
            for atom in special_radii[res]:
                res_str = "r. {} ".format(res) if res != "" else ""
                file.write("alter {} and {}and n. {}, vdw={}\n".format(sys.name, res_str, atom, special_radii[res][atom]))
        # Rebuild the system
        file.write("\nrebuild")
