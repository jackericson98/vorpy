import shutil
from os import path
from vorpy.src.chemistry import special_radii


def set_pymol_atoms(sys):

    """
    Generates a PyMOL script to configure atomic radii for visualization.

    This function creates a 'set_atoms.pml' script that sets the van der Waals radii for atoms
    in the PyMOL visualization. It handles both standard element radii and special cases for
    specific residues and atoms. The script can be used to ensure accurate representation of
    atomic sizes in PyMOL visualizations.

    Args:
        sys: System object containing atomic information including:
            - element_radii: Dictionary mapping element symbols to their radii
            - residues: List of residue objects with atom information
            - balls: DataFrame containing atom names and radii
            - type: System type ('foam' or 'coarse' for special handling)

    Returns:
        None: Creates a PyMOL script file in the system's output directory
    """
    # If we have special circumstances for the atoms in our base file, output the already created set pymol atoms
    if sys.type == 'foam' or sys.type == 'coarse':
        # Get the directory for the base_file and copy the set atoms file
        try:
            shutil.copyfile(path.dirname(sys.files['base_file']) + '/set_atoms.pml', sys.files['dir'] + '/sys/set_atoms.pml')
        except FileNotFoundError:
            # Create the file
            with open('set_atoms.pml', 'w') as file:
                for i, ball in sys.balls.iterrows():
                    file.write(
                        "alter r. {} and n. {}, vdw={}\n".format(ball['res_name'], ball['name'], ball['rad']))
                file.write("\nrebuild")
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