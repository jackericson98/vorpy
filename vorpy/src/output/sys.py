import os
from vorpy.src.output.system_info import export_sys_info
from vorpy.src.output.pdb import write_pdb
from vorpy.src.output.set_pymol_atoms import set_pymol_atoms


def export_sys(sys, all_=False, pdb=False, set_atoms=False, info=False, mol=False, cif=False, xyz=False,
               txt=False, print_output=False):
    """
    Manages system output operations and directory preparation.

    This function serves as the main entry point for exporting system data and ensuring proper output directory structure.
    It handles various output formats and configurations including:
    - PDB file generation
    - System information export
    - Network object export
    - PyMOL atom script generation
    
    The function ensures the output directory exists and is properly configured before performing any export operations.
    It maintains consistency across different output types and handles directory navigation appropriately.

    :param sys: System object containing the data to be exported
    :param all_: Boolean flag to enable all export options
    :param pdb: Boolean flag to enable PDB file export
    :param set_atoms: Boolean flag to enable PyMOL atom script generation
    :param info: Boolean flag to enable system information export
    :param mol: Boolean flag to enable MOL file export
    :param cif: Boolean flag to enable CIF file export
    :param xyz: Boolean flag to enable XYZ file export
    :param txt: Boolean flag to enable text file export
    :return: None
    """
    if sys.files['dir'] is None:
        sys.set_output_directory()
    # If the information is requested, export it
    if info or all_:
        os.chdir(sys.files['dir'])
        export_sys_info(sys)
    if pdb or all_:
        os.chdir(sys.files['dir'])
        # Export a pdb file for the system
        write_pdb([_ for i, _ in sys.balls.iterrows()], sys.name, sys)
        os.chdir(sys.files['dir'])
    # Write the alter atoms script
    if set_atoms or all_:
        os.chdir(sys.files['dir'])
        set_pymol_atoms(sys)
    # Print the output directory
    if print_output:
        print("\rOutput directory: {}".format(sys.files['dir']), end="")
