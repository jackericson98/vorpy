import multiprocessing as mp
import os

from Visualize.vorpy_gui import Vorpy, ErrorBox, LoadingBox
from System.system import System
from System.calcs import *
import sys as Sys


# If __name__ == '__main__' wrapper for parallel cpu and gpu processing
if __name__ == '__main__':


    #############################################  Functions  ##########################################################

    # Build surface function. Used for parallel computing at the main wrapper level
    def build_surf(mySurf):
        mySurf.build()

    # Analyze surface function. Used for analyzing each surface at the main wrapper level
    def analyze_surf(mySurf):
        mySurf.sa = calc_sa(mySurf)

    # Analyze cell function. Used for analyzing each atomic cell at the main wrapper level
    def analyze_cell(myAtom):
        myAtom.vol = calc_vol(myAtom)


    ##############################################  Settings GUI  ######################################################

    # Launch the GUI
    myVorpy = Vorpy()

    # Get the setting variables from the GUI
    sys_name = myVorpy.sys_file_name.get()
    sys_file = myVorpy.sys_file_address
    vert_file = myVorpy.vert_file_address
    user_atoms = myVorpy.sys_atom_list
    surf_res = myVorpy.sys_res_flt.get()
    box_multiplier = myVorpy.sys_box_x_flt.get()
    export_dir = myVorpy.sys_output_directory
    cpu_boost = myVorpy.parallelize.get()
    get_analysis = myVorpy.output_analysis.get()
    export_all = myVorpy.output_all.get()
    export_surfs = myVorpy.output_surfs.get()
    export_atoms = myVorpy.output_atoms.get()
    export_mols = myVorpy.output_mols.get()
    export_analysis = myVorpy.output_analysis.get()
    export_sys = myVorpy.output_sys.get()


    #############################################  Create the System  ##################################################

    # Check to see if the system has a file indicated. If any file is given the program builds that system
    if sys_file is not None:
        # Create the system
        mySys = System(sys_file, box_size=box_multiplier, min_dist=surf_res)
        mySys.name = sys_name
        # Check to see if a vertex file has been chosen
        if vert_file is not None:
            # Add the vertices to the system
            mySys.add_verts(vert_file)

    # If the user entered atoms build a list system
    elif len(user_atoms) > 0:
        mySys = System(user_atoms, box_size=box_multiplier, min_dist=surf_res)
        mySys.name = "User_Data"

    # Throw an error otherwise
    else:
        ErrorBox("Please select a file or enter atoms")
        mySys = None

    # If the cpu boost button was clicked set the recursion limit based on the size of the system *needs to be tweaked
    if cpu_boost:
        recursion_limit = max(1000, (len(mySys.atoms) // 1000) * 1000)
        Sys.setrecursionlimit(recursion_limit)


    #################################################  Find the vertices  ##############################################

    # Catch for if the verts have been loaded already.
    if vert_file is None:
        mySys.net.find_verts()

    # Export the vertices
    mySys.export_verts()


    ################################################  Connect network  #################################################

    # Connect the system
    mySys.net.connect()


    ################################################  Build the surfaces  ##############################################

    # If the cpu boost button was clicked run parallel processors
    if cpu_boost:
        # Set up the pool running with statement
        with mp.Pool(12) as calc_surf_pool:
            # Set the pool and get the progress
            for surf in calc_surf_pool.imap_unordered(build_surf, mySys.net.surfs):
                # Find the current percentage through the network and print it
                percentage = int((mySys.net.surfs.index(surf) + 1) / len(mySys.net.surfs) * 100)
                print("\rBuilding Surfaces: ",
                      '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')

    # If the cpu boost button wasn't clicked run through each surface in a for loop
    else:
        surfs = mySys.net.surfs
        for i in range(len(surfs)):
            # Calculate and print the running percentage for mesh calculations
            build_surf(surfs[i])
            # Calculate and print the running percentage for mesh calculations
            percentage = int((i + 1) / len(surfs) * 100)
            print("\rBuilding Surfaces: ",
                  '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')


    ##############################################  Analyze the network  ###############################################

    # Check to see if the user wants the network analyzed
    if get_analysis or export_all:

        # If the cpu boost button was clicked use the 'analyze surface' and 'analyze cell' functions
        if cpu_boost:
            # Analyze the surfaces in the system
            with mp.Pool() as anal_surf_pool:
                for surf in anal_surf_pool.imap(analyze_surf, mySys.net.surfs):
                    # Print the running percentage meter for analyzing the surfaces
                    percentage = int((mySys.net.surfs.index(surf) + 1) / len(mySys.net.surfs) * 100)
                    print("\rAnalyzing Surfaces: ",
                          '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
            # Run analysis on the cells of the system
            with mp.Pool() as anal_cell_pool:
                for surf in anal_cell_pool.imap(analyze_cell, mySys.net.surfs):
                    percentage = int((mySys.net.surfs.index(surf) + 1) / len(mySys.net.surfs) * 100)
                    print("\rAnalyzing Atomic Cells: ",
                          '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')

        # Without the cpu_boost button clicked loop through each of the surfaces and atoms and analyze them
        else:
            # Go through each of the surfaces in the network
            for surf in mySys.net.surfs:
                # Analyze the surface
                analyze_surf(surf)
                # Update the print statement
                percentage = int((mySys.net.surfs.index(surf) + 1) / len(mySys.net.surfs) * 100)
                print("\rAnalyzing Surfaces: ",
                      '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')
            # Go through each of the atoms in the network
            for atom in mySys.atoms:
                # Analyze each of the atom's cells
                calc_vol(atom)
                # Update the print statement
                percentage = int((mySys.net.surfs.index(surf) + 1) / len(mySys.net.surfs) * 100)
                print("\rAnalyzing Atomic Cells: ",
                      '#' * (percentage // 10) + ' ' * (10 - (percentage // 10)), percentage, "%", end='')

    ###################################################  Export the System  ############################################
    # Export the system
    mySys.export(export_all=export_all, export_sys=export_sys, export_mols=export_mols,
                 export_atoms=export_atoms, export_analysis=export_analysis, export_surfs=export_surfs)
