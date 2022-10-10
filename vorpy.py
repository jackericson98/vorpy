import multiprocessing as mp

from Visualize.vorpy_gui1 import Vorpy
from System.system import System, Network, Atom, verify_site
from System.calcs import *
import sys as Sys


# If __name__ == '__main__' wrapper for parallel cpu and gpu processing
if __name__ == '__main__':


    #############################################  Functions  ##########################################################

    #################################################  Find the vertices  ##############################################

    # Find vertices function. Used for parallel computing at the main wrapper level
    def find_vertices(net, gui=None, myCounter=None):
        net.find_verts(gui, myCounter)

    # Build surface function. Used for parallel computing at the main wrapper level
    def build_surf(mySurf):
        mySurf.build_surfs()

    # Analyze surface function. Used for analyzing each surface at the main wrapper level
    def analyze_surf(mySurf):
        mySurf.interface_sa = calc_sa(mySurf)

    # Analyze cell function. Used for analyzing each atomic cell at the main wrapper level
    def analyze_cell(myAtom):
        myAtom.cell_vol = calc_vol(myAtom)


    ##############################################  Settings GUI  ######################################################

    # Launch the GUI
    myVorpy = Vorpy()

    # Get the setting variables from the GUI
    sys_name = myVorpy.sys_name.get()
    sys_file = myVorpy.sys_file
    net_file = myVorpy.net_file
    surf_res = myVorpy.sys_res_flt.get()
    box_multiplier = myVorpy.sys_box_x_flt.get()
    beta_val = myVorpy.sys_alpha_value.get()
    export_dir = myVorpy.output_directory
    cpu_boost = myVorpy.parallelize.get()


    #############################################  Build the System  ##################################################

    # Check to see if the system has a file indicated. If any file is given the program builds that system
    if sys_file is not None:
        # Create the system
        mySys.base_file = sys_file
        mySys.box_size = box_multiplier
        mySys.min_dist = surf_res
        mySys.beta_val = beta_val
        mySys.name = sys_name
        # Check to see if a vertex file has been chosen
        if net_file is not None:
            # Add the vertices to the system
            mySys.load_net(net_file)


    # If the cpu boost button was clicked set the recursion limit based on the size of the system *needs to be tweaked
    if cpu_boost:
        recursion_limit = max(1000, (len(mySys.atoms) // 1000) * 1000)
        Sys.setrecursionlimit(recursion_limit)

    # Set the network object
    myNet = mySys.net


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
            percentage = float(np.round((i + 1) / len(surfs) * 100, 2))
            print("\rBuilding Surfaces: ",
                  '#' * (int(percentage) // 10) + ' ' * (10 - (int(percentage) // 10)), percentage, "%", end='')


    ##############################################  Analyze the network  ###############################################

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
            percentage = float(np.round((mySys.atoms.index(atom) + 1) / len(mySys.atoms) * 100, 2))
            print("\rAnalyzing Atomic Cells: ",
                  '#' * (int(percentage) // 10) + ' ' * (10 - (int(percentage) // 10)), percentage, "%", end='')

    ###################################################  Export the System  ############################################

    mySys.export(export_all=export_all, export_sys=export_sys, export_mols=export_mols,
                 export_atoms=export_atoms, export_surfs=export_surfs, export_reses=export_residues)
