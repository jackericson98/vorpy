import multiprocessing as mp
import os

from Visualize.vorpy_gui import Vorpy, ErrorBox, LoadingBox
from System.system import System, Network, Atom
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
        myAtom.cell_vol = calc_vol(myAtom)


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
    beta_val = myVorpy.sys_beta_value.get()
    export_dir = myVorpy.sys_output_directory
    cpu_boost = myVorpy.parallelize.get()
    export_all = myVorpy.output_all.get()
    export_surfs = myVorpy.output_surfs.get()
    export_atoms = myVorpy.output_atoms.get()
    export_mols = myVorpy.output_mols.get()
    export_sys = myVorpy.output_sys.get()
    export_residues = myVorpy.output_residues.get()


    #############################################  Create the System  ##################################################

    # Check to see if the system has a file indicated. If any file is given the program builds that system
    if sys_file is not None:
        # Create the system
        mySys = System(sys_file, box_size=box_multiplier, min_dist=surf_res, beta_val=beta_val)
        mySys.name = sys_name
        # Check to see if a vertex file has been chosen
        if vert_file is not None:
            # Add the vertices to the system
            mySys.add_verts(vert_file)

    # If the user entered atoms build a list system
    elif len(user_atoms) > 0:
        mySys = System(user_atoms=user_atoms, box_size=box_multiplier, min_dist=surf_res, beta_val=beta_val)
        mySys.name = "User_Data"

    # Throw an error otherwise
    else:
        ErrorBox("Please select a file or enter atoms")
        mySys = None

    # If the cpu boost button was clicked set the recursion limit based on the size of the system *needs to be tweaked
    if cpu_boost:
        recursion_limit = max(1000, (len(mySys.atoms) // 1000) * 1000)
        Sys.setrecursionlimit(recursion_limit)

    # Set the network object
    net = mySys.net


    #################################################  Find the vertices  ##############################################

    # Catch for if the verts have been loaded already.
    if vert_file is None:
        # For small systems (<= 200) run the normal algorithm
        if len(mySys.atoms) <= 2000:
            net.find_verts()
        # For large systems, split the atoms into separate smaller networks top search for vertices in
        else:
            # This gets us to average about 60 atoms per medium box
            n = max(int(np.cbrt(len(mySys.atoms) // 60)), 2)
            # The range to search for new vertices
            rnge = len(net.sub_boxes) // n
            # Get the atoms in the
            med_boxes = [[[[] for i in range(n)] for j in range(n)] for k in range(n)]
            # Set up the network list and the list of the atom indices from the main network for reference later
            test_nets = []
            test_nets_real_ndxs = []
            # Go through each of the boxes in the medium_boxes matrix to create networks
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        inc = 5
                        # Make sure there are enough atoms to do a valid search
                        while len(med_boxes[i][j][k]) < 100:
                            # Starting from the middle of the range search out the atoms within 0.75 + the increment
                            med_boxes[i][j][k] = net.get_atoms([[int((i + .5) * rnge), int((j + .5) * rnge),
                                                                 int((k + .5) * rnge)]], int(0.5 * rnge) + inc)
                            inc += 1
                        # Set up the indices' tracker list
                        new_atoms, old_ndxs = [], []
                        # Create a copy of each atom in the network
                        for atom in med_boxes[i][j][k]:
                            new_atoms.append(Atom(atom.loc, atom.rad))
                            old_ndxs.append(mySys.atoms.index(atom))
                        # Create the network
                        test_nets.append(Network(mySys, new_atoms, box_size=1.1, min_dist=surf_res, beta_val=beta_val))
                        test_nets_real_ndxs.append(old_ndxs)
                        rnge = len(net.sub_boxes) // n
            # Build the networks
            for i in range(len(test_nets)):
                num_verts = len(net.verts)
                net2 = test_nets[i]
                net2.find_verts([num_verts, num_verts +
                                 sum([len(test_nets[_].atoms) * 6 for _ in range(i, len(test_nets))])])
                # Sort through the vertices in each network
                for j in range(len(net2.verts)):
                    vert = net2.verts[j]
                    # Get the vertex's actual atoms
                    vert.atoms = [mySys.atoms[test_nets_real_ndxs[i][ndx]] for ndx in vert.ndx]
                    # Get the vertex's index
                    vert.ndx = [mySys.atoms.index(atom) for atom in vert.atoms]
                    vert.ndx.sort()
                    # Look for the vertex in the network
                    v_ndx = search_verts(net.vert_ndxs, vert.ndx)
                    # If found, verify it is not a doublet
                    if 0 == len(net.vert_ndxs) or len(net.vert_ndxs) <= v_ndx:
                        net.vert_ndxs.append(vert.ndx)
                        net.verts.append(vert)
                    elif vert.ndx == net.vert_ndxs[v_ndx]:
                        # If the location is different and we haven't indicated the vertex as a doublet already add it
                        if [round(vert.loc[k], 7) for k in range(3)] != \
                                [round(net.verts[v_ndx].loc[k], 7) for k in range(3)] \
                                and not net.verts[v_ndx].doublet:
                            # Add the doublet
                            net.verts[v_ndx].doublet = True
                            mySys.net.verts[v_ndx].loc2, mySys.net.verts[v_ndx].rad2 = vert.loc, vert.rad
                    else:
                        # Add the vertex to the system
                        mySys.net.vert_ndxs.insert(v_ndx, vert.ndx)
                        mySys.net.verts.insert(v_ndx, vert)

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
