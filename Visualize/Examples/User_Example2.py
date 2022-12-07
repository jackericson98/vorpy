from Visualize.mpl_visualize import *
import os
from System.system import System, Group
from System.output import *

os.chdir("../..")

"""
User Example. Follow the triple quotation comments in the code for instructions
"""


# Get group function. Interprets the strings from below
def get_group(selections, sys, name=None):
    print(selections)
    # Get the molecules from their names
    mol_atoms = []
    for mol in selections[0]:
        mol_atoms += sys.mols[sys.mol_names.index(mol)]
    res_atoms = []
    for res in selections[1]:
        res_atoms += sys.residues[sys.res_names.index(res)]
    my_atoms = [sys.atoms[int(atom)] for atom in selections[2]]
    ndx_atoms = []
    for ndx in selections[3]:
        ndx_atoms += sys.ndxs[sys.ndx_names.index(ndx)]

    return Group(net=sys.net, atoms=mol_atoms + res_atoms + my_atoms + ndx_atoms, name=name)


#################################### Create the System #################################################################

cube_atoms = [[[-1, 0, 0], 1], [[1, 0, 0], 1], [[0, 1, 0], 1], [[0, -1, 0], 1], [[0, 0, 1], 1], [[0, 0, -1], 1],
              [[0, 0, 0], 0.5]]

# Choose one of the following line for testing below (e.g. System(file=test_files[0]))  V V V
test_files = ["Na5", "Na7", "EDTA_Mg", "1BNA", "cambrin", "DB1976", "hairpin", "18L4_benzene"]
test_files = [os.getcwd() + "/Data/test_data/" + test_file + ".pdb" for test_file in test_files]

"""
1. Create the system and specify the different files. Either use the full file address or use a test file from above
    a. file: Main system file address. .pdb, .gro, .mol, .cif file types allowed. This gets priority over atoms
    b. atoms: List of either location and radii data or Atom objects (e.g. System(atoms=[[[1, 0, 0], 2], ... ])
    c. network_file: Network file address. Must be created with the same atoms as the current system  ...   VVV
    d. verts_file: Vertices file address ... 
    e. index_file: Indexes file address ...
    f. frames_file: List of frame file addresses ...
    g. output_directory: Directory address for desired output destination for network exports

"""
mySys = System(file=test_files[0],
               atoms=None,
               network_file=None,
               verts_file=None,
               index_file=None,
               frame_files=None,
               output_directory=None)

########################################## Build the Network ###########################################################


"""
 2. Build the network and specify the desired settings.
    a. surf_res: Resolution of surfaces created in network -|- From 0.01 to 1 A (Angstroms), recommended 0.1 A
    b. max_vert: Maximum allowed vertex radius for network -|- From 0.10 to 20 A, recommended 7 A 
    c. box_size: Allowed vertex retaining box size multiplier -|- From 1 to 10 A, recommended 1.5 A
"""
mySys.build_network(surf_res=0.2,
                    max_vert=10,
                    box_size=1.5,
                    sol_verts=False,
                    output=False)

# Print the network data
print("\r\n{} Network Built:\n\n    Time = {:.2f} seconds\n    Vertices Found = {}\n    Surfaces Built = {}\n"
      .format(mySys.name, mySys.net.my_time, len(mySys.net.verts), len(mySys.net.surfs)), end="")

######################################### Plot the network #############################################################

"""
3. Show (or dont show) the network objects and set the plot attributes
    a. net: Network object (use mySys.net)
    b. Show: Show the whole plot or not
    c. atoms: Show atoms in the plot
    d. verts: Show the vertices in the plot
    e. edges: Show the edges in the plot
    f. surfs: Show the surfaces in the plot
    g. bg_color: Set the background color for the plot
    h. grid: Show the grid for the plot

"""

plot_net(net=mySys.net,
         Show=True,
         atoms=True,
         verts=True,
         edges=True,
         surfs=True,
         bg_color='black',
         grid=False
         )

######################################## Create Groups for Analysis ####################################################


"""
4. Create Group 1:
    a. Using the names of the molecules, residues or atoms (e.g. mols = ["A", "C"], resids = ["DC 22"]) 
    b. From the index_file loaded above (e.g. g1_ndxs = ["PROTEIN"]). The index file can also be loaded after the fact 
       using mySys.load_ndx(file=path_str)
"""
g1_mols = []
g1_resids = []
g1_atoms = []
g1_ndxs = []
"""
4b. Create a Group 2 for interfacial analysis
"""
g2_mols = []
g2_resids = []
g2_atoms = []
g2_ndxs = []

# Get the groups
group1 = get_group([g1_mols, g1_resids, g1_atoms, g1_ndxs], mySys)
group2 = get_group([g2_mols, g2_resids, g2_atoms, g2_ndxs], mySys)

####################################### Export Selections ##############################################################


"""
5. Export the desired selections
    a. network: Export the reloadable network file created from a previously solved system
    b. pdb: Export a pdb file for the system (if a file was loaded it is a copy of that file)
    c. surfaces: Export all surfaces from the network
    d.no_sol_network_object: Export the outer surfaces for all of the 
"""

mySys.initial_export(network=True,
                     pdb=True,
                     surfaces=True,
                     full_network_object=True,
                     no_sol_network_object=True,
                     alter_atoms_script=True)

########################################## Open in pymol ###############################################################

"""
6. Check your output directory for the outputs (if you're not sure check vorpy/Data/User_data/) and open the outputs in
 a program that can handle molecule files and .off files like Pymol
 """


