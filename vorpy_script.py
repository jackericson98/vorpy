from Visualize.mpl_visualize import *
from System.system import System, Group
from System.sys_funcs.output import *


"""

User Example:
 
 To use this script follow the triple quotation comments in the code for instructions on set up. 
 
 Once everything is set up run this script. 
 
 Note: The variables work in their current state, so only change the ones needed


"""


############################################## Example system bank #####################################################


# Cube configuration of atoms
cube_atoms = [[[-1, 0, 0], 1], [[1, 0, 0], 1], [[0, 1, 0], 1], [[0, -1, 0], 1], [[0, 0, 1], 1], [[0, 0, -1], 1],
              [[0, 0, 0], 0.5]]

# Choose one of the following line for testing below (e.g. System(file=test_files[0]))  V V V
test_files = ["Na5", "Na7", "EDTA_Mg", "1BNA", "cambrin", "DB1976", "hairpin", "18L4_benzene"]
test_files = [os.getcwd() + "/Data/test_data/" + test_file + ".pdb" for test_file in test_files]


############################################## Create the System #######################################################


"""
1. Create the system and specify the different files. Use the full file address (or select test_files[i] for i in 0-7)
    a. file: Main system file address. .pdb, .gro, .mol, .cif file types allowed. This gets priority over atoms
    b. atoms: List of either location and radii data or Atom objects (e.g. System(atoms=[[[1, 0, 0], 2], ... ])
    c. network_file: Network file address. Must be created with the same atoms as the current system  ...   VVV
    d. verts_file: Vertices file address ... 
    e. index_file: Indexes file address ...
    f. frames_file: List of frame file addresses ...
    g. output_directory: Directory address for desired output destination for network exports

"""
mySys = System(file=test_files[7],
               atoms=None,
               network_file=None,
               verts_file="C:/Users/jacke/PycharmProjects/vorpy/Data/User_data/18L4_benzene2/18L4_benzene_verts.txt",
               index_file=None,
               frame_files=None,
               output_directory=None)


############################################### Build the Network ######################################################


"""
 2. Build the network and specify the desired settings.
    a. surf_res: Resolution of surfaces created in network -|- From 0.01 to 1 A (Angstroms), recommended 0.1 A
    b. max_vert: Maximum allowed vertex radius for network -|- From 0.10 to 20 A, recommended 7 A 
    c. box_size: Allowed vertex retaining box size multiplier -|- From 1 to 10 A, recommended 1.5 A
    d. sol_verts: Calculate the vertices between the solute atoms or just the molecule atoms
    e. output: 
"""
mySys.build_network(surf_res=0.1,
                    max_vert=7,
                    box_size=1.5,
                    sol_verts=True,
                    output=False)


############################################## Plot the Network ########################################################


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


################################################### Create Groups ######################################################


"""
4. Create Group 1:
    a. net: Network object (Use mySys.
"""
group1 = Group(net=mySys.net,
               mols=None,
               residues=None,
               atoms=[1, 5, 6],
               ndxs=None,
               name=mySys.name + "_g1")
"""
4b. (optional) Create a Group 2 for comparative and interfacial analysis between the two groups
"""
group2 = Group(net=mySys.net,
               mols=None,
               residues=None,
               atoms=[3, 4],
               ndxs=None,
               name=mySys.name + "_g2")


################################################ Export Selections #####################################################


"""
5. Export the desired selections
    a. network: Export the reloadable network file created from a previously solved system
    b. pdb: Export a pdb file for the system (if a file was loaded it is a copy of that file)
    c. surfaces: Export all surfaces from the network
    d. no_sol_network_object: Export the outer surfaces for all of the 
    e. alter_atoms_script: Export a script to alter the atoms in pymol to what they are set to here
    f. export_groups: Export the groups created above
    g. export_interface: Export the interface between the two groups above 
"""
mySys.exports(groups=[group1, group2],
              network=True,
              pdb=True,
              surfaces=True,
              full_network_object=True,
              no_sol_network_object=True,
              alter_atoms_script=True,
              export_groups=True,
              export_interface=True)


################################################# Run the program ######################################################


# Main running guy. Go get em guy!
if __name__ == '__main__':
    pass


################################################# Open in pymol ########################################################


"""
6. Check the output directory for the outputs (default is vorpy/Data/User_data/) and open the outputs in
 a program that can handle molecule files and .off files like Pymol
 """
