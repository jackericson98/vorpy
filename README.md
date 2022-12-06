# Vorpy

## Description


Vorpy is the first open-source comprehensive 3D Voronoi network generator for 3D spheres. It is designed with simulated molecular analysis in mind. Vorpy processes time-stamps of simulations of molecules by samples a network of points equidistant from the surfaces of all neighboring atoms in the system. This results in a partitioning of the system atom by atom and a provides a set of atom to atom partitioning surfaces for analysis. 

The general use case would be for analysis of simulated olecule files for volume, interface and void analysis. Once a frame is analyzed all of these data points can be deriven. For example, if a molecule is simulated in solution until it reaches a certain equilibrium, a frame can be processed through vorpy and tested for volume and voids. 

Created by Jack Ericson in collaboration with Georgia State University 

## Prerequisites
- Dependencies:
- Instalation:
Move to the vorpy directory in a shell or command prompt:
```
cd PATH/TO/vorpy
```
Install the requirements
```
pip install requirements.txt
```

## Usage

### Gui 

1. Download the repository and move to the main vorpy directory in a shell or command prompt
2. Install the requirements (see prerequisites section)
3. Run the following from the vorpy directory
   ```
   py vorpy.py
   ```
2. From the load screen of the gui load your pdb file (for testing go to the vorpy/Data/test_data folder)
3. From the build screen of the gui change the settings and then click build
4. From the analyze screen analyze the newly constructed network

### Jupyter Notebook

1. Download the repository and move to the main vorpy directory in a shell or command prompt
2. Install the requirements (see prerequisites section)
2. Import the system object and other helpful functions
   ```
   from System.system import *
   ```
2. Load your file:

   ```
   sys = System("PATH/TO/FILE.pdb")
   ```
2. Build the network:
   ```
   sys.build_network(sol_verts=False, surf_res=0.2, max_vert=10, box_size=1.5)
   ```
3. Plot the compnents
   ```
   from Visualize.mpl_visualize import *
   fig = plt.figure()
   ax = fig.add_subplot(projection="3d")
   plot_atoms(net.atoms[0:1], fig=fig, ax=ax)
   plot_verts(net.atoms[0].verts, fig=fig, ax=ax)
   plot_edges(net.atoms[0].edges, fig=fig, ax=ax)
   plot_surfs(net.atoms[0].surfs, fig=fig, ax=ax, Show=True)
   ```
4. Analyze the network
   The user can analyze atoms, surfaces, vertices or any other network object individually:
   ```
   my_vol = 0
   for atom in sys.mols[0]:
      my_vol += atom.vol
   print("Molecule {} has volume: {} Angstroms Cubed".format(sys.mol_names[0], my_vol)
   
   ```
   The user can also create groups for analysis
   ```
   g1 = Group(net=sys.net, atoms=sys.atoms[22:45], name="My_group1")
   g1.get_info()
   print(g1.body_vol)
   
   ```
 5. Export data
    Exporting the network will allow the network to be re-loaded without having to be calculated
    ```
    sys.export_net()
    ```
    Similarly the user can export vertices
    ```
    sys.export_verts()
    ```
    The user can export the surfaces associated with a group
    ```
    sys.export_selection(group1=g1, info=True)
    ```
    or the interface between two groups
    ```
    g2 = Group(net=sys.net, atoms=sys.atoms[:22], name="My_group2")
    sys.export_selection(group1=g1, group2=g2, info=True)
    ```

### Visualization

Once built the network can be viewed in a number of ways. The first being trough vorpy's built in visualization functions (plot_atoms, plot_verts, plot_edges, plot_surfs) or one of the following

#### Pymol

Pymol is currently the best way to view the data produced from vorpy and can be downloaded here:
https://pymol.org/2/ . Once downloaded run the software through one of the processes above and drag the output files (.off) and the system files (.pdb) into the pymol frame. Be sure to use the "set_pymol_atoms.pml" script to get accurately set atom radii.

#### VMD

I have never tried vmd
## Citation

## Contact
- Email: jericson1@gsu.edu
- Site: https://cas.gsu.edu/profile/greg-poon/
- Phone: (404)-413-5491
