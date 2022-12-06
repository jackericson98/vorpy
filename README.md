# Vorpy

## Description


Vorpy is the first open-source comprehensive 3D Voronoi network generator for 3D spheres. It is designed with simulated molecular analysis in mind. Vorpy time-stamps of simulations of molecules and samples a network of points equidistant from the surfaces of all neighboring atoms in the frame. The result is a partitioning of the frame atom by atom and a set of atom-atom partitioning surfaces. These surfaces can then be used for analysis and interpretation. 

The general use case would be for analysis of simulated olecule files for volume, interface and void analysis. Once a frame is analyzed all of these data points can be deriven. For example, if a molecule is simulated in solution until it reaches a certain equilibrium, a frame can be processed through vorpy and tested for volume and voids. 

Created by Jack Ericson in collaboration with Georgia State University 

## Usage (GUI)

1. Download the repository and move to the main vorpy directory in a shell or command prompt
2. Install the requirements
   ```
   pip install requirements.txt
   ```
3. Run the following from the vorpy directory
   ```
   py vorpy.py
   ```
2. From the load screen of the gui load your pdb file (for testing go to the vorpy/Data/test_data folder)
3. From the build screen of the gui change the settings and then click build
4. From the analyze screen analyze the newly constructed network

## Usage (Jupyter Notebook)

1. Import the system object and other helpful functions
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
3. Find your file in the /vorpy/Data/User_data/ folder


## Visualization

Once built the network can be viewed in a number of ways. The first being trough vorpy's built in visualization functions (plot_atoms, plot_verts, plot_edges, plot_surfs) or one of the following

### Pymol

Pymol is currently the best way to view the data produced from vorpy and can be downloaded here:
https://pymol.org/2/

## VMD

I have never tried vmd
