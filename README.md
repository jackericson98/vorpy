# Vorpy

## Description


Vorpy is the first open-source comprehensive 3D Voronoi network generator for 3D spheres. It is designed with simulated molecular analysis in mind. Vorpy processes time-stamps of simulations of molecules by samples a network of points equidistant from the surfaces of all neighboring atoms in the system. This results in a partitioning of the system atom by atom and a provides a set of atom to atom partitioning surfaces for analysis. 

The general use case would be for analysis of simulated olecule files for volume, interface and void analysis. Once a frame is analyzed all of these data points can be deriven. For example, if a molecule is simulated in solution until it reaches a certain equilibrium, a frame can be processed through vorpy and tested for volume and voids. 

Created by Jack Ericson in collaboration with Georgia State University 

## Prerequisites
- Dependencies: >= python 3.9
- Instalation:
Move to the vorpy directory in a shell or command prompt:
   ```
   cd PATH/TO/vorpy
   ```
   Install the requirements
   ```
   pip install requirements.txt
   ```
Note: If installing the requirements.txt fails, retry with just numpy and matplotlib (ex: pip install matplotlib) and only use the command prompt, jupyter notebook or script provided below

## Usage

### Gui 

1. Move to the main vorpy directory in a shell or command prompt
2. Run the either of following
   Basic GUI
   ```
   py vorpy_gui.py
   ```
   Fancy GUI
   ```
   py vorpy_gui1.py
   ```
   
3. From the load screen of the gui load your pdb file (for testing go to the vorpy/Data/test_data folder)
4. From the build screen of the gui change the settings and then click build
5. From the analyze screen analyze the newly constructed network
6. Use a software like pymol to view the outputs

### Command Line 

1. Move to the main vorpy directory in a shell or command prompt
2. Run the following
   ```
   py vorpy.py
   ```
3. Follow the prompts

### Script

1. Open vorpy_script.py in an ide
2. Change the load/build/output settings using the comments 
3. Run the script

### Jupyter Notebook

1. Make sure you have Anaconda downloaded
2. Move to the vorpy directory in an Anaconda shell and type:
   ```
   jupyter notebook
   ```
3. Select the vorpy_jupnot.ipynb file
4. Follow the instructions in the file

### Visualization

Once built the network can be viewed in a number of ways. The first being trough vorpy's built in visualization functions (plot_atoms, plot_verts, plot_edges, plot_surfs) or one of the following

#### Pymol

Pymol is currently the best way to view the data produced from vorpy and can be downloaded here:
https://pymol.org/2/ . Once downloaded run the software through one of the processes above and drag the output files (.off) and the system files (.pdb) into the pymol frame. Be sure to use the "set_pymol_atoms.pml" script to get accurately set atom radii.

#### VMD

I have never tried vmd
## Citation
[![DOI](https://zenodo.org/badge/502126698.svg)](https://zenodo.org/badge/latestdoi/502126698)


## Contact
- Email: jericson1@gsu.edu
- Site: https://cas.gsu.edu/profile/greg-poon/
- Phone: (404)-413-5491
