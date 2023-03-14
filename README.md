# Vorpy

## Description


Vorpy is a comprehensive Voronoi cell network generator for 3D spheres designed with simulated molecular dynamics (md) analysis in mind. It works by partitioning the space between simulated atoms and using these partitions to calculate volumes and surface areas of bodies and interfaces of interest (residue, protein/DNA, etc.). By probing inter-atomic partitions it is possible to better understand a number of features of simulated molecules ranging inter-atomic influence in different chemical bonds to the evolution of a protein-DNA complex. 

## Prerequisites

- Dependencies: >= python 3.9
- Basic Dependencies: Numpy, Matplotlib
- Instalation:
   Move to the vorpy directory in a shell or command prompt:
   ```
   cd PATH/TO/vorpy
   ```
   Install the requirements
   ```
   python3 -m pip install numpy
   python3 -m pip install matplotlib
   ```
   Note: If installing the requirements.txt fails, retry with just numpy and matplotlib (ex: pip install matplotlib) and only use the command prompt, jupyter notebook     or script provided below

## Usage

### Command Line 

1. Move to the main vorpy directory in a shell or command prompt
2. Run the following
   ```
   py vorpy.py
   ```
3. Follow the prompts

### Script

1. Open vorpy_script.py in an ide
2. Change the load/build/output settings following the provided comments for instruction
3. Run the script

### Jupyter Notebook

1. Make sure you have Anaconda downloaded
2. Move to the vorpy directory in an Anaconda shell and type:
   ```
   jupyter notebook
   ```
3. Select the vorpy_jupnot.ipynb file
4. Follow the instructions in the file

## Visualization

Once built the network can be viewed in a number of ways. The first being trough vorpy's built in visualization functions (plot_atoms, plot_verts, plot_edges, plot_surfs) or one of the following

### Pymol

Pymol is currently the best way to view the data produced from vorpy and can be downloaded here: https://pymol.org/2/ . Once downloaded run the software through one of the processes above and drag the output files (.off) and the system files (.pdb) into the pymol frame. Be sure to use the "set_pymol_atoms.pml" script to get accurately set atom radii.

### plot_net

For small systems, the plot_net function in Visualize.mpl_visualize can be used

## Documentation

Full documentation can be found in the docs.md file

## Citation

[![DOI](https://zenodo.org/badge/502126698.svg)](https://zenodo.org/badge/latestdoi/502126698)


## Contact
- Email: jericson1@gsu.edu
- Site: https://cas.gsu.edu/profile/greg-poon/
- Phone: (404)-413-5491
