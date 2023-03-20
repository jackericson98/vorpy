![image](https://user-images.githubusercontent.com/62311229/226451994-de2cd30f-4ee9-4d09-87f3-2fd2c1573b35.png)



# Vorpy

## Description


The first comprehensive Voronoi cell network generator for 3D spheres designed with simulated molecular dynamics (md) analysis in mind. It takes in text files containing the locations of atoms in simulated molecules and creates individually partitioned atomic Voronoi cells for analysis. Atoms are partitioned by sampling points along the Voronoi surfaces of neigboring atoms. Using these partitions, the volumes and surface areas of bodies and interfaces of interest (residue, protein/DNA, etc.) can be more accurately analyzed and visualized. 

# Usage

## Prerequisites

- Dependencies: >= python 3.9
- Basic Dependencies: numpy, scipy, matplotlib
- Instalation:
   Move to the vorpy directory in a shell or command prompt:
   ```
   cd PATH/TO/vorpy
   ```
   Install the requirements
   ```
   python3 -m pip install requirements.txt
   ```
   Note: If installing the requirements.txt fails, retry with just numpy and matplotlib (ex: pip install matplotlib) and only use the command prompt, jupyter notebook     or script provided below



## Command Line 

1. Move to the main vorpy directory in a shell or command prompt
2. Run the following
   ```
   py vorpy.py
   ```
3. Follow the prompts

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
