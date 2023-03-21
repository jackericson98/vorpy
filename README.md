![image](https://user-images.githubusercontent.com/62311229/226451994-de2cd30f-4ee9-4d09-87f3-2fd2c1573b35.png)



# Vorpy

## Description

A comprehensive Voronoi cell network generator for 3D spheres designed with simulated molecular dynamics (md) analysis in mind. It takes in text files containing the locations of atoms in simulated molecules and creates individually partitioned atomic Voronoi cells for analysis. Atoms are partitioned by sampling points along the Voronoi surfaces of neigboring atoms. Using these partitions, the volumes and surface areas of bodies and interfaces of interest (residue, protein/DNA, etc.) can be more accurately analyzed and visualized. 

why

When studying macro 


This allows chemistst to gather more accurate volume and area measurements. Visualitions of the range of influence of the Van der Waals radius can show the effect of 
Jack Ericson
Chemistry Department 
Georgia State University



## Usage

### Prerequisites

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



### Command Line 

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

### Pymol

[Pymol](https://pymol.org/2/) is one of the best ways to view the data produced from vorpy and can be downloaded here:  . Once downloaded run the software through one of the processes above and drag the output files (.off) and the system files (.pdb) into the pymol frame. Be sure to use the "set_pymol_atoms.pml" script to get accurately set atom radii.


## Theory

Full documentation can be found in the docs.md file

## Citation

[![DOI](https://zenodo.org/badge/502126698.svg)](https://zenodo.org/badge/latestdoi/502126698)


## Contact
- Email: jericson1@gsu.edu
- Site: https://cas.gsu.edu/profile/greg-poon/
- Phone: (404)-413-5491
