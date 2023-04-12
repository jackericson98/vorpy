![image](https://user-images.githubusercontent.com/62311229/226769798-5c1247aa-2bb6-4581-b3d5-ec334a92d28f.png)

# Vorpy

## Description

Volumetric, Area, and distance calculator for analysis of simulated molecules. Using the locations and radii of atoms in simulated molecules individually partitioned atomic cells are created for visualization and comparison. The atoms are partitioned by sampling points along the Voronoi, Power, or Delaunay surfaces of neigboring atoms. Using these partitions, data about atomic bodies and interfaces of interest (residue, protein/DNA, etc.) can be more accurately calculated.

Volume changes directly govern pressure-based technologies in microbiological sterilization and food processing. At the molecular level, the relative incompressibility of bulk water means that volume changes offer a sensitive window into the hydration and dynamic properties of macromolecular conformations in aqueous solution. In contrast, resolution of calorimetric observables such as entropy and heat capacity into hydration and dynamic components between solute and solvent remains theoretically challenging.

## Usage

### Prerequisites

Install the requirements
```
python3 -m pip install requirements.txt
```

### Basic Example

```
C:/.../vorpy> python3 vorpy.py EDTA_Mg.pdb
```
When just an atom file is provided, the program runs with the defaults for the flags.

### Change Group

```
C:/.../vorpy> python3 vorpy.py EDTA_Mg.pdb -g a 5
```
In this example the 5th atom is selected for construction

### Change Setting
```
C:/.../vorpy> python3 vorpy.py EDTA_Mg.pdb -s ar c 2
```
In this example all carbon atoms raddi are changed to 2 angstroms

### Specify exports
```
C:/.../vorpy> python3 vorpy.py EDTA_Mg.pdb -e v -e e -e i
```
In this example, the vertices, edges and information for the load file are exported

### Commands

Use these (separated by spaces) after the atom file

![image](https://user-images.githubusercontent.com/62311229/228362286-d19bb842-d1b5-40e9-9270-20b2ed666366.png)



### Command Line 

1. Move to the main vorpy directory in a shell or command prompt
2. Run the following
   ```
   py vorpy.py
   ```
3. Follow the prompts

### Jupyter Notebook

Run the User_Example.ipynb file from the ./Visualize/Notebooks folder

## Visualization

### Pymol

[Pymol](https://pymol.org/2/) is one of the best ways to view the data produced from vorpy and can be downloaded here:  . Once downloaded run the software through one of the processes above and drag the output files (.off) and the system files (.pdb) into the pymol frame. Be sure to use the "set_pymol_atoms.pml" script to get accurately set atom radii.


## Theory

In a simple Delaunay 3D partitioning the partitioned atomic cells would be constructed using a network of intersecting planes from neighboring atoms. For example, Atom 1 at (0, 0, 0) and Atom 2 at (2, 0, 0) would have a partitioning plane between them, normal to the x-axis including the point (1, 0, 0). 

In a weighted Delaunay 3D partitioning, the radii of the atoms are considered. For example, if Atoms 1 and 2 from the above example have radii 0.5 and 1, the partitioning plane would remain normal to the x-axis but would be closer to Atom 1’s location containing the point (0.75, 0, 0). 

In a Weighted 3D Voronoi S-Network (Medvedev et al, 2006) the partitioning surfaces are curved hyperboloidal surfaces whose curvature is determined by proximity and the relative radii of neighboring atoms. The program is designed to take in files containing atom location (centers) and radii and return a weighted 3D Voronoi S-Network partitioning of all points in space closest to each atom.

A 2D depiction of the differences between the partitioning schemes can be seen below:

![image](https://user-images.githubusercontent.com/62311229/228368891-b6af041c-b264-4be2-b233-547d3aa7cf30.png)


**Figure 1:**
Different partitioning schemes for two circles with radii 1 and 2. Delaunay (left, red) Power Distance (left middle, blue) Voronoi (right middle, purple) all three (right). Solid green lines represent vectors that are equal in size. Dashed lines represent the projections to the atoms centers. 

![image](https://user-images.githubusercontent.com/62311229/228373404-15fa0a1d-6e3a-456f-8e68-919291615b37.png)

**Figure 2:**
Different partitioning schemes for 5 circles with radii (top to bottom) 2, 1.5, 0.5, 1, 0.5. Delaunay (left, red) Power Distance (left middle, blue) Voronoi (right middle, purple) all three (right).

In the power distance and weighted voronoi digrams, smaller atoms bend partitionings toward themselves and create smaller cells. Due to the radial nature of atomic forces, these partitionings represent a better picture of reality for intera-atomnic interaction thresholds.   



### Process outline

The general overview of how the Voronoi S network is calculated is as follows:

![image](https://user-images.githubusercontent.com/62311229/226765763-623c310f-4ff0-4dee-98af-8e078c6c950c.png)

First, an atom file is loaded, the information is sorted (location, element, chain, residue, etc.), and the atom objects are created. Next, each vertex (point in space equidistant to the surfaces of 4 neighboring atoms) is found. Shared atoms between vertices are used to determine valid edges and valid surfaces. If two vertices share 3 atoms, we know there exists an edge (hyperbolic curve equidistant to the surfaces of 3 neighboring atoms). If the number of edges with a pair of atom matches the number of vertices with the same pair, there exists a surface (surface equidistant to the surfaces of 2 neighboring atoms). These possible surfaces are then constructed and used for determining surface area, volume, and curvature of these inter-atomic partitions. Lastly, different combinations of surfaces (in the form of .off files), selected atoms (in the form of .pdb files), a network checkpoint (in the form of .csv file), and network information (in the form of .txt files) are then exported at the user’s request.

### Vertex Calculation

Calculating the vertex made by four neighboring atoms is the foundation of the entire vertex finding process. The Voronoi-vertex is the center of a sphere which is tangential to the four input atoms. We want to determine (x, y, z, R) such that:

![image](https://user-images.githubusercontent.com/62311229/226765948-b57d697b-24de-4cd1-89f0-ce48db6dcdef.png)

### Surface Calcuation

The surface between two atoms can be represented as a hyperboloid, ranging from flat to elliptical in the extreme cases of equal radii atoms and enveloped atoms respectively. The calculation of the bisector function between two atoms followed the methods provided in (Hu et al, 2017) in which the solution to the equations for two spheres:

![image](https://user-images.githubusercontent.com/62311229/226766149-3a4a7fdd-db77-4921-a638-ac81b0d74334.png)

provides a general function for a hyperboloid of the form:

![image](https://user-images.githubusercontent.com/62311229/226766198-e5fc6201-504e-4bd6-8ddd-687c7f1be4ac.png)


## Citation

[![DOI](https://zenodo.org/badge/502126698.svg)](https://zenodo.org/badge/latestdoi/502126698)


## Contact
- Email: jericson1@gsu.edu
- Site: https://cas.gsu.edu/profile/greg-poon/
- Phone: (404)-413-5491
