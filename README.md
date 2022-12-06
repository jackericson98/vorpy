# Vorpy

## Description


Vorpy is a comprehensive 3D Voronoi cell generator designed for molecular analysis. It works by sampling a network of points equidistant from the surfaces of neighboring atoms. Then constructs the cells around a group of atoms or interfaces between two groups of atoms for analysis. 


## Installation



## Usage (GUI)

1. Run the following from the vorpy directory
   ```
   py vorpy.py
   ```
2. From the load screen of the gui load your pdb file
3. From the build screen of the gui change the settings and then click build
4. From the analyze screen analyze the newly constructed network

## Usage (Jupyter Notebook)

1. Load your file:

   ```
   sys = System("PATH/TO/FILE.pdb")
   ```
2. Build the network:
   ```
   sys.build_network(sol_verts=False, surf_res=0.2, max_vert=10, box_size=1.5)
   ```
3. Find your file in the /vorpy/Data/User_data/ folder

## Credits

...

## Outputs



