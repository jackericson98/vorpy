# Vorpy

## Description


Vorpy is a comprehensive 3D Voronoi cell generator designed for molecular analysis. Vorpy works by creating a network of 
vertices and edges vorpy creates a network it uses to build inter atomic interfacial surfaces that can be constructed 
into cells. These surfaces, volumes, curvatures, vertices, etc. can all be used to better understand molecular 
structures in a simulated environment.  


## Installation (Jupyter Notebook)

1. Make sure you have python 3.9 installed as well as jupyter notebook. I recommend getting the whole anaconda suite (it's free) and launching from there.
2. Clone the repository into a working directory
3. Import the system object from the vorpy directory:
    ```commandline
    from vorpy.System import System
    ```


## Usage

1. Load your file:

   ```
   sys = System("PATH/TO/FILE.pdb")
   ```
2. Build the network:
   ```
   sys.net.build(min_dist=0.5)
   ```
3. Run the analysis function
   ```
   sys.net.analyze()
   ```
4. 

## Credits

...

## Analysis

1. Determine the volume of closed cell molecules Voronoi S-networks.
2.  


## Tests


