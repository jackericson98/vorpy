![image](https://user-images.githubusercontent.com/62311229/226769798-5c1247aa-2bb6-4581-b3d5-ec334a92d28f.png)

# Vorpy

## Overview
Vorpy is a command line interface tool designed to process molecular files for computational chemistry and visualization tasks. It supports a variety of file formats and offers extensive options for loading files, setting simulation parameters, selecting specific molecular structures, conducting specialized calculations, and configuring output details.

## Installation
Download the Vorpy tool from the repository and ensure you have Python installed on your system.

## Usage

For first time users or virtual environment users check to see you have the requirements installed for python

    python -m pip install -r requirements.txt

The basic command structure for running Vorpy is:

    python vorpy.py <file> [options]

### File
- The first argument after `vorpy.py` should be the file address of the ball or atom file.
- If the file is located in the `Data/test_data` folder, specify the file name without the path or extension.
- Accepted file extensions include `.pdb`, `.mol`, `.gro`, `.cif`.

### Options
#### Load Flag `-l`
    -l <file>
Load additional files like vertex files from previous runs, log files, Voronota vertex files, or GROMACS index files.

#### Settings Flag `-s`
    -s <setting value>
Adjust various simulation parameters:
- `nt` - Network Type: Default = Additively Weighted `aw`, Power `pow`, Primitive `prm`, or Compare `com 'type1' 'type2'`
- `mv` - Maximum Vertex: Default = `40`
- `bm` - Box Multiplier: Default = `1.25`
- `sr` - Surface Resolution: Default = `0.2` 
- `sc` - Surface Color Map: Default = `viridis`, `plasma`, `rainbow`, or any other [matplotlib colormap](https://matplotlib.org/stable/gallery/color/colormap_reference.html) (note: '_r' inverts the scheme)
- `ss` - Surface Coloring Scheme: Default = curvature `curv`, inside vs outside spheres `nout`, distance from center `dist`
- `sf` - Surface Coloring Scale: Default = linear `lin`, log `log`, squared `square`, cube `cube`
- `ar` - Adjust Radii: `'element' 'value'` or `'atom name' 'value'` or `'residue' 'atom name' 'value'`. To see the current values for defaults for atomic radii go to the radii file (radii.py) or enter the radii flag`-r`

#### Group Flag `-g`
    -g <identifier>
Select specific balls or molecular elements:
- `b` - Ball Identifier. Used with a ball index `'index'` or range of indices `'index1'-'index2'`.
- `a` - Atom Identifier. Used with an atom element `'element'`, element name `'element name'`, index `'index'`, or range of indices `'index1'-'index2'`.
- `r` - Residue Identifier. Used with a residue name `'residue name'` and sequence number `'sequence number'` (optional), index `'index'`, or range of indices `'index1'-'index2'`.
- `c` - Chain Identifier. Used with a chain name `'chain name'`, index `'index'`, or range of indices `'index1'-'index2'`.

Note: If multiple of the above components are desired in the same group use the `and` qualifier between components. If multiple groups are desired use multiple group flags.  

#### Calculate Flag `-c`
    -c <calculation_type>
Identify additional calculations:
- Interfaces between all combinations of groups (iface, ifc, i)
- Calculate vertices up to the specified number of layers (layers)

#### Exports Flag `-e`
    -e <export_type>
Specify the intensity and type of exports:
- Groups of Exports: Default = `large`, `small`, `medium`, `all`
- Export choices : 

   Molecule File - `pdb`, `mol`, `cif`, `gro`, Set Atoms Radii PyMol Script - `set_atoms`, Group Information - `info`, Network Logs - `logs`, All Surfaces in One File - `surfs`, All Surfaces in Separate Files - `sep_surfs`, All Edges in One File - `edges`, All Edges in Separate Files - `sep_edges`, All Vertices in One File - `verts`, All Vertices in Separate Files - `sep_verts`, Surrounding Surfaces - `shell`, Surrounding Edges - `shell_edges`, Surrounding Vertices - `shell_verts`, Group Atoms - `atoms`, Atoms Surrounding Group - `surr_atoms`


## Notes
- Each option flag and its arguments must be separated by spaces.
- To use multiple commands for a single option use 'and' or repeat the flag (except for groups to avoid creating multiple groups).
- Any range can be set with a hyphen and no space (e.g. `-g a 0-100` is a group of the first 101 atoms) 

### Examples

Separately solve the tyrosine 2 and methionine 1 residues of the cambrin molecule, calculate their interface, and export the large export type of the results

    python vorpy.py cambrin -s sr 0.05 and mv 80 -g tyr 2 -g met 1 -c iface -e large

Calculate the primitive and power networks for the mg atom in the EDTA_Mg molecule and compare the difference

    python vorpy.py EDTA_Mg -s nt compare prm pow -g mg

Solve the network for hairpin and export the shell with the inside and outside parts of the surfaces highlighted at a high resolution

    python vorpy.py hairpin -s ss nout and sr 0.01 -e shell and pdb


## License
MIT License

Copyright (c) 2022 John Ericson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Citation

[![DOI](https://zenodo.org/badge/502126698.svg)](https://zenodo.org/badge/latestdoi/502126698)


## Contact
- Email: [jericson1@gsu.edu](mailto:jericson1@gsu.edu)
- Site: https://cas.gsu.edu/profile/greg-poon/
- Phone: +1 (404)-413-5491
