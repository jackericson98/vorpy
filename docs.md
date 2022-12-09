# Vorpy Documentation

# Table of contents

1. [Objects](#objects)
   1. [System](#system)
   2. [Network](#network)
   3. [Atom](#atom)
   4. [Group](#group)
   5. [Vertex](#vertex)
   6. [Edge](#edge)
   7. [Surface](#surface)
2. [Guis](#guis)
   1. [tkinter](#tk)
   2. [customtkinter](#ctk)
   3. [Flask](#flask)
3. [Functions](#functions)



## Objects <a name="objects"></a>

### System <a name="system"></a>
#### Attributes:
-	file – Address for the system file (.pdb, .gro, .mol, .cif)
-	atoms – List of locations and radii or Atom objects to be loaded into the system
-	vert_file – Address for the vertices file for the system that has been loaded
-	net_file - Address for the network file for the system that has been loaded
-	ndx_file - Address for the index file for the system that has been loaded. Gromacs format
-	frame_files - Address for the atom frame files for the system that has been loaded
-	output_directory – Folder for the export data 
-	gui - The GUI object (tkinter) associated with loading the system and loading/creating the network
#### Methods:
1. load_sys (self, file=None) - Loads the atoms and data from .pdb, .gro, .mol or .cif file types
   1. Parameters

Parameters:
        i.	file – Address for system file (.pdb, .gro, .mol, .cif)

2. load_verts (self, file=None)
Loads vertices from a previously saved vertex file
a.	Parameters:
i.	file – Address for the vertex file (vorpy style)

3. load_net (self, file=None, verts_only=False)
Loads a previously solved vorpy generated network into the system
a.	Parameters:
i.	file – Address for the network file (vorpy style)
4. load_ndx (self)
5. load_sys_atoms
6. random_system
7. sort_atoms
8. build_network
9. export_verts
10. export_net
11. export_selection
12. set_output_directory
13. initial_export
14. show_sys
15. show_net


### Network <a name="network"></a>
#### Attributes:
- sys
- atoms
- verts
- edges
- surfs
- doublets
- groups
- name
- box
- sub_boxes
- sub_box_size
- atoms_box
- max_atom_rad
- atom_ndxs
- min_dist
- max_vert
- box_size
- parallelize
- sol_verts
- curved_faces
- flat_faces
- verts_loaded
- cpu_time
- my_time

#### Methods:
1.	calc_box
2.	sort_atoms
3.	get_atoms
4.	connect
5.	find_verts
6.	build_edges
7.	build_surfs
8.	analyze
9.	build
10.	rebuild_net

### Atom <a name="atom"></a>
Atom object used to store the attributes and methods associated with each atom in a system
#### Attributes 
-	loc
-	rad
-	cell_vol
-	box
-	sys
-	verts
-	surfs
-	edges
-	load_ndxs
-	element
-	chain
-	res
-	res_seq
-	name
-	occupancy
-	t_fact
-	seg_id
Methods: 
1.	get_radius*

### Group <a name="group"></a>
#### Attributes: 
-	net
-	atoms
-	selects
-	select_strs
-	name
-	body_surfs
-	body_sa
-	body_vol
-	outer_body_atoms
-	surr_body_atoms
-	bff
-	iface_surfs
-	ifac_sa
#### Methods:
1.	get_info
2.	set_name
3.	add_sele
4.	undo_sele
### Vertex <a name="vertex"></a>
#### Attributes: 
-	net
-	atoms
-	edges
-	surfs
-	ndx
-	loc
-	rad
-	load_ndxs
-	doublet
-	d_type
-	loc2
-	rad2
-	flat_faced
-	ff_atoms
#### Methods:
1.	calc_vert
2.	make_ff_atoms
3.	calc_ff_vert

### Edge <a name="edge"></a>
#### Attributes: 
- ndx - Indices of the edges atoms in the net.atoms list used for identification 
- net - Network object from which the edge is constructed
- atoms - List of the edges atoms used for construction and identification
- verts - List of 2 vertices at the endpoints of the edge
- surfs - List of the (up to 3) surfaces attached to the edge
- load_ndxs - Lost of object load indices
- loc
- rad
- points
- pv0
- pv1
- pa
- doublet
- loc2
- 
#### Methods:
1.	get_loc
2.	find_pvals
3.	project
4.	build

### Surface <a name="surface"></a>
#### Attributes: 
-	ndx
-	net
-	atoms
-	verts
-	edges
-	load_ndxs
-	func
-	perimeter
-	points
-	flat_points
-	pflat_points
-	tris
-	sa
-	rn
-	center
-	com
-	doublet
-	flat
#### Methods:
1.	calc_func
2.	build
3.	build_vta

## GUIS <a name=guis></a>

### Tkinter <a name="tk"></a>
#### Attributes:
##### Methods:
1. load_sys_button
2. load_frames_button
3. build_atoms_button
4. load_verts
5. load_net
6. load_analyze_button
7. change_radius_button
8. load_ndx
9. build_network_button
10. change_output_directory
11. reset_all
12. undo_last
13. add_mol_button
14. add_res_button
15. add_atom_button
16. add_ndx_button
17. add_cell_group
18. add_interface_g1
19. add_interface_g2
20. export_cell
21. export_iface 
### CustomTkinter <a name="ctk"></a>
#### Attributes:
#### Methods:
1. load_sys_button
2. load_ndx_button
3. load_net_button
4. load_verts_button
5. change_atom_radius 
6. build_net_button	
7. set_out_dir_button	
8. set_show_list_mol, set_show_list_res, set_show_list_atom, set_show_list_ndx
9. get_current_selection_atoms	
10. add_g1_button 
11. add_g2_button 
12. undo_g1_button 
13. undo_g2_button 
14. reset_g1_button 
15. reset_g2_button 
16. export_g1_button 
17. export_g2_button 
18. export_iface_button
### Flask webapp <a name="flask"></a>



## Functions:
Load Functions: <a name="load_funcs"></a>
1.	Read atom functions:
a.	read_pdb(sys, file=None) 
Adds the atom data from pdb files and adds Atom objects to the system
i.	Parameters:
1.	sys – System object

b.	read_gro (sys, file=None)
Adds the atom data from gro files and adds Atom objects to the system
i.	Parameters:
1.	sys – System object

c.	read_mol (sys, file=None)
Adds the atom data from mol files and adds Atom objects to the system
i.	Parameters:
1.	sys – System object

d.	read_cif (sys, file=None)
Adds the atom data from cif files and adds Atom objects to the system
i.	Parameters:
1.	sys – System object

2.	Read Vorpy saved files
a.	read_verts (net, file)
Reads saved vertices files and adds the vertex data to the network provided
i.	Parameters:
1.	net – Network object
2.	file – file address for the system’s saved file address

b.	read_net (net, file, verts_only=False)
Adds network data from a previously saved network file for analysis
i.	Parameters:
1.	net – Network object
2.	file – File address for the network file
3.	verts_only – Boolean for whether to rebuild the surfaces or not

3.	Read other input files
a.	read_ndx (sys, file=None)
Takes in an index file and loads it into the list of indices in the given system.
i.	Parameters:
1.	sys – System object

b.	read_vta_data (sys, ball_file, vert_file)
Takes in Voronota data files and adds them to the system
i.	Parameters:
1.	sys – System object
2.	ball_file – Voronota ball file for the given system atoms
3.	vert_file – Voronota vertex file for the given system atoms

### Build functions: <a name="build_funcs"></a>
1. find_verts (net, a0=None) - Recursively searches vertices throughout the network
      1. Parameters:
         1. net – Network object
         2. a0 – The initial atom to start with for the vertex search
 <br></br>
2. find_site (net, edge_atoms, vn_1=None) - Follows the given edge and finds the other vertex. Returns vertex and it’s index 
   1. Parameters:
      1. net – Network object
      2. edge_atoms – The atoms used to find the next site
      3. vn_1 – The previous vertex on the edge

3. verify_site () -  	
4. find_v0 () - 
#### Connect network functions
1. build 
2. make_objects 
3. connect 
4. doublify
#### Build surface functions
1. make_mesh 
2. build_perimeter 
3. fill_mesh 
4. find_next_point 
5. calc_surf_point 
6. get_com 
7. find_simps 
8. filter_tris 
9. calc_tri_circ 
10. tri_within
### Analyze Functions: <a name="analyze_funcs"></a>
1. calc_sas	
2. calc_vols	
3. calc_curve	
4. find_sol_layers
### Export functions:<a name="export_funcs"></a>
1. export_verts	
2. export_net	
3. export_mySys	
4. export_iface	
5. export_body	
6. write_surfs	
7. write_pdb	
8. set_output_directory	
9. set_pymol_atoms
### Calculator functions:
1. calc_dist	
2. calc_angle	
3. calc_com	
4. calc_edges_com	
5. calc_circ	
6. calc_tri	
7. calc_bisector_val	
8. calc_tetra_vol	
9. calc_sa	
10. calc_vol 
11. inv_jac
12. rotate_points 
13. search_verts 
14. check_surf 
15. check_edge 
16. check_vert
### Visualization functions: <a name="vis_funcs"></a>
1. plot_atoms	
2. plot_verts	
3. plot_edges	
4. plot_surfs	
5. plot_simps	
6. setup_plot
