"""Explicit v1 object field registry. Additions require schema compatibility review.

Geometry tables carry named primitive scientific columns. Relationships in their
adjacency columns are network IDs, not Python objects. Disposable draw caches
are excluded by the table codec.
"""

GROUP_FIELDS = ('name', 'group_id', 'ball_ndxs', 'settings', 'interface_metadata', 'atms', 'mols', 'chns', 'rsds', 'sa', 'vol', 'density', 'mass', 'com', 'avg_mean_curv', 'avg_gauss_curv', 'int_mean_curv_face', 'int_mean_curv_edge', 'int_mean_curv', 'int_mean_curv_sq', 'int_gauss_curv_face', 'int_gauss_curv_edge', 'int_gauss_curv_vertex', 'int_gauss_curv', 'euler_characteristic', 'gauss_bonnet_expected', 'gauss_bonnet_error', 'gauss_bonnet_relative_error', 'boundary_is_complete', 'boundary_is_closed', 'boundary_is_manifold', 'boundary_is_orientable', 'boundary_component_count', 'boundary_face_count', 'boundary_edge_count', 'boundary_vertex_count', 'boundary_genus', 'boundary_missing_faces', 'boundary_missing_edges', 'boundary_missing_vertices', 'boundary_nonmanifold_edges', 'boundary_nonmanifold_vertices', 'vdw_vol', 'vdw_com', 'spatial_moment', 'moi', 'layer_atoms', 'layer_net_atoms', 'layer_verts', 'layer_edges', 'layer_surfs', 'layer_info', 'layer_body_topology', 'layer_time', 'oriented_int_mean_curv', 'oriented_avg_mean_curv', 'convex_sa', 'concave_sa', 'flat_sa', 'unoriented_sa', 'flat_surface_count', 'curved_unoriented_sa', 'curved_unoriented_surface_count', 'convex_sa_fraction', 'concave_sa_fraction', 'flat_sa_fraction', 'convex_int_mean_curv', 'concave_int_mean_curv', 'abs_int_mean_curv', 'convex_avg_mean_curv', 'concave_avg_mean_curv', 'oriented_surface_count', 'unoriented_surface_count', 'buried_water_metadata')

INTERFACE_FIELDS = ('group1_name', 'interface_id', 'settings', 'name', 'group1_indices', 'group2_indices', 'water_geometries', 'water_topology', 'group2_name')

NETWORK_FIELDS = ('group', 'iface_grps', 'group_name', 'settings', 'metrics', 'box',
    'network_mode', 'loaded_from_logs', 'completion_kind', 'vert_timing',
    'surface_timing', 'build_surf_timing', 'analysis_timing', 'max_curv', 'interface_vertex_state',
    'balls', 'verts', 'edges', 'surfs')
SYSTEM_FIELDS = ('name', 'type', 'foam_box', 'foam_data', 'ndxs', 'ndx_names',
    'atom_names', 'chn_names', 'res_names', 'group_names', 'elements',
    'element_radii', 'special_radii', 'round_to', 'export_type', 'file_type',
    'max_atom_rad', 'frame_index', 'frame_count', 'loaded_frame_count', 'data',
    'balls', 'segments')
RESIDUE_FIELDS = ('atoms', 'name', 'seq', 'id', 'print_name')
CHAIN_FIELDS = ('atoms', 'name', 'vol', 'sa')
