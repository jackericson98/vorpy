"""
Network module for Voronoi network construction and analysis.
"""

from vorpy.src.network.network import Network
from vorpy.src.network.analyze import analyze
from vorpy.src.network.build_edge import build_edge
from vorpy.src.network.build_net import build
from vorpy.src.network.build_surf import build_surf
from vorpy.src.network.build_surfs import build_surfs
from vorpy.src.network.edge_project import edge_project
from vorpy.src.network.fill import calc_surf_point
from vorpy.src.network.fill import calc_surf_point_from_plane
from vorpy.src.network.find_verts import find_verts
from vorpy.src.network.find_v0 import find_v0
from vorpy.src.network.fast import find_site_container
from vorpy.src.network.slow import find_site_container_slow
from vorpy.src.network.slow import find_site
from vorpy.src.network.split_net import split_net
from vorpy.src.network.split_net import split_net_slow
from vorpy.src.network.split_net import combine_nets

__all__ = [
    # Main Network Class
    'Network',
    
    # Network Analysis
    'analyze',
    
    # Network Construction
    'build',
    'build_edge',
    'build_surf',
    'build_surfs',
    
    # Vertex Finding
    'find_verts',
    'find_v0',
    'find_site_container',
    'find_site_container_slow',
    'find_site',
    
    # Surface Calculations
    'calc_surf_point',
    'calc_surf_point_from_plane',
    'edge_project',
    
    # Network Splitting
    'split_net',
    'split_net_slow',
    'combine_nets'
]
