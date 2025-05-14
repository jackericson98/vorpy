"""
API for calculations module - Core mathematical and geometric calculations
"""

from vorpy.src.calculations.calcs import (
    calc_dist,
    calc_angle,
    calc_tri,
    calc_tetra_vol,
    calc_cell_box,
    calc_cell_com,
    calc_cell_moi,
    calc_tetra_inertia,
    calc_spikes,
    calc_com,
    calc_length,
    calc_sphericity,
    calc_isoperimetric_quotient,
    calc_vol, 
    rotate_points
)

from vorpy.src.calculations.curvature import (
    gaussian_curvature,
    mean_curvature,
    calc_surf_tri_curvs
)

from vorpy.src.calculations.surf import (
    calc_surf_sa,
    calc_surf_func,
    calc_2d_surf_sa,
    calc_surf_tri_dists
)

from vorpy.src.calculations.vert import (
    calc_flat_vert,
    calc_vert,
    verify_site,
    verify_aw,
    verify_pow,
    verify_prm
)

from vorpy.src.calculations.sorting import (
    ndx_search
)

__all__ = [
    # Basic calculations
    'calc_dist',
    'calc_angle',
    'calc_tri',
    'calc_com',
    'calc_spikes',
    'calc_tetra_vol',
    'calc_cell_box',
    'calc_cell_com',
    'calc_cell_moi',
    'calc_tetra_inertia',
    'calc_length',
    'calc_sphericity',
    'calc_isoperimetric_quotient',
    'calc_vol',
    'rotate_points',
    # Curvature calculations
    'gaussian_curvature',
    'mean_curvature',
    'calc_surf_tri_curvs',
    
    # Surface calculations
    'calc_surf_sa',
    'calc_surf_func',
    'calc_2d_surf_sa',
    'calc_surf_tri_dists',
    
    # Vertex calculations
    'calc_flat_vert',
    'calc_vert',
    'verify_site',
    'verify_aw',
    'verify_pow',
    'verify_prm',

    # Sorting calculations
    'ndx_search'
] 