"""
API for output module - File writing and data export functionality
"""

from vorpy.src.output.atoms import (
    write_atom_cells,
    write_atom_radii
)

from vorpy.src.output.edges import (
    write_edges
)

from vorpy.src.output.surfs import (
    write_surfs
)

from vorpy.src.output.verts import (
    write_off_verts
)

from vorpy.src.output.draw import (
    draw_line, 
    draw_edge
)

from vorpy.src.output.output import (
    export_micro,
    export_tiny,
    export_med,
    export_large,
    export_all,
    other_exports
)

from vorpy.src.output.color_tris import color_tris

__all__ = [
    # Atom output
    'write_atom_cells',
    'write_atom_radii',
    
    # Edge output
    'write_edges'
    
    # Surface output
    'write_surfs'
    
    # Vertex output
    'write_off_verts'
    
    # Drawing functions
    'draw_edge',
    'draw_line',
    
    # Export functions
    'export_micro',
    'export_tiny',
    'export_med',
    'export_large',
    'export_all',
    'other_exports',
    
    # Coloring
    'color_tris'
]
