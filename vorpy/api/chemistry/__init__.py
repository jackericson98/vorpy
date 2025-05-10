from vorpy.src.chemistry.chemistry_interpreter import residue_names
from vorpy.src.chemistry.chemistry_interpreter import element_names
from vorpy.src.chemistry.chemistry_interpreter import my_masses
from vorpy.src.chemistry.chemistry_interpreter import residue_atoms
from vorpy.src.chemistry.Standard import special_radii as standard_special_radii
from vorpy.src.chemistry.Standard import element_radii as standard_element_radii
from vorpy.src.chemistry.Voronota import special_radii as voronota_special_radii
from vorpy.src.chemistry.Voronota import element_radii as voronota_element_radii

__all__ = [
    # Chemical name mappings
    'residue_names',  # Maps common names to standard residue codes
    'element_names',  # Maps common names to element symbols
    'residue_atoms',  # Maps residue codes to their constituent atoms
    
    # Atomic masses
    'my_masses',  # Dictionary of atomic masses
    
    # Standard radii (from Standard.py)
    'standard_special_radii',  # Special radii for amino acids
    'standard_element_radii',  # Standard atomic radii
    
    # Voronoi radii (from Voronota.py)
    'voronota_special_radii',  # Special radii for Voronoi calculations
    'voronota_element_radii'   # Element radii for Voronoi calculations
]
