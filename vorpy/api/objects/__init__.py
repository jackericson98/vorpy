"""
Objects module for molecular structure representation.
"""

from vorpy.src.objects.atom import Atom
from vorpy.src.objects.atom import make_atom
from vorpy.src.objects.atom import get_element
from vorpy.src.objects.atom import get_radius
from vorpy.src.objects.chain import Chain
from vorpy.src.objects.chain import Sol
from vorpy.src.objects.residue import Residue
from vorpy.src.objects.interface import Interface

__all__ = [
    # Atom-related
    'Atom',
    'make_atom',
    'get_element',
    'get_radius',
    
    # Chain-related
    'Chain',
    'Sol',
    
    # Residue-related
    'Residue',
    
    # Interface-related
    'Interface'
]

