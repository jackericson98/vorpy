"""
VORPY API - A Python package for Voronoi analysis of molecular structures
"""

from vorpy.api import calculations as calculations
from vorpy.api import inputs as inputs
from vorpy.api import interface as interface
from vorpy.api import output as output
from vorpy.api import system as system
from vorpy.api import chemistry as chemistry
from vorpy.api import command as command
from vorpy.api import group as group
from vorpy.api import network as network
from vorpy.api import GUI as GUI
from vorpy.api import objects as objects
from vorpy.api import visualize as visualize
from vorpy.api.GUI import VorPyGUI

__version__ = "1.0.3"

# Make everything available when importing from api
__all__ = [
    'calculations',
    'inputs',
    'interface',
    'output',
    'system',
    'chemistry',
    'command',
    'group',
    'network',
    'GUI',
    'objects',
    'visualize',
    '__version__',
    'VorPyGUI'
]
