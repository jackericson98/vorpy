"""
VORPY API - A Python package for Voronoi analysis of molecular structures
"""

from vorpy.api.calculations import *
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
from vorpy.src.version import __version__

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
