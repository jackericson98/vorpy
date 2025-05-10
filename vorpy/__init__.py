"""
VOROPY - Voronoi analysis of molecular structures
"""

import os
import sys

# Add the parent directory to sys.path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

__version__ = "1.0.3"

# Add everything in /api/ to the module search path.
__path__ = [os.path.dirname(__file__)]
__path__.append(os.path.join(os.path.dirname(__file__), "api"))

from vorpy.api import *
from vorpy.api import __version__

# Make VorPyGUI directly accessible from vorpy module
from vorpy.src.GUI.vorpy_gui import VorPyGUI as run

# Don't pollute namespace.
del os, sys
