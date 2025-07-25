"""
VORPY - Voronoi analysis of molecular structures
"""

import os
import sys

# Add the parent directory to sys.path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Add everything in /api/ to the module search path.
__path__ = [os.path.dirname(__file__), os.path.join(os.path.dirname(__file__), "api")]

from .api import *
from .src.version import __version__
from .src.run import Run

# Make VorPyGUI directly accessible from vorpy module
from .src.GUI.vorpy_gui import VorPyGUI


def run(file=None, load_files=None, settings=None, groups=None, exports=None):
    """Run the VorPy GUI application."""
    app = Run(file=file, load_files=load_files, settings=settings, groups=groups, exports=exports)


# Don't pollute namespace.
del os, sys
