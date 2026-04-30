import os
import sys
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog


vorpy_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..")
)
sys.path.append(vorpy_root)


from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2


"""
The basic goal of this script is to use existing logs files to calculate the surface energy of specified surfaces, so
that they may be compared to one another for active and non-active sites in complexes to determine if this may be a 
viable application of curvature for binding quality. 
"""


def surface_tension():

    """
    Since our solutes are in a water solvent, we can calculate the surface tension of a single surface by imagining it
    is the water that makes the surface. This means that by using water surface tension formulae and estimating the
    equivalent surface  
    """
    pass


def molar_free_surface_energy(surface_tension, molar_volume):
    """
    Meyers (2003) - Calculation of Molar Free Surface Energy.
    - Free energy cost to add more surface energy.
    - I would maybe compare the surface energy to a sphere of the same volume (similar to sphericity)
    """
    return 0.00845 * (molar_volume ** (2/3)) * surface_tension


def stokes_theorem():
    pass


def basic_surface_energy():
    pass

