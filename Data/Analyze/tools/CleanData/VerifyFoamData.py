import os
import sys

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))
sys.path.append(project_root)

import tkinter
from tkinter import filedialog

#
# def verify_foam_data(file=None):
#     if file is None:
