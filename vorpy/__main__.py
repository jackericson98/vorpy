"""
Main entry point for running vorpy as a module
"""

import os
import sys

# Add the parent directory to sys.path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)


from vorpy.src.GUI.vorpy_gui import VorPyGUI as run
from vorpy.src.command.argv import argv

if __name__ == "__main__":
    if len(sys.argv) == 1:
        app = run()
        app.mainloop()
    else:
        argv()        
