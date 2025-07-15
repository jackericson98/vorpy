import os
import sys

# Get the path to the root vorpy folder
vorpy_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# Add the root vorpy folder to the system path
sys.path.append(vorpy_root)

# Import the run function
from vorpy.src.command import argv, vorpy
from vorpy.src.system import System

if __name__ == "__main__":
    # Run the program with the command line arguments
    my_sys = System()
    # Run the program
    if len(sys.argv) > 1:
        argv(my_sys)
    else:
        my_sys.print_actions = True
        vorpy(my_sys=my_sys)
