from Visualize.cmnd.load import *
from Visualize.argv.vpy_argv import argv
from Visualize.cmnd.vpy_cmnd import vorpy
import sys

# Main run
if __name__ == '__main__':
    # Welcome introduction
    my_sys = System()
    # Check to see if the user input argvs
    if len(sys.argv) > 1:
        argv(my_sys)
    else:
        # Run vorpy
        my_sys.print_actions = True
        vorpy(my_sys=my_sys)
