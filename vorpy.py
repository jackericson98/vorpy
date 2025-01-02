import os
from Visualize.cmnd.load import *
from Visualize.cmnd.argv import argv
from Visualize.cmnd.vpy_cmnd import vorpy
import sys
from System.system import System
from System.sys_funcs.input import *
from System.sys_funcs.calcs import *
from System.sys_funcs.output import *
from System.Group.group import Group
from System.Network.network import Network



# Main run
if __name__ == '__main__':
    # Welcome introduction
    my_sys = System(root_dir=os.getcwd())
    # Check to see if the user input argvs
    if len(sys.argv) > 1:
        argv(my_sys)
    else:
        # Run vorpy
        my_sys.print_actions = True
        vorpy(my_sys=my_sys)
