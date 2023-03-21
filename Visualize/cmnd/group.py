from System.Group.group import Group
from Visualize.cmnd.interpret import *


def group(sys, usr_npt, bff=None):
    # Create the object and index variables
    my_obj, my_ndx = None, None
    if usr_npt[0] == 'ns':
        my_obj = 'm'
        my_ndx = [0, len(sys.chains) - 2]
    # User only input "export"
    elif len(usr_npt) == 0:
        # Tell the user to pick an object and an index
        my_obj = get_obj(sys=sys)
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # User entered "export obj" and needs an index
    elif len(usr_npt) == 1:
        # Add a check for system
        # Check the object provided by the user
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # If the user input an object and an index of their own
    elif len(usr_npt) >= 2:
        # Check the object
        my_obj = get_obj(sys=sys, obj=usr_npt[0])
        my_ndx = get_ndx(sys=sys, obj=my_obj, ndx_npt=usr_npt[1])
    # Get the group information
    obj_ndx = ['m', 'r', 'a', 'n'].index(my_obj)
    obj_list = [sys.chains, sys.residues, sys.atoms, sys.ndxs][obj_ndx]
    name_prfx = ['mol', 'resid', 'atom', 'ndx'][obj_ndx]
    my_list, name = None, None
    # Get the slice and name of the group
    if len(my_ndx) == 1 and my_ndx[0] < len(obj_list):
        my_list = [obj_list[my_ndx[0]]]
        name = name_prfx + '_' + str(my_ndx[0])
    elif len(my_ndx) <= 2:
        my_list = obj_list[max(0, my_ndx[0]):min(len(obj_list), my_ndx[1] + 1)]
        name = name_prfx + 's_' + str(my_ndx[0]) + '_' + str(my_ndx[1])
    # Create the group
    npt_list = [None] * 4
    npt_list[obj_ndx] = my_list
    my_group = Group(sys=sys, chains=npt_list[0], residues=npt_list[1], atoms=npt_list[2], indices=npt_list[3], name=name, bff=bff)
    return my_group
