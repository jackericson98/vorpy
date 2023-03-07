from System.sys_objs.group import Group
from Visualize.commands.interpret import *


def group(sys, usr_npt, group_2=None):
    # Create the object and index variables
    my_obj, my_ndx = None, None
    # User only input "export"
    if len(usr_npt) <= 1:
        # Tell the user to pick an object and an index
        my_obj = get_obj(sys=sys)
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # User entered "export obj" and needs an index
    elif len(usr_npt) == 2:
        # Add a check for system
        # Check the object provided by the user
        my_obj = get_obj(sys=sys, obj=usr_npt[1])
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # If the user input an object and an index of their own
    elif len(usr_npt) >= 3:
        # Check the object
        my_obj = get_obj(sys=sys, obj=usr_npt[1])
        my_ndx = get_ndx(sys=sys, obj=my_obj, ndx_npt=usr_npt[2])
    # Get the group information
    obj_ndx = ['m', 'r', 'a', 'n'].index(my_obj)
    obj_list = [sys.mols, sys.residues, sys.atoms, sys.ndxs][obj_ndx]
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
    my_group = Group(sys=sys, mols=npt_list[0], residues=npt_list[1], atoms=npt_list[2], indices=npt_list[3], name=name)
    if group_2 is not None:
        my_group.add_atoms(group_2.atoms)
    return my_group