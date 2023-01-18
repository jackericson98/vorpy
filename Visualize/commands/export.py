from Visualize.commands.interpret import *
from System.sys_objs.group import Group



def export(sys, usr_npt):
    """
    Takes in input strings and exports them based on their option choices
    :param sys:
    :param usr_npt:
    :return:
    """
    my_obj, my_ndx = None, None
    # User only input "export"
    if len(usr_npt) <= 1:
        # Tell the user to pick an object and an index
        my_obj = get_obj(sys=sys)
        my_ndx = get_ndx(sys=sys, obj=my_obj)
    # User entered "export obj" and needs an index
    elif len(usr_npt) == 2:
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
    if len(my_ndx) == 1:
        my_list = [obj_list[my_ndx[0]]]
        name = name_prfx + '_' + str(my_ndx[0])
    elif len(my_ndx) <= 2:
        my_list = obj_list[my_ndx[0]:my_ndx[1]]
        name = name_prfx + 's_' + str(my_ndx[0]) + '_' + str(my_ndx[1])
    # Create the group
    npt_list = [None]*4
    npt_list[obj_ndx] = my_list
    my_group = Group(sys=sys, mols=npt_list[0], residues=npt_list[1], atoms=npt_list[2], indices=npt_list[3], name=name)
    while True:
        # Export the group exports
        xpt_npt = input("choose one of the following to export (or type \'q\' to quit): 1. Shell, 2. Surfaces, 3. Layers, 4. Atoms, 5. Filled Body, 6. Info File\nexport >>>   ")
        # Check for a quit
        if xpt_npt.lower() in quits:
            return
        # Check for help request
        elif xpt_npt.lower() in helps:
            help_()
            continue
        # Export the shell:
        elif xpt_npt.lower() in ['1', '1.', 'shell', 'sh']:
            my_group.exports(shell=True)
            print("\r{} shell exported to {}".format(my_group.name, my_group.dir))
        # Export the Surfaces
        elif xpt_npt.lower() in ['2', '2.', 'surfs', 'surfaces']:
            my_group.exports(surfaces=True)
            print("\r{} surfaces exported to {}".format(my_group.name, my_group.dir + "/surfaces"))
        # Export the layers
        elif xpt_npt.lower() in ['3', '3.', 'layers', 'lyrs', 'l']:
            my_group.exports(layers=True)
            print("\r{} layers exported to {}".format(my_group.name, my_group.dir + "/layers"))
        # Export the atoms
        elif xpt_npt.lower() in ['4', '4.', 'atoms', 'a', 'atms']:
            my_group.exports(atoms=True)
            print("\r{} atoms exported to {}".format(my_group.name, my_group.dir))
        # Export the filled body
        elif xpt_npt.lower() in ['5', '5.', 'filled body', 'fb', 'f', 'filled_body', 'f_b', 'fld_bdy']:
            my_group.exports(fill=True)
            print("\r{} filled body exported to {}".format(my_group.name, my_group.dir))
        # Export the Info file
        elif xpt_npt.lower() in ['6', '6.', 'info', 'info_file', 'info file']:
            my_group.exports(info=True)
            print("\r{} info file exported to {}".format(my_group.name, my_group.dir))
        else:
            invalid_input(xpt_npt)
