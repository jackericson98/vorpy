from Visualize.commands.commands import *
from Visualize.commands.show import *
import os
from os import path


def get_ndx(sys, ndx=None, list_len=None, obj=None):
    """
    Asks the user for the index of the object they specified
    :return:
    """
    obj_name = None
    names = ["molecule", "residue", "atom", "index"]
    # Get the object and make sure it's real
    if obj is not None and obj in my_objects:
        for i in range(4):
            if obj in [mol_objs, res_objs, atom_objs, ndx_objs][i]:
                obj =  i + 1
                obj_name = names[i]
    elif type(obj) is int and obj <= 4:
        obj_name = names[obj - 1]
    else:
        obj = get_obj(sys=sys, obj=obj)
        obj_name = names[obj - 1]
    # Start asking for the index of the object
    asking = True
    while asking:
        max_ndx = len([sys.mols, sys.residues, sys.atoms, sys.groups][obj - 1])
        if max_ndx == 0:
            print("No {} to choose from. Choose another object or type \'h\' for help")
        elif ndx is None:
            extra = ""
            if list_len is not None:
                extra = " less than {}".format(list_len)
            prompt_str = "Enter {} index{} between 0 and {}\nindex >>>   ".format(obj_name, extra, max_ndx)
            ndx = input(prompt_str)
        if ndx.lower() in quits or ndx.lower() in dones:
            return
        elif ndx.lower() in helps:
            help_()
            continue
        elif ndx.lower() in show_cmds:
            show(obj_name)
            continue
        try:
            ndx = int(ndx)
            asking = False
        except ValueError:
            invalid_input(ndx)
            ndx = None
        if not asking and list_len is not None and ndx > list_len:
            invalid_input(ndx)
            asking = True
    return ndx


def get_obj(sys, obj=None, return_ndx=True):
    """
    Makes the user type a proper object
    :return: 1-4 based on if it is a 1. molecule 2. residue 3. atom or 4. index
    """
    my_input, choosing = obj, False
    # If obj not in my objects
    if obj is None or (type(obj) is str and obj.lower() not in my_objects) or (type(obj) is int and obj > 4):
        choosing = True
    # Keep asking the user to choose an object to export
    while choosing:
        # Prompt the user
        my_input = input("Enter an object type. (\'mol\', \'res\', \'atom\', or \'ndx\')\nobject >>>   ")
        # Check to see if the user gave a valid response or not
        if my_input.lower() in quits:
            return
        elif my_input.lower() in helps:
            help_()
        elif my_input.lower() not in my_objects:
            # Tell the user they suck and try again
            invalid_input(my_input)
            continue
        # Otherwise, we have a success
        else:
            choosing = False
    if return_ndx:
        # If the input is already an integer return it
        if type(my_input) is int and my_input <= 4:
            return my_input
        # Go through and find the type of object we are getting
        objs = [mol_objs, res_objs, atom_objs, ndx_objs]
        for i in range(4):
            if my_input.lower() in objs[i]:
                if len([sys.mols, sys.residues, sys.atoms, sys.ndxs][i]) > 0:
                    return i + 1
                else:
                    print("No {} in the system. Try again or typ \'h\' for help"
                          .format(["molecules", "residues", "atoms", "groups"][i]))

    # As a failsafe
    return my_input


def get_file(file=None):
    # Check if there is a file provided
    if file is None:
        print("Enter a file address. (Use \'./\' to load a file from the \'.../vorpy\' directory):")
    # Check the file
    checking_file = True
    while checking_file:
        # Get the file if None was specified
        if file is None:
            file = input("file >>>   ")
            if file.lower() in quits:
                return
            elif len(file) == 0:
                checking_file = True
                continue
            elif file.lower() in helps:
                help_()
            test_file = file.split()
            if test_file[0] in load_cmds:
                file = file[len(test_file[0]) + 1:]
        # Check if the initial file works
        if path.exists(file) and len(file) > 0:
            checking_file = False
        # Check if the file is in the ./Data/test_data folder
        elif path.exists("./Data/test_data/" + file) and len(file) > 0:
            file = os.getcwd() + "/Data/test_data/" + file
            checking_file = False
        # Check if it is just the raw name
        elif path.exists("./Data/test_data/" + file + ".pdb") and len(file) > 0:
            file = os.getcwd() + "/Data/test_data/" + file + ".pdb"
            checking_file = False
        # Otherwise, tell the user to try again
        else:
            invalid_input(file)
            file = None
            continue
    # Return the file
    return file


def get_set(setting=None, val=None):
    """
    Makes the user type a proper Value
    :return: 1-4 based on if it is a 1. molecule 2. residue 3. atom or 4. index
    """
    my_input, choosing = val, False
    # If obj not in my objects
    if val is None or val.lower() not in my_settings:
        choosing = True

    # Keep asking the user to choose an object to export
    while choosing:
        my_input = setting
        if setting is None:
            # Prompt the user
            my_input = input("Enter setting type. (\'surf_res\', \'max_vert\', \'box_size\', or \'sol_verts\')\nsetting >>>   ")
        # If they quit, then quit
        if my_input.lower() in quits:
            return
        elif my_input.lower() in helps:
            help_()
        # Check to see if the user gave a valid response or not
        if my_input.lower() not in my_settings:
            # Tell the user they suck and try again
            invalid_input(my_input)
            continue
        # Otherwise, we have a success
        else:
            choosing = False
    # As a failsafe
    return my_input


def get_val(setting=None, val=None):
    """
        Asks the user for the index of the object they specified
        :return:
        """
    # If the setting is not in the settingsd
    if setting is None or setting.lower() not in my_settings:
        setting = get_set()
    # Find the value for the setting
    asking = True
    while asking:
        # If no val has been provided
        if val is None:
            prompt_str = "Enter \'{}\' value \nvalue >>>   ".format(setting)
            val = input(prompt_str)

        if val.lower() in quits or val.lower() in dones:
            return
        elif val.lower() in helps:
            help_()
        if setting in sol_vertses:
            if val.lower() in ['t', 'true', 'tr'] + ys:
                val = True
            elif val.lower() in ['f', 'false', 'flse', 'fl', 'fa', 'fs', 'fls'] + ns:
                val = False
        else:
            try:
                val = float(val)
            except ValueError:
                val = None
        if val is not None:
            asking = False
    return val
