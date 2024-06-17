from Visualize.cmnd.commands import *
import os
from os import path


def get_ndx(sys, obj, ndx_npt=None):
    """
    Asks the user for the index of the object they specified
    :return:
    """
    # Create the naming dictionary
    name_dict = {'m': 'Molecule', 'r': 'Residue', 'a': 'Atom', 'n': 'Index'}
    # Get the name
    name = name_dict[obj]
    # Get the list
    obj_list = []
    if obj == 'm':
        obj_list, obj_num = sys.chains, 0
    elif obj == 'r':
        obj_list, obj_num = sys.residues, 1
    elif obj == 'a':
        obj_list, obj_num = sys.atoms, 2
    elif obj == 'n':
        obj_list, obj_num = sys.ndxs, 3
    # Start the ndx checking loop
    while True:
        # Get the index if the index doesn't exist
        if ndx_npt is None:
            ndx_npt = input("enter a {} index (range: 0 - {})\nindex >>>   ".format(name, len(obj_list) - 1))
        # Check for quits
        if ndx_npt.lower() in quits:
            return
        # Check for helps
        elif ndx_npt.lower in helps:
            help_()
            continue
        # Check the input
        ndx_npt = ndx_npt.split("-")
        # Get the list of index numbers
        try:
            return [int(_) for _ in ndx_npt]
        except ValueError:
            continue


def get_obj(sys, obj=None):
    """
    Makes the user type a proper object
    :return: 1-4 based on if it is a 1. molecule 2. residue 3. atom or 4. index
    """
    # Keep asking the user to choose an object to export
    while True:
        # If no input was given
        if obj is None:
            # Prompt the user
            obj = input("enter an object type. (\'mol\', \'res\', \'atom\', or \'ndx\')\nobject >>>   ")
            obj = obj.split()
        # Check to see if the user gave a valid response or not
        if obj[0].lower() in quits:
            return
        elif obj[0].lower() in helps:
            help_()
        elif obj[0].lower() not in my_objects:
            # Tell the user they suck and try again
            invalid_input(obj)
            obj = None
            continue
        # Otherwise, we have a success
        elif obj[0].lower() in chn_objs:
            return 'm'
        elif obj[0].lower() in res_objs:
            return 'r'
        elif obj[0].lower() in atom_objs:
            return 'a'
        elif obj[0].lower() in ndx_objs:
            return 'n'


def get_file(file=None):
    # Check if there is a file provided
    if file is None:
        print("enter a file address. (Use \'./\' to load a file from the \'.../vorpy\' directory):")
    # Check the file
    while True:
        # Get the file if None was specified
        if file is None:
            file = input("file >>>   ")
            if file.lower() in quits:
                return
            elif len(file) == 0:
                continue
            elif file.lower() in helps:
                help_()
            test_file = file.split()
            if test_file[0] in load_cmds:
                file = file[len(test_file[0]) + 1:]
        # Check if the initial file works
        if path.exists(file) and len(file) > 0:
            break
        # Check if the file is in the ./Data/test_data folder
        elif path.exists("./Data/test_data/" + file) and len(file) > 0:
            file = os.getcwd() + "/Data/test_data/" + file
            break
        # Check if it is just the raw name
        elif path.exists("./Data/test_data/" + file + ".pdb") and len(file) > 0:
            file = os.getcwd() + "/Data/test_data/" + file + ".pdb"
            break
        # Otherwise, tell the user to try again
        else:
            invalid_input(file)
            file = None
            continue
    # Return the file
    return file


def get_set(usr_npt=None):
    """
    Makes the user type a proper Value
    :return: base versions of the settings ('sr', 'mv', 'bm', 'bs', 'nt')
    """
    # Keep asking the user to choose an object to export
    while True:
        # Check if we need to start from scratch
        if usr_npt is None:
            # Prompt the user
            usr_npt = input("setting >>>   ")
        # If they quit, then quit
        if usr_npt.lower() in quits:
            return
        elif usr_npt.lower() in helps:
            help_()
        # Check to see if the user gave a valid response or not
        if usr_npt.lower() in surf_reses:
            # Return the base setting
            return 'sr'
        elif usr_npt.lower() in surf_colors:
            # Return the base setting
            return 'sc'
        elif usr_npt.lower() in surf_schemes:
            # Return the base setting
            return 'ss'
        elif usr_npt.lower() in max_verts:
            # Return the base setting
            return 'mv'
        elif usr_npt.lower() in box_sizes:
            # Return the base setting
            return 'bm'
        elif usr_npt.lower() in build_surfses:
            # Return the base setting
            return 'bs'
        elif usr_npt.lower() in net_types:
            # Return the base setting
            return 'nt'
        elif usr_npt.lower() in atom_radii:
            # Return the base setting
            return 'ar'
        else:
            # Tell the user they suck and try again
            print("\"{}\" is not a valid input. Enter a 181L value (\'surf_res\', \'max_vert\', \'box_size\', or "
                  "\'calc_surfs\')".format(usr_npt))
            usr_npt = None


def get_val(my_sys, setting=None, val=None):
    """
        Asks the user for the index of the object they specified
        :return:
        """
    # If the setting is not in the settingsd
    if setting is None or setting.lower() not in my_settings:
        setting = get_set()
    # Find the value for the setting
    while True:
        # If no val has been provided
        val1 = None
        if val is None:
            prompt_str = "Enter {} value \nvalue >>>   ".format(settings_dict[setting])
            val = input(prompt_str)
        if type(val) == list:
            if len(val) == 2:
                val, val1 = val
            else:
                val = val[0]
        # Quit if asked
        if val.lower() in quits or val.lower() in dones:
            return
        if setting in build_surfses:
            if val.lower() in ['t', 'true', 'tr'] + ys:
                val = True
            elif val.lower() in ['f', 'false', 'flse', 'fl', 'fa', 'fs', 'fls'] + ns:
                val = False
        # Test the validity of the user's true and false skills
        if setting in atom_radii:
            if val not in my_sys.radii:
                val = None
            else:
                try:
                    val1 = float(val1)
                except IndexError:
                    val.append(None)
                    while True:
                        val1 = input("enter new radius for {} (current radius = {})".format(val[0].upper(),
                                                                                            my_sys.radii[val[0]]))
                        try:
                            val1 = float(val1)
                            break
                        except ValueError:
                            val = None
                except ValueError:
                    val = None
            val = [val, val1]
        elif type(val) == str and setting in net_types:
            if val.lower() in trues + ys:
                val = 'del'
            elif val.lower() in falses + ns:
                val = 'vor'
            elif val.lower() in voronoi_vals:
                val = 'vor'
            elif val.lower() in delaunay_vals:
                val = 'del'
            elif val.lower() in power_vals:
                val = 'pow'
            elif val.lower() in compare_vals:
                val = 'com'
        # Test for a float value
        elif type(val) == str and setting in surf_reses + max_verts + box_sizes:
            try:
                val = float(val)
            except ValueError:
                val = None
        elif type(val) == str and setting in surf_colors:
            if val not in ["viridis", "plasma", "inferno", "cividis", "Greys", "Reds", "Greens", "Blues"]:
                val = None
        elif type(val) == str and setting in surf_schemes:
            if val not in ['dist', 'ins_out', 'curv']:
                val = None
        # Check if we cool
        if val is not None:
            break
    return val
