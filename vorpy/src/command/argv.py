import sys
import os
from copy import deepcopy
import tkinter as tk
from tkinter import filedialog
from vorpy.src.command.load import load
from vorpy.src.command.set2 import sett
from vorpy.src.command.group import ggroup
from vorpy.src.system import System
from vorpy.src.output import export_micro
from vorpy.src.output import export_tiny
from vorpy.src.output import export_med
from vorpy.src.output import export_large
from vorpy.src.output import export_all
from vorpy.src.output import other_exports
from vorpy.src.output import set_system_directory
from vorpy.src.command.commands import *



"""
Argv rules: 
1. Space delimited
2. flags (-l: load, -s: set, -g: group, -e: export)
3. For multiple inputs use &&
4. Defaults to no sol, default settings, export all
"""


def argv_export(my_sys, usr_npt, add_on=None):

    """
    Exports the specified elements from the system based on the provided export type.

    Parameters:
    -----------
    my_sys : System
        The system containing the data to export
    usr_npt : list
        List of export specifications and parameters
    add_on : str, optional
        Additional directory path to append to export directory

    Export Types:
    ------------
    'default' - Exports network, information, no Sol shell, pdb, set atoms from the system and atoms,
                surrounding atoms, shell, and info from each group
    'all' - Exports everything from the full system and each group (requires confirmation)
    'info' - Exports only the information files for the system and each group
    'surfs' - Exports the built surfaces individually from the system and all group surfaces + verts and edges
    'dir' or 'directory' - Sets the export directory path
    """

    # Go through each of the inputs in the exports
    if len(usr_npt) == 0:
        usr_npt.append(['default'])
    # Export the specified exports
    for npt in usr_npt:
        if npt[0].lower() in {'dir', 'directory'} and len(npt) == 2 and (os.path.isdir(npt[1]) or npt[1] == 'gsu_logs'):
            if npt[1] == 'gsu_logs':
                set_system_directory(my_sys, "C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/{}/logs".format(my_sys.name))
            else:
                if add_on is None:
                    my_sys.dir = npt[1]
                else:
                    my_sys.dir = npt[1] + add_on
    for npt in usr_npt:
        if npt[0].lower() in {'dir', 'directory'}:
            continue
        export_npt(my_sys, npt[0])


def export_npt(my_sys, usr_npt=None):
    """
    Handles the export of system data based on user-specified export type.

    Parameters:
    -----------
    my_sys : System
        The system object containing the data to be exported
    usr_npt : str, optional
        The export type specification. If None or 'default', performs a medium export.
        Valid options include:
        - 'default'/'2'/'medium'/'med': Medium export (default)
        - 'tiny'/'i'/'info'/'0'/'smallest': Small export
        - 'small'/'s'/'1': Medium-small export
        - 'large'/'l'/'3': Large export
        - 'all'/'a'/'everything': Full export
        - Other custom export types

    Returns:
    --------
    None
        The function performs exports but does not return any values.
    """
    # If nothing is specified export the defaults
    if usr_npt is None or usr_npt.lower() in {'default', '2', 'medium', '', 'med'}:
        export_med(sys=my_sys)

    # Small export
    elif usr_npt.lower() in {"tiny", "i", "info", "0", "smallest"}:
        export_micro(my_sys)

    # Medium small export
    elif usr_npt.lower() in {"small", "s", "1"}:
        export_tiny(my_sys)

    # Large Export
    elif usr_npt.lower() in {"large", "l", "3"}:
        export_large(my_sys)

    # Export all
    elif usr_npt.lower() in {'all', 'a', 'everything'}:
        export_all(my_sys)

    else:
        other_exports(my_sys, usr_npt)


def interpret_argvs(counter=0):
    """
    Interprets command line arguments and organizes them into a structured dictionary.

    This function processes command line arguments starting from a specified counter position,
    organizing them into different command categories based on their flags. It handles various
    command types including loading (-l), settings (-s), grouping (-g), building (-b),
    exporting (-e), and interface (-i) commands.

    Parameters:
    -----------
    counter : int, optional
        The starting position in sys.argv to begin processing arguments. Default is 0.

    Returns:
    --------
    dict
        A dictionary containing organized commands with the following structure:
        {
            'npt': list of load commands,
            'set': list of setting commands,
            'grp': dict of group commands (indexed by group number),
            'bld': list of build commands,
            'xpt': list of export commands,
            'ifc': list of interface commands
        }
    """
    # Separate the rest of the argv args
    my_args = sys.argv[2 + counter:]
    # Set up the commands dictionary
    cmnds = {'npt': [], 'set': [], 'grp': {}, 'bld': [], 'xpt': [], 'ifc': []}
    # Set the arg to load as a default
    arg = '-l'
    group_counter = -1
    # Go through the arguments
    while my_args:
        # Remove the first argument flag
        if my_args[0] in ands:
            # If the argument is a flag, remove it
            my_args.pop(0)
        else:
            # If the argument is not a flag, set it as the current argument
            arg = my_args.pop(0)
            # If the argument is a group flag, increment the group counter
            if arg == '-g':
                group_counter += 1
                # Initialize the group command list
                cmnds['grp'][group_counter] = []
        # Gather the cmnd and the flag
        arg_cmnds = []
        # Keep gathering the commands for the flag
        while True:
            # If the argument is a flag or the end of the list, break
            if len(my_args) == 0 or my_args[0][0] == '-' or my_args[0] in ands:
                break
            else:
                # Keep gathering the commands for the flag
                arg_cmnds.append(my_args.pop(0))
        # Add the command to the list
        if arg.lower() == '-l':
            # Add the load command to the list
            cmnds['npt'].append(arg_cmnds)
        elif arg.lower() == '-s':
            # Add the setting command to the list
            cmnds['set'].append(arg_cmnds)
        elif arg.lower() == '-g':
            # Add the group command to the list
            cmnds['grp'][group_counter].append(arg_cmnds)
        elif arg.lower() == '-b':
            # Add the build command to the list
            cmnds['bld'].append(arg_cmnds)
        elif arg.lower() == '-e':
            # If the argument is 'logs', add the build type and logs command
            if arg_cmnds == 'logs':
                cmnds['set'].append(['bt', 'logs'])
            # If the argument is a directory command, format it
            if arg_cmnds[0] == 'dir':
                arg_cmnds = ['dir', " ".join(arg_cmnds[1:])]
            # Add the export command to the list
            cmnds['xpt'].append(arg_cmnds)
        elif arg.lower() == '-i':
            # Add the interface command to the list
            cmnds['ifc'].append(arg_cmnds)
    # Return the lists
    return cmnds


def argv(my_sys=None):
    """
    Processes command line arguments for the VORPY program.

    This function handles the command line interface for VORPY, allowing users to:
    - Load molecular structure files
    - Set various parameters and settings
    - Build and export surfaces
    - Group and manipulate molecular components
    - Configure interface calculations

    The function supports multiple command flags:
    - -l: Load files
    - -s: Set parameters
    - -g: Group operations
    - -b: Build surfaces
    - -e: Export results
    - -i: Interface calculations

    Parameters:
    -----------
    my_sys : System
        The main system object that stores molecular data and settings

    Returns:
    --------
    None
    """
    if my_sys is None:
        my_sys = System()
    # First check if the argv value is in helps
    if sys.argv[1].lower() in helps:
        print_help()
        return
    # Load the atom file
    if len(sys.argv) == 2:
        # Load the atom file
        load(my_sys, [["", sys.argv[1]]], balls_file=True)
        counter = 0
    # If the second argument is a flag, load the atom file
    elif sys.argv[2][0] == '-':
        load(my_sys, [["", sys.argv[1]]], balls_file=True)
        counter = 0
    else:
        # If the second argument is not a flag, load the atom file
        counter = 1
        while '-' not in sys.argv[1 + counter]:
            counter += 1
        # Join the arguments into a single string
        my_file = " ".join(sys.argv[1:1 + counter])
        # Load the atom file
        load(my_sys, [["", my_file]], balls_file=True)
    # Interpret the commands
    cmnds = interpret_argvs(counter)
    # Set the system commands
    my_sys.cmnds = cmnds
    # Go through each of the ls
    load(my_sys, cmnds['npt'])
    # Go through each of the export commands
    for commandaroonski in cmnds['xpt']:
        # If the command is a directory, check if it exists
        if commandaroonski[0] == 'dir':
            if os.path.isdir(commandaroonski[1]):
                my_sys.files['dir'] = commandaroonski[1]
            # If the command is a browse, open a dialog to choose a folder
            elif commandaroonski[1].lower() in browse_names:
                # Create a new Tkinter root window
                my_root = tk.Tk()
                # Hide the root window
                my_root.withdraw()
                # Make the root window the topmost window
                my_root.wm_attributes('-topmost', 1)
                # Open a dialog to choose a folder
                folder = filedialog.askdirectory(title='Choose Output Folder')
                # If the folder exists, set the directory
                if os.path.exists(folder):
                    my_sys.files['dir'] = folder
                else:
                    print(f"{folder} is not a valid folder")
                # Remove the command from the list
            cmnds['xpt'].pop(cmnds['xpt'].index(commandaroonski))

    # Declare the settings variable
    settings = None
    # Go through the user inputs loading files
    for my_set in cmnds['set']:
        # Alter the settings
        settings = sett(my_set[0], my_set[1:], settings)
    # Update the sphere radii in the system
    if settings is not None and settings['atom_rad'] is not None:
        my_sys.set_radii(settings['atom_rad']['element'], settings['atom_rad']['special'])

    # Compare the groups
    if my_sys.groups is None or len(my_sys.groups) == 0:
        ggroup(my_sys, cmnds['grp'], settings)
    else:
        # Compare the groups
        verts = my_sys.groups[0].verts
        # Delete the groups 
        my_sys.groups = None
        # Recalculate the groups
        ggroup(my_sys, cmnds['grp'], settings)
        # Reset the vertices
        my_sys.groups[0].verts = verts
    # If the groups are not valid, calculate the whole molecule
    if my_sys.groups is None or len(my_sys.groups) == 0:
        # Print the error message
        print('{} not a valid group command. Calculating whole molecule'.format(cmnds['grp']))
        # Calculate the whole molecule
        ggroup(my_sys, [['ns']])
    # Build the groups
    if settings is not None and len(settings['net_type']) > 1 and settings['net_type'][0] == 'com':
        # Create a new list of groups
        new_groups = []
        # Go through each of the groups
        for grp in my_sys.groups:
            # Copy the group
            copy_group = deepcopy(grp)
            # Change the name of the group
            copy_group.name = copy_group.name + '_' + settings['net_type'][1]
            # Change the net type of the group
            copy_group.settings['net_type'] = settings['net_type'][1]
            # Change the name of the group
            grp.settings['net_type'] = settings['net_type'][2]
            grp.name = grp.name + '_' + settings['net_type'][2]
            # Delete the vertices from the grp group because they are only available in the original group
            grp.verts = None
            # Build the groups
            copy_group.build()
            grp.build()
            # Add the new groups to the list
            new_groups.append(copy_group)
            # Compare the two networks
            my_sys.compare_networks(group1=copy_group, group2=grp)
        # Add the new groups to the list
        my_sys.groups += new_groups

    else:
        # Go through each of the groups
        for grp in my_sys.groups:
            # Build the groups
            grp.build()
    # Make the system's interfaces
    my_sys.make_interfaces()
    # Export everything
    argv_export(my_sys, cmnds['xpt'])


if __name__ == '__main__':
    mySys = System()
    argv(mySys)
