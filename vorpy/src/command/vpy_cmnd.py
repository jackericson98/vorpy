import sys
import tkinter as tk
from tkinter import filedialog
import os
from copy import deepcopy
from vorpy.src.inputs.net import read_net
from vorpy.src.command.commands import *
from vorpy.src.system.system import System
from vorpy.src.command.interpret import get_file
from vorpy.src.command.set import sett
from vorpy.src.command.group import ggroup
from vorpy.src.command.command_export import argv_export
from vorpy.src.command.interface import build_interfaces


class Command:
    def __init__(self, sys=None, settings=None):
        self.sys = sys
        self.base_file = None
        self.load_commands = []
        self.groups = {}
        self.builds = []
        self.exports = []
        self.interfaces = []
        self.settings_cmnds = []
        self.settings_dict = settings
        self.logs_files = []

    def run(self):
        self._run_pipeline()

    def _run_pipeline(self):
        """
        Starts the reading process. Pulls out the base file if not already loaded.
        """
        
        # Check for the first argv and whether it is a 
        if sys.argv[1][-3:] in {"pdb", "gro", "mol", "cif", 'txt'}:
            self.base_file = sys.argv[1]
        else:
            if get_file(sys.argv[1]) is not None:
                self.base_file = get_file(sys.argv[1])
            else:
                self.base_file = None
                # Print the invalid input and return us from the commands
                print("{} is not a valid input file".format(sys.argv[1]))
                return
        # Check to see if a system has been loaded or not
        if self.sys is None:
            # Define the system
            self.sys = System(file=self.base_file)
        
        # Interpret the rest of the argv
        self.parse_commands()

        # Load the other files
        self.load_files()

        # Check for whether the 
        self.apply_settings()

        # Compare the groups
        self.create_groups()
        # Build interface networks instead of complete group networks
        if len(self.interfaces) > 0:
            build_interfaces(
                sys=self.sys,
                num_requested_groups=len(self.groups),
            )

        # Compare network types
        elif (
                self.settings_dict is not None
                and len(self.settings_dict["net_type"]) > 1
                and self.settings_dict["net_type"][0] == "com"
        ):
            new_groups = []

            for grp in self.sys.groups:
                copy_group = deepcopy(grp)

                copy_group.name = (
                        copy_group.name + "_" + self.settings_dict["net_type"][1]
                )
                copy_group.settings["net_type"] = self.settings_dict["net_type"][1]

                grp.settings["net_type"] = self.settings_dict["net_type"][2]
                grp.name = grp.name + "_" + self.settings_dict["net_type"][2]

                grp.verts = None

                copy_group.build()
                grp.build()

                new_groups.append(copy_group)

                self.sys.compare_networks(
                    group1=copy_group,
                    group2=grp,
                )

            self.sys.groups += new_groups

        else:
            # Build the groups
            if self.settings_dict is not None and len(self.settings_dict['net_type']) > 1 and self.settings_dict['net_type'][0] == 'com':
                # Create a new list of groups
                new_groups = []
                # Go through each of the groups
                for grp in self.sys.groups:
                    # Copy the group
                    copy_group = deepcopy(grp)
                    # Change the name of the group
                    copy_group.name = copy_group.name + '_' + self.settings_dict['net_type'][1]
                    # Change the net type of the group
                    copy_group.settings['net_type'] = self.settings_dict['net_type'][1]
                    # Change the name of the group
                    grp.settings['net_type'] = self.settings_dict['net_type'][2]
                    grp.name = grp.name + '_' + self.settings_dict['net_type'][2]
                    # Delete the vertices from the grp group because they are only available in the original group
                    grp.verts = None
                    # Build the groups
                    copy_group.build()
                    grp.build()
                    # Add the new groups to the list
                    new_groups.append(copy_group)
                    # Compare the two networks
                    self.sys.compare_networks(group1=copy_group, group2=grp)
                # Add the new groups to the list
                self.sys.groups += new_groups

            else:

                build_type = None

                if self.settings_dict is not None:
                    build_type = self.settings_dict.get("bld_type", None)

                if build_type == "logs":

                    self.build_groups_from_logs()
                else:
                    for grp in self.sys.groups:
                        grp.build()

        # Make the system's interfaces
        if len(self.interfaces) == 0:
            self.sys.make_interfaces()

        # Export everything
        self.run_exports()

    def build_groups_from_logs(self):
        """
        Rebuilds command-line groups from existing vorpy logs instead of recomputing vertices.

        Expected command style:
            py vorpy system.pdb -b logs path/to/aw_logs.csv -e med

        If one logs file is provided, it is used for the first/default group.
        If multiple groups are provided, pass one logs file per group in the same order.
        """

        if len(self.logs_files) == 0:
            raise ValueError(
                "Build type was set to logs, but no logs file was provided.\n"
                "Example:\n"
                "  py vorpy system.pdb -b logs path/to/aw_logs.csv -e med"
            )

        if len(self.logs_files) == 1 and len(self.sys.groups) >= 1:
            logs_by_group = [self.logs_files[0]]

        elif len(self.logs_files) == len(self.sys.groups):
            logs_by_group = self.logs_files

        else:
            raise ValueError(
                "Number of logs files does not match number of groups.\n"
                f"groups = {len(self.sys.groups)}\n"
                f"logs files = {len(self.logs_files)}"
            )

        for grp, logs_file in zip(self.sys.groups, logs_by_group):
            if not os.path.exists(logs_file):
                raise FileNotFoundError(f"Logs file not found: {logs_file}")

            print(f"\nLoading group from logs:")
            print(f"  group = {grp.name}")
            print(f"  logs  = {logs_file}")

            read_net(
                group=grp,
                net=grp.net,
                file_name=logs_file,
                rebuild_edges=True,
                rebuild_surfs=True,
                analyze=True,
                store_points=True
            )

    def parse_commands(self, counter=0):
        """
        Splits the user inputs into the different commands and flags
        """
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
                    self.groups[group_counter] = []
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
                self.load_commands.append(arg_cmnds)
            elif arg.lower() == '-s':
                # Add the setting command to the list
                self.settings_cmnds.append(arg_cmnds)
            elif arg.lower() == '-g':
                # Add the group command to the list
                self.groups[group_counter].append(arg_cmnds)
            elif arg.lower() == '-b':
                self.builds.append(arg_cmnds)

            if len(arg_cmnds) > 0 and arg_cmnds[0].lower() in {"logs", "log"}:
                self.settings_cmnds.append(["bt", "logs"])

                if len(arg_cmnds) > 1:
                    self.logs_files.append(" ".join(arg_cmnds[1:]))
            elif arg.lower() == '-e':
                # If the argument is 'logs', add the build type and logs command
                if arg_cmnds == 'logs':
                    self.settings_cmnds.append(['bt', 'logs'])
                # If the argument is a directory command, format it
                if arg_cmnds[0] == 'dir':
                    # Check if the direcory is in the browse names
                    if arg_cmnds[1] in browse_names:
                        # Launch the browse window
                        my_root = tk.Tk()
                        my_root.withdraw()
                        my_root.wm_attributes('-topmost', 1)
                        folder = filedialog.askdirectory(title='Choose Output Folder')
                        if os.path.exists(folder):
                            self.sys.files['dir'] = folder
                        else:
                            print(f"{folder} is not a valid folder")
                    else:
                        base_dir = arg_cmnds[1]
                        system_dir = os.path.join(base_dir, self.sys.name)

                        print("Directory set to: {}".format(system_dir))
                        self.sys.set_output_directory(system_dir)
                else:
                    # Add the export command to the list
                    self.exports.append(arg_cmnds)
            elif arg.lower() == '-i':
                # Add the interface command to the list
                self.interfaces.append(arg_cmnds)

    def load_files(self):
        """
        Loads molecular structure files and associated data into the system.

        This function handles the loading of various file types:
        - Molecular structure files (.pdb, .mol, .gro, .cif)
        - Vertex files (.txt with 'verts' or 'vertices' in name)
        - Network files (.txt with 'net' in name)

        The function provides interactive confirmation prompts when:
        - Replacing an existing system
        - Replacing existing vertex files
        - Replacing existing network files

        Parameters:
        -----------
        sys : System
            The system object to load data into
        usr_npt : list
            List of file specifications to load
        balls_file : bool, optional
            Flag indicating if the file should be treated as a molecular structure file
            regardless of extension. Default is False.

        Returns:
        --------
        System or None
            Returns the updated system object if successful, None if loading is cancelled
        """

        # Process each file in the list
        for my_file in self.load_commands:
            # Interpret the file
            file = get_file(my_file)
            # Check to see what type of file it is
            if file[-3:] == 'pdb' or file[-3:] == 'mol' or file[-3:] == 'gro' or file[-3:] == 'cif':
                # If the system already exists, prompt the user to confirm replacement
                if self.sys.name is not None and \
                        (self.sys.atoms is not None or self.sys.files['verts_file'] is not None or self.sys.files['net_file'] is not None):
                    reset_sys = input("replacing {} with {}\nconfirm >>>   "
                                    .format(self.sys.name, file))
                    # If the user confirms the replacement, create a new system
                    if reset_sys.lower() in ys:
                        self.sys = System(file)
                        print(self.sys.name + " loaded - {} atoms, {} molecules, solute: {}"
                              .format(len(self.sys.atoms), len(self.sys.chains), self.sys.sol.name))
                        return self.sys
                    # If the user requests help, print the help message
                    elif reset_sys.lower() in helps:
                        print_help()
                    # If the user quits, return None
                    elif reset_sys.lower() in quits:
                        return
                # If the system does not exist, load the new system
                else:
                    self.sys.load_sys(file=file)
                    # noinspection PyUnresolvedReferences
                    self.sys.print_info()
                    return self.sys
            # If the loaded file is a vertex or network file load them accordingly
            elif file[-3:] == 'txt':
                # If the new file is a vertex file load it
                if file[-9:-4].lower() == 'verts' or file[-12:-4].lower() == 'vertices':
                    # If a vertex file has already been loaded make sure the user wants to load it if not load it
                    if self.sys.files['verts_file'] is not None and self.sys.vert_file != "":
                        replace_vert_file = input("replacing {} with {}\n "
                                                "confirm >>>   ".format(self.sys.files['verts_file'], file))
                        # If the user confirms the replacement, load the vertices
                        if replace_vert_file.lower() in ys or replace_vert_file.lower() in dones:
                            self.sys.load_verts(file, vta_ball_file=self.sys.ball_file)
                            print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                                .format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'], self.sys.net.settings['box_size']))
                        # If the user requests help, print the help message
                        elif replace_vert_file.lower() in helps:
                            print_help()
                        # If the user quits, return None
                        elif replace_vert_file.lower() in quits:
                            return
                    # If the vertex file has not been loaded, load it
                    else:
                        self.sys.load_verts(file, vta_ball_file=self.sys.files['ball_file'])
                        # print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                        #       .format(sys.name, len(sys.net.vta_verts), sys.net.max_vert, sys.net.box_size))
                elif file[-9:-4].lower() == 'balls':
                    self.sys.ball_file = file
                # If the new file is a network file load it
                elif file[-11:-4].lower() in 'network':
                    # If a vertex file has already been loaded make sure the user wants to load it if not load it
                    if self.sys.net_file is not None or self.sys.net_file != "":
                        replace_net_file = input("replacing {} with {}\n "
                                                "confirm >>>   ".format(self.sys.net_file, file))
                        # If the user confirms the replacement, load the network
                        if replace_net_file in ys:
                            self.sys.load_net(file)
                            print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box"
                                " size: {} x\n".format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'],
                                                        self.sys.net.settings['box_size']))
                        # If the user requests help, print the help message
                        elif replace_net_file in helps:
                            print_help()
                        # If the user quits, return None
                        else:
                            return
                    else:
                        # Load the file
                        self.sys.load_net(file)
                        if len(sys.net.surfs) > 0:
                            print("{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box size: {} x\n"
                                .format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'], self.sys.net.settings['box_size']))
                        else:
                            print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                                .format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'], self.sys.net.settings['box_size']))
            # Check to see if it is a new network file
            elif file[-3:] == 'csv':
                # Check to see that this is a network file
                if file[-7:-4].lower() == 'net':

                    self.sys.load_net(file=file)

            # If the file is an index file load it accordingly
            elif file[-3:] == 'ndx':
                self.sys.load_ndx(file)
                print(self.sys.ndx_file + "loaded -  {}".format(self.sys.ndx_names[:min(len(self.sys.ndx_names) - 1, 10)]))
            # In all other case print an error and give the user a chance to try again
            else:
                print("\'{}\' is not a valid input. allowed file types: .pdb, .mol, .cif, .gro, .txt, .ndx. type "
                    "\'h\' for help".format(file))
                return

    def apply_settings(self):
        # Go through the user inputs loading files
        for my_set in self.settings_cmnds:
            # Alter the settings
            self.settings_dict = sett(my_set[0], my_set[1:], self.settings_dict)
        # Update the sphere radii in the system
        if self.settings_dict is not None and self.settings_dict['atom_rad'] is not None:
            self.sys.set_radii(self.settings_dict['atom_rad']['element'], self.settings_dict['atom_rad']['special'])

    def create_groups(self):
        # Set the groups up
        ggroup(self.sys, self.groups, self.settings_dict)
    
    def run_exports(self):

        # Export everything
        argv_export(self.sys, self.exports)
    
