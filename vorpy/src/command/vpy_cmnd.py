import os
import sys
from pathlib import Path
import tkinter as tk
from tkinter import filedialog
from itertools import combinations
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
        self.interface_mode = False
        self.settings_cmnds = []
        self.settings_dict = settings
        self.logs_files = []

        # Global diagnostic-output flag.
        # Timers/metrics may still be collected when False; this controls
        # whether verbose diagnostic information is printed.
        self.verbose = False

    def run(self):
        self._run_pipeline()

    def _run_pipeline(self):
        """
        Run the complete command-line workflow.

        Pipeline order:
            1. Resolve the base input file.
            2. Create the System.
            3. Parse command-line arguments.
            4. Load additional files.
            5. Apply settings.
            6. Create groups.
            7. Build interface, comparison, logs, or normal networks.
            8. Create standard interfaces when not using interface mode.
            9. Run exports.
        """

        # Resolve the base input file
        input_arg = sys.argv[1]

        if input_arg[-3:].lower() in {"pdb", "gro", "mol", "cif", "txt"}:
            self.base_file = input_arg
        else:
            resolved_file = get_file(input_arg)

            if resolved_file is None:
                self.base_file = None
                print(f"{input_arg} is not a valid input file")
                return

            self.base_file = resolved_file

        # Create the system if one was not supplied
        if self.sys is None:
            self.sys = System(file=self.base_file)

        # Parse the remaining command-line arguments
        self.parse_commands()

        # Expose the CLI verbosity state at the system level so downstream
        # build/analyze/export code has one common place to query it.
        self.sys.verbose = self.verbose

        # Load any additional files
        self.load_files()

        # Apply command-line settings
        self.apply_settings()

        # Create the requested groups
        self.create_groups()

        comparison_mode = (
                self.settings_dict is not None
                and isinstance(self.settings_dict.get("net_type"), list)
                and len(self.settings_dict["net_type"]) > 1
                and self.settings_dict["net_type"][0] == "com"
        )

        build_type = None

        if self.settings_dict is not None:
            build_type = self.settings_dict.get("bld_type")

        # Comparison mode:
        # Build two network types for each group and compare them.
        if comparison_mode:
            new_groups = []

            first_net_type = self.settings_dict["net_type"][1]
            second_net_type = self.settings_dict["net_type"][2]

            for grp in self.sys.groups:
                copy_group = deepcopy(grp)

                copy_group.name = f"{copy_group.name}_{first_net_type}"
                copy_group.settings["net_type"] = first_net_type

                grp.name = f"{grp.name}_{second_net_type}"
                grp.settings["net_type"] = second_net_type

                # Vertices from the original group should not be reused blindly
                # for the second network type.
                grp.verts = None

                copy_group.build()
                grp.build()

                new_groups.append(copy_group)

                self.sys.compare_networks(
                    group1=copy_group,
                    group2=grp,
                )

            self.sys.groups += new_groups

        # Reconstruct networks from logs
        elif build_type == "logs":
            self.build_groups_from_logs()

        # Preserve the existing interface workflow during normal builds.
        # Interface mode handles its own interface construction.
        elif self.interface_mode:
            num_requested_groups = len(self.groups)
            num_created_groups = len(self.sys.groups or [])

            if num_requested_groups >= 2 > num_created_groups:
                raise ValueError(
                    "Interface mode requested multiple groups, but fewer than "
                    "two non-empty groups were created.\n"
                    f"Requested group commands: {self.groups}\n"
                    f"Created groups: "
                    f"{[group.name for group in (self.sys.groups or [])]}"
                )

            if num_created_groups != num_requested_groups:
                print("\nWARNING: Requested and created group counts differ.")
                print(f"  requested groups: {num_requested_groups}")
                print(f"  created groups: {num_created_groups}")

            if num_created_groups >= 2:

                interface_pairs = list(combinations(self.sys.groups, 2))
            else:
                interface_pairs = build_interfaces(sys=self.sys, num_requested_groups=num_requested_groups)

            self.sys.make_interfaces(interface_pairs)

        # Standard full-group build
        else:
            for grp in self.sys.groups:
                grp.build()

        # Export requested outputs
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
        my_args = list(sys.argv[2 + counter:])

        # -v / --verbose is a universal argumentless flag. Remove it before
        # the existing command parser processes argument/value groups.
        if '-v' in my_args or '--verbose' in my_args:
            self.verbose = True
            my_args = [
                arg for arg in my_args
                if arg not in {'-v', '--verbose'}
            ]

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
            elif arg.lower() == '-i':
                self.interface_mode = True

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
                        (self.sys.atoms is not None or self.sys.files['verts_file'] is not None or self.sys.files[
                            'net_file'] is not None):
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
                                  .format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'],
                                          self.sys.net.settings['box_size']))
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
                            print(
                                "{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box"
                                " size: {} x\n".format(self.sys.name, len(self.sys.net.verts),
                                                       self.sys.net.settings['max_vert'],
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
                            print(
                                "{} network loaded - surface resolution: {}\u208B, maximum vertex radius: {} \u208B, box size: {} x\n"
                                .format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'],
                                        self.sys.net.settings['box_size']))
                        else:
                            print("{} vertices loaded - {} vertices, maximum vertex radius: {} \u208B, box size: {} x\n"
                                  .format(self.sys.name, len(self.sys.net.verts), self.sys.net.settings['max_vert'],
                                          self.sys.net.settings['box_size']))
            # Check to see if it is a new network file
            elif file[-3:] == 'csv':
                # Check to see that this is a network file
                if file[-7:-4].lower() == 'net':
                    self.sys.load_net(file=file)

            # If the file is an index file load it accordingly
            elif file[-3:] == 'ndx':
                self.sys.load_ndx(file)
                print(self.sys.ndx_file + "loaded -  {}".format(
                    self.sys.ndx_names[:min(len(self.sys.ndx_names) - 1, 10)]))
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
        # Interface mode only needs group definitions; the Interface creates
        # the single network that will actually be solved.
        ggroup(
            self.sys,
            self.groups,
            self.settings_dict,
            make_net=not self.interface_mode
        )

        # Propagate the universal verbose flag into every created group and
        # network. This lets existing timing code use:
        #
        #     if net.settings.get('verbose', False):
        #         ...
        #
        # without making verbose output a calculation setting.
        for group in (self.sys.groups or []):
            if getattr(group, 'settings', None) is not None:
                group.settings['verbose'] = self.verbose

            net = getattr(group, 'net', None)
            if net is not None:
                if getattr(net, 'settings', None) is None:
                    net.settings = {}
                net.settings['verbose'] = self.verbose

    def run_exports(self):
        # Export everything
        argv_export(self.sys, self.exports)

