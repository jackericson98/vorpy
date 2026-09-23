import os
import sys
import gc
from pathlib import Path
import tkinter as tk
from tkinter import filedialog
from itertools import combinations
from copy import deepcopy
from contextlib import closing
from vorpy.src.inputs.frames import iter_pdb_frames, count_pdb_frames
from vorpy.src.inputs.net import read_net
from vorpy.src.command.commands import *
from vorpy.src.system.system import System
from vorpy.src.command.interpret import get_file
from vorpy.src.command.set import sett
from vorpy.src.command.group import ggroup
from vorpy.src.command.command_export import argv_export
from vorpy.src.command.interface import build_interfaces


FRAME_OPTIONS = ('load_commands', 'groups', 'builds', 'exports',
                 'settings_cmnds', 'logs_files', 'interface_mode',
                 'verbose', 'diagnose_edges')


def available_cpu_count():
    """Count logical CPUs available to this process where supported."""
    count = getattr(os, 'process_cpu_count', os.cpu_count)() or 1
    if hasattr(os, 'sched_getaffinity'):
        count = min(count, len(os.sched_getaffinity(0)))
    return max(1, count)


def prompt_frame_workers(total_frames):
    """Return a worker count, or None for sequential frame processing."""
    while True:
        try:
            answer = input('Parallelize frame runs? [Y/n] ').strip().lower()
        except EOFError:
            raise SystemExit('No response received. Use --all-frames for sequential runs '
                             'or --parallel-frames N for parallel runs.') from None
        if answer in {'n', 'no'}:
            return None
        if answer in {'', 'y', 'yes'}:
            break
        print('Please enter y (parallel) or n (sequential).')

    cpus = available_cpu_count()
    print(f'{cpus} logical CPUs (hardware threads) available.')
    while True:
        try:
            answer = input('Enter CPU usage percentage (1-100) [50]: ').strip().rstrip('%')
        except EOFError:
            raise SystemExit('No response received. Use --parallel-frames N to choose '
                             'the worker count without a prompt.') from None
        try:
            percentage = float(answer or '50')
        except ValueError:
            percentage = 0
        if 1 <= percentage <= 100:
            workers = min(total_frames, max(1, int(cpus * percentage / 100)))
            print(f'Using {workers} frame worker(s) for {total_frames} frames '
                  f'({percentage:g}% of {cpus} logical CPUs, rounded down; minimum 1).')
            print('Frames are queued; each free worker starts the next frame. '
                  'This selects concurrency, not a strict CPU utilization limit.')
            return workers
        print('Enter a percentage from 1 to 100.')


def frame_cli_options(args):
    """Remove the parallel-frame option before the legacy grouped parser."""
    remaining = []
    workers = None
    args = iter(args)
    for arg in args:
        if arg == '--parallel-frames':
            if workers is not None:
                raise SystemExit('--parallel-frames may only be specified once.')
            try:
                workers = int(next(args))
            except (StopIteration, ValueError):
                raise SystemExit('--parallel-frames requires a positive integer worker count.') from None
            if workers < 1:
                raise SystemExit('--parallel-frames requires a positive integer worker count.')
        else:
            remaining.append(arg)
    return remaining, workers


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
        self.diagnose_edges = False
        self.save_network_path = None

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
        self._arguments = list(sys.argv[1:])
        if self._arguments and self._arguments[0] == '--load-network':
            self._arguments.pop(0)
        if not self._arguments:
            raise SystemExit('Provide a structure or .vpy network archive path')
        input_arg = self._arguments[0]

        self.base_file = get_file(input_arg, interactive=False)
        if self.base_file is None:
            raise SystemExit(f"Input file not found: {input_arg}. "
                             "Provide an existing path or a filename from vorpy/data.")

        _, workers = frame_cli_options(sys.argv[2:])
        if workers is not None:
            self.run_all_frames(workers=workers)
            return

        if '--all-frames' in sys.argv[2:]:
            self.run_all_frames()
            return

        if Path(self.base_file).suffix.lower() == '.pdb':
            total_frames = count_pdb_frames(self.base_file)
            if total_frames > 1:
                while True:
                    try:
                        answer = input(f'{total_frames} frames found. Run all? [Y/n] ').strip().lower()
                    except EOFError:
                        raise SystemExit('No response received. Use --all-frames to run all frames without a prompt.')
                    if answer in {'', 'y', 'yes'}:
                        workers = prompt_frame_workers(total_frames)
                        if workers is None:
                            self.run_all_frames(total_frames=total_frames)
                        else:
                            self.run_all_frames(total_frames=total_frames, workers=workers)
                        return
                    if answer in {'n', 'no'}:
                        print('Running frame 1 only.')
                        break
                    print('Please enter y (all frames) or n (frame 1 only).')

        # Create the system if one was not supplied
        if self.sys is None:
            self.sys = System(file=self.base_file, make_dir=False)

        # Parse the remaining command-line arguments
        self.parse_commands()

        self._process_system()

    def run_all_frames(self, total_frames=None, workers=None):
        """Run the normal workflow independently for every input PDB frame."""
        if Path(self.base_file).suffix.lower() != '.pdb':
            raise ValueError('Frame processing currently supports PDB input files only.')
        if workers is not None:
            from vorpy.src.command.parallel_frames import run_parallel_frames
            return run_parallel_frames(self, workers, total_frames)
        initial_settings = deepcopy(self.settings_dict)
        output_root = None
        if total_frames is None:
            total_frames = count_pdb_frames(self.base_file)
        print(f'{Path(self.base_file).name}: {total_frames} frames found; processing all frames sequentially.')
        option_names = ('load_commands', 'groups', 'builds', 'exports',
                        'settings_cmnds', 'logs_files', 'interface_mode',
                        'verbose', 'diagnose_edges', 'save_network_path')
        with closing(iter_pdb_frames(self.base_file)) as frames:
            for ordinal, frame_file in frames:
                frame_system = System(file=frame_file, make_dir=False)
                frame_system.frame_index = ordinal
                frame_system.frame_count = total_frames
                command = Command(sys=frame_system, settings=deepcopy(initial_settings))
                command.base_file = frame_file
                if output_root is None:
                    command.parse_commands()
                    if frame_system.files['dir'] is None:
                        frame_system.set_output_directory()
                    output_root = Path(frame_system.files['dir']) / 'frames'
                    options = {name: deepcopy(getattr(command, name)) for name in FRAME_OPTIONS}
                else:
                    for name, value in options.items():
                        setattr(command, name, deepcopy(value))
                frame_system.set_output_directory(str(output_root / f'frame_{ordinal:04d}'))
                print(f'Processing frame {ordinal}: {frame_system.files["dir"]}')
                command._process_system()
                print(f'\nFrame {ordinal}/{total_frames} complete. Releasing frame memory.')
                # Systems, groups, networks and atom tables reference each other.
                # Drop all batch-owned references and collect those cycles before
                # allocating the next frame's geometry. Results remain on disk.
                self._release_frame_system(frame_system)
                self.sys = None
                del command, frame_system
                gc.collect()

    @staticmethod
    def _release_frame_system(system):
        """Break references through NumPy object tables that GC cannot trace.

        Only call after a batch frame has finished exporting: this consumes its
        in-memory state. In particular, atom tables hold references back to the
        System, so merely deleting local variables does not release them.
        """
        owners = list(system.groups or [])
        owners.extend(getattr(system, 'interfaces', None) or [])
        owners.extend(getattr(system, 'ifaces', None) or [])
        for owner in owners:
            network = getattr(owner, 'net', None)
            if network is not None:
                network.__dict__.clear()
            owner.__dict__.clear()
        system.__dict__.clear()

    def _process_system(self):
        """Build and export a loaded system using the parsed commands."""
        # Select the default only after export-directory commands have been applied.
        if self.sys.files['dir'] is None:
            self.sys.set_output_directory()

        # Expose the CLI verbosity state at the system level so downstream
        # build/analyze/export code has one common place to query it.
        self.sys.verbose = self.verbose

        # Load any additional files
        self.load_files()

        if getattr(self.sys, 'frame_count', 0):
            print(f'{self.sys.name}: {self.sys.frame_count} frames found; '
                  f'{self.sys.loaded_frame_count} frame loaded '
                  f'(frame {self.sys.frame_index} only), {len(self.sys.balls):,} atoms.')

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

        # Run optional read-only diagnostics after construction and before
        # export so the report describes the completed in-memory networks.
        if self.diagnose_edges:
            self.run_edge_diagnostics()

        # Export requested outputs
        self.run_exports()
        self._save_archive()

    def _save_archive(self):
        if self.save_network_path is not None:
            from vorpy.src.io import save_network
            destination = self.save_network_path
            if getattr(self.sys, 'frame_count', 1) > 1:
                destination = destination.with_name(f'{destination.stem}_frame_{self.sys.frame_index:04d}.vpy')
            save_network(self.sys, destination)
            print(f'Solved network saved: {destination}')

    def _run_archive(self):
        from vorpy.src.io import load_network, group_from_network
        network = load_network(self.base_file)
        self.sys = network.sys
        self.sys.set_output_directory(self.sys.files['dir'])
        self.parse_commands()
        if self.builds or self.load_commands or self.settings_cmnds or self.interface_mode:
            raise ValueError('Archive input reuses solved geometry. Use -e export options or -g selections; rebuild/settings commands require a structure input.')
        if self.groups:
            previous = list(self.sys.groups)
            ggroup(self.sys, self.groups, dict(network.settings), make_net=False)
            definitions = [group for group in self.sys.groups if group not in previous]
            # Some command paths replace the group list rather than append it.
            self.sys.groups = previous
            for definition in definitions:
                group_from_network(network, definition.ball_ndxs, definition.name)
        if self.diagnose_edges:
            self.run_edge_diagnostics()
        self.run_exports()
        self._save_archive()

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
        my_args = list(getattr(self, '_arguments', sys.argv[1:])[1 + counter:])
        if '--save-network' in my_args:
            index = my_args.index('--save-network')
            if index + 1 >= len(my_args) or my_args[index + 1].startswith('-'):
                raise ValueError('--save-network requires an output .vpy path')
            self.save_network_path = Path(my_args[index + 1]).expanduser().resolve()
            del my_args[index:index + 2]
        my_args = list(sys.argv[2 + counter:])
        my_args, _ = frame_cli_options(my_args)
        my_args = [arg for arg in my_args if arg != '--all-frames']

        # -v / --verbose is a universal argumentless flag. Remove it before
        # the existing command parser processes argument/value groups.
        if '-v' in my_args or '--verbose' in my_args:
            self.verbose = True
            my_args = [
                arg for arg in my_args
                if arg not in {'-v', '--verbose'}
            ]

        # Analytic edge auditing is also a universal argumentless flag.
        # Remove it before the legacy grouped parser handles option values.
        if '--diagnose-edges' in my_args:
            self.diagnose_edges = True
            my_args = [
                arg for arg in my_args
                if arg != '--diagnose-edges'
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
                if arg_cmnds[0] in {'dir', 'directory'}:
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
        # Establish CLI defaults before applying explicit command-line overrides.
        if self.settings_dict is None:
            self.settings_dict = sett('mv', ['5'])
        elif self.settings_dict.get('max_vert') is None:
            self.settings_dict['max_vert'] = 5.0
        # Go through the user inputs loading files
        for my_set in self.settings_cmnds:
            # Alter the settings
            self.settings_dict = sett(my_set[0], my_set[1:], self.settings_dict)
        # Update the sphere radii in the system
        if self.settings_dict is not None and self.settings_dict['atom_rad'] is not None:
            self.sys.set_radii(self.settings_dict['atom_rad']['element'], self.settings_dict['atom_rad']['special'])
        if self.settings_dict is not None:
            self.sys.round_to = self.settings_dict.get('round_to', self.sys.round_to)
            self.sys.file_type = self.settings_dict.get('file_type', 'off')

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

    def _diagnostic_networks(self):
        """Yield each unique named network created by this command."""
        seen = set()
        owners = list(self.sys.groups or [])
        owners.extend(list(getattr(self.sys, 'interfaces', None) or []))

        for owner in owners:
            net = getattr(owner, 'net', None)
            if net is None or id(net) in seen:
                continue
            seen.add(id(net))
            name = (
                getattr(owner, 'name', None)
                or getattr(net, 'group_name', None)
                or 'Network'
            )
            yield str(name), net

    @staticmethod
    def _diagnostic_number(value):
        """Format a finite diagnostic value while preserving unavailable data."""
        try:
            if value is None or not float('-inf') < float(value) < float('inf'):
                return 'n/a'
            return f"{float(value):.6g}"
        except (TypeError, ValueError):
            return 'n/a'

    def run_edge_diagnostics(self, max_problem_edges=10):
        """Print analytic AW edge audits for the networks built by the CLI."""
        networks = list(self._diagnostic_networks())
        if not networks:
            print("\nNo constructed networks are available for edge diagnostics.")
            return

        for name, net in networks:
            settings = getattr(net, 'settings', None) or {}
            net_type = settings.get('net_type', 'aw')
            if net_type != 'aw':
                print(
                    f"\nSkipping analytic edge diagnostic for {name}: "
                    f"network type '{net_type}' is not AW."
                )
                continue

            try:
                report = net.diagnose_aw_edges()
                summary = report.summary()
            except (AttributeError, TypeError, ValueError) as error:
                print(f"\nAnalytic edge diagnostic failed for {name}: {error}")
                continue

            counts = summary['status_counts']
            failed = summary['total_edges'] - summary['matched_edges']
            print("\n" + "="*70)
            print(f"ANALYTIC AW EDGE DIAGNOSTIC: {name}")
            print("="*70)
            print(f"Total edges:                  {summary['total_edges']:,}")
            print(f"Matched analytic curves:      {counts.get('matched_curve', 0):,}")
            print(f"Matched turning conics:       {counts.get('matched_nonsingular', 0):,}")
            print(f"Matched straight edges:       {counts.get('matched_line', 0):,}")
            print(f"Unsupported / failed:         {failed:,}")
            print(f"Match coverage:               {100*summary['matched_fraction']:.2f} %")
            print(
                "Maximum sample residual:      "
                f"{self._diagnostic_number(summary['maximum_sample_error'])} A"
            )
            print(
                "RMS sample residual:          "
                f"{self._diagnostic_number(summary['rms_sample_error'])} A"
            )
            print(
                "Maximum |length difference|: "
                f"{self._diagnostic_number(summary['maximum_absolute_length_difference'])} A"
            )
            print(f"Status counts:                {counts}")

            problems = [
                record for record in report.records
                if not record.status.startswith('matched')
            ]
            if problems:
                print(f"Problem edges (first {min(len(problems), max_problem_edges)}):")
                for record in problems[:max_problem_edges]:
                    print(
                        f"  edge {record.edge_index}: {record.status}; "
                        f"balls={record.ball_indices}; verts={record.vertex_indices}; "
                        f"{record.reason}"
                    )
                if len(problems) > max_problem_edges:
                    print(
                        f"  ... {len(problems) - max_problem_edges:,} additional "
                        "problem edges omitted"
                    )
            print("="*70)
