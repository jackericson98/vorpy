import gc
import os
import re
import sys
import tkinter as tk
from contextlib import closing
from copy import deepcopy
from itertools import combinations
from pathlib import Path
from tkinter import filedialog

from vorpy.src.boundary import BoundaryConfig, resolve_boundary
from vorpy.src.command.command_export import argv_export
from vorpy.src.command.commands import *
from vorpy.src.command.group import ggroup
from vorpy.src.command.interface import build_interfaces
from vorpy.src.command.interpret import get_file
from vorpy.src.command.set import sett
from vorpy.src.inputs.frames import count_pdb_frames, iter_pdb_frames
from vorpy.src.inputs.net import read_net
from vorpy.src.system.system import System

FRAME_OPTIONS = ('load_commands', 'groups', 'builds', 'exports',
                 'settings_cmnds', 'logs_files', 'interface_mode',
                 'verbose', 'diagnose_edges', 'boundary_config',
                 'boundary_explicitly_requested', 'boundary_warning_printed',
                 'apollonius_requested', 'alpha_value', 'power_settings',
                 'power_vs_aw_requested', 'power_cli_arguments',
                 'aw_alpha_values', 'aw_alpha_pair_overlap_only',
                 'visualize_dual_requested', 'visualize_alpha_value')


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


def boundary_cli_options(args):
    """Extract boundary options before the legacy grouped command parser."""
    remaining = []
    values = {}
    args = iter(args)
    for arg in args:
        if arg in {'--boundary', '--probe-radius', '--shell-spacing', '--shell-offset'}:
            try:
                value = next(args)
            except StopIteration:
                raise SystemExit(f'{arg} requires a value') from None
            key = {'--boundary': 'mode', '--probe-radius': 'probe_radius',
                   '--shell-spacing': 'shell_spacing', '--shell-offset': 'shell_offset'}[arg]
            if key in values:
                raise SystemExit(f'{arg} may only be specified once')
            try:
                values[key] = value if key == 'mode' else float(value)
            except ValueError:
                raise SystemExit(f'{arg} requires a number') from None
        else:
            remaining.append(arg)
    try:
        return remaining, BoundaryConfig(**values), 'mode' in values
    except (ValueError, TypeError) as error:
        raise SystemExit(str(error)) from None


def visualization_cli_options(args):
    """Extract read-only dual visualization flags before legacy parsing."""
    remaining = []
    requested = False
    seen_dual_flag = False
    alpha = None
    args = iter(args)
    for arg in args:
        if arg == '--export-dual-visualization':
            if seen_dual_flag:
                raise SystemExit('--export-dual-visualization may only be specified once')
            seen_dual_flag = True
            requested = True
        elif arg == '--export-alpha-visualization':
            if alpha is not None:
                raise SystemExit('--export-alpha-visualization may only be specified once')
            try:
                alpha = float(next(args))
            except (StopIteration, ValueError):
                raise SystemExit('--export-alpha-visualization requires a finite alpha value') from None
            if not (-float('inf') < alpha < float('inf')):
                raise SystemExit('--export-alpha-visualization requires a finite alpha value')
            requested = True
        else:
            remaining.append(arg)
    return remaining, requested, alpha


def apollonius_cli_options(args):
    """Extract opt-in dual and alpha-complex flags from legacy CLI arguments."""
    remaining = []
    requested = False
    seen_flag = False
    alpha = None
    args = iter(args)
    for arg in args:
        if arg == '--apollonius':
            if seen_flag:
                raise SystemExit('--apollonius may only be specified once')
            seen_flag = True
            requested = True
        elif arg == '--alpha-complex':
            if alpha is not None:
                raise SystemExit('--alpha-complex may only be specified once')
            try:
                alpha = float(next(args))
            except (StopIteration, ValueError):
                raise SystemExit('--alpha-complex requires a finite numeric value in Å') from None
            if not (-float('inf') < alpha < float('inf')):
                raise SystemExit('--alpha-complex requires a finite numeric value in Å')
            requested = True
        else:
            remaining.append(arg)
    return remaining, requested, alpha


def prompt_boundary(config):
    """Choose an optional boundary on a terminal; batch runs use no boundary."""
    if not sys.stdin.isatty():
        return config
    print("\nNo explicit solvent was detected. Exterior cells may be unbounded.\n"
          "Boundary methods:\n"
          "    [0] No boundary (default; use existing generators only)\n"
          "    [1] Virtual solvent shell\n"
          "SES, SAS, and AW probe boundaries are planned and unavailable.")
    while True:
        try:
            answer = input("Choose boundary [0]: ").strip().lower()
        except EOFError:
            answer = ""
        modes = {"": "none", "0": "none", "none": "none", "1": "shell", "shell": "shell"}
        if answer in modes:
            return BoundaryConfig(modes[answer], config.probe_radius,
                                  config.shell_spacing, config.shell_offset)
        print("Enter 0 (none) or 1 (shell).")


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
        self.boundary_config = BoundaryConfig()
        self.boundary_explicitly_requested = False
        self.boundary_warning_printed = False
        self.apollonius_requested = False
        self.alpha_value = None
        self.power_settings = {'preset': None, 'probe_radius': 1.4, 'alpha': 0.0,
                               'condition_beta_m': 5.0, 'compare_aw': False}
        self.power_vs_aw_requested = False
        self.power_cli_arguments = []
        self.power_aw_water_records_omitted = 0
        self.aw_alpha_values = []
        self.aw_alpha_pair_overlap_only = False
        self.aw_alpha_solvent_chains_normalized = 0
        self.visualize_dual_requested = False
        self.visualize_alpha_value = None

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
        from vorpy.src.analyze.power_interface_cli import (
            extract_power_interface_options,
        )
        self._arguments, self.power_settings = extract_power_interface_options(self._arguments)
        self.power_vs_aw_requested = self.power_settings['compare_aw']
        if not self._arguments:
            raise SystemExit('Provide a structure or .vpy network archive path')
        input_arg = self._arguments[0]

        self._arguments, self.boundary_config, self.boundary_explicitly_requested = boundary_cli_options(self._arguments)
        (self._arguments, self.visualize_dual_requested,
         self.visualize_alpha_value) = visualization_cli_options(self._arguments)
        if len(self._arguments) > 1:
            from vorpy.src.geometry.aw_alpha.cli import (
                extract_aw_alpha_options,
                extract_aw_alpha_pair_only,
            )
            tail, self.aw_alpha_values = extract_aw_alpha_options(self._arguments[1:])
            tail, self.aw_alpha_pair_overlap_only = extract_aw_alpha_pair_only(tail)
            self._arguments = [self._arguments[0], *tail]
            if self.aw_alpha_pair_overlap_only and not any(
                value == 0.0 for value in self.aw_alpha_values
            ):
                raise SystemExit(
                    "--aw-alpha-pairs-only-experimental requires "
                    "--aw-alpha-experimental 0."
                )
        if not self._arguments:
            raise SystemExit('Provide a structure or .vpy network archive path')
        input_arg = self._arguments[0]

        self.base_file = get_file(input_arg, interactive=False)
        if self.base_file is None:
            raise SystemExit(f"Input file not found: {input_arg}. "
                             "Provide an existing path or a filename from vorpy/data.")

        if self.power_settings['preset'] and not self.power_vs_aw_requested:
            from vorpy.src.analyze.power_interface_cli import (
                chain_groups_from_legacy_args,
                power_output_directory,
                run_power_interface_pdb,
            )
            chains_a, chains_b = chain_groups_from_legacy_args(self._arguments[1:])
            destination = power_output_directory(self._arguments[1:], self.base_file)
            try:
                run = run_power_interface_pdb(
                    self.base_file, chains_a, chains_b, destination,
                    probe_radius=self.power_settings['probe_radius'],
                    alpha=self.power_settings['alpha'],
                    condition_beta_m=self.power_settings['condition_beta_m'],
                )
            except (ValueError, RuntimeError) as error:
                raise SystemExit(str(error)) from None
            print(run['summary'], end='')
            print(f"Power-interface exports written to: {destination}")
            return

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

        # The legacy PDB reader dereferences its optional solvent container
        # for crystallographic waters. For the two-protein AW comparison only,
        # use a temporary copy without water HETATM records; neither the source
        # PDB nor the independent Cazals/Power preparation is altered.
        comparison_temp = None
        system_input = self.base_file
        if self.power_vs_aw_requested and Path(self.base_file).suffix.lower() == '.pdb':
            from vorpy.src.analyze.power_interface_cli import (
                aw_comparison_structure_without_waters,
            )
            comparison_temp, self.power_aw_water_records_omitted = aw_comparison_structure_without_waters(self.base_file)
            system_input = comparison_temp
        elif self.aw_alpha_values and Path(self.base_file).suffix.lower() == '.pdb':
            from vorpy.src.geometry.aw_alpha.cli import parser_safe_aw_alpha_structure
            comparison_temp, self.aw_alpha_solvent_chains_normalized = parser_safe_aw_alpha_structure(self.base_file)
            if comparison_temp is not None:
                system_input = comparison_temp
        try:
            if self.sys is None:
                self.sys = System(file=str(system_input), make_dir=False, boundary_config=self.boundary_config)
                if comparison_temp is not None:
                    self.sys.name = Path(self.base_file).stem
            else:
                self.sys.boundary_config = self.boundary_config
        finally:
            if comparison_temp is not None:
                comparison_temp.unlink(missing_ok=True)

        # Parse the remaining command-line arguments
        self.parse_commands()
        if self.power_vs_aw_requested:
            self.interface_mode = True
            self.power_cli_arguments = list(self._arguments[1:])

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
        option_names = (*FRAME_OPTIONS, 'save_network_path')
        with closing(iter_pdb_frames(self.base_file)) as frames:
            for ordinal, frame_file in frames:
                frame_input = frame_file
                parser_temp = None
                normalized_solvent = 0
                if self.aw_alpha_values:
                    from vorpy.src.geometry.aw_alpha.cli import (
                        parser_safe_aw_alpha_structure,
                    )
                    parser_temp, normalized_solvent = parser_safe_aw_alpha_structure(frame_file)
                    if parser_temp is not None:
                        frame_input = parser_temp
                try:
                    frame_system = System(
                        file=frame_input, make_dir=False,
                        boundary_config=self.boundary_config,
                    )
                finally:
                    if parser_temp is not None:
                        parser_temp.unlink(missing_ok=True)
                frame_system.frame_index = ordinal
                frame_system.frame_count = total_frames
                command = Command(sys=frame_system, settings=deepcopy(initial_settings))
                command.base_file = frame_file
                command._arguments = list(self._arguments)
                command.aw_alpha_values = list(self.aw_alpha_values)
                command.aw_alpha_pair_overlap_only = self.aw_alpha_pair_overlap_only
                command.aw_alpha_solvent_chains_normalized = normalized_solvent
                command.visualize_dual_requested = self.visualize_dual_requested
                command.visualize_alpha_value = self.visualize_alpha_value
                command.boundary_config = deepcopy(self.boundary_config)
                command.boundary_explicitly_requested = self.boundary_explicitly_requested
                command.boundary_warning_printed = self.boundary_warning_printed
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
                self.boundary_warning_printed = command.boundary_warning_printed
                self.boundary_config = command.boundary_config
                if output_root is not None:
                    options['boundary_warning_printed'] = self.boundary_warning_printed
                    options['boundary_config'] = deepcopy(self.boundary_config)
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

        molecular_suffixes = {'.pdb', '.cif', '.gro', '.mol', '.sdf', '.mol2', '.txt'}
        if Path(self.base_file).suffix.lower() in molecular_suffixes:
            if (not self.boundary_explicitly_requested and not self.boundary_warning_printed
                    and not self.sys.has_explicit_solvent):
                self.boundary_config = prompt_boundary(self.boundary_config)
                self.boundary_warning_printed = True
            try:
                self.sys.boundary_config = resolve_boundary(
                    self.boundary_config, self.sys.has_explicit_solvent
                )
            except (ValueError, NotImplementedError) as error:
                raise SystemExit(str(error)) from None
            self.sys.boundary_mode = self.sys.boundary_config.mode.value


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

        if self.power_vs_aw_requested:
            self.run_power_aw_comparison()
            # This analysis mode writes its own compact Power/AW artifacts.
            # Avoid the legacy default bulk exporter, which is unrelated to
            # the analysis output and currently does not support Interface's
            # selected export keyword set.
            return

        if self.apollonius_requested:
            self.run_apollonius()

        if self.visualize_dual_requested:
            self.run_dual_visualization()

        if self.aw_alpha_values:
            self.run_aw_alpha_experimental()
            self._save_archive()
            return

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
        from vorpy.src.io import group_from_network, load_network
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
        if self.apollonius_requested:
            self.run_apollonius()
        if self.visualize_dual_requested:
            self.run_dual_visualization()
        if self.aw_alpha_values:
            self.run_aw_alpha_experimental()
            self._save_archive()
            return
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
        my_args, boundary_config, boundary_requested = boundary_cli_options(my_args)
        if boundary_requested:
            self.boundary_config = boundary_config
            self.boundary_explicitly_requested = True
            if self.sys is not None:
                self.sys.boundary_config = boundary_config
        my_args, apollonius_requested, alpha_value = apollonius_cli_options(my_args)
        self.apollonius_requested = self.apollonius_requested or apollonius_requested
        if alpha_value is not None:
            self.alpha_value = alpha_value
        from vorpy.src.geometry.aw_alpha.cli import extract_aw_alpha_options
        my_args, aw_alpha_values = extract_aw_alpha_options(my_args)
        if aw_alpha_values:
            self.aw_alpha_values = aw_alpha_values
        from vorpy.src.geometry.aw_alpha.cli import extract_aw_alpha_pair_only
        my_args, pair_overlap_only = extract_aw_alpha_pair_only(my_args)
        if pair_overlap_only:
            self.aw_alpha_pair_overlap_only = True
        if '--save-network' in my_args:
            index = my_args.index('--save-network')
            if index + 1 >= len(my_args) or my_args[index + 1].startswith('-'):
                raise ValueError('--save-network requires an output .vpy path')
            self.save_network_path = Path(my_args[index + 1]).expanduser().resolve()
            del my_args[index:index + 2]
        # Use the normalized argv captured by _run_pipeline so custom options
        # such as --boundary never reach the legacy grouped parser.
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
                        base_dir = os.path.abspath(os.path.expanduser(arg_cmnds[1]))
                        system_dir = os.path.join(base_dir, self.sys.name)

                        os.makedirs(system_dir, exist_ok=True)

                        self.sys.dir = system_dir
                        self.sys.files['dir'] = system_dir

                        print(f"Directory set to: {system_dir}")
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
                    reset_sys = input(f"replacing {self.sys.name} with {file}\nconfirm >>>   "
                                      )
                    # If the user confirms the replacement, create a new system
                    if reset_sys.lower() in ys:
                        self.sys = System(file)
                        print(self.sys.name + f" loaded - {len(self.sys.atoms)} atoms, {len(self.sys.chains)} molecules, solute: {self.sys.sol.name}"
                              )
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
                        replace_net_file = input(f"replacing {self.sys.net_file} with {file}\n "
                                                 "confirm >>>   ")
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
                print(self.sys.ndx_file + f"loaded -  {self.sys.ndx_names[:min(len(self.sys.ndx_names) - 1, 10)]}")
            # In all other case print an error and give the user a chance to try again
            else:
                print(f"\'{file}\' is not a valid input. allowed file types: .pdb, .mol, .cif, .gro, .txt, .ndx. type "
                      "\'h\' for help")
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

    def run_apollonius(self):
        """Build/export derived complexes from completed in-memory networks."""
        from vorpy.src.analyze.apollonius import ApolloniusComplex
        from vorpy.src.analyze.apollonius_export import export_apollonius

        networks = list(self._diagnostic_networks())
        if not networks:
            raise SystemExit('Apollonius analysis requested, but no solved network is available.')
        root = Path((getattr(self.sys, 'files', {}) or {}).get('dir') or Path.cwd())
        for name, network in networks:
            settings = getattr(network, 'settings', None) or {}
            if settings.get('net_type', 'aw') != 'aw':
                raise SystemExit(
                    f"Apollonius analysis requires an AW network; '{name}' uses "
                    f"'{settings.get('net_type', 'aw')}'."
                )
            complex_ = ApolloniusComplex(network, name=name)
            alpha_complex = None
            if self.alpha_value is not None:
                complex_.calculate_alpha_births()
                alpha_complex = complex_.alpha_complex(self.alpha_value)
                closure_errors = alpha_complex.validate()
                if closure_errors:
                    raise SystemExit(
                        'Extracted alpha complex is not closed: '
                        + '; '.join(closure_errors[:5])
                    )
            slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('._') or 'network'
            destination = root / 'apollonius' / slug
            print(complex_.summary(), end='')
            if alpha_complex is not None:
                print(alpha_complex.summary(), end='')
            export_apollonius(complex_, destination, alpha_complex=alpha_complex)
            print(f"Apollonius exports written to: {destination}")

    def run_aw_alpha_experimental(self):
        """Analyze already-solved AW networks with the experimental filtration."""
        from vorpy.src.geometry.aw_alpha.cli import run_aw_alpha_experimental

        networks = list(self._diagnostic_networks())
        root = Path((getattr(self.sys, 'files', {}) or {}).get('dir') or Path.cwd())
        try:
            results = run_aw_alpha_experimental(
                networks, root / 'aw_alpha', self.aw_alpha_values,
                pair_overlap_only=self.aw_alpha_pair_overlap_only,
            )
        except (TypeError, ValueError, RuntimeError) as error:
            raise SystemExit(str(error)) from None
        for name, result in results:
            print(f"Experimental AW alpha filtration: {name}; "
                  f"unresolved births={len(result['unresolved'])}; "
                  f"pair-overlap-only={self.aw_alpha_pair_overlap_only}; "
                  f"water/ion records temporarily normalized for PDB loading="
                  f"{self.aw_alpha_solvent_chains_normalized}; "
                  f"output={root / 'aw_alpha'}")

    def run_dual_visualization(self):
        """Export solved-network dual layers and optional native-alpha layers."""
        from vorpy.src.geometry.duals import build_dual
        from vorpy.src.geometry.filtrations import build_alpha_filtration
        from vorpy.src.geometry.visualization import export_dual_visualization

        networks = list(self._diagnostic_networks())
        if not networks:
            raise SystemExit('Dual visualization requested, but no solved network is available.')
        root = Path((getattr(self.sys, 'files', {}) or {}).get('dir') or Path.cwd())
        for name, network in networks:
            dual = build_dual(network)
            filtration = None
            if self.visualize_alpha_value is not None:
                scheme = str((getattr(network, 'settings', None) or {}).get('net_type', 'aw'))
                max_dimension = 1 if scheme == 'aw' else 3
                filtration = build_alpha_filtration(network, max_dimension=max_dimension)
            slug = re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('._') or 'network'
            destination = root / 'visualization' / slug
            metadata = export_dual_visualization(
                network, dual, destination, filtration=filtration,
                alpha=self.visualize_alpha_value,
                structure_path=self.base_file if Path(self.base_file).suffix.lower() == '.pdb' else None,
            )
            print(f"Dual visualization written: {destination}; "
                  f"full={metadata['full_dual_counts']}; alpha={metadata['alpha_counts']}")

    def run_power_aw_comparison(self):
        """Run Power beside a full-system AW network and audit the interface net."""
        from copy import deepcopy

        import numpy as np

        from vorpy.src.analyze.aw_interface_curvature import (
            analyze_aw_interface_curvature,
        )
        from vorpy.src.analyze.aw_pipeline_audit import (
            audit_interface_surfaces,
            compare_surface_tables,
            write_rows,
        )
        from vorpy.src.analyze.power_interface_cli import (
            chain_groups_from_legacy_args,
            power_output_directory,
            run_power_interface_pdb,
        )
        from vorpy.src.analyze.power_vs_aw import compare_power_and_aw
        from vorpy.src.group import Group

        if Path(self.base_file).suffix.lower() != '.pdb':
            raise SystemExit('--power-vs-aw currently requires a PDB structure')
        if not getattr(self.sys, 'ifaces', None):
            raise SystemExit('--power-vs-aw requires two non-empty -g groups')
        interface = self.sys.ifaces[0]
        if interface.net is None:
            raise SystemExit('The requested AW interface network was not built')
        if (getattr(interface.net, 'settings', None) or {}).get('net_type', 'aw') != 'aw':
            raise SystemExit('--power-vs-aw requires the AW network type')
        chains_a, chains_b = chain_groups_from_legacy_args(self.power_cli_arguments)
        root = Path(self.sys.files.get('dir') or power_output_directory(self.power_cli_arguments, self.base_file))
        power_run = run_power_interface_pdb(
            self.base_file, chains_a, chains_b, root / 'power_interface',
            probe_radius=self.power_settings['probe_radius'],
            alpha=self.power_settings['alpha'],
            condition_beta_m=self.power_settings['condition_beta_m'],
        )

        # An Interface Network intentionally retains only vertices/edges/
        # surfaces spanning the two requested groups. Its cells are therefore
        # not whole-cell complexes, and its Apollonius subcomplex need not have
        # all same-group faces. Build a separate normal, full-system network for
        # cell completeness and supported AW analysis. Both networks use the
        # same System coordinates/radii and the same settings (including mv).
        full_group = Group(
            sys=self.sys, name='power_vs_aw_full_system',
            atoms=list(range(len(self.sys.balls))),
            settings=deepcopy(interface.settings), make_net=True,
            build_net=False, print_metrics=False,
        )
        if full_group in (self.sys.groups or []):
            self.sys.groups.remove(full_group)
        full_group.build()
        full_net = full_group.net
        if not full_net.balls.index.equals(interface.net.balls.index):
            raise RuntimeError('Full and interface AW networks do not share generator IDs.')
        coordinate_identity = all(
            np.array_equal(np.asarray(full_net.balls.at[index, 'loc'], dtype=float),
                           np.asarray(interface.net.balls.at[index, 'loc'], dtype=float))
            and float(full_net.balls.at[index, 'rad']) == float(interface.net.balls.at[index, 'rad'])
            for index in full_net.balls.index
        )
        if not coordinate_identity:
            raise RuntimeError('Full and interface AW networks have different generator coordinates/radii.')

        filtered_audit = audit_interface_surfaces(
            interface.net, interface.group1_indices, interface.group2_indices,
            full_reference=full_net,
        )
        old_surface_csv = (Path(__file__).resolve().parents[3] / 'output' / '2KAI_aw'
                           / 'data' / 'aw_interface_2KAI_supported' / 'aw_interface_surfaces.csv')
        surface_comparison = compare_surface_tables(
            filtered_audit, full_net, interface.group1_indices,
            interface.group2_indices,
            previous_csv=old_surface_csv if Path(self.base_file).stem.upper() == '2KAI' else None,
        )
        write_rows(root / 'aw_surface_completeness_audit.csv', filtered_audit)
        write_rows(root / 'aw_pipeline_surface_comparison.csv', surface_comparison)
        (root / 'aw_pipeline_audit.txt').write_text(
            'AW Pipeline Reconciliation\n'
            f'Source PDB: {self.base_file}\n'
            f'AW parser atoms retained: {len(self.sys.balls)}; waters omitted for parser compatibility: '
            f'{self.power_aw_water_records_omitted}\n'
            f'Interface network: {len(interface.net.balls)} generators, {len(interface.net.verts)} vertices, '
            f'{len(interface.net.edges)} edges, {len(interface.net.surfs)} surfaces; '
            f'complete cells={int(interface.net.balls.complete.astype(bool).sum())}\n'
            f'Full network: {len(full_net.balls)} generators, {len(full_net.verts)} vertices, '
            f'{len(full_net.edges)} edges, {len(full_net.surfs)} surfaces; '
            f'complete cells={int(full_net.balls.complete.astype(bool).sum())}\n'
            f'Network settings (full/interface): {full_net.settings} / {interface.net.settings}\n'
            f'Full/interface coordinate and radius identity: {coordinate_identity}\n'
            'The interface network is intentionally topology-filtered to cross-group features. Its cell.complete '
            'flags cannot establish completeness of whole molecular cells. The full network is used for supported '
            'surface/edge analysis; no completeness criterion was relaxed.\n'
            f'Filtered interface candidates: {len(filtered_audit)}; feature-complete by closed-boundary audit: '
            f'{sum(bool(row["feature_complete"]) for row in filtered_audit)}\n'
            f'Surface identity comparison rows: {len(surface_comparison)}; prior standalone artifact: '
            f'{old_surface_csv if old_surface_csv.exists() and Path(self.base_file).stem.upper() == "2KAI" else "unavailable"}\n',
            encoding='utf-8',
        )

        aw = analyze_aw_interface_curvature(
            full_net, interface.group1_indices, interface.group2_indices,
            selection_mode='supported',
        )
        aw.export_csv(root / 'aw_interface')
        lys15 = set()
        # Use the loader's explicit author-facing residue columns. Residue and
        # Chain object labels differ across readers (e.g. sequence vs seq), so
        # attribute probing can silently miss a valid author-chain selection.
        atoms = self.sys.balls
        lys15.update(int(index) for index, row in atoms.iterrows()
                     if str(row.get('chain_name', '')).strip() == 'I'
                     and str(row.get('res_name', '')).strip().upper() == 'LYS'
                     and str(row.get('res_seq', '')).strip() == '15')
        comparison = compare_power_and_aw(power_run, aw, root,
                                          aw_lys15_ids=lys15,
                                          input_notes=(f"AW parser input excluded {self.power_aw_water_records_omitted} crystallographic water HETATM records using a temporary copy; original PDB unchanged.",),
                                          group_labels=(chains_a, chains_b))
        print(comparison['summary'], end='')
        print(f"Power/AW comparison written to: {root}")

    def _diagnostic_networks(self):
        """Yield each unique named network created by this command."""
        seen = set()
        owners = list(self.sys.groups or [])
        owners.extend(list(getattr(self.sys, 'ifaces', None) or []))
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
