import gc
import json
import os
import re
import shutil
import sys
import tempfile
import tkinter as tk
import zipfile
from contextlib import closing
from copy import deepcopy
from itertools import combinations
from pathlib import Path
from tkinter import filedialog
from time import perf_counter, process_time

from vorpy.src.boundary import (
    BoundaryConfig,
    ION_RESIDUES,
    SolventMode,
    WATER_RESIDUES,
    normalize_solvent_mode,
    resolve_boundary,
)
from vorpy.src.command.command_export import argv_export
from vorpy.src.command.calculation_plan import (
    build_calculation_plan,
    expand_group_tokens,
)
from vorpy.src.command.commands import *
from vorpy.src.command.group import ggroup
from vorpy.src.command.interface import build_interfaces
from vorpy.src.command.interpret import get_file
from vorpy.src.command.set import sett
from vorpy.src.inputs.frames import count_pdb_frames, iter_pdb_frames
from vorpy.src.inputs.source_context import describe_source, source_display_label
from vorpy.src.inputs.net import read_net
from vorpy.src.output.layout import slug, unique_output_directory
from vorpy.src.system.system import System

FRAME_OPTIONS = ('load_commands', 'groups', 'builds', 'exports',
                 'settings_cmnds', 'logs_files', 'interface_mode',
                 'verbose', 'diagnose_edges', 'boundary_config',
                 'boundary_explicitly_requested', 'boundary_warning_printed',
                 'solvent_mode',
                 '_solvent_mode_explicit',
                 'apollonius_requested', 'alpha_value', 'power_settings',
                 'power_vs_aw_requested', 'power_cli_arguments',
                 'aw_alpha_values', 'aw_alpha_pair_overlap_only',
                 'visualize_dual_requested', 'visualize_alpha_value',
                 'execution_profile', 'calculation_plan', 'timing_report_path',
                 '_build_group_networks', '_full_system_separate',
                 '_planned_interface_operands')

EXECUTION_PROFILES = frozenset({'geometry', 'interface', 'analysis', 'full'})
BOUNDARY_PROMPT_CHOICES = (
    ('0', 'No boundary (default; use existing generators only)', 'none'),
    ('1', 'Virtual solvent shell', 'shell'),
)
_BULK_COUNTERION_RESIDUES = frozenset({'NA', 'NA+', 'SOD', 'CL', 'CL-'})

_LEGACY_INTERFACE_RESULT_REASON = (
    'CLI interface geometry did not produce a validated typed Workbench '
    'InterfaceResult; producer provenance and Results dependencies are unresolved.'
)


def _source_composition(system):
    """Return compact chemistry counts from the parsed source table."""
    balls = getattr(system, 'balls', None)
    if balls is None or not hasattr(balls, 'columns'):
        return {
            'atoms': 0, 'molecular_atoms': 0, 'water_atoms': 0, 'ion_atoms': 0,
            'residues': 0, 'waters': 0, 'molecular_chains': 0,
        }

    columns = set(balls.columns)
    names = (
        balls['res_name'].fillna('').astype(str).str.strip().str.upper()
        if 'res_name' in columns else None
    )
    water_mask = names.isin(WATER_RESIDUES) if names is not None else None
    ion_mask = names.isin(ION_RESIDUES) if names is not None else None
    if water_mask is None:
        water_mask = balls.index.to_series().map(lambda _: False)
        ion_mask = balls.index.to_series().map(lambda _: False)
    molecular_mask = ~(water_mask | ion_mask)

    residues = list(getattr(system, 'residues', None) or ())
    solvent = getattr(system, 'sol', None)
    residues.extend(list(getattr(solvent, 'residues', None) or ()))
    water_residues = sum(
        str(getattr(residue, 'name', '')).strip().upper() in WATER_RESIDUES
        for residue in residues
    )

    chain_column = next(
        (
            name for name in ('chain_id', 'chain', 'auth_asym_id', 'chain_name')
            if name in columns
            and balls.loc[molecular_mask, name].fillna('').astype(str).str.strip().ne('').any()
        ),
        None,
    )
    if chain_column is None:
        molecular_chains = len(getattr(system, 'chains', None) or ())
    else:
        chain_values = (
            balls.loc[molecular_mask, chain_column]
            .fillna('')
            .astype(str)
            .str.strip()
        )
        molecular_chains = len({value for value in chain_values.tolist() if value})

    return {
        'atoms': len(balls),
        'molecular_atoms': int(molecular_mask.sum()),
        'water_atoms': int(water_mask.sum()),
        'ion_atoms': int(ion_mask.sum()),
        'residues': len(residues),
        'waters': int(water_residues),
        'molecular_chains': molecular_chains,
    }


def _print_source_summary(system):
    counts = _source_composition(system)
    print()
    print(f'{getattr(system, "name", "Structure")} loaded')
    print(
        f'  {counts["atoms"]:,} atoms, {counts["residues"]:,} residues, '
        f'{counts["molecular_chains"]:,} chains, {counts["waters"]:,} waters'
    )


def _legacy_interface_result_state():
    """Describe the legacy CLI boundary without certifying its geometry."""
    return {
        'status': 'UNRESOLVED',
        'typed_result': False,
        'provenance': False,
        'dependency': False,
        'reason': _LEGACY_INTERFACE_RESULT_REASON,
    }


def _typed_cli_interface_result_state(handoffs, unresolved_count=0):
    """Summarize live CLI handoffs without promoting scientific status."""
    statuses = tuple(
        getattr(getattr(handoff, 'result', None), 'status', 'UNRESOLVED')
        for handoff in handoffs
    )
    status_values = tuple(
        getattr(status, 'value', str(status)) for status in statuses
    )
    status = status_values[0] if len(set(status_values)) == 1 else 'PARTIAL'
    reason = (
        'Typed CLI ordered A/B InterfaceResult handoff; scientific '
        'completeness and curvature certification remain explicit.'
    )
    if unresolved_count:
        reason += f' {unresolved_count} unsupported route(s) remain UNRESOLVED.'
    return {
        'status': status,
        'typed_result': True,
        'provenance': all(
            getattr(getattr(handoff, 'result', None), 'provenance', None) is not None
            for handoff in handoffs
        ),
        'dependency': all(getattr(handoff, 'dependency', None) is not None for handoff in handoffs),
        'reason': reason,
        'statuses': status_values,
    }


def execution_profile_cli_options(args):
    """Extract the opt-in execution profile before legacy parsing."""
    remaining = []
    profile = None
    args = iter(args)
    for arg in args:
        if arg == '--profile':
            if profile is not None:
                raise SystemExit('--profile may only be specified once')
            try:
                profile = str(next(args)).strip().lower()
            except StopIteration:
                raise SystemExit('--profile requires geometry, interface, analysis, or full') from None
            if profile not in EXECUTION_PROFILES:
                raise SystemExit(
                    '--profile requires geometry, interface, analysis, or full'
                )
        else:
            remaining.append(arg)
    return remaining, profile


def cache_export_profile_cli_options(args):
    """Extract the opt-in archive-load/export profile before legacy parsing."""
    remaining = []
    enabled = False
    for arg in args:
        if arg == '--profile-cache-export':
            if enabled:
                raise SystemExit('--profile-cache-export may only be specified once')
            enabled = True
        else:
            remaining.append(arg)
    return remaining, enabled


def fail_if_exists_cli_options(args):
    """Extract the opt-in exclusive output/archive policy."""
    remaining = []
    enabled = False
    for arg in args:
        if arg == '--fail-if-exists':
            if enabled:
                raise SystemExit('--fail-if-exists may only be specified once')
            enabled = True
        else:
            remaining.append(arg)
    return remaining, enabled


def timing_report_cli_options(args):
    """Extract the opt-in machine-readable CLI timing report destination."""
    remaining = []
    report = None
    args = iter(args)
    for arg in args:
        if arg != '--timing-report':
            remaining.append(arg)
            continue
        if report is not None:
            raise SystemExit('--timing-report may only be specified once')
        try:
            value = next(args)
        except StopIteration:
            raise SystemExit('--timing-report requires a JSON output path') from None
        if value.startswith('-'):
            raise SystemExit('--timing-report requires a JSON output path')
        report = Path(value).expanduser().resolve()
    return remaining, report


def _prepare_output_directory(directory, *, fail_if_exists=False):
    """Create an output directory, optionally reserving it atomically."""
    directory = Path(directory).expanduser().resolve()
    try:
        directory.mkdir(parents=True, exist_ok=not fail_if_exists)
    except FileExistsError as error:
        if fail_if_exists:
            raise FileExistsError(
                f'Refusing to reuse existing output directory: {directory}'
            ) from error
        raise
    return directory


def _check_archive_target_available(path):
    """Fail before a solve when an opt-in archive target is already present."""
    path = Path(path).expanduser().resolve()
    timing = path.with_name(path.name + '.timing.json')
    if path.exists() or timing.exists():
        raise SystemExit(
            f'Refusing to overwrite existing archive or timing sidecar: {path}'
        )
    return path


def _archive_member_counts(path):
    """Inspect ZIP metadata only; do not read archive payloads."""
    if path is None or not Path(path).is_file():
        return None
    with zipfile.ZipFile(path, 'r') as archive:
        infos = archive.infolist()
    return {
        'member_count': len(infos),
        'array_member_count': sum(info.filename.startswith('arrays/') for info in infos),
        'table_member_count': sum(info.filename.startswith('tables/') for info in infos),
        'json_member_count': sum(info.filename.endswith('.json') for info in infos),
        'uncompressed_bytes': sum(info.file_size for info in infos),
        'compressed_bytes': sum(info.compress_size for info in infos),
        'archive_bytes': Path(path).stat().st_size,
    }


def _cache_export_state(system):
    """Return a read-only summary of archived interface caches."""
    interfaces = []
    state_counts = {}
    for index, interface in enumerate(getattr(system, 'ifaces', None) or []):
        state = str(getattr(interface, 'archive_cache_state', None) or 'UNKNOWN')
        state_counts[state] = state_counts.get(state, 0) + 1
        caches = getattr(interface, 'representation_caches', None) or {}
        interfaces.append({
            'index': index,
            'interface_id': getattr(interface, 'archive_id', None),
            'network_id': getattr(getattr(interface, 'net', None), 'archive_id', None),
            'archive_cache_state': state,
            'geometry_analysis': getattr(interface, 'geometry_analysis', None) is not None,
            'dual': getattr(interface, '_geometry_dual', None) is not None,
            'filtration': getattr(interface, '_geometry_filtration', None) is not None,
            'representation': getattr(interface, '_geometry_representation', None) is not None,
            'representation_cache_keys': sorted(str(key) for key in caches),
        })
    return {
        'interface_count': len(interfaces),
        'state_counts': state_counts,
        'interfaces': interfaces,
    }


def _atomic_copy(source, destination, *, fail_if_exists=False):
    """Copy one completed artifact without exposing a partial destination."""
    source = Path(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == destination.resolve():
        return
    fd, temporary = tempfile.mkstemp(prefix='.vorpy-copy-', suffix='.tmp', dir=destination.parent)
    os.close(fd)
    try:
        shutil.copyfile(source, temporary)
        if fail_if_exists:
            try:
                os.link(temporary, destination)
            except FileExistsError as error:
                raise FileExistsError(
                    f'Refusing to overwrite existing archive: {destination}'
                ) from error
        else:
            os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _copy_archive_artifacts(source, destination, *, fail_if_exists=False):
    """Reuse an automatic archive and its timing sidecar for explicit output."""
    source = Path(source)
    destination = Path(destination)
    _atomic_copy(source, destination, fail_if_exists=fail_if_exists)
    source_timing = source.with_name(source.name + '.timing.json')
    destination_timing = destination.with_name(destination.name + '.timing.json')
    if not source_timing.is_file():
        return
    try:
        manifest = json.loads(source_timing.read_text(encoding='utf-8'))
        manifest['archive'] = destination.name
        encoded = json.dumps(manifest, allow_nan=False, separators=(',', ':')) + '\n'
        fd, temporary = tempfile.mkstemp(
            prefix='.vorpy-timing-', suffix='.tmp', dir=destination_timing.parent
        )
        try:
            with os.fdopen(fd, 'w', encoding='utf-8', newline='') as stream:
                stream.write(encoded)
            if fail_if_exists:
                if destination_timing.exists():
                    raise FileExistsError(
                        f'Refusing to overwrite existing archive timing sidecar: {destination_timing}'
                    )
                os.link(temporary, destination_timing)
            else:
                os.replace(temporary, destination_timing)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    except (OSError, TypeError, ValueError):
        _atomic_copy(
            source_timing, destination_timing,
            fail_if_exists=fail_if_exists,
        )


def _profile_interval(started):
    return {
        'wall_seconds': perf_counter() - started[0],
        'cpu_seconds': process_time() - started[1],
    }


def _write_cache_export_profile(system, input_archive, profile):
    destination = Path(system.files['dir']) / 'cache_export_profile.json'
    fail_if_exists = getattr(system, '_fail_if_exists', False)
    payload = {
        'format': 'vorpy.cache_export_profile',
        'profile_version': 1,
        'input_archive': Path(input_archive).name,
        'output_directory': str(Path(system.files['dir']).resolve()),
        **profile,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.vorpy-profile-', suffix='.tmp', dir=destination.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', newline='') as stream:
            json.dump(payload, stream, allow_nan=False, indent=2, sort_keys=True)
            stream.write('\n')
        if fail_if_exists:
            try:
                os.link(temporary, destination)
            except FileExistsError as error:
                raise FileExistsError(
                    f'Refusing to overwrite existing cache export profile: {destination}'
                ) from error
        else:
            os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


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


def solvent_cli_options(args):
    """Extract explicit solvent-generator-universe options."""
    remaining = []
    mode = None
    args = iter(args)
    for arg in args:
        if arg == '--solvent-mode':
            if mode is not None:
                raise SystemExit('--solvent-mode may only be specified once')
            try:
                mode = normalize_solvent_mode(next(args))
            except StopIteration:
                raise SystemExit(
                    '--solvent-mode requires MOLECULE_ONLY, SOLVENT_CONTEXT, or FULL_SOLVENT'
                ) from None
            except ValueError as error:
                raise SystemExit(str(error)) from None
        elif arg in {'--solvent-context', '--full-solvent'}:
            if mode is not None:
                raise SystemExit('solvent mode may only be specified once')
            mode = (
                SolventMode.SOLVENT_CONTEXT
                if arg == '--solvent-context'
                else SolventMode.FULL_SOLVENT
            )
        else:
            remaining.append(arg)
    return remaining, mode


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
          "Boundary methods:")
    for key, label, _ in BOUNDARY_PROMPT_CHOICES:
        print(f"    [{key}] {label}")
    print("SES, SAS, and AW probe boundaries are planned and unavailable.")
    modes = {"": "none"}
    for key, _, mode in BOUNDARY_PROMPT_CHOICES:
        modes[key] = mode
        modes[mode] = mode
    while True:
        try:
            answer = input("Choose boundary [0]: ").strip().lower()
        except EOFError:
            answer = ""
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
        self.solvent_mode = SolventMode.MOLECULE_ONLY
        self._solvent_mode_explicit = False
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
        self.execution_profile = None
        self.cache_export_profile = False
        self.fail_if_exists = False
        self.timing_report_path = None
        self._cli_timing_started = None
        self._cli_timing = None
        self.calculation_plan = None
        self._build_group_networks = False
        self._full_system_separate = False
        self._full_system_group = None
        self._planned_interface_operands = ()
        self.interface_result_state = None
        self.interface_results = ()
        self.interface_result_handoffs = ()
        self.interface_result_errors = ()

    def run(self):
        self._cli_timing_started = (perf_counter(), process_time())
        self._run_pipeline()

    def _timing_phase_start(self):
        if self._cli_timing is None:
            return None
        return perf_counter(), process_time()

    def _record_timing_phase(self, name, started, *, source):
        if self._cli_timing is None or started is None:
            return
        interval = _profile_interval(started)
        interval['source'] = source
        self._cli_timing['phases'][name] = interval

    @staticmethod
    def _timing_number(value):
        try:
            number = float(value)
        except (TypeError, ValueError):
            return 0.0
        return number if number >= 0.0 else 0.0

    def _record_network_timings(self):
        """Use FORGE's completed network metrics without a second build timer."""
        if self._cli_timing is None:
            return
        rows = []
        geometry_seconds = 0.0
        analysis_seconds = 0.0
        for name, network in self._diagnostic_networks():
            metrics = getattr(network, 'metrics', {}) or {}
            total = self._timing_number(metrics.get('tot'))
            analysis = min(self._timing_number(metrics.get('anal')), total)
            geometry = max(total - analysis, 0.0)
            geometry_seconds += geometry
            analysis_seconds += analysis
            rows.append({
                'name': str(name),
                'network_seconds': total,
                'geometry_seconds': geometry,
                'analysis_seconds': analysis,
            })
        self._cli_timing['networks'] = rows
        self._cli_timing['phases']['geometry_network'] = {
            'wall_seconds': geometry_seconds,
            'cpu_seconds': None,
            'source': "FORGE Network.metrics['tot'] - ['anal']",
        }
        self._cli_timing['phases']['analysis'] = {
            'wall_seconds': analysis_seconds,
            'cpu_seconds': None,
            'source': "FORGE Network.metrics['anal']",
        }

    def _record_export_archive_timings(self, export_started):
        """Consume VECTOR/ATLAS timings and keep archive time separate."""
        if self._cli_timing is None:
            return
        export_timing = dict(getattr(self.sys, 'export_timing', {}) or {})
        archive_timing = dict(getattr(self.sys, 'archive_timing', {}) or {})
        export_total = self._timing_number(export_timing.get('total'))
        archive_total = self._timing_number(archive_timing.get('total'))
        if export_total:
            archive_in_export = 'network archive' in export_timing
            export_seconds = (
                max(export_total - archive_total, 0.0)
                if archive_in_export else export_total
            )
            self._cli_timing['phases']['export'] = {
                'wall_seconds': export_seconds,
                'cpu_seconds': None,
                'source': 'VECTOR ExportProgress.total',
                'run_exports_seconds': export_total,
            }
            self._cli_timing['export_operations'] = export_timing
        elif export_started is not None:
            self._record_timing_phase(
                'export', export_started, source='PRISM CLI fallback boundary'
            )
        if archive_total:
            self._cli_timing['phases']['archive_save'] = {
                'wall_seconds': archive_total,
                'cpu_seconds': None,
                'source': 'ATLAS save_network.total',
            }
            self._cli_timing['archive_operations'] = archive_timing
        else:
            self._cli_timing['phases']['archive_save'] = {
                'wall_seconds': 0.0,
                'cpu_seconds': None,
                'source': 'ATLAS save_network not requested',
            }

    def _write_timing_report(self, output_directory):
        if self._cli_timing is None or self.timing_report_path is None:
            return
        if self._cli_timing_started is not None:
            total = _profile_interval(self._cli_timing_started)
            total['source'] = 'PRISM CLI run boundary'
            self._cli_timing['phases']['total_wall'] = total
        policy = getattr(self, '_source_generator_policy', None) or {}
        settings = getattr(self, 'settings_dict', None) or {}
        self._cli_timing.update({
            'format': 'vorpy.cli_timing',
            'output_directory': str(Path(output_directory).resolve()),
            'input': str(Path(self.base_file).resolve()) if self.base_file else None,
            'structure': getattr(self.sys, 'name', None),
            'frame': {
                'index': getattr(self.sys, 'frame_index', None),
                'count': getattr(self.sys, 'frame_count', None),
            },
            'network_type': settings.get('net_type'),
            'solvent_mode': getattr(
                getattr(self, 'solvent_mode', None), 'value', self.solvent_mode
            ),
            'source_generator_policy': {
                key: policy.get(key)
                for key in (
                    'source_generator_count', 'active_generator_count',
                    'excluded_generator_count', 'excluded_water_generator_count',
                    'excluded_bulk_counterion_generator_count',
                )
                if key in policy
            },
            'interfaces': [
                getattr(interface, 'name', None)
                for interface in (getattr(self.sys, 'ifaces', None) or [])
            ],
        })
        destination = self.timing_report_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(
            prefix='.vorpy-cli-timing-', suffix='.tmp', dir=destination.parent
        )
        try:
            with os.fdopen(fd, 'w', encoding='utf-8', newline='') as stream:
                json.dump(self._cli_timing, stream, indent=2, sort_keys=True)
                stream.write('\n')
            if self.fail_if_exists:
                try:
                    os.link(temporary, destination)
                except FileExistsError as error:
                    raise SystemExit(
                        f'Refusing to overwrite timing report: {destination}'
                    ) from error
            else:
                os.replace(temporary, destination)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

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
        self._arguments, self.execution_profile = execution_profile_cli_options(
            self._arguments
        )
        self._arguments, self.cache_export_profile = cache_export_profile_cli_options(
            self._arguments
        )
        self._arguments, self.fail_if_exists = fail_if_exists_cli_options(
            self._arguments
        )
        self._arguments, self.timing_report_path = timing_report_cli_options(
            self._arguments
        )
        if self.timing_report_path is not None:
            self._cli_timing = {
                'schema_version': 1,
                'unit': 'seconds',
                'phases': {},
            }
            if self._cli_timing_started is None:
                self._cli_timing_started = (perf_counter(), process_time())
            if self.fail_if_exists and self.timing_report_path.exists():
                raise SystemExit(
                    f'Refusing to overwrite timing report: {self.timing_report_path}'
                )
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
        self._arguments, requested_solvent_mode = solvent_cli_options(self._arguments)
        if requested_solvent_mode is not None:
            self.solvent_mode = requested_solvent_mode
            self._solvent_mode_explicit = True
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
        if self.cache_export_profile and Path(self.base_file).suffix.lower() != '.vpy':
            raise SystemExit('--profile-cache-export requires a .vpy archive input')
        if Path(self.base_file).suffix.lower() == '.vpy':
            self._run_archive()
            return

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
                self.sys = System(
                    file=str(system_input), make_dir=False,
                    boundary_config=self.boundary_config,
                    solvent_mode=self.solvent_mode,
                )
                if comparison_temp is not None:
                    self.sys.name = Path(self.base_file).stem
            else:
                self.sys.boundary_config = self.boundary_config
                self.sys.solvent_mode = self.solvent_mode
            self.sys._fail_if_exists = self.fail_if_exists
        finally:
            if comparison_temp is not None:
                comparison_temp.unlink(missing_ok=True)

        # Parse the remaining command-line arguments
        self.parse_commands()
        self._record_timing_phase(
            'loading_initialization',
            self._cli_timing_started,
            source='PRISM CLI input resolution, System load, and command parsing',
        )
        if self.power_vs_aw_requested:
            self.interface_mode = True
            self.power_cli_arguments = list(self._arguments[1:])

        try:
            self._process_system()
        except BaseException:
            finish_progress = getattr(self.sys, 'finish_progress_line', None)
            if callable(finish_progress):
                finish_progress()
            raise

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
                        solvent_mode=self.solvent_mode,
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
                command.execution_profile = self.execution_profile
                command.fail_if_exists = self.fail_if_exists
                command.visualize_dual_requested = self.visualize_dual_requested
                command.visualize_alpha_value = self.visualize_alpha_value
                command.boundary_config = deepcopy(self.boundary_config)
                command.boundary_explicitly_requested = self.boundary_explicitly_requested
                command.boundary_warning_printed = self.boundary_warning_printed
                command.solvent_mode = self.solvent_mode
                command._solvent_mode_explicit = self._solvent_mode_explicit
                if output_root is None:
                    command.parse_commands()
                    if frame_system.files['dir'] is None:
                        frame_system.set_output_directory(
                            fail_if_exists=self.fail_if_exists
                        )
                    output_root = Path(frame_system.files['dir']) / 'frames'
                    options = {name: deepcopy(getattr(command, name)) for name in FRAME_OPTIONS}
                else:
                    for name, value in options.items():
                        setattr(command, name, deepcopy(value))
                frame_system.set_output_directory(
                    str(output_root / f'frame_{ordinal:04d}'),
                    fail_if_exists=self.fail_if_exists,
                )
                print(f'Processing frame {ordinal}: {frame_system.files["dir"]}')
                try:
                    command._process_system()
                except BaseException:
                    finish_progress = getattr(frame_system, 'finish_progress_line', None)
                    if callable(finish_progress):
                        finish_progress()
                    raise
                self.solvent_mode = command.solvent_mode
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

    @staticmethod
    def _interface_selector_key(selector):
        return tuple(token.lower() for token in selector)

    def _planned_interface_pairs(self):
        """Resolve planned selectors to groups without changing geometry."""
        plan = self.calculation_plan
        groups = list(self.sys.groups or [])
        if plan is None or not plan.interface_requests:
            if len(groups) >= 2:
                return list(combinations(groups, 2))
            return build_interfaces(
                sys=self.sys,
                num_requested_groups=len(self.groups),
            )

        # Group definitions are kept in the same order as the legacy parser.
        by_selector = {}
        for index, commands in self.groups.items():
            key = tuple(
                tuple(str(token).lower() for token in command)
                for command in commands
            )
            if index < len(groups):
                by_selector[key] = groups[index]

        pairs = []
        seen = set()
        for request in plan.interface_requests:
            if not request.operands:
                candidates = groups
            else:
                candidates = []
                for operand in request.operands:
                    key = (self._interface_selector_key(operand),)
                    group = by_selector.get(key)
                    if group is None:
                        raise ValueError(
                            'Interface selector did not resolve to a declared group: '
                            + ' '.join(operand)
                        )
                    if group not in candidates:
                        candidates.append(group)

            if len(request.operands) == 1:
                requested_pairs = [(candidates[0], None)]
            elif len(request.operands) >= 2:
                requested_pairs = list(combinations(candidates, 2))
            elif len(candidates) >= 2:
                requested_pairs = list(combinations(candidates, 2))
            else:
                requested_pairs = [(candidates[0], None)]

            for group1, group2 in requested_pairs:
                marker = (id(group1), id(group2) if group2 is not None else None)
                if marker not in seen:
                    seen.add(marker)
                    pairs.append((group1, group2))
        return pairs

    def _build_full_system_target(self):
        """Build the one full-system source network requested by the plan."""
        plan = self.calculation_plan
        full_system_requested = bool(
            self._full_system_separate
            or (plan is not None and plan.full_system_request)
            or (plan is not None and plan.interface_requests)
        )
        if not full_system_requested or self._full_system_group is not None:
            return
        self._resolve_full_system_solvent_mode()
        from vorpy.src.group import Group

        settings = deepcopy(self.settings_dict) or {}
        if self.interface_mode or (plan is not None and plan.interface_requests):
            settings['_progress_network_name'] = self.sys.name
            # Interface extraction consumes the same completed source Network
            # as a regular group solve; keep the public phase terminology
            # identical to the normal CLI path.
            settings['_completion_kind'] = 'network'
        source_indices = getattr(self, '_source_system_indices', None)
        if source_indices is None:
            source_indices = tuple(range(len(self.sys.balls)))
        self._full_system_group = Group(
            sys=self.sys,
            name=f'{self.sys.name}_full_system',
            atoms=list(source_indices),
            settings=settings,
            make_net=True,
            build_net=False,
            print_metrics=False,
        )
        calculate_curvature = getattr(
            self.sys, '_network_calculate_curvature', True
        )
        self._full_system_group.build(calculate_curvature=calculate_curvature)
        source_network = getattr(self._full_system_group, 'net', None)
        if source_network is not None:
            self._record_generator_policy(source_network)
            self._persist_geometry_checkpoint(
                source_network,
                scope='full_system_source',
            )

    def _full_system_reuse_enabled(self):
        """Whether selected targets must be read-only views of the source solve."""
        plan = self.calculation_plan
        source = getattr(self._full_system_group, 'net', None)
        return bool(
            (self._full_system_separate
             or (plan is not None and plan.full_system_request)
             or (plan is not None and plan.interface_requests))
            and source is not None
            and getattr(source, 'verts', None) is not None
            and getattr(source, 'edges', None) is not None
            and getattr(source, 'surfs', None) is not None
        )

    def _persist_geometry_checkpoint(self, network, *, scope):
        """Persist solved geometry before optional analysis can fail."""
        output_dir = (getattr(self.sys, 'files', None) or {}).get('dir')
        if output_dir is None:
            return None
        from vorpy.src.io import save_network

        internal_interface_checkpoint = bool(
            self.interface_mode
            and self.execution_profile not in {'geometry', 'interface'}
        )
        checkpoint_dir = Path(output_dir)
        if internal_interface_checkpoint:
            checkpoint_dir = checkpoint_dir / '.vorpy' / 'checkpoints'
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        provenance = dict(getattr(network, 'archive_provenance', {}) or {})
        provenance.update({
            'checkpoint_kind': 'solved_geometry',
            'checkpoint_scope': scope,
            'geometry_status': 'GEOMETRY_COMPLETE',
            'analysis_status': 'NOT_CALCULATED',
            'analysis_complete': False,
        })
        network.archive_provenance = provenance
        destination = checkpoint_dir / (
            f'{getattr(network, "group_name", None) or getattr(self.sys, "name", "network")}'
            '.geometry.vpy'
        )
        archive = save_network(
            network, destination, fail_if_exists=getattr(self, 'fail_if_exists', False)
        )
        self.sys._geometry_checkpoint_path = str(archive)
        self.sys._geometry_checkpoint_status = 'GEOMETRY_COMPLETE'
        self.sys._analysis_status = 'NOT_CALCULATED'
        if not internal_interface_checkpoint or self.verbose:
            print(f'Solved geometry checkpoint saved: {archive}')
        return archive

    def _full_system_source_network_id(self):
        """Assign one producer namespace shared by every source-derived view."""
        source = self._full_system_group.net
        source_network_id = getattr(source, 'source_network_id', None)
        if isinstance(source_network_id, str) and source_network_id.strip():
            return source_network_id
        from vorpy.src.results.model import stable_id

        context = source.get_solve_context_key()
        source_network_id = stable_id(
            'producer:cli-full-system',
            {
                'source': str(self.base_file or getattr(self.sys, 'name', 'system')),
                'frame': int(getattr(self.sys, 'frame_index', 1) or 1),
                'solve_universe_id': context.solve_universe_id,
            },
        )
        source.source_network_id = source_network_id
        return source_network_id

    def _attach_full_system_group_views(self):
        """Attach read-only group views without creating subgroup Networks."""
        plan = self.calculation_plan
        # Interface-only commands use groups as interface sides.  Validating
        # whole group views here would pull unrelated incomplete cells into the
        # dependency check before the interface facets are selected.
        if plan is None or not plan.full_system_request:
            return
        if not self._full_system_reuse_enabled():
            return
        from vorpy.src.group.shared_view import analyze_group_view, make_group_view

        source = self._full_system_group.net
        source_network_id = self._full_system_source_network_id()
        for group in self.sys.groups or ():
            view = make_group_view(
                group, source, source_network_id=source_network_id
            )
            group.source_view = view
            group.source_network = source
            group.source_view_analysis = analyze_group_view(view)
            # Keep existing scalar consumers useful while leaving curvature
            # and subgroup topology explicitly unavailable on the Group.
            group.sa = group.source_view_analysis.area
            group.vol = group.source_view_analysis.volume

    def _source_interface_analysis(self, view, interface_id):
        """Summarize source facets; do not run alpha/dual/interface geometry."""
        from vorpy.src.calculations import calc_surf_sa
        from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis

        group_a = set(view.group_a_ids)
        group_b = set(view.group_b_ids)
        facet_pairs = set()
        area = 0.0
        for surface in view.surfaces:
            area_value = surface.data.get('sa')
            if area_value is None:
                try:
                    area_value = calc_surf_sa(
                        tris=surface.data['tris'], points=surface.data['points']
                    )
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError(
                        f'Source surface {surface.source_id!r} has no analytic area'
                    ) from error
            area += float(area_value)
            side_a = group_a.intersection(surface.generator_system_ids)
            side_b = group_b.intersection(surface.generator_system_ids)
            facet_pairs.update((int(a), int(b)) for a in side_a for b in side_b)

        reason = (
            'Source-derived A/B facets from the solved full-system network; '
            'no independent interface solve or alpha selection was performed.'
        )
        unresolved_reasons = [reason, *view.completeness.reasons]
        area_scope = (
            'supported full-system source A/B facets; unresolved source facets '
            'were omitted'
            if view.completeness.reasons
            else 'complete full-system source A/B facets'
        )
        side = {
            'interface_atoms': len({
                index for surface in view.surfaces
                for index in surface.generator_system_ids
                if index is not None
            }),
            'area': area,
            'area_scope': area_scope,
            'area_selection_complete': bool(view.completeness.proven),
            'selected_pairs': len(view.surfaces),
            'selected_physical_pairs': len(view.surfaces),
            'topology_faces': len(view.surfaces),
            'topology_edges': len(view.edges),
            'topology_vertices': len(view.vertices),
            'mean_certified': False,
            'gaussian_total_certified': False,
            'mean_curvature_certified': False,
            'gaussian_curvature_certified': False,
            'manifold_incidence_certified': False,
            'unresolved_reasons': unresolved_reasons,
            'unresolved_surface_ids': view.completeness.unresolved_surface_ids,
            'unresolved_surface_area': view.completeness.unresolved_surface_area,
            'source_surface_ids': tuple(view.source_ids('surface')),
            'source_edge_ids': tuple(view.source_ids('edge')),
            'source_vertex_ids': tuple(view.source_ids('vertex')),
            'source_generator_pairs': tuple(sorted(facet_pairs)),
        }
        metadata = {
            'schema_version': 2,
            'interface_id': interface_id,
            'scheme': view.scheme,
            'alpha_convention': 'not_applicable_source_facets',
            'alpha_value': None,
            'alpha_units': None,
            'curvature_requested': False,
            'additional_network_solves': 0,
            'source_geometry': 'full_system_network_view',
            'source_network_id': view.source_network_id,
            'source_view_id': view.selection_id,
            'solve_universe_id': view.solve_universe_id,
            'geometry_status': 'GEOMETRY_COMPLETE',
            'analysis_status': 'NOT_CALCULATED',
            'analysis_complete': False,
            'alpha_status': 'NOT_CALCULATED',
            'raw_physical_geometry_available': True,
            'result_status': 'PARTIAL',
            'unresolved_surface_ids': view.completeness.unresolved_surface_ids,
            'unresolved_surface_area': view.completeness.unresolved_surface_area,
            'solvent_mode': getattr(
                self._full_system_group.net.solvent_mode,
                'value',
                self._full_system_group.net.solvent_mode,
            ),
        }
        return InterfaceGeometryAnalysis(
            metadata=metadata,
            alpha_selection={
                'certified': False,
                'status': 'not_calculated',
                'coverage': None,
                'method': 'not_applicable_source_facets',
                'source_surface_ids': tuple(view.source_ids('surface')),
                'source_generator_pairs': tuple(sorted(facet_pairs)),
            },
            voronoi_side_1=side,
            voronoi_side_2=dict(side),
            unresolved=unresolved_reasons,
        )

    def _attach_full_system_interface_views(self, interface_pairs):
        """Create interface records backed by the full source network only."""
        if not self._full_system_reuse_enabled():
            return
        from vorpy.src.interface import Interface
        from vorpy.src.network.views import NetworkViewError

        source = self._full_system_group.net
        source_network_id = self._full_system_source_network_id()
        self.sys.ifaces = []
        for group_a, group_b in interface_pairs:
            interface = Interface(
                sys=self.sys,
                group1=group_a,
                group2=group_b,
            )
            actual_group_b = interface.group2
            try:
                view = source.make_view(
                    iface_grps=(
                        tuple(sorted(int(index) for index in group_a.ball_ndxs)),
                        tuple(sorted(int(index) for index in actual_group_b.ball_ndxs)),
                    ),
                    source_network_id=source_network_id,
                )
            except NetworkViewError as error:
                # The solved source network is still a valid completed
                # checkpoint; only this requested derived interface failed.
                self.interface_result_errors += (
                    f'{type(error).__name__}: {error}',
                )
                continue
            # Interface.net remains the source owner for archive/link identity;
            # all selected geometry is carried by source_view, never rebuilt.
            interface.net = source
            interface.source_view = view
            interface.source_network = source
            interface.geometry_source = 'FULL_SYSTEM_NETWORK_VIEW'
            interface.group1_indices = set(group_a.ball_ndxs)
            interface.group2_indices = set(actual_group_b.ball_ndxs)
            interface.group2_name = actual_group_b.name
            interface.geometry_analysis = self._source_interface_analysis(
                view, interface.interface_id
            )
            interface._update_group_metadata(network_created=False, built=True)
            for group in (group_a, actual_group_b):
                metadata = group.interface_metadata[interface.interface_id]
                metadata.update({
                    'geometry_source': interface.geometry_source,
                    'source_network_id': view.source_network_id,
                    'source_view_id': view.selection_id,
                    'completeness_proven': bool(view.completeness.proven),
                })
            self.sys.ifaces.append(interface)

    @staticmethod
    def _interface_group_label(group):
        label = getattr(group, 'name', None) or getattr(group, 'group_id', None)
        if isinstance(label, str) and label.strip():
            return label
        return type(group).__name__

    def _prepare_interface_output_directories(self):
        root = Path(self.sys.files['dir']).resolve()
        occupied = []
        for interface in getattr(self.sys, 'ifaces', ()) or ():
            name = slug(getattr(interface, 'name', 'interface'))
            directory = root / name
            if directory.exists():
                if self.fail_if_exists:
                    raise SystemExit(
                        f'Refusing to reuse existing interface output directory: '
                        f'{directory}'
                    )
                directory = unique_output_directory(root, name, occupied)
            directory.mkdir(parents=True, exist_ok=False)
            interface.dir = str(directory)
            occupied.append(directory)

    @staticmethod
    def _interface_arrow():
        arrow = '↔'
        encoding = getattr(sys.stdout, 'encoding', None)
        if encoding:
            try:
                arrow.encode(encoding)
            except UnicodeEncodeError:
                return '<->'
        return arrow

    def _print_interface_statistics(self):
        """Print compact source/interface counts without object reprs."""
        for interface in getattr(self.sys, 'ifaces', ()) or ():
            view = getattr(interface, 'source_view', None)
            analysis = getattr(interface, 'geometry_analysis', None)
            side = getattr(analysis, 'voronoi_side_1', {}) or {}
            surfaces = side.get('topology_faces')
            edges = side.get('topology_edges')
            vertices = side.get('topology_vertices')
            area = side.get('area')
            if view is not None:
                surfaces = len(view.surfaces) if surfaces is None else surfaces
                edges = len(view.edges) if edges is None else edges
                vertices = len(view.vertices) if vertices is None else vertices
            if surfaces is None or edges is None or vertices is None:
                continue
            try:
                unit = 'A²'
                encoding = getattr(sys.stdout, 'encoding', None)
                if encoding:
                    try:
                        unit.encode(encoding)
                    except UnicodeEncodeError:
                        unit = 'A^2'
                area_text = f', area {float(area):,.3f} {unit}' if area is not None else ''
            except (TypeError, ValueError):
                area_text = ''
            completeness = getattr(view, 'completeness', None)
            if completeness is not None:
                completeness_text = f' - geometry {completeness.status.lower()}'
            else:
                network = getattr(interface, 'net', None)
                balls = getattr(network, 'balls', None)
                complete = None
                if balls is not None and hasattr(balls, 'columns') and 'complete' in balls:
                    complete = int(balls['complete'].astype(bool).sum())
                    completeness_text = f' - geometry {complete:,}/{len(balls):,} cells'
                else:
                    completeness_text = ''
            print(
                f'"{self._interface_group_label(interface)}" built - '
                f'{int(surfaces):,} surfaces, {int(edges):,} edges, '
                f'{int(vertices):,} vertices{area_text}{completeness_text}'
            )

    def _resolve_generator_policy(self):
        """Resolve source membership without treating a wet PDB as a box."""
        if self.sys is None:
            return ()

        previous = getattr(self, '_source_generator_policy', None)
        if previous is not None:
            return tuple(previous.get('excluded_system_indices', ()))

        from vorpy.src.boundary import prepare_network_geometry

        mode = normalize_solvent_mode(self.solvent_mode)
        balls = getattr(self.sys, 'balls', None)
        if balls is None:
            return ()
        rows = [
            balls.iloc[index] if hasattr(balls, 'iloc') else balls[index]
            for index in range(len(balls))
        ]
        def residue_name(row):
            return str(
                row.get('res_name', '') if hasattr(row, 'get')
                else getattr(row, 'res_name', '')
            ).strip().upper()

        retained = set(getattr(self.sys, 'retained_system_indices', ()) or ())
        structural_ions = set()
        if mode is SolventMode.MOLECULE_ONLY:
            for index, row in enumerate(rows):
                residue = residue_name(row)
                if residue in ION_RESIDUES and residue not in _BULK_COUNTERION_RESIDUES:
                    structural_ions.add(index)
            retained.update(structural_ions)

        self.sys.retained_system_indices = tuple(sorted(retained))
        geometry = prepare_network_geometry(
            self.sys,
            solvent_mode=mode,
            retained_system_indices=self.sys.retained_system_indices,
        )
        active = tuple(
            int(index) for index in geometry.system_indices if index is not None
        )
        active_set = set(active)
        excluded = tuple(index for index in range(len(rows)) if index not in active_set)
        excluded_water = tuple(
            index for index in excluded
            if residue_name(rows[index]) in WATER_RESIDUES
        )
        excluded_bulk_ions = tuple(
            index for index in excluded
            if residue_name(rows[index]) in _BULK_COUNTERION_RESIDUES
        )
        source = describe_source(self.base_file)
        source_provenance = {
            'path': str(Path(self.base_file).resolve()) if self.base_file else None,
            'dataset_id': getattr(source, 'dataset_id', None),
            'role': getattr(source, 'role', None),
            'evidence': getattr(source, 'evidence', None),
        }
        stable_by_system = {
            int(system_index): str(stable_id)
            for stable_id, system_index in zip(
                geometry.stable_ids, geometry.system_indices
            )
            if system_index is not None
        }
        interface_source = bool(
            self.interface_mode
            or (self.calculation_plan is not None and self.calculation_plan.interface_requests)
        )
        resolution = (
            'explicit' if self._solvent_mode_explicit
            else 'implicit_interface_full_solvent'
            if interface_source and mode is SolventMode.FULL_SOLVENT
            else 'implicit_default'
        )
        policy = {
            'name': 'molecular_only' if mode is SolventMode.MOLECULE_ONLY else mode.value.lower(),
            'solvent_mode': mode.value,
            'solvent_mode_resolution': resolution,
            'source_provenance': source_provenance,
            'source_generator_count': len(rows),
            'active_generator_count': len(active),
            'excluded_generator_count': len(excluded),
            'excluded_water_generator_count': len(excluded_water),
            'excluded_bulk_counterion_generator_count': len(excluded_bulk_ions),
            'excluded_system_indices': excluded,
            'retained_system_indices': tuple(sorted(retained)),
            'retained_structural_ion_indices': tuple(sorted(structural_ions)),
            'retained_structural_ion_ids': tuple(
                stable_by_system[index] for index in sorted(structural_ions)
                if index in stable_by_system
            ),
            'active_system_indices': active,
        }
        self._source_generator_policy = policy
        self._source_system_indices = active
        self.sys._source_generator_policy = policy
        self.sys.solvent_mode = mode
        self.sys.solvent_mode_resolution = resolution
        if self.settings_dict is None:
            self.settings_dict = {}
        self.settings_dict['solvent_mode'] = mode.value
        self.settings_dict['solvent_mode_resolution'] = resolution
        self.settings_dict['generator_policy'] = policy
        return excluded

    def _record_generator_policy(self, network):
        """Persist source membership evidence in settings and archive provenance."""
        policy = getattr(self, '_source_generator_policy', None)
        if policy is None or network is None:
            return
        if not isinstance(getattr(network, 'settings', None), dict):
            network.settings = {}
        network.settings['generator_policy'] = policy
        network.settings['solvent_mode_resolution'] = policy['solvent_mode_resolution']
        for key in (
            'source_generator_count', 'active_generator_count',
            'excluded_generator_count', 'excluded_water_generator_count',
            'excluded_bulk_counterion_generator_count',
        ):
            network.settings[key] = policy[key]
        provenance = dict(getattr(network, 'archive_provenance', {}) or {})
        provenance['source_generator_policy'] = policy
        network.archive_provenance = provenance

    def _resolve_full_system_solvent_mode(self):
        """Resolve the full-system source membership before Group creation."""
        plan = self.calculation_plan
        if not (
            plan is not None
            and (self._full_system_separate or plan.full_system_request or plan.interface_requests)
            and self.sys is not None
        ):
            return ()
        return self._resolve_generator_policy()

    def _resolve_default_solvent_mode(self):
        """Use all wet-input generators for implicit interface source solves."""
        if self.sys is None:
            return False
        plan = self.calculation_plan
        interface_source = bool(
            self.interface_mode
            or (plan is not None and plan.interface_requests)
        )
        if interface_source and not self._solvent_mode_explicit:
            self.solvent_mode = SolventMode.FULL_SOLVENT
            self.sys.solvent_mode = self.solvent_mode
        self._resolve_generator_policy()
        return True

    def _make_cli_interface_results(self, interface_pairs):
        """Adapt explicit live CLI interfaces through the Workbench seam."""
        self.interface_results = ()
        self.interface_result_handoffs = ()
        # Preserve extraction failures recorded before typed-result adaptation.
        interface_result_errors = tuple(
            getattr(self, 'interface_result_errors', ()) or ()
        )
        self.interface_result_errors = interface_result_errors
        self.sys.cli_interface_results = ()
        self.sys.cli_interface_result_handoffs = ()
        interfaces = list(getattr(self.sys, 'ifaces', ()) or ())
        if not interface_pairs or len(interfaces) < len(interface_pairs):
            return ()

        from vorpy.src.network.views import NetworkViewError
        from vorpy.src.results.model import stable_id
        from vorpy.workbench.services.vorpy_backend import (
            InterfaceSolveSelection,
            VorPyBackend,
        )

        source = Path(self.base_file) if self.base_file else None
        system_id = (
            source.stem
            if source is not None and source.stem
            else str(getattr(self.sys, 'name', None) or 'cli-system')
        )
        source_path = str(source.resolve()) if source is not None else system_id
        frame_index = int(getattr(self.sys, 'frame_index', 1) or 1)
        handoffs = []

        # System.make_interfaces appends in request order; only pair entries
        # with an explicit second group are eligible for ordered A/B Results.
        for (group_a, group_b), interface in zip(
            interface_pairs, interfaces[-len(interface_pairs):]
        ):
            if group_b is None:
                continue
            try:
                group_a_id = str(
                    getattr(group_a, 'group_id', None)
                    or getattr(group_a, 'name', '')
                )
                group_b_id = str(
                    getattr(group_b, 'group_id', None)
                    or getattr(group_b, 'name', '')
                )
                selection = InterfaceSolveSelection(
                    interface_id=str(interface.interface_id),
                    group_a_id=group_a_id,
                    group_a_system_indices=tuple(sorted(int(index) for index in group_a.ball_ndxs)),
                    group_b_id=group_b_id,
                    group_b_system_indices=tuple(sorted(int(index) for index in group_b.ball_ndxs)),
                )
                network = interface.net
                context = network.get_solve_context_key()
                if context.scheme not in {'aw', 'pow'}:
                    continue

                source_view = getattr(interface, 'source_view', None)
                source_network_id = (
                    getattr(source_view, 'source_network_id', None)
                    if source_view is not None
                    else getattr(network, 'source_network_id', None)
                )
                if not isinstance(source_network_id, str) or not source_network_id.strip():
                    source_network_id = stable_id(
                        'producer:cli-interface',
                        {
                            'source': source_path,
                            'frame': frame_index,
                            'interface_id': selection.interface_id,
                            'group_A_system_indices': selection.group_a_system_indices,
                            'group_B_system_indices': selection.group_b_system_indices,
                            'solve_universe_id': context.solve_universe_id,
                        },
                    )
                    network.source_network_id = source_network_id

                source_group = self._full_system_group
                source_network = getattr(source_group, 'net', None)
                if source_network is None:
                    raise ValueError(
                        'Validated full-system source network is required for '
                        'ordered CLI InterfaceResult handoff'
                    )
                if getattr(source_network, 'source_network_id', None) is None:
                    source_network.source_network_id = source_network_id
                backend = VorPyBackend()
                if source_view is not None:
                    view = source_view
                else:
                    view = backend.make_validated_interface_view(
                        interface,
                        selection,
                        source_network,
                        source_network_id=source_network_id,
                    )
                result_network_id = getattr(network, 'archive_id', None) or stable_id(
                    'archive:cli-interface',
                    {
                        'source': source_path,
                        'frame': frame_index,
                        'interface_id': selection.interface_id,
                        'solve_universe_id': context.solve_universe_id,
                    },
                )
                if result_network_id == source_network_id:
                    result_network_id = stable_id(
                        'archive:cli-interface',
                        {'source_network_id': source_network_id, 'result': True},
                    )
                solvent_mode = getattr(network.solvent_mode, 'value', network.solvent_mode)
                environment = (
                    'solvent_competing'
                    if solvent_mode != SolventMode.MOLECULE_ONLY.value
                    else 'dry'
                )
                handoff = backend.make_interface_result_handoff(
                    interface,
                    selection,
                    interface.geometry_analysis,
                    system_id=system_id,
                    result_network_id=result_network_id,
                    radius_configuration_id='vorpy-backend-radii-v1',
                    selection_id=view.selection_id,
                    group_a_selection_id=stable_id(
                        'view',
                        {
                            'group_id': selection.group_a_id,
                            'system_ids': selection.group_a_system_indices,
                        },
                    ),
                    group_b_selection_id=stable_id(
                        'view',
                        {
                            'group_id': selection.group_b_id,
                            'system_ids': selection.group_b_system_indices,
                        },
                    ),
                    orientation='A -> B',
                    source_identity_status='AUTHORITATIVE',
                    source_revision=f'cli-interface:{context.solve_universe_id}',
                    environment=environment,
                    # A source-derived view is intentionally distinct from an
                    # independently solved interface; its lineage stays
                    # MATERIALIZED_VIEW and its curvature remains uncertified.
                    completeness_proven=bool(view.completeness.proven),
                )
            except (
                AttributeError,
                IndexError,
                KeyError,
                NetworkViewError,
                TypeError,
                ValueError,
            ) as error:
                # Legacy and unsupported routes remain visible as unresolved.
                self.interface_result_errors += (
                    f'{type(error).__name__}: {error}',
                )
                continue
            interface.typed_results = (handoff.result,)
            interface.results_dependency = handoff.dependency
            interface.geometry_lineage = handoff.geometry_lineage
            handoffs.append(handoff)

        self.interface_result_handoffs = tuple(handoffs)
        self.interface_results = tuple(handoff.result for handoff in handoffs)
        self.sys.cli_interface_results = self.interface_results
        self.sys.cli_interface_result_handoffs = self.interface_result_handoffs
        return self.interface_result_handoffs

    def _process_system(self):
        """Build and export a loaded system using the parsed commands."""
        # Select the default only after export-directory commands have been applied.
        if self.sys.files['dir'] is None:
            self.sys.set_output_directory(fail_if_exists=self.fail_if_exists)

        self.sys._fail_if_exists = self.fail_if_exists

        # Expose the CLI verbosity state at the system level so downstream
        # build/analyze/export code has one common place to query it.
        self.sys.verbose = self.verbose
        self.sys._write_archive_timing_manifest = bool(
            self.verbose or self.diagnose_edges or self.execution_profile is not None
            or self.cache_export_profile or self.timing_report_path is not None
        )
        interface_export_names = {
            str(command[0]).strip().lower()
            for command in self.exports if command
        }
        self.sys._interface_water_requested = bool(
            interface_export_names & {
                'water', 'waters', 'interface_water',
                'tiny', 'small', 'med', 'medium', 'large', 'all',
            }
        )
        self.sys._interface_curvature_requested = bool(
            interface_export_names & {
                'curvature', 'curv', 'diagnostic', 'diagnostics',
                'logs', 'dual', 'apollonius', 'small', 'med', 'medium',
                'large', 'all',
            }
        )
        interface_export_modifiers = {
            'dir', 'directory', 'ft', 'file_type', 'file_format',
            'format', 'mesh_format', 'no_archive', 'no-archive',
            'no_network', 'no-network',
        }
        self.sys._interface_default_export = bool(
            self.interface_mode
            and not any(name not in interface_export_modifiers
                        for name in interface_export_names)
        )
        self.sys._compact_interface_workflow = bool(
            self.interface_mode and (
                not self.exports or 'visualize' in interface_export_names
            )
        )
        self.sys._execution_profile = self.execution_profile
        if self.execution_profile in {'geometry', 'interface'}:
            self.sys._interface_curvature_requested = False
            self.sys._interface_water_requested = False
            self.sys._network_calculate_curvature = False
        elif self.execution_profile == 'analysis':
            self.sys._interface_curvature_requested = True
            self.sys._network_calculate_curvature = True
        elif self.execution_profile == 'full':
            self.sys._interface_curvature_requested = True
            self.sys._interface_water_requested = True
            self.sys._network_calculate_curvature = True
        else:
            self.sys._network_calculate_curvature = True

        # Load any additional files
        self.load_files()

        _print_source_summary(self.sys)

        if getattr(self.sys, 'frame_count', 0) > 1:
            if self.verbose:
                if Path(self.base_file).suffix.lower() == '.pdb':
                    print(
                        'Input context: '
                        + source_display_label(
                            self.base_file,
                            frame_index=getattr(self.sys, 'frame_index', 1),
                            frame_count=getattr(self.sys, 'frame_count', 1),
                        )
                    )
                print(f'{self.sys.name}: {self.sys.frame_count} frames found; '
                      f'{self.sys.loaded_frame_count} frame loaded '
                      f'(frame {self.sys.frame_index} only), {len(self.sys.balls):,} atoms.')
            else:
                print(f'  Frame: {self.sys.frame_index}/{self.sys.frame_count} loaded')

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

        self._resolve_default_solvent_mode()
        self._resolve_full_system_solvent_mode()

        if not self.verbose:
            net_type = (self.settings_dict or {}).get('net_type', 'aw')
            if isinstance(net_type, str):
                print(f'\nSolving {net_type.upper()} network...', flush=True)

        self._build_full_system_target()

        # Create the requested groups after the source solve exists.  In a
        # full-system plan these are selection containers, not new Networks.
        self.create_groups()
        # Interface-only requests need direct A/B source facets, not whole
        # selected-group cell snapshots. Avoid rejecting a valid interface
        # because an unrelated selected cell is incomplete.
        plan = self.calculation_plan
        interface_only = bool(
            self.interface_mode
            or (plan is not None and plan.interface_requests)
        )
        if not interface_only:
            self._attach_full_system_group_views()

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

                calculate_curvature = getattr(
                    self.sys, '_network_calculate_curvature', True
                )
                copy_group.build(calculate_curvature=calculate_curvature)
                grp.build(calculate_curvature=calculate_curvature)

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

            interface_pairs = self._planned_interface_pairs()
            interface_timing_started = None

            if self._build_group_networks and not self._full_system_reuse_enabled():
                calculate_curvature = getattr(
                    self.sys, '_network_calculate_curvature', True
                )
                for group in self.sys.groups or []:
                    if getattr(group, 'net', None) is None:
                        group.make_net()
                    group.build(calculate_curvature=calculate_curvature)

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

            if self._full_system_reuse_enabled():
                interface_timing_started = self._timing_phase_start()
                self._attach_full_system_interface_views(interface_pairs)
            else:
                self.sys.make_interfaces(
                    interface_pairs,
                    analyze_waters=getattr(self.sys, '_interface_water_requested', False),
                )
                interface_timing_started = self._timing_phase_start()
            if not self.verbose:
                print()
            self._print_interface_statistics()
            handoffs = self._make_cli_interface_results(interface_pairs)
            if handoffs:
                self.interface_result_state = _typed_cli_interface_result_state(
                    handoffs,
                    unresolved_count=len(interface_pairs) - len(handoffs),
                )
            else:
                self.interface_result_state = _legacy_interface_result_state()
                if self.interface_result_errors:
                    self.interface_result_state['reason'] += (
                        ' Handoff diagnostics: '
                        + ' | '.join(self.interface_result_errors)
                    )
                    print(
                        'Warning: interface extraction failed: '
                        + self.interface_result_errors[0]
                        + '; see interface_info.txt'
                    )
                else:
                    print(
                        'Warning: interface Results unresolved; '
                        'see interface_info.txt'
                    )
            self.sys.cli_interface_result_state = self.interface_result_state
            self._record_timing_phase(
                'interface_extraction',
                interface_timing_started,
                source='PRISM source-view selection and typed handoff',
            )

        # Standard full-group build
        else:
            for grp in self.sys.groups:
                if self._full_system_reuse_enabled():
                    continue
                if self.execution_profile is None:
                    grp.build()
                else:
                    grp.build(
                        calculate_curvature=getattr(
                            self.sys, '_network_calculate_curvature', True
                        )
                    )

        if self._full_system_group is not None and self._full_system_group not in self.sys.groups:
            self.sys.groups.append(self._full_system_group)

        for _, network in self._diagnostic_networks():
            self._record_generator_policy(network)
        self._record_network_timings()

        # Run optional read-only diagnostics after construction and before
        # export so the report describes the completed in-memory networks.
        optional_analysis = []
        if self.diagnose_edges:
            optional_analysis.append(('AW edge diagnostics', self.run_edge_diagnostics))
        if self.apollonius_requested:
            optional_analysis.append(('Apollonius analysis', self.run_apollonius))
        if self.visualize_dual_requested:
            optional_analysis.append(('dual visualization', self.run_dual_visualization))
        if self.aw_alpha_values:
            optional_analysis.append(('AW alpha analysis', self.run_aw_alpha_experimental))

        for name, action in optional_analysis:
            try:
                action()
            except (Exception, SystemExit) as error:
                self._record_optional_analysis_failure(name, error)

        if self.power_vs_aw_requested:
            self.run_power_aw_comparison()
            # This analysis mode writes its own compact Power/AW artifacts.
            # Avoid the legacy default bulk exporter, which is unrelated to
            # the analysis output and currently does not support Interface's
            # selected export keyword set.
            return

        if self.interface_mode:
            self._prepare_interface_output_directories()

        # Export requested outputs
        self.sys._last_export_archive_path = None
        self.sys._last_export_archive_output_path = None
        export_timing_started = self._timing_phase_start()
        if not self.verbose:
            print('\nExporting results...', flush=True)
        self.run_exports()
        if self.verbose:
            print(flush=True)
        self._save_archive()
        output_directory = self.sys.files['dir']
        if self.interface_mode and len(getattr(self.sys, 'ifaces', ()) or ()) == 1:
            output_directory = self.sys.ifaces[0].dir or output_directory
        output_directory = os.path.normpath(os.path.abspath(str(output_directory)))
        self._record_export_archive_timings(export_timing_started)
        self._write_timing_report(output_directory)
        print(f'Files exported to: {output_directory}', flush=True)
        print(flush=True)

    def _save_archive(self):
        if self.save_network_path is not None:
            destination = self.save_network_path
            if getattr(self.sys, 'frame_count', 1) > 1:
                destination = destination.with_name(f'{destination.stem}_frame_{self.sys.frame_index:04d}.vpy')
            automatic = getattr(self.sys, '_last_export_archive_path', None)
            if automatic and Path(automatic).is_file():
                _copy_archive_artifacts(
                    automatic, destination,
                    fail_if_exists=self.fail_if_exists,
                )
            else:
                from vorpy.src.io import save_network
                if self.fail_if_exists:
                    save_network(self.sys, destination, fail_if_exists=True)
                else:
                    save_network(self.sys, destination)
            self.sys._last_export_archive_output_path = str(destination)
            self.sys._last_export_archive_path = None
            print(f'Solved network saved: {destination}')

    def _run_archive_interface_recovery(self, network, archive_path):
        """Recover an ordered interface from an already solved source network."""
        from types import SimpleNamespace

        from vorpy.src.io import validate_network_source

        metadata = dict(getattr(network, 'archive_metadata', {}) or {})
        source_path = metadata.get('original_structure_path')
        validation = validate_network_source(
            network,
            source_path,
            expected_scheme='pow',
            expected_solvent_mode='MOLECULE_ONLY',
        )
        self.sys.archive_recovery_validation = validation
        self.sys.archive_recovery_checkpoint = str(Path(archive_path).resolve())
        # Interface Result identity is the original source, while the loaded
        # checkpoint remains the geometry producer and is never overwritten.
        self.base_file = Path(validation['source_path'])
        self.sys.files['base_file'] = str(self.base_file)
        self._full_system_group = SimpleNamespace(net=network)
        self._full_system_separate = True
        self._build_group_networks = False

        ggroup(self.sys, self.groups, dict(network.settings), make_net=False)
        interface_pairs = self._planned_interface_pairs()
        if not interface_pairs or any(group_b is None for _, group_b in interface_pairs):
            raise ValueError(
                'Checkpoint recovery requires two ordered interface selections: '
                '-g <group A> -g <group B> -i'
            )

        # This reuses the existing source-view path. It only selects rows from
        # loaded geometry; no Group or Network build is reachable here.
        self._attach_full_system_interface_views(interface_pairs)
        self._print_interface_statistics()
        handoffs = self._make_cli_interface_results(interface_pairs)
        if handoffs:
            self.interface_result_state = _typed_cli_interface_result_state(
                handoffs,
                unresolved_count=len(interface_pairs) - len(handoffs),
            )
        else:
            self.interface_result_state = _legacy_interface_result_state()
            if self.interface_result_errors:
                self.interface_result_state['reason'] += (
                    ' Handoff diagnostics: '
                    + ' | '.join(self.interface_result_errors)
                )
        self.sys.cli_interface_result_state = self.interface_result_state

        self.sys._interface_default_export = not self.exports
        self.sys._compact_interface_workflow = not self.exports or any(
            command and str(command[0]).strip().lower() == 'visualize'
            for command in self.exports
        )
        self.sys._interface_water_requested = False
        self.sys._interface_curvature_requested = False
        self.sys._network_calculate_curvature = False
        self._prepare_interface_output_directories()
        self.sys._last_export_archive_path = None
        self.sys._last_export_archive_output_path = None
        interfaces = list(getattr(self.sys, 'ifaces', ()) or ())
        export_interfaces = self._partial_interface_export_views(interfaces)
        self.sys.ifaces = export_interfaces
        try:
            self.run_exports()
        finally:
            self.sys.ifaces = interfaces
        self._annotate_partial_interface_exports(interfaces)
        self._save_archive()

    @staticmethod
    def _partial_interface_export_views(interfaces):
        """Give legacy writers only locally valid rows from an incomplete view."""
        from copy import copy

        exported = []
        for interface in interfaces:
            view = getattr(interface, 'source_view', None)
            if view is None or view.completeness.proven:
                exported.append(interface)
                continue
            source = interface.net
            if source is None:
                exported.append(interface)
                continue
            detached_network = copy(source)
            detached_network.balls = source.balls.copy()
            surface_positions = [int(feature.topology_index) for feature in view.surfaces]
            edge_positions = [int(feature.topology_index) for feature in view.edges]
            vertex_positions = [int(feature.topology_index) for feature in view.vertices]
            edge_map = {source_index: local_index for local_index, source_index in enumerate(edge_positions)}
            vertex_map = {source_index: local_index for local_index, source_index in enumerate(vertex_positions)}
            surface_map = {source_index: local_index for local_index, source_index in enumerate(surface_positions)}

            def local_rows(table_name, positions):
                table = getattr(source, table_name).iloc[positions].copy().reset_index(drop=True)
                maps = {
                    'surfs': {'edges': edge_map, 'verts': vertex_map},
                    'edges': {'verts': vertex_map, 'surfs': surface_map},
                    'verts': {'edges': edge_map, 'surfs': surface_map},
                }.get(table_name, {})
                for column, mapping in maps.items():
                    if column not in table:
                        continue
                    table[column] = [
                        [mapping[int(value)] for value in values if int(value) in mapping]
                        for values in table[column]
                    ]
                return table

            detached_network.surfs = local_rows('surfs', surface_positions)
            detached_network.edges = local_rows('edges', edge_positions)
            detached_network.verts = local_rows('verts', vertex_positions)
            group_a_ids = set(int(value) for value in interface.group1_indices)
            group_b_ids = set(int(value) for value in interface.group2_indices)
            def direct_first(values):
                direct = []
                remainder = []
                for value in values:
                    system_id = int(source.balls.iloc[int(value)].get('system_num', value))
                    (direct if system_id in group_a_ids or system_id in group_b_ids else remainder).append(value)
                ordered = [
                    value for value in direct
                    if int(source.balls.iloc[int(value)].get('system_num', value)) in group_a_ids
                ][:1]
                ordered += [
                    value for value in direct
                    if int(source.balls.iloc[int(value)].get('system_num', value)) in group_b_ids
                ][:1]
                return ordered + [value for value in values if value not in ordered]
            detached_network.surfs['balls'] = [
                direct_first(values) for values in detached_network.surfs['balls']
            ]
            detached = copy(interface)
            detached.net = detached_network
            detached.source_network = source
            detached.source_view = None
            detached._recovery_partial_completeness = view.completeness
            exported.append(detached)
        return exported

    @staticmethod
    def _annotate_partial_interface_exports(interfaces):
        """Keep view-level incompleteness visible in human-readable info files."""
        for interface in interfaces:
            view = getattr(interface, 'source_view', None)
            if view is None or view.completeness.proven or not getattr(interface, 'dir', None):
                continue
            info = Path(interface.dir) / 'interface_info.txt'
            if not info.is_file():
                continue
            analysis = getattr(interface, 'geometry_analysis', None)
            side = getattr(analysis, 'voronoi_side_1', {}) or {}
            area = side.get('area', 'not available')
            with info.open('a', encoding='utf-8') as stream:
                stream.write(
                    '\nRecovery source-view status:\n'
                    '  Completeness: INCOMPLETE\n'
                    f'  Unresolved source facets: {len(view.completeness.unresolved_surface_ids)}\n'
                    f'  Source-view physical facets: {len(view.surfaces)} surfaces, '
                    f'{len(view.edges)} edges, {len(view.vertices)} vertices\n'
                    f'  Source-view analytic area (Å²): {area}\n'
                    '  Local physical facets were exported; global closure and curvature remain unresolved.\n'
                )

    def _run_archive(self):
        from vorpy.src.io import group_from_network, load_network
        profile = None
        profile_started = None
        if self.cache_export_profile:
            profile = {
                'timing': {},
                'export_phases': {},
            }
            profile_started = (perf_counter(), process_time())
            started = (perf_counter(), process_time())
        network = load_network(self.base_file)
        if profile is not None:
            profile['timing']['load'] = _profile_interval(started)
        self.sys = network.sys
        self.sys._loaded_from_archive = True
        self.sys._fail_if_exists = self.fail_if_exists
        self.sys.set_output_directory(self.sys.files['dir'])
        if profile is not None:
            started = (perf_counter(), process_time())
        self.parse_commands()
        self.sys.verbose = self.verbose
        self.sys._write_archive_timing_manifest = bool(
            self.verbose or self.diagnose_edges or self.execution_profile is not None
            or self.cache_export_profile
        )
        self.sys._compact_interface_workflow = any(
            command and str(command[0]).strip().lower() == 'visualize'
            for command in self.exports
        )
        if self.interface_mode:
            self._run_archive_interface_recovery(network, self.base_file)
            if profile is not None:
                profile['timing']['command_setup'] = _profile_interval(started)
                profile['timing']['export'] = 0.0
                self._finish_cache_export_profile(profile, profile_started)
            return
        if self.builds or self.load_commands or self.settings_cmnds:
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
        if profile is not None:
            profile['timing']['command_setup'] = _profile_interval(started)
            self.sys._last_export_archive_path = None
            self.sys._last_export_archive_output_path = None
        if self.aw_alpha_values:
            started = (perf_counter(), process_time())
            self.run_aw_alpha_experimental()
            if profile is not None:
                profile['timing']['export'] = _profile_interval(started)
            started = (perf_counter(), process_time())
            self._save_archive()
            if profile is not None:
                profile['timing']['explicit_archive'] = _profile_interval(started)
                self._finish_cache_export_profile(profile, profile_started)
            return
        started = (perf_counter(), process_time())
        self.run_exports()
        if profile is not None:
            profile['timing']['export'] = _profile_interval(started)
        started = (perf_counter(), process_time())
        self._save_archive()
        if profile is not None:
            profile['timing']['explicit_archive'] = _profile_interval(started)
            self._finish_cache_export_profile(profile, profile_started)

    def _finish_cache_export_profile(self, profile, profile_started):
        """Persist the opt-in profile after cache-only export work completes."""
        from vorpy.src.io.network_archive import _process_memory

        profile['timing']['total'] = _profile_interval(profile_started)
        output_archive = getattr(self.sys, '_last_export_archive_output_path', None)
        output_archive = output_archive or getattr(self.sys, '_last_export_archive_path', None)
        profile['archive_members'] = {
            'input': _archive_member_counts(self.base_file),
            'output': _archive_member_counts(output_archive),
        }
        profile['cache_state'] = _cache_export_state(self.sys)
        profile['export_phases'] = {
            'run_exports': profile['timing'].get('export'),
            'operations': getattr(self.sys, 'export_timing', {}),
            'mesh': getattr(self.sys, 'export_mesh_timing', {}),
        }
        rss, source = _process_memory()
        profile['peak_rss_bytes'] = rss if source == 'resource_ru_maxrss' else None
        profile['peak_rss_supported'] = source == 'resource_ru_maxrss'
        profile['rss_bytes'] = rss
        profile['memory_source'] = source
        profile['output_archive'] = None if output_archive is None else Path(output_archive).name
        _write_cache_export_profile(self.sys, self.base_file, profile)

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

    def _legacy_parse_commands(self, counter=0):
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
        my_args, requested_solvent_mode = solvent_cli_options(my_args)
        if requested_solvent_mode is not None:
            self.solvent_mode = requested_solvent_mode
            self._solvent_mode_explicit = True
            if self.sys is not None:
                self.sys.solvent_mode = requested_solvent_mode
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
            if self.fail_if_exists:
                _check_archive_target_available(self.save_network_path)
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
                if arg_cmnds and arg_cmnds[0].lower() in {'logs', 'log'}:
                    self.settings_cmnds.append(['bt', 'logs'])
                    if len(arg_cmnds) > 1:
                        self.logs_files.append(' '.join(arg_cmnds[1:]))
            elif arg.lower() == '-i':
                self.interface_mode = True

            if arg.lower() == '-e':
                # If the argument is a directory command, format it
                if arg_cmnds and arg_cmnds[0] in {'dir', 'directory'}:
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

                        try:
                            _prepare_output_directory(
                                system_dir,
                                fail_if_exists=self.fail_if_exists,
                            )
                        except FileExistsError as error:
                            raise SystemExit(str(error)) from None

                        self.sys.dir = system_dir
                        self.sys.files['dir'] = system_dir

                        if self.verbose:
                            print(f"Directory set to: {system_dir}")
                else:
                    # Add the export command to the list
                    self.exports.append(arg_cmnds)

    def parse_commands(self, counter=0):
        """Parse CLI blocks into a calculation plan and legacy command state."""
        my_args = list(getattr(self, '_arguments', sys.argv[1:])[1 + counter:])
        my_args, timing_report = timing_report_cli_options(my_args)
        if timing_report is not None and self.timing_report_path is None:
            self.timing_report_path = timing_report
            self._cli_timing = {
                'schema_version': 1,
                'unit': 'seconds',
                'phases': {},
            }
            if self._cli_timing_started is None:
                self._cli_timing_started = (perf_counter(), process_time())
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
        from vorpy.src.geometry.aw_alpha.cli import (
            extract_aw_alpha_options,
            extract_aw_alpha_pair_only,
        )
        my_args, aw_alpha_values = extract_aw_alpha_options(my_args)
        if aw_alpha_values:
            self.aw_alpha_values = aw_alpha_values
        my_args, pair_overlap_only = extract_aw_alpha_pair_only(my_args)
        if pair_overlap_only:
            self.aw_alpha_pair_overlap_only = True
        planning_args = list(my_args)
        if '--save-network' in my_args:
            index = my_args.index('--save-network')
            if index + 1 >= len(my_args) or my_args[index + 1].startswith('-'):
                raise ValueError('--save-network requires an output .vpy path')
            self.save_network_path = Path(my_args[index + 1]).expanduser().resolve()
            if self.fail_if_exists:
                _check_archive_target_available(self.save_network_path)
            del my_args[index:index + 2]
        my_args, _ = frame_cli_options(my_args)
        my_args = [arg for arg in my_args if arg != '--all-frames']

        if '-v' in my_args or '--verbose' in my_args:
            self.verbose = True
            my_args = [arg for arg in my_args if arg not in {'-v', '--verbose'}]
        if '--diagnose-edges' in my_args:
            self.diagnose_edges = True
            my_args = [arg for arg in my_args if arg != '--diagnose-edges']

        self.calculation_plan = build_calculation_plan(
            [self.base_file or 'input', *planning_args],
            input_source=self.base_file,
            analysis_profile=self.execution_profile,
        )
        plan = self.calculation_plan
        if self.execution_profile is None and plan.analysis_profile is not None:
            self.execution_profile = plan.analysis_profile
        self.interface_mode = bool(plan.interface_requests)
        self._planned_interface_operands = tuple(
            request.operands for request in plan.interface_requests
        )
        self._build_group_networks = bool(
            plan.all_groups_request
            or plan.independent_group_requests
        )
        self._full_system_separate = bool(
            plan.interface_requests
            or (
                plan.full_system_request
                and plan.declared_group_requests
            )
        )

        # Preserve the old execution containers while making block parsing
        # exact and ordered.  ``and`` remains a selector separator inside a
        # block, rather than accidentally ending that block.
        group_counter = -1
        declared_keys = set()
        for block in plan.blocks:
            if block.flag == 'load':
                self.load_commands.append(list(block.args))
            elif block.flag == 'settings':
                # ``solver_settings`` already contains all settings blocks;
                # apply each setting once even when several -s blocks occur.
                for setting in plan.solver_settings:
                    setting = list(setting)
                    if setting not in self.settings_cmnds:
                        self.settings_cmnds.append(setting)
            elif block.flag == 'group':
                group_counter += 1
                commands = []
                current = []
                for token in expand_group_tokens(block.args):
                    if token.lower() in {'&', 'and', 'nd', 'also', '+', '&&'}:
                        if current:
                            commands.append(current)
                            current = []
                    else:
                        current.append(token)
                if current:
                    commands.append(current)
                self.groups[group_counter] = commands
                declared_keys.add(tuple(tuple(token.lower() for token in command)
                                        for command in commands))
            elif block.flag == 'build':
                build = list(block.args)
                self.builds.append(build)
                if build and build[0].lower() in {'logs', 'log'}:
                    self.settings_cmnds.append(['bt', 'logs'])
                    if len(build) > 1:
                        self.logs_files.append(' '.join(build[1:]))
            elif block.flag == 'interface':
                self.interface_mode = True
            elif block.flag == 'export':
                # Export options are flattened so ``-e a and b`` retains the
                # intended two export commands.
                exports = []
                current = []
                for token in block.args:
                    if token.lower() in {'&', 'and', 'nd', 'also', '+', '&&'}:
                        if current:
                            exports.append(tuple(current))
                            current = []
                    else:
                        current.append(token)
                if current:
                    exports.append(tuple(current))
                for export in exports:
                    if not export:
                        continue
                    command = list(export)
                    if command[0].lower() in {'dir', 'directory'}:
                        if len(command) < 2:
                            raise ValueError('Export directory requires a path.')
                        if command[1] in browse_names:
                            my_root = tk.Tk()
                            my_root.withdraw()
                            my_root.wm_attributes('-topmost', 1)
                            folder = filedialog.askdirectory(title='Choose Output Folder')
                            if os.path.exists(folder):
                                self.sys.files['dir'] = folder
                            else:
                                print(f"{folder} is not a valid folder")
                        else:
                            base_dir = os.path.abspath(os.path.expanduser(command[1]))
                            system_dir = os.path.join(base_dir, self.sys.name)
                            try:
                                _prepare_output_directory(
                                    system_dir,
                                    fail_if_exists=self.fail_if_exists,
                                )
                            except FileExistsError as error:
                                raise SystemExit(str(error)) from None
                            self.sys.dir = system_dir
                            self.sys.files['dir'] = system_dir
                            if self.verbose:
                                print(f"Directory set to: {system_dir}")
                    else:
                        self.exports.append(command)

        # Explicit interface operands and calculate-chain targets define
        # groups without turning them into independent solves.
        for request in plan.declared_group_requests:
            if request.key in declared_keys or not request.commands:
                continue
            group_counter += 1
            commands = [list(command) for command in request.commands]
            self.groups[group_counter] = commands
            declared_keys.add(request.key)

        if plan.full_system_request and not self.groups and not plan.interface_requests:
            self.groups[0] = [['f']]

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
        default_selection = not self.groups
        make_net = (
            not self._full_system_reuse_enabled()
            and ((not self.interface_mode) or self._build_group_networks)
        )
        # The molecular default group is assembled from System.residues, while
        # structural ions are held by the loader's solvent container.
        retained_indices = tuple(getattr(self.sys, 'retained_system_indices', ()) or ())
        if default_selection and retained_indices:
            make_net = False
        ggroup(
            self.sys,
            self.groups,
            self.settings_dict,
            make_net=make_net,
        )

        plan = self.calculation_plan
        implicit_full_default = bool(
            self._full_system_reuse_enabled()
            and plan is not None
            and plan.full_system_request
            and not plan.declared_group_requests
            and not plan.interface_requests
        )
        for group in (self.sys.groups or []):
            if default_selection and retained_indices:
                group.add_balls(retained_indices)
            elif implicit_full_default:
                active = set(getattr(self, '_source_system_indices', ()))
                group.ball_ndxs = [
                    int(index) for index in group.ball_ndxs if int(index) in active
                ]
                group.atms = [
                    int(index) for index in getattr(group, 'atms', ())
                    if int(index) in active
                ]

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

    def _record_optional_analysis_failure(self, name, error):
        """Keep optional analysis failures visible while preserving exports."""
        failure = {
            'analysis': str(name),
            'status': 'UNRESOLVED',
            'error_type': type(error).__name__,
            'message': str(error),
        }
        failures = list(getattr(self.sys, '_optional_analysis_failures', ()) or ())
        failures.append(failure)
        self.sys._optional_analysis_failures = failures
        for _, network in self._diagnostic_networks():
            provenance = dict(getattr(network, 'archive_provenance', {}) or {})
            recorded = list(provenance.get('optional_analysis_failures', ()) or ())
            recorded.append(failure)
            provenance['optional_analysis_failures'] = recorded
            provenance['analysis_status'] = 'PARTIAL'
            network.archive_provenance = provenance
        print(
            f'Optional analysis {name} failed ({type(error).__name__}: {error}); '
            'solved geometry exports continue.',
            file=sys.stderr,
        )

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
