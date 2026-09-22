"""Bounded process workers for independent PDB trajectory frames."""

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import closing, redirect_stderr, redirect_stdout
from copy import deepcopy
import gc
import multiprocessing
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import traceback

from vorpy.src.inputs.frames import count_pdb_frames, iter_pdb_frames


def run_frame(frame_file, directory, ordinal, total_frames, options, settings):
    """Run one frame in a child process; keep geometry and cwd process-local."""
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.system.system import System

    original_directory = os.getcwd()
    system = None
    command = None
    with open(Path(directory) / 'run.log', 'w', encoding='utf-8', buffering=1) as log:
        with redirect_stdout(log), redirect_stderr(log):
            try:
                system = System(file=frame_file, make_dir=False)
                system.files['dir'] = directory
                system.frame_index = ordinal
                system.frame_count = total_frames
                command = Command(sys=system, settings=deepcopy(settings))
                command.base_file = frame_file
                for name, value in options.items():
                    setattr(command, name, deepcopy(value))
                command._process_system()
                return None
            except (Exception, SystemExit) as error:
                traceback.print_exc()
                return f'{type(error).__name__}: {error}'
            finally:
                os.chdir(original_directory)
                if system is not None:
                    Command._release_frame_system(system)
                del command, system
                gc.collect()


def run_parallel_frames(parent, workers, total_frames=None):
    """Stream at most one temporary frame per worker; report every failure."""
    from vorpy.src.command.vpy_cmnd import Command, FRAME_OPTIONS
    from vorpy.src.system.system import System

    if workers < 1:
        raise ValueError('The frame worker count must be positive.')
    if total_frames is None:
        total_frames = count_pdb_frames(parent.base_file)
    if not total_frames:
        raise ValueError(f'No atom frames found in PDB file: {parent.base_file}')
    workers = min(workers, total_frames)
    settings = deepcopy(parent.settings_dict)
    failures = []
    with closing(iter_pdb_frames(parent.base_file)) as frames, TemporaryDirectory(prefix='vorpy_parallel_') as staging:
        first_ordinal, first_file = next(frames)
        # Resolve export options once in the parent (including folder browsing).
        system = System(file=first_file, make_dir=False)
        setup = Command(sys=system, settings=deepcopy(settings))
        setup.base_file = first_file
        setup.parse_commands()
        if system.files['dir'] is None:
            system.set_output_directory()
        output_root = Path(system.files['dir']).resolve() / 'frames'
        options = {name: deepcopy(getattr(setup, name)) for name in FRAME_OPTIONS}
        Command._release_frame_system(system)
        del setup, system
        gc.collect()
        parent.sys = None
        print(f'{Path(parent.base_file).name}: {total_frames} frames found; '
              f'processing with {workers} worker processes.', flush=True)

        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
            pending = {}

            def submit(ordinal, source):
                # iter_pdb_frames reuses its path: copy before advancing it.
                frame_dir = Path(staging) / f'frame_{ordinal:04d}'
                frame_dir.mkdir()
                frame_file = frame_dir / Path(source).name
                shutil.copyfile(source, frame_file)
                directory = output_root / f'frame_{ordinal:04d}'
                directory.mkdir(parents=True, exist_ok=False)
                future = pool.submit(run_frame, str(frame_file), str(directory),
                                     ordinal, total_frames, options, settings)
                pending[future] = (ordinal, frame_file, directory)
                print(f'Frame {ordinal}/{total_frames} started: {directory}', flush=True)

            submit(first_ordinal, first_file)
            exhausted = False
            while pending or not exhausted:
                while not exhausted and len(pending) < workers:
                    frame = next(frames, None)
                    if frame is None:
                        exhausted = True
                    else:
                        submit(*frame)
                if not pending:
                    break
                completed, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in completed:
                    ordinal, frame_file, directory = pending.pop(future)
                    try:
                        error = future.result()
                    except Exception as exception:
                        error = f'{type(exception).__name__}: {exception}'
                    frame_file.unlink()
                    frame_file.parent.rmdir()
                    if error is None:
                        print(f'Frame {ordinal}/{total_frames} complete: {directory}', flush=True)
                    else:
                        failures.append((ordinal, error))
                        print(f'Frame {ordinal}/{total_frames} FAILED: {error}\n'
                              f'  Log: {directory / "run.log"}', flush=True)

    if failures:
        details = '; '.join(f'frame {ordinal}: {error}' for ordinal, error in sorted(failures))
        raise SystemExit(f'{len(failures)} of {total_frames} frames failed. {details}')
    print(f'All {total_frames} frames complete: {output_root}', flush=True)
