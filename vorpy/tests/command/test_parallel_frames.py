import os
from pathlib import Path
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Lock

import pytest

from vorpy.src.command.vpy_cmnd import Command, frame_cli_options, prompt_frame_workers
from vorpy.src.command.parallel_frames import run_frame
from vorpy.src.command import parallel_frames


@pytest.mark.parametrize('cpus, frames, answers, expected', [
    (12, 10, ['y', '50'], 6),
    (12, 10, ['', ''], 6),
    (12, 10, ['yes', '50%'], 6),
    (12, 10, ['y', '100'], 10),
    (12, 10, ['y', '1'], 1),
    (3, 10, ['y', '50'], 1),
    (12, 10, ['n'], None),
    (12, 10, ['maybe', 'y', 'bad', '0', '101', 'nan', '50'], 6),
])
def test_interactive_worker_selection(monkeypatch, cpus, frames, answers, expected):
    monkeypatch.setattr('vorpy.src.command.vpy_cmnd.available_cpu_count', lambda: cpus)
    replies = iter(answers)
    monkeypatch.setattr('builtins.input', lambda prompt: next(replies))
    assert prompt_frame_workers(frames) == expected


@pytest.mark.parametrize('answers', [[], ['y']])
def test_interactive_worker_selection_requires_response(monkeypatch, answers):
    replies = iter(answers)

    def reply(prompt):
        try:
            return next(replies)
        except StopIteration:
            raise EOFError

    monkeypatch.setattr('builtins.input', reply)
    with pytest.raises(SystemExit, match='No response received'):
        prompt_frame_workers(10)


def test_all_frames_prompt_routes_percentage_to_parallel_runner(tmp_path, monkeypatch):
    source = tmp_path / 'trajectory.pdb'
    source.write_text(('ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\nEND\n') * 10)
    monkeypatch.setattr(sys, 'argv', ['vorpy', str(source)])
    monkeypatch.setattr('vorpy.src.command.vpy_cmnd.available_cpu_count', lambda: 12)
    replies = iter(['y', 'y', '50'])
    prompts = []

    def reply(prompt):
        prompts.append(prompt)
        return next(replies)

    monkeypatch.setattr('builtins.input', reply)
    command = Command()
    calls = []
    monkeypatch.setattr(command, 'run_all_frames', lambda **kwargs: calls.append(kwargs))
    command.run()
    assert calls == [{'total_frames': 10, 'workers': 6}]
    assert prompts == ['10 frames found. Run all? [Y/n] ',
                       'Parallelize frame runs? [Y/n] ',
                       'Enter CPU usage percentage (1-100) [50]: ']


@pytest.mark.parametrize('args', [
    ['--parallel-frames'], ['--parallel-frames', '0'],
    ['--parallel-frames', '-2'], ['--parallel-frames', 'two'],
    ['--parallel-frames', '2', '--parallel-frames', '3'],
])
def test_invalid_worker_counts(args):
    with pytest.raises(SystemExit, match='--parallel-frames'):
        frame_cli_options(args)


def test_parallel_flag_preserves_grouped_options(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['vorpy', 'sample.pdb', '-s', 'mv', '5',
                                    '--parallel-frames', '2', '--all-frames',
                                    '-e', 'only', 'pdb'])
    command = Command()
    command.parse_commands()
    assert command.settings_cmnds == [['mv', '5']]
    assert command.exports == [['only', 'pdb']]


def test_worker_failure_records_traceback_and_restores_cwd(tmp_path, monkeypatch):
    source = tmp_path / 'input.pdb'
    source.write_text('ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\nEND\n')
    output = tmp_path / 'output'
    output.mkdir()
    original_cwd = Path.cwd()

    def fail(command):
        os.chdir(output)
        raise ValueError('test frame failure')

    monkeypatch.setattr(Command, '_process_system', fail)
    error = run_frame(str(source), str(output), 2, 3, {}, None)
    assert error == 'ValueError: test frame failure'
    assert Path.cwd() == original_cwd
    assert 'Traceback' in (output / 'run.log').read_text()
    assert 'test frame failure' in (output / 'run.log').read_text()


@pytest.mark.parametrize('workers, total', [(2, 3), (6, 10)])
def test_scheduler_bounds_active_frames_and_continues_after_failure(tmp_path, monkeypatch, capsys, workers, total):
    source = tmp_path / 'trajectory.pdb'
    source.write_text(''.join(
        f'ATOM      1  CA  ALA A   1    {float(i):8.3f}   0.000   0.000  1.00  0.00           C\nEND\n'
        for i in range(1, total + 1)))
    monkeypatch.setattr(sys, 'argv', ['vorpy', str(source), '--parallel-frames', str(workers),
                                    '-e', 'dir', str(tmp_path / 'results')])
    barrier = Barrier(workers)
    lock = Lock()
    active = 0
    maximum = 0
    seen = {}

    def worker(frame_file, directory, ordinal, total, options, settings):
        nonlocal active, maximum
        with lock:
            active += 1
            maximum = max(maximum, active)
        if ordinal <= workers:
            barrier.wait(timeout=10)
        # Both workers see their own frame even after the generator advances.
        assert float(Path(frame_file).read_text().splitlines()[0][30:38]) == ordinal
        seen[ordinal] = frame_file
        Path(directory, 'result.txt').write_text(str(ordinal))
        with lock:
            active -= 1
        return 'intentional failure' if ordinal == 2 else None

    monkeypatch.setattr(parallel_frames, 'run_frame', worker)
    monkeypatch.setattr(parallel_frames, 'ProcessPoolExecutor',
                        lambda max_workers, mp_context: ThreadPoolExecutor(max_workers=max_workers))
    with pytest.raises(SystemExit, match=f'1 of {total} frames failed.*frame 2: intentional failure'):
        Command().run()
    assert maximum == workers
    assert sorted(seen) == list(range(1, total + 1))
    assert all(not Path(file).exists() for file in seen.values())
    assert f'Frame {total}/{total} complete' in capsys.readouterr().out


@pytest.mark.parametrize('workers, entry', [(1, 'directory'), (2, 'module'), (None, 'module')])
def test_spawned_cli_builds_and_exports_every_frame(tmp_path, workers, entry):
    root = Path(__file__).resolve().parents[2]
    atoms = [line for line in (root / 'data' / 'DB1976.pdb').read_text().splitlines(True)
             if line.startswith(('ATOM  ', 'HETATM'))]
    atoms = atoms[:12]
    source = tmp_path / 'trajectory.pdb'
    source.write_text(''.join(
        'MODEL        %d\n%sENDMDL\n' % (i, ''.join(
            line[:30] + f'{float(line[30:38]) + i:8.3f}' + line[38:] for line in atoms
        )) for i in (1, 2, 3)))
    output = tmp_path / 'results with spaces'
    invocation = [str(root)] if entry == 'directory' else ['-m', 'vorpy']
    parallel_args = ['--parallel-frames', str(workers)] if workers is not None else []
    result = subprocess.run(
        [sys.executable, *invocation, str(source), *parallel_args,
         '-s', 'nt', 'pow', '-s', 'mv', '5', '-s', 'sr', '1',
         '-e', 'dir', str(output), '-e', 'only', 'pdb'],
        cwd=root.parent, capture_output=True, text=True, timeout=90,
        input='y\ny\n50\n' if workers is None else None,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'All 3 frames complete' in result.stdout
    frames = output / 'trajectory' / 'frames'
    assert sorted(path.name for path in frames.iterdir()) == [f'frame_{i:04d}' for i in (1, 2, 3)]
    for i in (1, 2, 3):
        directory = frames / f'frame_{i:04d}'
        assert (directory / 'trajectory.pdb').is_file()
        exported = (directory / 'trajectory.pdb').read_text().splitlines()
        first_atom = next(line for line in exported if line.startswith(('ATOM  ', 'HETATM')))
        assert float(first_atom[30:38]) == pytest.approx(float(atoms[0][30:38]) + i)
        log = (directory / 'run.log').read_text()
        assert f'frame {i} only' in log
        assert 'network built' in log
