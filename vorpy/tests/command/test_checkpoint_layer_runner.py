import csv
import importlib.util
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location(
    'checkpoint_layer_runner', Path(__file__).resolve().parents[3] / 'scripts/run_checkpoint_layers.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def test_manifest_orders_by_size_and_rejects_colliding_outputs(tmp_path):
    manifest = tmp_path / 'jobs.csv'
    (tmp_path / 'small.vpy').write_bytes(b'1')
    (tmp_path / 'large.vpy').write_bytes(b'123')
    manifest.write_text('checkpoint,output\nsmall.vpy,out/small\nlarge.vpy,out/large\n')
    jobs = runner.load_jobs(manifest, tmp_path)
    assert [job[0].name for job in jobs] == ['large.vpy', 'small.vpy']
    assert runner.build_command(*jobs[0][:2])[4:] == [
        '-c', 'l', '4', '-e', 'small', '-e', 'dir', str(tmp_path / 'out/large'), '-v']
    manifest.write_text('checkpoint,output\nsmall.vpy,out\nlarge.vpy,out\n')
    with pytest.raises(ValueError, match='Duplicate output'):
        runner.load_jobs(manifest, tmp_path)


def test_queue_bounds_workers_and_records_failed_jobs(tmp_path, monkeypatch):
    manifest = tmp_path / 'jobs.csv'
    rows = ['checkpoint,output']
    for index in range(3):
        (tmp_path / f'{index}.vpy').write_bytes(b'x' * (index + 1))
        rows.append(f'{index}.vpy,out/{index}')
    manifest.write_text('\n'.join(rows))
    active = []
    started = []

    class Process:
        def __init__(self, command, **kwargs):
            assert len(active) < 2
            assert kwargs['env']['OPENBLAS_NUM_THREADS'] == '1'
            assert kwargs['env']['OMP_NUM_THREADS'] == '1'
            assert kwargs['env']['VORPY_PROGRESS_STREAM'] == '1'
            self.code = 1 if Path(command[3]).name == '1.vpy' else 0
            self.polls = 0
            started.append(Path(command[3]).name)
            active.append(self)

        def poll(self):
            self.polls += 1
            if self.polls == 1:
                return None
            active.remove(self)
            return self.code

    monkeypatch.setattr(runner.subprocess, 'Popen', Process)
    monkeypatch.setattr(runner.time, 'sleep', lambda _: None)
    assert runner.main([str(manifest), '--workers', '2']) == 1
    assert started == ['2.vpy', '1.vpy', '0.vpy']
    with next(tmp_path.glob('vorpy_layer4_runlogs/*/summary.csv')).open() as stream:
        results = list(csv.DictReader(stream))
    assert len(results) == 3
    assert [row['exit_code'] for row in results] == ['0', '1', '0']
    assert not active


def test_live_display_updates_each_worker_without_interleaving(tmp_path, monkeypatch):
    import io
    import os
    log = tmp_path / '0001.log'
    log.write_bytes(b'old stage\rVertices - 25%\r\x1b[32mSurfaces - 50%\x1b[0m\n')
    assert runner.latest_output(log) == 'Surfaces - 50%'
    stream = io.StringIO()
    monkeypatch.setattr(runner.sys, 'stdout', stream)
    monkeypatch.setattr(runner, 'enable_live_console', lambda: True)
    monkeypatch.setattr(runner.shutil, 'get_terminal_size', lambda _: os.terminal_size((140, 30)))
    display = runner.ProgressDisplay(tmp_path, 2, 3)
    source = tmp_path / 'PEG/Checkpoints/aw/EG/frame_0001/EG.vpy'
    active = [(None, None, source, None, log, runner.time.monotonic())]
    display.render(active, 1, 0, 'Job finished', force=True)
    first = stream.getvalue()
    assert '1/3 completed | 1 running | 0 failed' in first
    assert 'Job 1/3 | Surfaces - 50% | 0s | PEG/aw/EG/frame_0001' in first
    assert 'Idle' in first
    log.write_text('Analysis - 75%\n')
    display.render(active, 1, 0, 'Job finished', force=True)
    assert '\x1b[4F\x1b[J' in stream.getvalue()[len(first):]
    assert 'Analysis - 75%' in stream.getvalue()[len(first):]


def test_progress_percentage_stays_visible_in_a_narrow_console(tmp_path, monkeypatch):
    import io
    import os
    log = tmp_path / '0007.log'
    log.write_text('0:01:23.45 - water_frontier_2 - Vertices: 120 / ~1,000 - 12.00 %\n')
    stream = io.StringIO()
    monkeypatch.setattr(runner.sys, 'stdout', stream)
    monkeypatch.setattr(runner, 'enable_live_console', lambda: True)
    monkeypatch.setattr(runner.shutil, 'get_terminal_size', lambda _: os.terminal_size((80, 30)))
    display = runner.ProgressDisplay(tmp_path, 1, 429)
    source = tmp_path / 'Heparin/Checkpoints/aw/dp36/frame_0001/dp36_checkpoint.vpy'
    display.render([(None, None, source, None, log, runner.time.monotonic())], 0, 0, '', force=True)
    row = stream.getvalue().splitlines()[1]
    assert row.startswith('Job 7/429 | Vertices 12.00% (120 / ~1,000) | water_frontier_2')
    assert len(row) <= 79


def test_real_workers_keep_full_logs_and_show_progress(tmp_path, monkeypatch, capsys):
    import sys
    manifest = tmp_path / 'jobs.csv'
    (tmp_path / 'one.vpy').write_bytes(b'x')
    (tmp_path / 'two.vpy').write_bytes(b'xx')
    manifest.write_text('checkpoint,output\none.vpy,out/one\ntwo.vpy,out/two\n')
    monkeypatch.setattr(runner, 'enable_live_console', lambda: False)
    code = ('import time\nfrom vorpy.src.system import System\n'
            'system = System(simple=True, make_dir=False)\n'
            'system.update_progress("Connections", 0, network="water_frontier_1")\n'
            'time.sleep(0.6)\n'
            'system.update_progress("Connections", 12, network="water_frontier_1")\n'
            'time.sleep(0.7)\nprint("Finished", flush=True)\n')
    monkeypatch.setattr(runner, 'build_command', lambda *_: [sys.executable, '-u', '-c', code])
    assert runner.main([str(manifest), '--workers', '2']) == 0
    output = capsys.readouterr().out
    assert 'Connections 12.00%' in output
    logs = list(tmp_path.glob('vorpy_layer4_runlogs/*/*.log'))
    assert len(logs) == 2
    assert all('Connections - 12.00 %' in log.read_text() and 'Finished' in log.read_text() for log in logs)
