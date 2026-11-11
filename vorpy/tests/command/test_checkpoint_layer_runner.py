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
        '-c', 'l', '4', '-e', 'small', '-e', 'dir', str(tmp_path / 'out/large')]
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
