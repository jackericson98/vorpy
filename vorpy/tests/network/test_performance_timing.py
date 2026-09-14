"""Boundary accounting must survive nested calls and exceptions."""
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    'performance_timing', Path(__file__).resolve().parents[3] / 'scripts/performance_timing.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_nested_spans_reconcile_without_double_counting():
    ticks = iter([0, 1, 3, 5])
    recorder = module.Recorder(wall=lambda: next(ticks), cpu=lambda: 0)
    with recorder.span('parent', 'build'):
        with recorder.span('child', 'export'):
            pass
    result = recorder.report()
    parent, child = result['events']
    assert parent['wall_seconds'] == 5
    assert child['wall_seconds'] == 2
    assert parent['wall_children_seconds'] == 2
    assert parent['wall_unattributed_seconds'] == 3
    assert sum(item['wall_seconds'] for item in result['exclusive_categories'].values()) == 5


def test_wrappers_restore_after_exception_and_keep_call_tree():
    class Example:
        def fail(self):
            raise ValueError('expected')
    original = Example.fail
    recorder = module.Recorder()
    with pytest.raises(ValueError, match='expected'):
        with module.instrument(recorder, [(Example, 'fail', 'fail', 'setup')]):
            with recorder.span('total', 'other'):
                Example().fail()
    assert Example.fail is original
    assert recorder.stack == []
    assert recorder.report()['events'][1]['parent'] == 0


def test_cpu_accounting_and_repeated_sibling_calls():
    wall = iter([0, 1, 3, 4, 7, 10])
    cpu = iter([0, 1, 2, 3, 5, 8])
    recorder = module.Recorder(wall=lambda: next(wall), cpu=lambda: next(cpu))
    with recorder.span('total', 'other'):
        for _ in range(2):
            with recorder.span('work', 'build'):
                pass
    report = recorder.report()
    assert report['exclusive_categories']['build']['calls'] == 2
    assert report['exclusive_categories']['build']['wall_seconds'] == 5
    assert report['exclusive_categories']['build']['cpu_seconds'] == 3
    assert report['events'][0]['cpu_unattributed_seconds'] == 5


def test_output_normalization_preserves_numerical_changes(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[3] / 'scripts'))
    import benchmark_pipeline as benchmark
    roots = [tmp_path / name / 'artifacts' for name in ('control', 'trial')]
    for index, root in enumerate(roots):
        root.mkdir(parents=True)
        (root / 'EDTA_info.txt').write_text(
            f'Output Directory: {root}\nVertex Time: {index + 1}.0 s\nVolume: 123.456\n')
        (root / 'EDTA_logs.csv').write_text(
            f'Name,Completion Date,Vertex Time,Volume\nEDTA,date-{index},{index},123.456\n')
    assert benchmark.output_manifest(roots[0]) == benchmark.output_manifest(roots[1])
    path = roots[1] / 'EDTA_logs.csv'
    path.write_text(path.read_text().replace('123.456', '123.457'))
    assert benchmark.output_manifest(roots[0]) != benchmark.output_manifest(roots[1])


def test_snapshot_digest_preserves_sequence_order_and_exact_values(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[3] / 'scripts'))
    import benchmark_pipeline as benchmark
    assert benchmark.digest(benchmark.canonical([1, 2])) != benchmark.digest(benchmark.canonical([2, 1]))
    assert benchmark.digest(benchmark.canonical([1.0])) != benchmark.digest(benchmark.canonical([1.0 + 1e-12]))
    assert benchmark.canonical({1, 2}) == benchmark.canonical({2, 1})
