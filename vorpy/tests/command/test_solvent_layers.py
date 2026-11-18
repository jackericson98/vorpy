import csv
import hashlib
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.command.calculation_plan import build_calculation_plan
from vorpy.src.group.solvent_layers import water_neighbors, boundary_indices, complete_layer_count
from vorpy.src.io import load_network


@pytest.mark.parametrize('alias', ['l', 'layer', 'layers'])
def test_layer_count_aliases(alias):
    plan = build_calculation_plan(['EDTA', '-c', alias, '4', '-e', 'large'])
    assert plan.solvent_layers == 4
    assert not plan.full_system_request
    assert plan.export_settings == (('large',),)


@pytest.mark.parametrize('args', [[], ['0'], ['-1'], ['1.5'], ['four'], ['g']])
def test_layer_count_requires_positive_integer(args):
    with pytest.raises(ValueError, match='positive integer'):
        build_calculation_plan(['EDTA', '-c', 'l', *args])


def test_layer_requests_combine_with_other_targets():
    plan = build_calculation_plan(['EDTA', '-c', 'g,l', '2', '-c', 'layers', '4'])
    assert plan.solvent_layers == 4
    assert plan.all_groups_request


def graph():
    residues = [SimpleNamespace(name='MOL', atoms=[0]), SimpleNamespace(name='MG', atoms=[1]),
                SimpleNamespace(name='SOL', atoms=[2, 3, 4]),
                SimpleNamespace(name='SOL', atoms=[5, 6, 7]),
                SimpleNamespace(name='SOL', atoms=[8, 9, 10])]
    system = SimpleNamespace(balls=pd.DataFrame([
        {'res': residue} for residue in residues for _ in residue.atoms]))
    balls = pd.DataFrame({'num': range(100, 111), 'system_num': range(11)})
    faces = pd.DataFrame({'balls': [[100, 102], [100, 101], [101, 105], [103, 106], [104, 108]],
                          'sa': [1., 1., 1., 1., 0.], 'edges': [[0], [1], [2], [3], [4]],
                          'verts': [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]})
    net = SimpleNamespace(balls=balls, surfs=faces,
                          edges=pd.DataFrame({'verts': [[i, i + 1] for i in range(5)]}))
    return system, net


def test_water_graph_uses_whole_molecules_and_excludes_ions_and_zero_faces():
    system, net = graph()
    assert water_neighbors(net, system, {0}, {0}) == {2, 3, 4}
    assert water_neighbors(net, system, {2, 3, 4}, {0, 2, 3, 4}) == {5, 6, 7}
    assert water_neighbors(net, system, {5, 6, 7}, {0, 2, 3, 4, 5, 6, 7}) == set()
    surfaces, edges, verts = boundary_indices(net, [0], [[2, 3, 4], [5, 6, 7], [], []])
    assert surfaces == [[0], [3], [], []]
    assert edges == [[0], [3], [], []]
    assert verts == [[0, 1], [3, 4], [], []]


@pytest.mark.parametrize('complete_atoms, expected', [([], 0), ([2, 3, 4], 1),
    ([2, 3, 4, 5, 6, 7], 2), (list(range(2, 11)), 3), ([2, 3, 5, 6, 7, 8, 9, 10], 0)])
def test_complete_layers_require_every_cell_and_stop_at_first_failure(complete_atoms, expected):
    system, net = graph()
    net.balls['complete'] = pd.Series([i in complete_atoms for i in range(11)], dtype=object)
    net.balls.loc[0, 'complete'] = None
    group = SimpleNamespace(net=net, solvent_layers={'atoms': [[2, 3, 4], [5, 6, 7], [8, 9, 10], []]})
    assert complete_layer_count(group) == expected


def test_incomplete_layer_warning_and_summary(tmp_path, capsys):
    from vorpy.src.output.solvent import export_solvent_layers
    system, net = graph()
    net.balls['complete'] = [i in range(2, 8) for i in range(11)]
    group = SimpleNamespace(sys=system, net=net, dir=str(tmp_path), solvent_layers={
        'atoms': [[2, 3, 4], [5, 6, 7], [8, 9, 10], []],
        'surfs': [[0], [3], [], []], 'edges': [[0], [3], [], []]})
    export_solvent_layers(group, [])
    assert capsys.readouterr().out.strip() == (
        'Warning: 4 layers were requested; only 2 complete layers exist. '
        'Layers 3-4 will not be included in logs.')
    summary = pd.read_csv(tmp_path/'layers.csv')
    assert summary.complete.tolist() == [True, True, False, False]
    assert summary.included_in_logs.tolist() == [True, True, False, False]


@pytest.mark.parametrize('complete', [True, False])
def test_solvent_info_never_reports_a_partial_or_zero_filled_layer(tmp_path, complete):
    from vorpy.src.output.solvent import export_solvent_layers
    system, net = graph()
    net.balls['complete'] = [complete and i in [2, 3, 4] for i in range(11)]
    calls = []

    class TestGroup:
        def exports(self, **options):
            calls.append(list(self.ball_ndxs))
            (tmp_path / 'group_info.txt').write_text(str(self.ball_ndxs))

    group = TestGroup()
    group.sys, group.net, group.dir = system, net, str(tmp_path)
    group.ball_ndxs = [2, 3, 4, 5, 6, 7]
    group.solvent_layers = dict(atoms=[[2, 3, 4], [5, 6, 7]], surfs=[[], []], edges=[[], []])
    export_solvent_layers(group, [('info', {'info': True})])
    assert group.ball_ndxs == [2, 3, 4, 5, 6, 7]
    if complete:
        assert calls == [[2, 3, 4]]
    else:
        assert not calls
        assert 'summary unavailable' in (tmp_path / 'group_info.txt').read_text()


def test_streamed_incomplete_layer_exports_partial_totals(tmp_path):
    from vorpy.src.output.solvent import export_solvent_layers
    system, net = graph()
    net.balls['complete'] = [False] * len(net.balls)
    calls = []

    class TestGroup:
        def exports(self, **options):
            calls.append(options)

    group = TestGroup()
    group.sys, group.net, group.dir = system, net, str(tmp_path)
    group.ball_ndxs = [2, 3, 4]
    group.vol, group.sa = 12.5, 8.25
    group.int_mean_curv, group.int_gauss_curv = 3.5, 1.25
    group.contained_volume, group.metrics_status = 112.5, 'partial'
    group.solvent_layers = dict(atoms=[[2, 3, 4]], surfs=[[]], edges=[[]],
                                numbers=[1], allow_partial=True)
    export_solvent_layers(group, [('logs', {'logs': True})])
    summary = pd.read_csv(tmp_path/'layers.csv').iloc[0]
    assert calls == [{'logs': True}]
    assert summary['included_in_logs']
    assert summary['metrics_status'] == 'partial'
    assert summary['contained_volume'] == 112.5
    assert summary['layer_volume'] == 12.5
    assert summary['layer_surface_area'] == 8.25


def atom_records(path):
    with path.open(newline='') as stream:
        rows = csv.reader(stream)
        for row in rows:
            if row == ['Atoms']:
                header = next(rows)
                result = []
                for values in rows:
                    if values == ['Surfaces']:
                        return result
                    result.append(dict(zip(header, values)))
    raise AssertionError('Missing Atoms/Surfaces sections')


def test_cli_layer_exports_match_single_sequential_and_parallel_frames(tmp_path, monkeypatch, capsys):
    root = Path(__file__).resolve().parents[3]
    # Real molecular coordinates, bounded fixture size; complete input remains context.
    records = [line for line in (root/'vorpy/data/EDTA.pdb').read_text().splitlines(True)
               if line.startswith(('ATOM  ', 'HETATM'))][:69]
    reference_tree = None
    reference_membership = None
    for mode in ['single', 'sequential', 'parallel']:
        source = tmp_path/mode/'EDTA.pdb'
        source.parent.mkdir()
        source.write_text(''.join(records) if mode == 'single' else
            ''.join('MODEL        %d\n%sENDMDL\n' % (i, ''.join(records)) for i in [1, 2]))
        output = tmp_path/(mode+'_results')
        checkpoint = tmp_path/(mode+'_checkpoint.vpy')
        batch = [] if mode == 'single' else (['--all-frames'] if mode == 'sequential' else ['--parallel-frames', '2'])
        entry = [str(root/'vorpy')] if mode == 'single' else ['-m', 'vorpy']
        result = subprocess.run([sys.executable, *entry, str(source), *batch,
            '--boundary', 'none', '-c', 'l', '4', '-s', 'nt', 'pow', 'and', 'mv', '5', 'and', 'sr', '1',
            '-e', 'large', 'and', 'dir', str(output), '--save-network', str(checkpoint)],
            cwd=root, capture_output=True, text=True, timeout=180)
        assert result.returncode == 0, result.stdout[-3000:] + result.stderr
        directories = [output/'EDTA'] if mode == 'single' else sorted((output/'EDTA/frames').iterdir())
        assert len(directories) == (1 if mode == 'single' else 2)
        for index, directory in enumerate(directories, 1):
            assert (directory/'EDTA_checkpoint.vpy').is_file()
            assert (directory/'groups/EDTA/pow_logs.csv').is_file()
            solvent = directory/'solvent'
            assert len(list(solvent.glob('layer_*/*_logs.csv'))) == 4
            assert len(list(solvent.glob('layer_*/layer_*.pdb'))) == 4
            assert len(list(solvent.glob('layer_*/layer_*_surfs.off'))) == 4
            assert len(list(solvent.glob('layer_*/layer_*_edges.off'))) == 4
            rows = [row for path in sorted(solvent.glob('layer_*/*_logs.csv'))
                    for row in atom_records(path)]
            assert all(row['Residue'] == 'SOL' for row in rows)
            assert all(1 <= int(row['Solvent Layer']) <= 4 for row in rows)
            assert len({row['Index'] for row in rows}) == len(rows)
            from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
            parsed = pd.concat([
                read_logs2(path, all_=False, balls=True)['atoms']
                for path in sorted(solvent.glob('layer_*/*_logs.csv'))
            ], ignore_index=True)
            if rows:
                assert parsed['Solvent Layer'].tolist() == [int(row['Solvent Layer']) for row in rows]
            from vorpy.src.command.vpy_cmnd import Command
            loader = Command()
            archived = loader._load_canonical_checkpoint(directory/'EDTA_checkpoint.vpy')
            assert f'Vertices - {len(archived.verts)} vertices' in (directory/'pow_verts.txt').read_text().splitlines()[0]
            assert not any(getattr(group, 'solvent_layers', None) for group in archived.sys.groups)
            summary = pd.read_csv(solvent/'layers.csv')
            if not summary.complete.all():
                assert 'will be included with explicitly partial totals' in result.stdout
            assert summary.metrics_status.isin(['complete', 'partial']).all()
            for layer, atom_count in zip(summary.layer, summary.atoms):
                pdb = solvent/f'layer_{layer:02d}'/f'layer_{layer:02d}.pdb'
                assert sum(line.startswith(('ATOM  ', 'HETATM')) for line in pdb.read_text().splitlines()) == atom_count
            tree = {str(path.relative_to(directory)) for path in directory.rglob('*')
                    if path.is_file() and path.name != 'run.log'}
            if reference_tree is None:
                reference_tree, reference_membership = tree, summary.atoms.tolist()
            else:
                assert tree == reference_tree
                assert summary.atoms.tolist() == reference_membership
            expected = checkpoint if mode == 'single' else checkpoint.with_name(f'{checkpoint.stem}_frame_{index:04d}.vpy')
            assert expected.is_file()

            loader._cleanup_checkpoint_source()

    # The recovery checkpoint contains the main solve, not released layer topology.
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.network.network import Network
    original_build = Network.build
    def forbidden(*args, **kwargs):
        pytest.fail('Checkpoint export attempted a network solve')
    monkeypatch.setattr(Network, 'build', forbidden)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['vorpy', str(tmp_path/'single_checkpoint.vpy'),
                                    '-e', 'large', 'and', 'dir', str(tmp_path/'reloaded')])
    Command().run()
    assert not (tmp_path/'reloaded/EDTA/solvent').exists()
    monkeypatch.setattr(Network, 'build', original_build)

    checkpoint = tmp_path/'single_checkpoint.vpy'
    original_digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    capsys.readouterr()
    monkeypatch.setattr(sys, 'argv', ['vorpy', str(checkpoint), '-c', 'l', '4',
                                    '-e', 'small', 'and', 'dir', str(tmp_path/'reused_layers')])
    Command().run()
    assert 'Calculating water layers from checkpoint' in capsys.readouterr().out
    assert len(list((tmp_path/'reused_layers/EDTA/solvent').glob('layer_*/layer_*.pdb'))) == 4
    assert hashlib.sha256(checkpoint.read_bytes()).hexdigest() == original_digest


@pytest.mark.parametrize('source_mode, target', [
    ('SOLVENT_CONTEXT', 'MOLECULAR'), ('MOLECULE_ONLY', 'MOLECULAR'), ('SOLVENT_CONTEXT', 'ALL')])
def test_checkpoint_layer_calculation_preserves_source_geometry(tmp_path, monkeypatch, source_mode, target):
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.group import Group
    from vorpy.src.network import Network
    from vorpy.src.output.checkpoint import write_checkpoint
    from vorpy.src.system import System
    root = Path(__file__).resolve().parents[3]
    source = tmp_path/'EDTA.pdb'
    records = [line for line in (root/'vorpy/data/EDTA.pdb').read_text().splitlines(True)
               if line.startswith(('ATOM  ', 'HETATM'))][:69]
    source.write_text(''.join(records))
    saved_dir = tmp_path/'saved'
    saved_dir.mkdir()
    system = System(file=str(source), make_dir=False, solvent_mode=source_mode)
    system.files['dir'] = str(saved_dir)
    body = Group(system, name='EDTA', atoms=list(range(32)), solvent_mode=source_mode, cell_target=target,
        settings=dict(net_type='pow', max_vert=5., surf_res=1., sys_dir=str(saved_dir)))
    system.groups = [body]
    body.build()
    checkpoint = write_checkpoint(system)
    original_digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
    original_build = Network.build
    builds = []

    def tracked_build(network, *args, **kwargs):
        assert network.group_name != 'EDTA', 'The saved solute network must never be rebuilt'
        builds.append(network.group_name)
        return original_build(network, *args, **kwargs)

    monkeypatch.setattr(Network, 'build', tracked_build)
    monkeypatch.setattr(sys, 'argv', ['vorpy', str(checkpoint), '-e', 'small', '-c', 'l', '4',
                                    '-e', 'dir', str(tmp_path/'rerun')])
    command = Command()
    try:
        command.run()
        saved_body = next(group for group in command.sys.groups if group.name == 'EDTA')
        assert saved_body.net is None, 'The exported main solve should be released before layer solves'
        if target == 'ALL':
            assert builds == [], 'A fully cached solvent network needs no new solve'
        else:
            assert builds
            assert ('water_seed_context' in builds) == (source_mode == 'MOLECULE_ONLY')
        directory = tmp_path/'rerun/EDTA/solvent'
        assert len(list(directory.glob('layer_*/*_logs.csv'))) == 4
        assert len(list(directory.glob('layer_*/layer_*.pdb'))) == 4
        summary = pd.read_csv(directory/'layers.csv')
        assert {'contained_volume', 'layer_volume', 'layer_surface_area',
                'integrated_mean_curvature', 'integrated_gaussian_curvature',
                'metrics_status'}.issubset(summary.columns)
        assert hashlib.sha256(checkpoint.read_bytes()).hexdigest() == original_digest
    finally:
        command._cleanup_checkpoint_source()


def test_virtual_boundary_is_not_a_water_layer():
    system, net = graph()
    net.balls = pd.concat([net.balls, pd.DataFrame({'num': [111], 'system_num': [-1]})], ignore_index=True)
    net.surfs.loc[5] = [[100, 111], 1., [], []]
    assert water_neighbors(net, system, {0}, {0}) == {2, 3, 4}
    assert boundary_indices(net, [0], [[2, 3, 4]])[0] == [[0]]


def test_geometry_only_frontier_is_not_reused_as_analyzed_solvent(tmp_path, monkeypatch):
    from copy import copy
    import vorpy.src.group.solvent_layers as module
    system, source = graph()
    source.surfs = source.surfs.loc[[0]].copy()
    source.settings = {'net_type': 'pow'}
    source.group = [100]
    system.files = {'dir': str(tmp_path)}
    system.groups = [SimpleNamespace(net=source, ball_ndxs=[0], settings=source.settings,
                                    name='solute', group_id='solute', dir=None)]
    builds = []

    class TestGroup:
        def __init__(self, system, name, atoms, settings, **kwargs):
            self.ball_ndxs, self.settings, self.name = list(atoms), settings, name
            self.net = None

        def build(self, *, calculate_curvature):
            builds.append((self.name, calculate_curvature))
            self.net = copy(source)
            self.net.balls = source.balls.copy()
            self.net.group = source.balls.iloc[self.ball_ndxs]['num'].tolist()
            if calculate_curvature:
                for column in ['vol', 'vdw_vol', 'com', 'moi', 'int_mean_curv', 'int_gauss_curv']:
                    self.net.balls[column] = 0.

    monkeypatch.setattr(module, 'Group', TestGroup)
    result = module.build_solvent_layers(system, 4)
    assert builds == [('water_frontier_1', False), ('solvent', True)]
    assert result.solvent_layers['atoms'] == [[2, 3, 4], [], [], []]
