from pathlib import Path
from types import SimpleNamespace
import weakref

import pytest

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.inputs.frames import iter_pdb_frames
from vorpy.src.inputs.pdb import read_pdb
from vorpy.src.output.pdb import write_pdb
from vorpy.src.network.network import Network


def atom(x):
    return (f'ATOM      1  CA  ALA A   1    {x:8.3f}{0.:8.3f}{0.:8.3f}'
            '  1.00  0.00           C\n')


@pytest.mark.parametrize('delimiter', ['ENDMDL', 'END'])
def test_stream_frames_and_cleanup(tmp_path, delimiter):
    source = tmp_path / 'trajectory.pdb'
    source.write_text('TITLE first\nMODEL        1\n' + atom(1) + delimiter + '\n'
                      'TITLE second\nMODEL        2\n' + atom(2))
    frames = iter_pdb_frames(source)
    number, path = next(frames)
    assert number == 1
    assert atom(1) in Path(path).read_text()
    assert atom(2) not in Path(path).read_text()
    number, path = next(frames)
    assert number == 2
    assert 'TITLE second' in Path(path).read_text()
    assert atom(1) not in Path(path).read_text()
    assert atom(2) in Path(path).read_text()
    frames.close()
    assert not Path(path).exists()


def test_empty_trajectory(tmp_path):
    source = tmp_path / 'empty.pdb'
    source.write_text('HEADER empty\nEND\n')
    with pytest.raises(ValueError, match='No atom frames'):
        list(iter_pdb_frames(source))


@pytest.mark.parametrize('delimiter', ['ENDMDL', 'END'])
def test_default_reader_and_export_use_only_frame_one(tmp_path, monkeypatch, delimiter):
    source = tmp_path / 'trajectory.pdb'
    source.write_text('MODEL        1\n' + atom(1) + delimiter + '\n'
                      'MODEL        2\n' + atom(2) + delimiter + '\n')
    system = SimpleNamespace(files={'base_file': str(source), 'vpy_dir': None, 'dir': None},
                             type='mol', sol=None, print_actions=False)
    read_pdb(system)
    assert system.frame_count == 2
    assert system.loaded_frame_count == 1
    assert system.frame_index == 1
    assert len(system.balls) == 1
    assert system.balls.iloc[0]['loc'][0] == 1.
    monkeypatch.chdir(tmp_path)
    write_pdb([0], 'result', system)
    output = (tmp_path / 'result.pdb').read_text()
    assert atom(1) in output
    assert atom(2) not in output


def test_network_progress_identifies_frame():
    events = []
    system = SimpleNamespace(frame_count=11, frame_index=1,
                             update_progress=lambda **kwargs: events.append(kwargs))
    net = SimpleNamespace(sys=system, progress_process_prefix=None,
                          progress_network_name=None, group_name='trajectory')
    Network.update_progress(net, 'Vertices: 27', 0.4)
    assert events[0]['network'] == 'trajectory - Frame 1/11'


@pytest.mark.parametrize('explicit', [True, False])
def test_batch_uses_independent_systems_and_standard_exports(tmp_path, monkeypatch, capsys, explicit):
    source = tmp_path / 'trajectory.pdb'
    source.write_text('MODEL        1\n' + atom(1) + 'ENDMDL\n'
                      'MODEL        2\n' + atom(2) + 'ENDMDL\n')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('builtins.input', lambda prompt: pytest.fail('Unexpected prompt') if explicit else '')
    monkeypatch.setattr('sys.argv', ['vorpy', str(source)] + (['--all-frames'] if explicit else []) + [
                                   '-e', 'dir', str(tmp_path / 'results'),
                                   '-e', 'only', 'pdb', '-s', 'nt', 'pow'])
    seen = []
    monkeypatch.setattr(Command, 'create_groups', lambda self: None)
    original_export = Command.run_exports

    def export(command):
        seen.append((command.sys.balls.iloc[0]['loc'][0],
                     command.settings_dict['net_type'], command.sys.files['dir']))
        original_export(command)

    monkeypatch.setattr(Command, 'run_exports', export)
    Command().run()
    output = capsys.readouterr().out
    assert '2 frames found' in output
    assert '1 frame loaded (frame 1 only)' in output
    assert '1 frame loaded (frame 2 only)' in output
    assert [item[:2] for item in seen] == [(1., 'pow'), (2., 'pow')]
    for index, (_, _, directory) in enumerate(seen, 1):
        assert Path(directory) == tmp_path / 'results' / 'trajectory' / 'frames' / f'frame_{index:04d}'
        assert (Path(directory) / 'trajectory.pdb').exists()


@pytest.mark.parametrize('answers, expected', [(['y'], 'all'), (['YES'], 'all'),
                                             ([''], 'all'), (['n'], 'first'),
                                             (['no'], 'first'), (['maybe', 'y'], 'all')])
def test_frame_prompt_routes_answer(tmp_path, monkeypatch, answers, expected):
    source = tmp_path / 'trajectory.pdb'
    source.write_text(atom(1) + 'END\n' + atom(2) + 'END\n')
    monkeypatch.setattr('sys.argv', ['vorpy', str(source)])
    replies = iter(answers)
    prompts = []
    def respond(prompt):
        prompts.append(prompt)
        return next(replies)
    monkeypatch.setattr('builtins.input', respond)
    command = Command(sys=SimpleNamespace())
    actions = []
    monkeypatch.setattr(command, 'run_all_frames', lambda total_frames: actions.append(('all', total_frames)))
    monkeypatch.setattr(command, 'parse_commands', lambda: None)
    monkeypatch.setattr(command, '_process_system', lambda: actions.append(('first', 1)))
    command.run()
    assert actions == [(expected, 2 if expected == 'all' else 1)]
    assert prompts == ['2 frames found. Run all? [Y/n] '] * len(answers)


def test_batch_collects_previous_frame_before_loading_next(tmp_path, monkeypatch):
    from vorpy.src.command import vpy_cmnd

    source = tmp_path / 'trajectory.pdb'
    source.write_text(atom(1) + 'END\n' + atom(2) + 'END\n')
    monkeypatch.setattr('sys.argv', ['vorpy', str(source), '--all-frames',
                                   '-e', 'dir', str(tmp_path / 'results')])
    real_system = vpy_cmnd.System
    references = []

    def load_frame(*args, **kwargs):
        # Check before constructing frame 2, not merely after the whole run.
        assert all(ref() is None for ref in references)
        system = real_system(*args, **kwargs)
        references.append(weakref.ref(system))
        return system

    def process(command):
        # Model the cycles retained by groups, networks and their caches.
        command.sys.test_network = SimpleNamespace(system=command.sys, geometry=bytearray(1024))

    monkeypatch.setattr(vpy_cmnd, 'System', load_frame)
    monkeypatch.setattr(Command, '_process_system', process)
    command = Command()
    command.run()
    assert len(references) == 2
    assert all(ref() is None for ref in references)
    assert command.sys is None
