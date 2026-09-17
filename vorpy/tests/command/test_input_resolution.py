from pathlib import Path
from types import SimpleNamespace

import pytest

from vorpy.src.command.interpret import get_file
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.inputs.pdb import read_pdb


@pytest.mark.parametrize('name', ['EDTA.pdb', 'edta'])
@pytest.mark.parametrize('all_frames', [False, True])
def test_cli_resolves_bundled_input(tmp_path, monkeypatch, name, all_frames):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('sys.argv', ['vorpy', name] + (['--all-frames'] if all_frames else []))
    command = Command(sys=SimpleNamespace())
    seen = []
    monkeypatch.setattr(command, 'parse_commands', lambda: None)
    monkeypatch.setattr(command, '_process_system', lambda: seen.append(command.base_file))
    monkeypatch.setattr(command, 'run_all_frames', lambda: seen.append(command.base_file))
    command.run()
    assert len(seen) == 1
    assert Path(seen[0]).is_file()
    assert Path(seen[0]).name == 'EDTA.pdb'


def test_local_file_takes_precedence(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    local = tmp_path / 'EDTA.pdb'
    local.write_text('END\n')
    assert get_file('EDTA.pdb', interactive=False) == str(local.resolve())


def test_missing_cli_file_exits_without_prompt_or_build(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('sys.argv', ['vorpy', 'missing_trajectory_123.pdb'])
    monkeypatch.setattr('builtins.input', lambda *args: pytest.fail('Unexpected prompt'))
    with pytest.raises(SystemExit, match='Input file not found'):
        Command().run()


def test_missing_pdb_raises_in_reader(tmp_path):
    system = SimpleNamespace(files={'vpy_dir': None, 'dir': None})
    with pytest.raises(FileNotFoundError, match='PDB input file not found'):
        read_pdb(system, str(tmp_path / 'missing.pdb'))
