import os
from types import SimpleNamespace

import pytest

from vorpy.workbench.domain import AnalysisResult
from vorpy.workbench.services.export import export_result


def test_export_presets_restore_state_and_keep_previous_exports(tmp_path):
    calls = []
    group = SimpleNamespace(name='sample', dir=None, sys=SimpleNamespace(files={'dir': 'original'}))
    def write(**options):
        calls.append(options)
        os.chdir(group.dir)
        (tmp_path / 'called').touch()
    group.exports = write
    result = AnalysisResult(None, 'sample', export_group=group)
    cwd = os.getcwd()
    first = export_result(result, tmp_path, 'Small', {}, 4, 'ply')
    second = export_result(result, tmp_path, 'Medium', {}, 3, 'off')
    assert first != second and first.is_dir() and second.is_dir()
    assert calls[0] == dict(info=True, shell_surfs=True, logs=True, round_to=4, file_type='ply')
    assert calls[1]['surr_resids'] and calls[1]['verts']
    assert os.getcwd() == cwd
    assert group.dir is None and group.sys.files['dir'] == 'original'


def test_export_failure_restores_state(tmp_path):
    group = SimpleNamespace(name='sample', dir=None, sys=SimpleNamespace(files={'dir': 'original'}))
    def fail(**options):
        os.chdir(group.dir)
        raise RuntimeError('write failed')
    group.exports = fail
    result = AnalysisResult(None, 'sample', export_group=group)
    cwd = os.getcwd()
    with pytest.raises(RuntimeError, match='write failed'):
        export_result(result, tmp_path, 'Custom', {'edges': True}, 3, 'vtp')
    assert os.getcwd() == cwd
    assert group.dir is None and group.sys.files['dir'] == 'original'


def test_imported_result_requires_solve(tmp_path):
    with pytest.raises(ValueError, match='Solve'):
        export_result(AnalysisResult(None, 'imported'), tmp_path, 'Small', {}, 3, 'off')


def test_atomic_formats_full_molecule_and_subset(tmp_path):
    from vorpy.workbench.domain import Atom
    from vorpy.workbench.services.atomic_export import export_atomic_selections
    atoms = [Atom(0, 1, 'C', 'C', (1., 2., 3.)), Atom(1, 2, 'O', 'O', (4., 5., 6.))]
    result = AnalysisResult(None, 'example', atoms=atoms,
                            export_group=SimpleNamespace(ball_ndxs=[1]))
    for fmt in ('pdb', 'xyz'):
        export_atomic_selections(result, {'atoms': True, 'full_molecule': True}, tmp_path, fmt)
    assert (tmp_path / 'group_atoms.xyz').read_text().splitlines()[0] == '1'
    assert (tmp_path / 'full_molecule.xyz').read_text().splitlines()[0] == '2'
    assert sum(line.startswith('ATOM') for line in (tmp_path / 'full_molecule.pdb').read_text().splitlines()) == 2
    assert 'O' in (tmp_path / 'group_atoms.pdb').read_text()


def test_atomic_pdb_uses_current_frame(tmp_path):
    from vorpy.workbench.domain import Atom
    from vorpy.src.output.pdb import make_pdb_line
    from vorpy.workbench.services.atomic_export import write_atomic_structure
    source = tmp_path / 'frames.pdb'
    source.write_text('MODEL        1\n' + make_pdb_line(ser_num=1, name='C', x=1.) +
                      'ENDMDL\nMODEL        2\n' + make_pdb_line(ser_num=1, name='C', x=9.) + 'ENDMDL\n')
    result = AnalysisResult(source, 'frames', atoms=[Atom(0, 1, 'C', 'C', (9., 0., 0.))], frame_index=2)
    write_atomic_structure(result, [0], tmp_path, 'current', 'pdb')
    atom_line = next(line for line in (tmp_path / 'current.pdb').read_text().splitlines() if line.startswith('ATOM'))
    assert float(atom_line[30:38]) == 9.


def test_export_color_provider_restored_after_failure(tmp_path, monkeypatch):
    from vorpy.workbench.services import export_colors
    provider = object()
    monkeypatch.setattr(export_colors, 'make_color_provider', lambda *args: provider)
    net = SimpleNamespace()
    def fail(**kwargs):
        assert net._export_color_provider is provider
        assert 'full_molecule' not in kwargs
        raise RuntimeError('failed')
    group = SimpleNamespace(name='sample', net=net, dir=None,
                            sys=SimpleNamespace(files={'dir': 'original'}), exports=fail)
    result = AnalysisResult(None, 'sample', export_group=group)
    with pytest.raises(RuntimeError):
        export_result(result, tmp_path, 'Custom', {'full_molecule': True}, 3, 'off',
                      atom_format='xyz', color_overrides={})
    assert not hasattr(net, '_export_color_provider')
    assert group.dir is None
