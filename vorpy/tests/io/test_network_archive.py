"""Real solved-network persistence and hostile/malformed archive regression tests."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import zipfile

import numpy as np
import pytest

from vorpy.src.io import ArchiveError, group_from_network, load_network, save_network
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


@pytest.fixture(scope='module')
def solved(tmp_path_factory):
    directory = tmp_path_factory.mktemp('archive-solve')
    source = directory / 'molecule.pdb'
    shutil.copy(Path(__file__).resolve().parents[2] / 'data' / 'EDTA.pdb', source)
    cwd = Path.cwd()
    try:
        result = VorPyBackend(VorPySolveSettings(network_type='pow')).solve(source, lambda *a: None, lambda: False)
        yield result
    finally:
        os.chdir(cwd)


@pytest.fixture
def archive(solved, tmp_path):
    path = tmp_path / 'network.vpy'
    save_network(solved.export_group.net, path)
    return path


def same(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            same(a[key], b[key])
    elif isinstance(a, (list, tuple, set, np.ndarray)):
        if isinstance(a, set):
            assert a == b
        else:
            assert len(a) == len(b)
            for left, right in zip(a, b):
                same(left, right)
    elif isinstance(a, (float, np.floating)):
        assert float(a) == pytest.approx(float(b), nan_ok=True)
    else:
        assert a == b


def test_round_trip_all_scientific_columns(solved, archive, monkeypatch):
    from vorpy.src.network import Network
    from vorpy.src.group import Group
    def forbidden(*a, **kw):
        pytest.fail('Load attempted to build or analyze')
    for name in ('build', 'build_edges', 'build_surfaces', 'find_verts'):
        monkeypatch.setattr(Network, name, forbidden)
    monkeypatch.setattr(Group, 'get_info', forbidden)
    loaded = load_network(archive)
    original = solved.export_group.net
    for key in ('balls', 'surfs', 'edges', 'verts'):
        left, right = getattr(original, key), getattr(loaded, key)
        assert list(left.columns) == list(right.columns)
        assert left.index.tolist() == right.index.tolist()
        for column in left.columns:
            same(left[column].tolist(), right[column].tolist())
    same(original.settings, loaded.settings)
    same(original.box, loaded.box)
    same(original.sys.balls['loc'].tolist(), loaded.sys.balls['loc'].tolist())
    same(original.sys.balls['rad'].tolist(), loaded.sys.balls['rad'].tolist())
    assert loaded.sys.groups[0].net is loaded
    assert loaded.sys.balls.iloc[0]['res'] is loaded.sys.residues[0]
    assert loaded.sys.residues[0].sys is loaded.sys
    assert loaded.sys.balls.iloc[0]['chn'] in loaded.sys.chains
    assert loaded.sys.files['base_file'] is None


def test_export_equivalence_and_no_source(solved, archive, tmp_path):
    from vorpy.src.output.surfs import write_surfs
    from vorpy.src.output.edges import write_edges
    from vorpy.src.output.verts import write_off_verts
    from vorpy.src.output.pdb import write_pdb
    original = solved.export_group.net
    loaded = load_network(archive)
    old_source = solved.source
    hidden_source = old_source.with_suffix('.hidden')
    old_source.rename(hidden_source)
    cwd = Path.cwd()
    try:
        # Export all supported mesh formats from the same scientific values.
        for fmt in ('off', 'ply', 'vtp'):
            dirs = []
            for index, network in enumerate((original, loaded)):
                directory = tmp_path / f'{fmt}-{index}'
                directory.mkdir()
                dirs.append(directory)
                for name, function in (('surfs', write_surfs), ('edges', write_edges), ('verts', write_off_verts)):
                    function(network, [0, 1], name, directory=str(directory), file_type=fmt)
            for file in dirs[0].iterdir():
                assert file.read_bytes() == (dirs[1] / file.name).read_bytes()
        group = loaded.sys.groups[0]
        loaded.sys.files['dir'] = str(tmp_path)
        group.dir = str(tmp_path / 'group-export')
        Path(group.dir).mkdir()
        group.exports(atoms=True, info=True, logs=True, shell_surfs=True, shell_edges=True,
                      shell_verts=True, surr_atoms=True, surr_resids=True, atom_surfs=True,
                      atom_edges=True, atom_verts=True, sep_verts=True, file_type='ply')
        assert (Path(group.dir) / 'group_atoms.pdb').stat().st_size > 0
        assert (Path(group.dir) / 'info.txt').is_file()
        assert list(Path(group.dir).glob('*logs.csv'))
        assert list((Path(group.dir) / 'atoms').glob('*.ply'))
        assert list((Path(group.dir) / 'verts').glob('*.ply'))
        # A new complete-cell selection reuses existing geometry.
        selection = group_from_network(loaded, [0, 1], 'new_selection')
        assert selection.net.surfs is loaded.surfs
        selection.exports(atoms=True, shell_surfs=True, file_type='off')
        assert (Path(selection.dir) / 'shell_surfs.off').is_file()
    finally:
        os.chdir(cwd)
        hidden_source.rename(old_source)


def rewrite(path, member, transform):
    with zipfile.ZipFile(path) as z:
        files = {name: z.read(name) for name in z.namelist()}
    files[member] = json.dumps(transform(json.loads(files[member]))).encode()
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, data in files.items():
            z.writestr(name, data)


@pytest.mark.parametrize('field,value,match', [('format','wrong','Not a VorPy'),
                                               ('format_version',999,'Unsupported'),
                                               ('root_network',999,'root network')])
def test_invalid_metadata(archive, field, value, match):
    rewrite(archive, 'metadata.json', lambda data: dict(data, **{field:value}))
    with pytest.raises(ArchiveError, match=match):
        load_network(archive)


def test_count_mismatch(archive):
    rewrite(archive, 'metadata.json', lambda data: dict(data, counts=dict(data['counts'], verts=99)))
    with pytest.raises(ArchiveError, match='count mismatch'):
        load_network(archive)


def test_dangling_object_reference(archive):
    def corrupt(state):
        state['objects'][0]['fields']['sys'] = {'ref': 999999}
        return state
    rewrite(archive, 'state.json', corrupt)
    with pytest.raises(ArchiveError, match='Dangling object'):
        load_network(archive)


def test_missing_array(archive):
    with zipfile.ZipFile(archive) as z:
        files = {name:z.read(name) for name in z.namelist()}
    del files[next(name for name in files if name.startswith('arrays/'))]
    with zipfile.ZipFile(archive,'w') as z:
        for name,data in files.items():
            z.writestr(name,data)
    with pytest.raises(ArchiveError, match='array'):
        load_network(archive)


def test_process_exit_load_export(archive, tmp_path):
    import sys
    code = '''from vorpy.src.io import load_network
from pathlib import Path
import sys
n=load_network(sys.argv[1]); g=n.sys.groups[0]
g.exports(atoms=True,shell_surfs=True,verts=True,file_type='off')
assert (Path(g.dir)/'group_atoms.pdb').is_file()
'''
    completed = subprocess.run([sys.executable, '-c', code, str(archive)], capture_output=True, text=True, timeout=30)
    assert completed.returncode == 0, completed.stderr


@pytest.mark.timeout(120)
def test_aw_curvature_and_interface_records_round_trip(tmp_path):
    from vorpy.src.interface import Interface
    source = Path(__file__).resolve().parents[2] / 'data' / 'EDTA.pdb'
    cwd = Path.cwd()
    try:
        result = VorPyBackend().solve(source, lambda *a:None, lambda:False)
        group = result.export_group
        group.get_info()
        # Definitions and contact relationships on the existing solved topology.
        left = group_from_network(group.net, [0, 1], 'left')
        right = group_from_network(group.net, [2, 3], 'right')
        interface = Interface(group.sys, left, right)
        interface.net = group.net
        interface.group1_indices = {0, 1}
        interface.group2_indices = {2, 3}
        interface.water_topology = {'waters': {}, 'test_reference': group.sys.residues[0]}
        group.sys.ifaces.append(interface)
        archive = save_network(group.net, tmp_path / 'aw.vpy')
        loaded = load_network(archive)
        for table_name in ('balls','surfs','edges','verts'):
            original = getattr(group.net, table_name)
            restored = getattr(loaded, table_name)
            for column in original.columns:
                if 'curv' in column or 'energy' in column or 'angle' in column:
                    same(original[column].tolist(), restored[column].tolist())
        loaded_interface = loaded.sys.ifaces[0]
        assert loaded_interface.net is loaded
        assert loaded_interface.group1.interfaces[interface.interface_id] is loaded_interface
        assert loaded_interface.water_topology['test_reference'] is loaded.sys.residues[0]
        same(group.vol, loaded.sys.groups[0].vol)
        same(group.int_mean_curv, loaded.sys.groups[0].int_mean_curv)
        loaded.sys.set_output_directory(str(tmp_path / 'interface_output'))
        loaded_interface.export(atoms=True, surfs=True, edges=True, verts=True, info=True)
        assert (Path(loaded_interface.dir) / 'interface_atoms.pdb').is_file()
    finally:
        os.chdir(cwd)


def test_cli_archive_input_without_build(archive, monkeypatch, tmp_path):
    import sys
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.network import Network
    monkeypatch.setattr(Network, 'build', lambda *a, **kw: pytest.fail('CLI tried to rebuild'))
    monkeypatch.setattr(sys, 'argv', ['vorpy', '--load-network', str(archive), '-e', 'only', 'logs',
                                    '--save-network', str(tmp_path / 'resaved.vpy')])
    cwd = Path.cwd()
    try:
        command = Command()
        command.run()
        assert (tmp_path / 'resaved.vpy').is_file()
        assert list(Path(command.sys.files['dir']).rglob('*logs.csv'))
    finally:
        os.chdir(cwd)
