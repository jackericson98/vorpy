from types import SimpleNamespace
from contextlib import closing
from pathlib import Path

import numpy as np
import pytest

from vorpy.src.input_progress import iter_input_lines
from vorpy.src.inputs.pdb import read_pdb
from vorpy.src.inputs.frames import iter_pdb_frames
from vorpy.src.inputs.gro import read_gro
from vorpy.src.inputs.mol import read_mol
from vorpy.src.inputs.mol2 import read_mol2


def system(path):
    return SimpleNamespace(files={'base_file': str(path), 'vpy_dir': None, 'dir': None},
                           type='mol', sol=None, print_actions=False)


def test_progress_is_bounded_and_tracks_bytes(tmp_path):
    path = tmp_path / 'lines.txt'
    path.write_text('abc\n' * 10000)
    events = []
    assert len(list(iter_input_lines(path, lambda label, value: events.append(value)))) == 10000
    assert events[0] == 0 and events[-1] == 100
    assert events == sorted(events)
    assert len(events) < 10


def test_pdb_streaming_preserves_hierarchy_and_silent_mode(tmp_path, capsys):
    path = tmp_path / 'tiny.pdb'
    path.write_text('REMARK\n' + ''.join(
        f'ATOM  {i:5d} {name:^4s} ILE A   1    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {element:>2s}\n'
        for i, name, x, element in [(1, 'N', 0., 'N'), (2, 'CA', 1., 'C')]))
    sys = system(path)
    read_pdb(sys)
    assert len(sys.balls) == 2
    assert sys.chains[0].atoms == [0, 1]
    assert sys.residues[0].atoms == [0, 1]
    assert capsys.readouterr().out == ''


def test_pdb_skips_blank_element_water_virtual_sites(tmp_path, monkeypatch):
    path = tmp_path / 'water.pdb'
    path.write_text(''.join(
        f'ATOM  {i:5d} {name:^4s} SOL A   1    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00\n'
        for i, name, x in [(1, 'OW', 0.), (2, 'MW', .1),
                           (3, 'HW1', .9), (4, 'HW2', -.9)]))
    sys = system(path)
    read_pdb(sys)
    assert sys.balls['name'].tolist() == ['OW', 'HW1', 'HW2']
    assert sys.balls['num'].tolist() == [0, 1, 2]
    assert sys.balls['element'].tolist() == ['O', 'H', 'H']
    assert all(sys.balls['rad'] > 0)
    from vorpy.src.output.pdb import write_pdb
    monkeypatch.chdir(tmp_path)
    sys.frame_count = 2
    write_pdb([0, 1, 2], 'exported', sys)
    names = [line[12:16].strip() for line in (tmp_path / 'exported.pdb').read_text().splitlines()
             if line.startswith('ATOM')]
    assert names == ['OW', 'HW1', 'HW2']


def test_failed_frame_load_closes_input_before_temporary_cleanup(tmp_path):
    source = tmp_path / 'unknown.pdb'
    source.write_text('ATOM      1  XX  UNK A   1       1.000   0.000   0.000  1.00  0.00          XX\nEND\n')
    with pytest.raises(ValueError, match="No radius available for atom 'XX'") as error:
        with closing(iter_pdb_frames(source)) as frames:
            _, frame_path = next(frames)
            read_pdb(system(frame_path))
    # Keep the exception/traceback alive to catch handles retained by frames.
    assert error.value.__traceback__ is not None
    assert not Path(frame_path).parent.exists()


def test_gro_fixed_fields_units_and_indices(tmp_path):
    path = tmp_path / 'tiny.gro'
    path.write_text('Title\n2\n' + ''.join(
        f'{42:5d}{"ABCDE":<5s}{name:>5s}{serial:5d}{x:8.3f}{0.:8.3f}{0.:8.3f}\n'
        for name, serial, x in [('C1', 99999, 0.1), ('H1', 0, 0.2)]) + '1.0 1.0 1.0\n')
    sys = system(path)
    read_gro(sys)
    assert sys.balls['res_name'].tolist() == ['ABCDE', 'ABCDE']
    assert sys.balls['num'].tolist() == [0, 1]
    assert sys.balls['gro_id'].tolist() == [99999, 0]
    assert sys.balls['res_seq'].tolist() == [42, 42]
    np.testing.assert_allclose(np.stack(sys.balls['loc']), [[1, 0, 0], [2, 0, 0]])
    path.write_text('Title\n2\n')
    with pytest.raises(ValueError, match='atom record 1'):
        read_gro(system(path))


def test_mol_v3000_continuations_and_bonds(tmp_path):
    path = tmp_path / 'tiny.mol'
    path.write_text('Title\nProgram\n\n  0  0  0     0  0            999 V3000\n'
                    'M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\n'
                    'M  V30 BEGIN ATOM\nM  V30 10 C 0 0 0 0 -\nM  V30 CHG=1\n'
                    'M  V30 20 O 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n'
                    'M  V30 1 1 10 20\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n')
    sys = system(path)
    read_mol(sys)
    assert sys.balls['element'].tolist() == ['C', 'O']
    assert sys.balls.iloc[0]['charge'] == '1'
    assert sys.bonds == [(0, 1, '1')]


def test_mol_v2000_streaming_dispatch(tmp_path):
    path = tmp_path / 'tiny.mol'
    path.write_text('Title\nProgram\n\n  2  1  0  0  0  0            999 V2000\n' + ''.join(
        f'{x:10.4f}{0.:10.4f}{0.:10.4f} {element:<3s} 0  0  0  0\n'
        for x, element in [(0., 'C'), (1., 'O')]) + '  1  2  1  0\nM  END\n')
    sys = system(path)
    read_mol(sys)
    assert sys.balls['element'].tolist() == ['C', 'O']
    assert sys.bonds == [(0, 1, '1')]


def test_mol2_two_pass_metadata_and_progress(tmp_path):
    path = tmp_path / 'tiny.mol2'
    path.write_text('@<TRIPOS>MOLECULE\nTiny\n2 1 1\nSMALL\nUSER_CHARGES\n'
                    '@<TRIPOS>ATOM\n10 C1 0 0 0 C.3 1 ILE7 0.2\n20 O1 1 0 0 O.2 1 ILE7 -0.2\n'
                    '@<TRIPOS>BOND\n1 10 20 1\n'
                    '@<TRIPOS>SUBSTRUCTURE\n1 ILE7 10 RESIDUE 1 A ILE\n')
    sys = system(path)
    events = []
    read_mol2(sys, progress=lambda label, value: events.append(value))
    assert sys.bonds == [(0, 1, '1')]
    assert sys.residues[0].atoms == [0, 1]
    assert sys.chains[0].name == 'A'
    assert events == sorted(events) and events[-1] == 100
