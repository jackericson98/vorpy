from pathlib import Path
from types import SimpleNamespace

import pytest

from vorpy.src.boundary import SolventMode
from vorpy.src.command.vpy_cmnd import Command, _source_composition
from vorpy.src.system.system import System


SOURCE = Path(__file__).resolve().parents[2] / 'data' / '1bna_DB1976.pdb'


def _command(arguments):
    system = System(file=str(SOURCE), make_dir=False)
    command = Command(sys=system)
    command.base_file = str(SOURCE)
    command._arguments = [str(SOURCE), *arguments]
    command.parse_commands()
    return command, system


def test_exporting_logs_does_not_select_log_reconstruction():
    command, _ = _command([
        '-g', 'c A and c B', '-g', 'c C', '-i',
        '-s', 'nt', 'aw', 'and', 'sr', '0.2', 'and', 'mv', '5',
        '-e', 'surfaces', '-e', 'edges', '-e', 'vertices', '-e', 'dual',
        '-e', 'logs', '-e', 'info',
    ])

    command.apply_settings()

    assert command.exports[-2:] == [['logs'], ['info']]
    assert command.logs_files == []
    assert command.settings_dict['bld_type'] is None


def test_explicit_log_build_keeps_input_validation(tmp_path):
    log_file = tmp_path / 'missing.csv'
    command, system = _command(['-b', 'logs', str(log_file)])
    command.apply_settings()
    system.groups = [SimpleNamespace(net=object())]

    assert command.settings_dict['bld_type'] == 'logs'
    assert command.logs_files == [str(log_file)]
    with pytest.raises(FileNotFoundError, match='Logs file not found'):
        command.build_groups_from_logs()


def test_1bna_model1_interface_uses_full_source_default():
    command, system = _command(['-g', 'c A and c B', '-g', 'c C', '-i'])
    command.apply_settings()
    command._resolve_default_solvent_mode()
    command._build_group_networks = False
    command.create_groups()

    # The bounded fixture is currently a single-frame MODEL 1 input.  The
    # parser must preserve that selected frame without assuming the historical
    # 11-frame source file.
    assert system.frame_count >= 1
    assert system.loaded_frame_count == 1
    assert system.frame_index == 1

    composition = _source_composition(system)
    assert composition['molecular_atoms'] == 805
    assert composition['molecular_chains'] == 3
    assert composition['water_atoms'] == 9669
    assert composition['ion_atoms'] == 20

    groups = system.groups
    assert [len(group.ball_ndxs) for group in groups] == [758, 47]
    assert {
        tuple(sorted(set(system.balls.loc[group.ball_ndxs, 'chain_name'])))
        for group in groups
    } == {('A', 'B'), ('C',)}

    assert command.solvent_mode.value == 'FULL_SOLVENT'
    assert system.solvent_mode.value == 'FULL_SOLVENT'
    assert command.settings_dict['solvent_mode_resolution'] == (
        'implicit_interface_full_solvent'
    )
    source = system.balls.loc[list(command._source_system_indices)]
    assert len(source) == 10494
    assert command._source_generator_policy['excluded_generator_count'] == 0
    assert command._source_generator_policy['active_generator_count'] == 10494
    ligand = system.balls[system.balls['chain_name'].astype(str).str.strip() == 'C']
    selenium = ligand[ligand['element'].astype(str).str.strip().str.upper() == 'SE']
    assert len(selenium) == 1
    assert selenium['rad'].iloc[0] == pytest.approx(1.9)


def test_1bna_interface_explicit_molecule_only_remains_dry():
    command, system = _command(['-g', 'c A and c B', '-g', 'c C', '-i'])
    command.solvent_mode = SolventMode.MOLECULE_ONLY
    system.solvent_mode = SolventMode.MOLECULE_ONLY
    command._solvent_mode_explicit = True

    command.apply_settings()
    command._resolve_default_solvent_mode()

    assert command.solvent_mode is SolventMode.MOLECULE_ONLY
    assert len(command._source_system_indices) == 805
    assert command.settings_dict['solvent_mode_resolution'] == 'explicit'


def test_normal_interface_build_path_is_used_for_logs_export(tmp_path, capsys):
    command, system = _command([
        '-g', 'c A and c B', '-g', 'c C', '-i', '-e', 'logs',
    ])
    command.apply_settings()
    system.files['dir'] = str(tmp_path)
    built = []

    class StubGroup:
        name = 'stub'
        net = None

        def make_net(self):
            built.append('make_net')

        def build(self, **kwargs):
            built.append('build')

    command.load_files = lambda: None
    command._build_full_system_target = lambda: None
    command.create_groups = lambda: setattr(system, 'groups', [StubGroup(), StubGroup()])
    command.run_exports = lambda: None
    command._save_archive = lambda: None
    command.build_groups_from_logs = lambda: pytest.fail('log reconstruction selected')
    system.make_interfaces = lambda pairs, *, analyze_waters: None

    command._process_system()
    output = capsys.readouterr().out

    assert built == ['make_net', 'build', 'make_net', 'build']
    assert 'Solving AW network...' in output
    assert 'Extracting ' not in output
    assert 'Source generators:' not in output
    assert 'Active groups:' not in output
    assert 'Interface Result status:' not in output
    if system.frame_count > 1:
        assert f'Frame: 1/{system.frame_count} loaded' in output
    assert 'Input context:' not in output
    assert 'frames found' not in output
