import sys
from pathlib import Path

import pytest

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.inputs.frames import count_pdb_frames


ROOT = Path(__file__).resolve().parents[2]
WET = ROOT / 'data' / '1bna_DB1976.pdb'
DRY = ROOT / 'data' / '1bna_DB1976_no_sol.pdb'


@pytest.mark.parametrize(
    ('run_name', 'source', 'network_type', 'solvent_mode', 'active_generators'),
    (
        ('aw_wet', WET, 'aw', 'FULL_SOLVENT', 10494),
        ('pow_wet', WET, 'pow', 'FULL_SOLVENT', 10494),
        ('aw_dry', DRY, 'aw', 'MOLECULE_ONLY', 805),
        ('pow_dry', DRY, 'pow', 'MOLECULE_ONLY', 805),
    ),
)
def test_benchmark_cli_matrix_is_parser_valid_and_output_isolated(
    tmp_path, monkeypatch, capsys, run_name, source, network_type, solvent_mode,
    active_generators,
):
    if not source.is_file():
        pytest.skip(f'prepared benchmark fixture is absent: {source.name}')
    if count_pdb_frames(source) != 1:
        pytest.skip(
            f'prepared benchmark fixture is not single-frame: {source.name}'
        )
    destination = tmp_path / 'output' / 'benchmarks' / '1bna_DB1976' / run_name
    captured = []

    def stop_before_geometry(self):
        # Exercise settings, source-universe resolution, and group selection;
        # the real solver is deliberately not entered in this parser contract test.
        self.apply_settings()
        self._resolve_full_system_solvent_mode()
        self.create_groups()
        captured.append(self)

    monkeypatch.setattr(Command, '_process_system', stop_before_geometry)
    monkeypatch.setattr(
        sys,
        'argv',
        [
            'vorpy', str(source),
            '-g', 'c A and c B', '-g', 'c C', '-i',
            '-s', 'nt', network_type, 'and', 'sr', '0.2', 'and', 'mv', '5',
            '--boundary', 'none',
            '--solvent-mode', solvent_mode,
            '--profile', 'analysis',
            '--fail-if-exists',
            '--timing-report', str(destination / 'timing.json'),
            '-e', 'dir', str(destination),
        ],
    )

    Command().run()
    command = captured[0]
    system = command.sys

    assert command.execution_profile == 'analysis'
    assert command.boundary_config.mode.value == 'none'
    assert command.solvent_mode.value == solvent_mode
    assert command.fail_if_exists is True
    assert command.timing_report_path == (destination / 'timing.json').resolve()
    assert command.settings_dict['net_type'] == network_type
    assert command.settings_dict['surf_res'] == pytest.approx(0.2)
    assert command.settings_dict['max_vert'] == pytest.approx(5.0)
    assert command.calculation_plan.solver_settings == (
        ('nt', network_type), ('sr', '0.2'), ('mv', '5')
    )
    assert command.interface_mode is True
    assert command.exports == []  # -e dir configures the destination only.

    expected_system_dir = destination / source.stem
    assert Path(system.files['dir']) == expected_system_dir
    assert expected_system_dir.is_dir()
    assert command._source_generator_policy['active_generator_count'] == active_generators
    assert [len(group.ball_ndxs) for group in system.groups] == [758, 47]
    assert [group.name for group in system.groups] == ['c_A_and_c_B', 'c_C']

    # Generator policy and group membership remain in the plan/archive; they
    # are intentionally not repeated in the normal terminal summary.
    assert capsys.readouterr().out == ''

    # The actual interface child is created after the solve from the interface
    # name; it remains below the already-resolved structure directory.
    system.ifaces = [type('InterfaceStub', (), {
        'name': 'c_A_and_c_B_c_C_interface',
    })()]
    command._prepare_interface_output_directories()
    assert Path(system.ifaces[0].dir) == (
        expected_system_dir / 'c_A_and_c_B_c_C_interface'
    )
    assert not list(destination.glob('*.geometry.vpy'))
    assert not list(destination.rglob('visualization_manifest.json'))


def test_benchmark_cli_fail_if_exists_refuses_reuse(tmp_path, monkeypatch):
    if not WET.is_file():
        pytest.skip('prepared wet 1BNA fixture is absent')
    if count_pdb_frames(WET) != 1:
        pytest.skip('prepared wet 1BNA fixture is not single-frame')
    destination = tmp_path / 'output' / 'benchmarks' / '1bna_DB1976' / 'aw_wet'

    monkeypatch.setattr(Command, '_process_system', lambda self: None)
    argv = [
        'vorpy', str(WET), '-g', 'c A and c B', '-g', 'c C', '-i',
        '--boundary', 'none', '--solvent-mode', 'FULL_SOLVENT',
        '--profile', 'analysis', '--fail-if-exists',
        '-e', 'dir', str(destination),
    ]
    monkeypatch.setattr(sys, 'argv', argv)
    Command().run()

    with pytest.raises(SystemExit, match='existing output directory'):
        Command().run()
