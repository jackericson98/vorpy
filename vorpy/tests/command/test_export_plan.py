import os
import os
import subprocess
import sys
from types import SimpleNamespace
from pathlib import Path

import pytest

from vorpy.src.command.calculation_plan import (
    build_calculation_plan,
    build_export_plan,
)
from vorpy.src.command.vpy_cmnd import Command


def test_canonical_preset_aliases_and_modifier_examples():
    tiny = build_export_plan([('micro',)], scope='interface')
    assert tiny.preset == 'tiny'
    assert tiny.included_components == ('info', 'edges', 'interface_atoms')

    small = build_export_plan([('small', 'add', 'surfs')], scope='interface')
    assert small.preset == 'small'
    assert 'surfaces' in small.included_components
    assert small.operations == (('preset', 'small'), ('add', 'surfaces'))

    no_shell = build_export_plan([('large', 'no', 'shell')], scope='group')
    assert no_shell.preset == 'large'
    assert not {'shell', 'shell_edges', 'shell_vertices'} & set(
        no_shell.included_components
    )
    assert no_shell.excluded_components == (
        'shell', 'shell_edges', 'shell_vertices'
    )

    no_atom_surfaces = build_export_plan(
        [('all', 'no', 'atom_surfs')], scope='group'
    )
    assert not {
        'atom_surfaces', 'atom_edges', 'atom_vertices'
    } & set(no_atom_surfaces.included_components)

    tiny_logs = build_export_plan([('tiny', 'add', 'logs')], scope='interface')
    assert tiny_logs.included_components == (
        'info', 'edges', 'interface_atoms', 'logs'
    )


def test_none_is_a_plan_sentinel_and_later_components_can_be_added():
    none = build_export_plan([('none',)], scope='group')
    assert none.included_components == ()
    assert none.checkpoint_policy == 'shared'
    assert none.operations == (('none',),)

    later = build_export_plan([('none',), ('logs',)], scope='group')
    assert later.included_components == ('logs',)
    assert later.checkpoint_policy == 'shared'


def test_none_runtime_dispatch_keeps_the_mandatory_checkpoint(monkeypatch):
    from vorpy.src.command.command_export import argv_export
    from vorpy.src.output import output

    checkpoints = []
    system = SimpleNamespace(
        name='synthetic', balls=[object()], groups=[], ifaces=[], files={'dir': 'unused'},
        update_progress=lambda **kwargs: None,
    )
    monkeypatch.setattr(
        output, 'prepare_canonical_root',
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        output, '_export_network_archive',
        lambda **kwargs: checkpoints.append(kwargs),
    )

    argv_export(system, [['none']])
    assert system._export_plan.checkpoint_policy == 'shared'
    assert len(checkpoints) == 1


def test_export_suppression_controls_keep_the_mandatory_checkpoint():
    assert build_export_plan([('no', 'archive')]).checkpoint_policy == 'shared'
    assert build_export_plan([('no_checkpoint',)]).checkpoint_policy == 'shared'


@pytest.mark.parametrize('export_tokens', [('none',), ('no', 'archive')])
def test_public_cli_export_suppression_still_writes_checkpoint(tmp_path, export_tokens):
    source = tmp_path / 'checkpoint-source.pdb'
    source.write_text(
        'ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        'ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C\n'
        'ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C\n'
        'ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C\n'
        'END\n',
        encoding='ascii',
    )
    destination = tmp_path / 'output'
    root = Path(__file__).resolve().parents[3]
    command = [
        sys.executable, '-m', 'vorpy', '--fail-if-exists', str(source),
        '-s', 'nt', 'pow', 'and', 'mv', '10', '--boundary', 'shell',
        '-e', *export_tokens, '-e', 'dir', str(destination),
    ]
    completed = subprocess.run(
        command, cwd=root, capture_output=True, text=True, timeout=120,
        env={**os.environ, 'PYTHONPATH': str(root)},
    )
    assert completed.returncode == 0, completed.stdout + '\n' + completed.stderr
    checkpoints = list(destination.rglob('*_checkpoint.vpy'))
    assert len(checkpoints) == 1


def test_checkpoint_is_staged_before_optional_export_failure(tmp_path, monkeypatch):
    from vorpy.src.command.command_export import _run_component_plan
    from vorpy.src.output import output

    events = []
    checkpoint = tmp_path / 'synthetic_checkpoint.vpy'
    system = SimpleNamespace(
        name='synthetic',
        balls=[object()],
        groups=[],
        ifaces=[],
        files={'dir': str(tmp_path)},
        update_progress=lambda **kwargs: None,
        exports=lambda **kwargs: (_ for _ in ()).throw(
            RuntimeError('optional export failure')
        ),
    )
    plan = build_export_plan([('tiny',)], scope='system')

    def save_checkpoint(**kwargs):
        events.append('checkpoint')
        checkpoint.write_bytes(b'validated checkpoint')

    monkeypatch.setattr(output, '_export_network_archive', save_checkpoint)

    with pytest.raises(RuntimeError, match='optional export failure'):
        _run_component_plan(system, plan)

    assert events == ['checkpoint']
    assert checkpoint.read_bytes() == b'validated checkpoint'


def test_canonical_checkpoint_path_is_reused_for_explicit_archive(tmp_path, monkeypatch):
    from vorpy.src.output import output

    checkpoint = tmp_path / 'synthetic_checkpoint.vpy'
    checkpoint.write_bytes(b'validated checkpoint')
    system = SimpleNamespace(
        files={'dir': str(tmp_path)},
        _canonical_output_layout=True,
        _export_archive_written=False,
    )

    monkeypatch.setattr(
        'vorpy.src.output.checkpoint.write_checkpoint',
        lambda current: checkpoint,
    )

    output._export_network_archive(system)

    assert system._last_export_archive_path == str(checkpoint.resolve())
    assert system._export_archive_written is True


def test_interface_tiers_write_scoped_artifacts_and_one_checkpoint(tmp_path):
    from vorpy.src.command.command_export import _run_component_plan
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    expected = {
        'tiny': {'interface_info.txt', 'interface_edges.off', 'interface_atoms.pdb'},
        'small': {'interface_info.txt', 'interface_edges.off', 'interface_surfs.off', 'interface_atoms.pdb'},
        'medium': {
            'interface_info.txt', 'interface_edges.off', 'interface_surfs.off',
            'interface_verts.off', 'interface_atoms.pdb', 'interface_surfaces.csv',
            'interface_edges.csv', 'interface_vertices.csv', 'surr_atoms.pdb',
            'surr_residues.pdb',
        },
        'large': {
            'interface_info.txt', 'interface_edges.off', 'interface_surfs.off',
            'interface_verts.off', 'interface_atoms.pdb', 'interface_surfaces.csv',
            'interface_edges.csv', 'interface_vertices.csv', 'surr_atoms.pdb',
            'surr_residues.pdb',
        },
        'all': {
            'interface_info.txt', 'interface_edges.off', 'interface_surfs.off',
            'interface_verts.off', 'interface_atoms.pdb', 'interface_surfaces.csv',
            'interface_edges.csv', 'interface_vertices.csv', 'surr_atoms.pdb',
            'surr_residues.pdb',
        },
    }

    for tier, required in expected.items():
        root = tmp_path / tier
        root.mkdir()
        source, interface = _contract_case(root)
        system = interface.sys
        source.sys = system
        system.update_progress = lambda **kwargs: None
        system._export_plan = build_export_plan([(tier,)], scope='interface')

        _run_component_plan(system, system._export_plan)

        export_root = Path(system.files['dir'])
        interface_dir = next((export_root / 'interfaces').iterdir())
        assert required <= {path.name for path in interface_dir.iterdir()}
        assert (export_root / f'{system.name}_checkpoint.vpy').is_file()
        assert not list((export_root / 'system').glob('*.off'))


def test_runtime_dispatch_uses_components_without_rebuilding_networks(tmp_path, monkeypatch):
    from vorpy.src.command.command_export import argv_export
    from vorpy.src.output import output

    calls = []
    builds = []
    group = SimpleNamespace(
        name='A', net=SimpleNamespace(), settings={},
        build=lambda: builds.append('build'),
        exports=lambda **kwargs: calls.append(kwargs),
    )
    system = SimpleNamespace(
        name='synthetic', groups=[group], ifaces=[], files={'dir': str(tmp_path)},
        update_progress=lambda **kwargs: None,
    )
    monkeypatch.setattr(output, '_export_network_archive', lambda **kwargs: None)

    argv_export(system, [['large', 'no', 'shell']])

    assert builds == []
    assert any(call.get('surfs') for call in calls)
    assert not any(call.get('shell_surfs') or call.get('shell_edges')
                   or call.get('shell_verts') for call in calls)


def test_repeated_export_blocks_use_last_conflicting_operation():
    plan = build_export_plan(
        [
            ('large', 'no', 'shell'),
            ('small', 'add', 'shell_edges'),
            ('large', 'no', 'shell_edges'),
        ],
        scope='group',
    )

    assert plan.preset == 'large'
    assert 'shell_edges' not in plan.included_components
    assert 'shell_edges' in plan.excluded_components

    readded = build_export_plan(
        [('large', 'no', 'shell'), ('add', 'shell_edges')],
        scope='group',
    )
    assert 'shell_edges' in readded.included_components
    assert 'shell_edges' not in readded.excluded_components
    assert 'shell' in readded.excluded_components


def test_export_plan_keeps_analysis_profile_independent():
    plan = build_calculation_plan(
        ['input.pdb', '--profile', 'geometry', '-e', 'all', 'no', 'shell']
    )

    assert plan.analysis_profile == 'geometry'
    assert plan.export_plan.preset == 'all'
    assert 'shell' not in plan.export_plan.included_components

    simple = build_calculation_plan(['input.pdb', '--profile', 'simple'])
    assert simple.analysis_profile == 'geometry'


def test_scope_routing_matches_interface_contract():
    direct = build_calculation_plan(
        ['input.pdb', '-g', 'c', 'A', '-g', 'c', 'B', '-i']
    )
    assert direct.interface_only is True
    assert direct.independent_group_requests == ()
    assert direct.full_system_request is False
    assert direct.export_plan.scope == 'interface'

    groups = build_calculation_plan(
        ['input.pdb', '-g', 'c', 'A', '-g', 'c', 'B', '-i', '-c', 'groups']
    )
    assert groups.all_groups_request is True
    assert groups.interface_only is False
    assert groups.network_requests[0] == ('groups',)
    assert groups.export_plan.scope == 'mixed'

    system = build_calculation_plan(
        ['input.pdb', '-g', 'c', 'A', '-g', 'c', 'B', '-i', '-c', 'system']
    )
    assert system.full_system_request is True
    assert system.independent_group_requests == ()
    assert system.export_plan.scope == 'mixed'

    independent = build_calculation_plan(
        ['input.pdb', '-g', 'c', 'A', '-g', 'c', 'B']
    )
    assert len(independent.independent_group_requests) == 2
    assert independent.interface_requests == ()
    assert independent.export_plan.scope == 'group'


def test_reversed_interface_identity_does_not_add_a_second_network_request():
    plan = build_calculation_plan(
        [
            'input.pdb', '-i', 'g', 'c', 'A', 'g', 'c', 'B',
            '-i', 'g', 'c', 'B', 'g', 'c', 'A',
        ]
    )
    assert plan.interface_requests[0].operands == (('c', 'A'), ('c', 'B'))
    assert plan.interface_requests[1].operands == (('c', 'B'), ('c', 'A'))
    assert plan.network_requests == (
        ('interface', (('c', 'a'), ('c', 'b'))),
    )

    command = Command(sys=SimpleNamespace(groups=[object(), object()]))
    command.groups = {0: [['c', 'A']], 1: [['c', 'B']]}
    command.calculation_plan = plan
    assert len(command._planned_interface_pairs()) == 1


def test_empty_interface_request_remains_a_valid_single_target():
    plan = build_calculation_plan(['input.pdb', '-i', '-g', 'c', 'A', '-g', 'c', 'B'])

    assert plan.interface_requests[0].pair_count() == 1
    assert plan.network_requests == (('interface', ()),)
    assert plan.export_plan.scope == 'interface'


def test_public_cli_loads_canonical_checkpoint_without_source(tmp_path, monkeypatch):
    from vorpy.src.network import Network
    from vorpy.src.output.checkpoint import _load_checkpoint
    from vorpy.workbench.services.vorpy_backend import (
        VorPyBackend,
        VorPySolveSettings,
    )

    source = tmp_path / 'synthetic.pdb'
    source.write_text(
        'ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        'ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C\n'
        'ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C\n'
        'ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C\n'
        'END\n',
        encoding='ascii',
    )
    solved = VorPyBackend(
        VorPySolveSettings(
            network_type='aw',
            boundary_mode='shell',
            max_vertices=5,
            typed_group_selections={'synthetic': (0,)},
        )
    ).solve(source, lambda *args: None, lambda: False)
    system = solved.export_group.net.sys
    producer_root = tmp_path / 'producer'
    producer_root.mkdir()
    output = producer_root / 'exports'
    system.files['dir'] = str(output)
    system.update_progress = lambda **kwargs: None
    system._export_plan = build_export_plan([('tiny',)], scope='group')
    from vorpy.src.command.command_export import _run_component_plan
    _run_component_plan(system, system._export_plan)

    checkpoint = output / f'{system.name}_checkpoint.vpy'
    assert _load_checkpoint(checkpoint).payloads
    source.unlink()
    cache_output = tmp_path / 'cache-output'

    monkeypatch.setattr(
        Network, 'build',
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError('canonical checkpoint reload rebuilt geometry')
        ),
    )
    monkeypatch.setattr(
        sys, 'argv',
        [
            'vorpy', str(checkpoint), '-e', 'dir', str(cache_output),
            '-e', 'only', 'logs',
        ],
    )
    command = Command()
    command.run()

    assert command.sys._loaded_from_canonical_checkpoint is True
    assert command.sys._checkpoint_document.source.filename == 'synthetic.pdb'
    assert list(cache_output.rglob('*logs.csv'))


def test_interface_fail_if_exists_reserves_before_nested_creation(tmp_path):
    from vorpy.src.output.layout import canonical_interface_directory
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    fresh_root = tmp_path / 'fresh'
    fresh_root.mkdir()
    _, interface = _contract_case(fresh_root)
    system = interface.sys
    system.files['dir'] = str(tmp_path / 'fresh-output')
    system._export_plan = build_export_plan([('tiny',)], scope='interface')
    command = Command(sys=system)
    command.fail_if_exists = True

    command._prepare_interface_output_directories()
    command._prepare_interface_output_directories()
    assert canonical_interface_directory(system, interface).is_dir()

    existing_root = tmp_path / 'existing'
    existing_root.mkdir()
    _, existing_interface = _contract_case(existing_root)
    existing_system = existing_interface.sys
    existing_system.files['dir'] = str(tmp_path / 'existing-output')
    existing_system._export_plan = build_export_plan([('tiny',)], scope='interface')
    existing_dir = canonical_interface_directory(existing_system, existing_interface)
    existing_dir.mkdir(parents=True)
    guarded = Command(sys=existing_system)
    guarded.fail_if_exists = True

    with pytest.raises(SystemExit, match='existing interface output directory'):
        guarded._prepare_interface_output_directories()


def test_same_run_interface_reservation_is_reusable_but_not_inherited(
    tmp_path,
):
    from vorpy.src.output.layout import canonical_interface_directory
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    _, interface = _contract_case(tmp_path)
    system = interface.sys
    system.files['dir'] = str(tmp_path / 'reserved-output')
    system._export_plan = build_export_plan([('tiny',)], scope='interface')
    command = Command(sys=system)
    command.fail_if_exists = True
    pair = [(interface.group1, interface.group2)]

    command._reserve_interface_output_destinations(pair)
    reserved = canonical_interface_directory(system, interface).resolve()
    assert reserved.is_dir()

    # This is the original failure boundary: reservation happens before the
    # Interface exists, then preparation runs after it has been materialized.
    command._prepare_interface_output_directories()
    command._prepare_interface_output_directories()
    assert Path(interface.dir).resolve() == reserved
    assert Path(interface.dir).name == reserved.name

    # A new Command represents a new run, even when a caller reuses System.
    foreign = Command(sys=system)
    foreign.fail_if_exists = True
    with pytest.raises(SystemExit, match='existing interface output directory'):
        foreign._prepare_interface_output_directories()


def test_public_cli_replays_canonical_interface_without_rebuilding(
    tmp_path, monkeypatch,
):
    from vorpy.src.network import Network
    from vorpy.src.output.checkpoint import write_checkpoint
    from vorpy.workbench.services.vorpy_backend import (
        VorPyBackend,
        VorPySolveSettings,
    )
    source = tmp_path / 'synthetic.pdb'
    source.write_text(
        'ATOM      1  C1  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        'ATOM      2  C2  ALA A   2       2.000   0.000   0.000  1.00  0.00           C\n'
        'ATOM      3  C3  ALA A   3       0.000   2.000   0.000  1.00  0.00           C\n'
        'ATOM      4  C4  ALA A   4       0.000   0.000   2.000  1.00  0.00           C\n'
        'END\n',
        encoding='ascii',
    )
    producer = tmp_path / 'producer'
    solved = VorPyBackend(
        VorPySolveSettings(
            network_type='pow',
            boundary_mode='shell',
            max_vertices=30,
        )
    ).solve(source, lambda *args: None, lambda: False)
    system = solved.export_group.net.sys
    system.files['dir'] = str(producer)
    producer.mkdir()
    checkpoint = write_checkpoint(system, producer / 'synthetic_checkpoint.vpy')
    destination = tmp_path / 'replayed'

    monkeypatch.setattr(
        Network,
        'build',
        lambda *args, **kwargs: pytest.fail('canonical interface replay rebuilt geometry'),
    )
    monkeypatch.setattr(
        sys,
        'argv',
        [
            'vorpy', str(checkpoint),
            '-g', 'a', '0', '-g', 'a', '1', '-i',
            '-e', 'dir', str(destination), '-e', 'interface',
        ],
    )

    command = Command()
    command.run()

    assert command.sys._loaded_from_canonical_checkpoint is True
    assert len(command.sys.ifaces) == 1
    assert command.sys.ifaces[0].source_view is not None
    assert command.interface_results


def test_source_less_legacy_interface_recovery_fails_closed(tmp_path):
    from vorpy.src.io import ArchiveError

    system = SimpleNamespace(
        files={'base_file': None},
        _loaded_from_canonical_checkpoint=False,
    )
    command = Command(sys=system)
    network = SimpleNamespace(
        archive_metadata={},
        settings={'net_type': 'aw', 'solvent_mode': 'MOLECULE_ONLY'},
    )

    with pytest.raises(ArchiveError, match='source-less legacy archive'):
        command._run_archive_interface_recovery(network, tmp_path / 'legacy.vpy')
    assert system.archive_recovery_source_kind == 'legacy_source_less'


def test_subprocess_fresh_interface_fail_if_exists_reserves_once(tmp_path):
    """Exercise the live CLI reservation path with a bounded synthetic PDB."""
    source = tmp_path / 'synthetic_ab.pdb'
    atoms = (
        ('A', 1, 0.0, 0.0, 0.0), ('A', 2, 1.6, 0.0, 0.0),
        ('A', 3, 0.0, 1.6, 0.0), ('A', 4, 0.0, 0.0, 1.6),
        ('B', 5, 2.0, 1.0, 0.0), ('B', 6, 3.6, 1.0, 0.0),
        ('B', 7, 2.0, 2.6, 0.0), ('B', 8, 2.0, 1.0, 1.6),
    )
    source.write_text(
        ''.join(
            f'ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}'
            f'    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n'
            for serial, (chain, residue, x, y, z) in enumerate(atoms, 1)
        ) + 'END\n',
        encoding='ascii',
    )
    output = tmp_path / 'cli-output'
    root = Path(__file__).resolve().parents[3]
    command = [
        sys.executable, '-m', 'vorpy', '--fail-if-exists', str(source),
        '-g', 'c', 'A', '-g', 'c', 'B', '-i',
        '-s', 'nt', 'pow', 'and', 'mv', '30',
        '--profile', 'geometry', '--boundary', 'shell',
        '-e', 'dir', str(output),
    ]
    completed = subprocess.run(
        command, cwd=root, capture_output=True, text=True, timeout=120,
        env={**os.environ, 'PYTHONPATH': str(root)},
    )
    assert completed.returncode == 0, completed.stdout + '\n' + completed.stderr
    assert 'Solving POW interface network...' in completed.stdout
    assert 'Export preparation...' in completed.stdout
    assert 'Final summary...' in completed.stdout
    interface_infos = list(output.rglob('interface_info.txt'))
    assert len(interface_infos) == 1
    assert interface_infos[0].parent.name != 'interfaces'
