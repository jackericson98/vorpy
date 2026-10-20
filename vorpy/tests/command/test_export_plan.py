from types import SimpleNamespace
import sys

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
    assert none.checkpoint_policy == 'none'
    assert none.operations == (('none',),)

    later = build_export_plan([('none',), ('logs',)], scope='group')
    assert later.included_components == ('logs',)
    assert later.checkpoint_policy == 'none'


def test_none_runtime_dispatch_does_not_prepare_or_write_exports(monkeypatch):
    from vorpy.src.command.command_export import argv_export
    from vorpy.src.output import output

    system = SimpleNamespace(
        name='synthetic', groups=[], ifaces=[], files={'dir': 'unused'},
    )
    monkeypatch.setattr(
        output, 'prepare_canonical_root',
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError('none must not prepare export directories')
        ),
    )
    monkeypatch.setattr(
        output, '_export_network_archive',
        lambda **kwargs: (_ for _ in ()).throw(
            AssertionError('none must disable checkpoint export')
        ),
    )

    argv_export(system, [['none']])
    assert system._export_plan.checkpoint_policy == 'none'


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
