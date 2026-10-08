import pytest
from types import SimpleNamespace

from vorpy.src.command.calculation_plan import (
    build_calculation_plan,
    parse_command_blocks,
)
from vorpy.src.command.vpy_cmnd import Command


def test_blocks_are_ordered_and_negative_values_are_not_flags():
    source, blocks = parse_command_blocks(
        ['p53tet', '-s', 'mv', '-5', '-g', 'c', '0', '-i']
    )
    assert source == 'p53tet'
    assert [(block.flag, block.args) for block in blocks] == [
        ('settings', ('mv', '-5')),
        ('group', ('c', '0')),
        ('interface', ()),
    ]


def test_group_declarations_are_independent_without_interface():
    plan = build_calculation_plan(['p53tet', '-g', 'c', '0'])
    assert [request.selector for request in plan.independent_group_requests] == [
        ('c', '0')
    ]
    assert plan.interface_requests == ()


def test_interface_keeps_top_level_group_solves():
    plan = build_calculation_plan(
        ['p53tet', '-g', 'c', '0', '-g', 'c', '1', '-i']
    )
    assert [request.selector for request in plan.independent_group_requests] == [
        ('c', '0'), ('c', '1')
    ]
    assert plan.interface_requests[0].operands == ()
    assert len(plan.declared_group_requests) == 2


@pytest.mark.parametrize(
    ('argv', 'expected'),
    [
        (['p53tet', '-i'], ()),
        (['p53tet', '-i', 'g', 'c', '0'], (('c', '0'),)),
        (['p53tet', '-i', 'g', 'c', '0', 'g', 'c', '1'],
         (('c', '0'), ('c', '1'))),
        (['p53tet', '-i', 'g', 'c', '0', 'g', 'c', '1', 'g', 'c', '2'],
         (('c', '0'), ('c', '1'), ('c', '2'))),
    ],
)
def test_interface_operand_cardinality(argv, expected):
    request = build_calculation_plan(argv).interface_requests[0]
    assert request.operands == expected
    assert request.pair_count() == (1 if len(expected) < 2 else len(expected) * (len(expected) - 1) // 2)


def test_calculate_targets_aliases_and_repetition_are_normalized():
    plan = build_calculation_plan(
        ['p53tet', '--calculate', 'i,s', '-c', 'g', '-c', 'chain', '0',
         '-c', 'c', '0']
    )
    assert plan.full_system_request is True
    assert plan.all_groups_request is True
    assert plan.all_interfaces_request is True
    assert len(plan.interface_requests) == 1
    assert plan.independent_group_requests == ()
    assert len(plan.calculation_requests) == 4


def test_duplicate_network_requests_are_deduplicated():
    plan = build_calculation_plan(
        ['p53tet', '-c', 'c', '0', '-c', 'c', '0', '-i',
         'g', 'c', '0', 'g', 'c', '1', '-i',
         'g', 'c', '0', 'g', 'c', '1']
    )
    assert len(plan.independent_group_requests) == 1
    assert len(plan.interface_requests) == 1
    assert len(plan.network_requests) == 2


def test_all_group_target_covers_declared_groups_once():
    plan = build_calculation_plan(
        ['p53tet', '-g', 'c', '0', '-g', 'c', '1', '-c', 'g']
    )
    assert plan.all_groups_request is True
    assert plan.independent_group_requests == ()
    assert plan.network_requests == (('groups',),)


def test_system_and_interface_requests_can_coexist():
    plan = build_calculation_plan(['p53tet', '-c', 's', '-c', 'i'])
    assert plan.full_system_request is True
    assert len(plan.interface_requests) == 1
    assert plan.network_requests[0] == ('system',)


def test_invalid_calculation_syntax_is_rejected_before_execution():
    with pytest.raises(ValueError, match='requires a selector'):
        build_calculation_plan(['p53tet', '-c', 'c'])
    with pytest.raises(ValueError, match='Unknown calculation target'):
        build_calculation_plan(['p53tet', '-c', 'not-a-target'])


def test_command_adapter_accepts_long_aliases_and_negative_settings():
    command = Command()
    command._arguments = [
        'p53tet', '--settings', 'mv', '-5', '--group', 'c', '0',
        '--interface',
    ]
    command.parse_commands()
    assert command.settings_cmnds == [['mv', '-5']]
    assert command.groups == {0: [['c', '0']]}
    assert command.interface_mode is True
    assert command.calculation_plan.interface_only is False


def test_plan_keeps_archive_settings():
    plan = build_calculation_plan(['p53tet', '--save-network', 'result.vpy'])
    assert [(block.flag, block.args) for block in plan.archive_settings] == [
        ('save_network', ('result.vpy',))
    ]


def test_command_adapter_materializes_specific_chain_calculation():
    command = Command()
    command._arguments = ['p53tet', '-c', 'c', '0']
    command.parse_commands()
    assert command.groups == {0: [['c', '0']]}
    assert command.interface_mode is False
    assert command._build_group_networks is True


def test_interface_plan_resolves_selected_groups_and_pairwise_requests():
    groups = [SimpleNamespace(name=f'g{index}') for index in range(3)]
    command = Command(sys=SimpleNamespace(groups=groups))
    command.groups = {
        0: [['c', '0']],
        1: [['c', '1']],
        2: [['c', '2']],
    }
    command.calculation_plan = build_calculation_plan(
        ['p53tet', '-i', 'g', 'c', '0', 'g', 'c', '1', 'g', 'c', '2']
    )
    assert command._planned_interface_pairs() == [
        (groups[0], groups[1]),
        (groups[0], groups[2]),
        (groups[1], groups[2]),
    ]


@pytest.mark.parametrize(
    ('argv', 'independent', 'full_system', 'interface_only'),
    [
        (
            ['p53tet', '-g', 'c', '0', '-g', 'c', '1', '-i'],
            [('c', '0'), ('c', '1')], False, False,
        ),
        (
            ['p53tet', '-i', 'g', 'c', '0', 'g', 'c', '1'],
            [], False, True,
        ),
        (
            ['p53tet', '-g', 'c', '0', '-g', 'c', '1', '-i', '-c', 's'],
            [('c', '0'), ('c', '1')], True, False,
        ),
        (
            ['p53tet', '-i', 'g', 'c', '0', 'g', 'c', '1', '-c', 's'],
            [], True, True,
        ),
    ],
)
def test_final_calculation_contract_resolves_execution_targets(
    argv, independent, full_system, interface_only,
):
    command = Command()
    command._arguments = argv
    command.parse_commands()
    plan = command.calculation_plan

    assert [request.selector for request in plan.independent_group_requests] == independent
    assert plan.full_system_request is full_system
    assert bool(plan.interface_requests) is True
    assert command._build_group_networks is bool(independent)
    assert command._full_system_separate is full_system
    assert plan.interface_only is interface_only
    assert len(plan.network_requests) == (
        len(independent) + 1 + (1 if full_system else 0)
    )

    # These are the concrete group definitions handed to the existing
    # executor; operand-only groups remain distinct from top-level -g groups.
    assert list(command.groups.values()) == [[['c', '0']], [['c', '1']]]


def test_group_backend_build_flag_matches_independent_targets(monkeypatch):
    calls = []
    monkeypatch.setattr(
        'vorpy.src.command.vpy_cmnd.ggroup',
        lambda sys, groups, settings, make_net: calls.append(make_net),
    )

    command = Command(sys=SimpleNamespace(groups=[]))
    command.groups = {0: [['c', '0']], 1: [['c', '1']]}
    command.interface_mode = True
    command._build_group_networks = False
    command.create_groups()
    command._build_group_networks = True
    command.create_groups()

    assert calls == [False, True]


def test_full_system_target_is_built_once_separately(monkeypatch):
    builds = []

    class FakeFullGroup:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def build(self, **kwargs):
            builds.append((self.kwargs['atoms'], kwargs))

    monkeypatch.setattr('vorpy.src.group.Group', FakeFullGroup)
    command = Command(
        sys=SimpleNamespace(name='synthetic', balls=[0, 1, 2], groups=[]),
    )
    command.settings_dict = {}
    command._full_system_separate = True
    command._build_full_system_target()
    command._build_full_system_target()

    assert builds == [([0, 1, 2], {'calculate_curvature': True})]
