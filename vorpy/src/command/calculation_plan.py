"""Parsing and planning for the calculation-oriented VorPy CLI.

The original command runner is intentionally retained as the execution
backend.  This module only turns argv into an ordered, validated description
of requested work.
"""

from dataclasses import dataclass
from itertools import combinations
from typing import Optional


_CONJUNCTIONS = frozenset({'&', 'and', 'nd', 'also', '+', '&&'})
_TARGET_ALIASES = {
    's': 'system',
    'system': 'system',
    'g': 'groups',
    'group': 'groups',
    'groups': 'groups',
    'i': 'interfaces',
    'interface': 'interfaces',
    'interfaces': 'interfaces',
    'c': 'chain',
    'chain': 'chain',
}

# These are exact tokens.  In particular, ``-5`` is deliberately not a
# flag, so numeric settings retain their historical meaning.
_FLAG_ALIASES = {
    '-l': 'load', '--load': 'load',
    '-s': 'settings', '--settings': 'settings', '--set': 'settings',
    '-g': 'group', '--group': 'group',
    '-b': 'build', '--build': 'build',
    '-i': 'interface', '--interface': 'interface',
    '-c': 'calculate', '--calculate': 'calculate',
    '-e': 'export', '--export': 'export',
    '--profile': 'profile',
    '--save-network': 'save_network',
    '--load-network': 'load_network',
    '--all-frames': 'all_frames',
    '--parallel-frames': 'parallel_frames',
    '--verbose': 'verbose',
    '--diagnose-edges': 'diagnose_edges',
    '--boundary': 'boundary',
    '--probe-radius': 'probe_radius',
    '--shell-spacing': 'shell_spacing',
    '--shell-offset': 'shell_offset',
    '--apollonius': 'apollonius',
    '--alpha-complex': 'alpha_complex',
}


def _canonical_flag(token):
    """Return a canonical flag name, or ``None`` for an ordinary token."""
    return _FLAG_ALIASES.get(str(token).strip().lower())


@dataclass(frozen=True)
class CommandBlock:
    """One top-level flag and all tokens up to the next recognized flag."""

    flag: str
    args: tuple[str, ...] = ()


def parse_command_blocks(argv):
    """Parse ordered top-level command blocks without interpreting values.

    ``argv`` may include the input source as its first positional token.  The
    returned input source is separate so planning and execution can use the
    same representation.  Unknown tokens are retained in the current block;
    only exact recognized flags begin a new block.
    """
    tokens = [str(token) for token in argv]
    input_source = None
    index = 0
    if tokens and _canonical_flag(tokens[0]) is None and not tokens[0].startswith('-'):
        input_source = tokens[0]
        index = 1

    blocks = []
    while index < len(tokens):
        canonical = _canonical_flag(tokens[index])
        if canonical is None:
            if input_source is None:
                input_source = tokens[index]
                index += 1
                continue
            raise ValueError(f"Unexpected token outside a command block: {tokens[index]!r}")

        index += 1
        start = index
        while index < len(tokens) and _canonical_flag(tokens[index]) is None:
            index += 1
        blocks.append(CommandBlock(canonical, tuple(tokens[start:index])))

    return input_source, tuple(blocks)


def _without_conjunctions(tokens):
    return tuple(token for token in tokens if token.lower() not in _CONJUNCTIONS)


def _commands(tokens):
    """Split a legacy command block into ``and``-separated commands."""
    result = []
    current = []
    for token in tokens:
        if token.lower() in _CONJUNCTIONS:
            if current:
                result.append(tuple(current))
                current = []
        else:
            current.append(token)
    if current:
        result.append(tuple(current))
    return tuple(result)


@dataclass(frozen=True)
class GroupRequest:
    """A molecular group described by one or more legacy selectors."""

    commands: tuple[tuple[str, ...], ...]
    source: str = 'group'

    @property
    def selector(self):
        """Compatibility-friendly flat selector for single-command groups."""
        if len(self.commands) == 1:
            return self.commands[0]
        return tuple(token for command in self.commands for token in command)

    @property
    def key(self):
        return tuple(tuple(token.lower() for token in command) for command in self.commands)


@dataclass(frozen=True)
class InterfaceRequest:
    """An interface request before it is resolved to System Group objects."""

    operands: tuple[tuple[str, ...], ...] = ()

    @property
    def key(self):
        return tuple(tuple(token.lower() for token in operand) for operand in self.operands)

    def pair_count(self):
        if len(self.operands) < 2:
            return 1
        return len(tuple(combinations(self.operands, 2)))


@dataclass(frozen=True)
class CalculationRequest:
    """One ordered calculation target from a ``-c`` block."""

    target: str
    selector: tuple[str, ...] = ()


@dataclass(frozen=True)
class CalculationPlan:
    """Validated calculation plan consumed by the existing command runner."""

    input_source: Optional[str]
    independent_group_requests: tuple[GroupRequest, ...] = ()
    interface_requests: tuple[InterfaceRequest, ...] = ()
    full_system_request: bool = False
    all_groups_request: bool = False
    all_interfaces_request: bool = False
    calculation_requests: tuple[CalculationRequest, ...] = ()
    declared_group_requests: tuple[GroupRequest, ...] = ()
    solver_settings: tuple[tuple[str, ...], ...] = ()
    analysis_profile: Optional[str] = None
    export_settings: tuple[tuple[str, ...], ...] = ()
    archive_settings: tuple[CommandBlock, ...] = ()
    blocks: tuple[CommandBlock, ...] = ()

    @property
    def interface_only(self):
        """Whether all group definitions are interface operands only."""
        return bool(self.interface_requests) and not (
            self.all_groups_request or self.independent_group_requests
        )

    @property
    def network_requests(self):
        """Return stable, de-duplicated request keys in execution order."""
        requests = []
        seen = set()

        def add(value):
            if value not in seen:
                seen.add(value)
                requests.append(value)

        if self.full_system_request:
            add(('system',))
        if self.all_groups_request:
            add(('groups',))
        for request in self.independent_group_requests:
            add(('group', request.key))
        for request in self.interface_requests:
            operands = request.operands
            if len(operands) >= 3:
                for first, second in combinations(operands, 2):
                    add(('interface', (first, second)))
            else:
                add(('interface', request.key))
        return tuple(requests)


def _dedupe(items, key=lambda item: item):
    result = []
    seen = set()
    for item in items:
        marker = key(item)
        if marker not in seen:
            seen.add(marker)
            result.append(item)
    return tuple(result)


def _interface_operands(tokens):
    tokens = [token for token in tokens if token.lower() not in _CONJUNCTIONS]
    if not tokens:
        return ()

    # The explicit interface grammar is ``g <selector> g <selector>``.
    if any(token.lower() in {'g', 'group'} for token in tokens):
        operands = []
        current = []
        for token in tokens:
            if token.lower() in {'g', 'group'}:
                if current:
                    operands.append(tuple(current))
                current = []
            else:
                current.append(token)
        if current:
            operands.append(tuple(current))
    else:
        # Also accept repeated object selectors such as ``c 0 c 1``.
        object_markers = {
            'c', 'chain', 'r', 'res', 'residue', 'a', 'atom', 'n', 'index',
        }
        operands = []
        current = []
        for token in tokens:
            if current and token.lower() in object_markers:
                operands.append(tuple(current))
                current = []
            current.append(token)
        if current:
            operands.append(tuple(current))

    if any(not operand for operand in operands):
        raise ValueError('Interface operands must contain a group selector.')
    return tuple(operands)


def _target_requests(args):
    """Parse one or more comma-separated calculation target names."""
    tokens = []
    for token in args:
        if token.lower() in _CONJUNCTIONS:
            continue
        tokens.extend(part for part in token.split(',') if part)

    requests = []
    index = 0
    while index < len(tokens):
        target = tokens[index].lower()
        canonical = _TARGET_ALIASES.get(target)
        if canonical is None:
            raise ValueError(f'Unknown calculation target: {tokens[index]}')
        index += 1
        if canonical == 'chain':
            if index >= len(tokens):
                raise ValueError('Calculation target c/chain requires a selector.')
            selector = ('c', tokens[index])
            index += 1
            requests.append(CalculationRequest('group', selector))
        else:
            requests.append(CalculationRequest(canonical))
    return tuple(requests)


def build_calculation_plan(argv, *, input_source=None, analysis_profile=None):
    """Build and validate a :class:`CalculationPlan` from command-line argv."""
    parsed_input, blocks = parse_command_blocks(argv)
    if input_source is None:
        input_source = parsed_input

    declared = []
    interface_requests = []
    calculations = []
    solver_settings = []
    exports = []
    archive = []
    full_system = False
    all_groups = False
    all_interfaces = False

    for block in blocks:
        if block.flag == 'group':
            declared.append(GroupRequest(_commands(block.args), source='group'))
        elif block.flag == 'interface':
            interface_requests.append(InterfaceRequest(_interface_operands(block.args)))
        elif block.flag == 'calculate':
            for request in _target_requests(block.args):
                calculations.append(request)
                if request.target == 'system':
                    full_system = True
                elif request.target == 'groups':
                    all_groups = True
                elif request.target == 'interfaces':
                    all_interfaces = True
                elif request.target == 'group':
                    declared.append(GroupRequest((request.selector,), source='calculate'))
        elif block.flag == 'settings':
            solver_settings.extend(_commands(block.args))
        elif block.flag == 'export':
            exports.extend(_commands(block.args))
        elif block.flag in {'load_network', 'save_network'}:
            archive.append(block)
        elif block.flag == 'profile':
            if len(block.args) != 1:
                raise ValueError('--profile requires one profile name.')
            analysis_profile = block.args[0].lower()

    if analysis_profile is not None:
        analysis_profile = str(analysis_profile).lower()
        if analysis_profile not in {'geometry', 'interface', 'analysis', 'full'}:
            raise ValueError(
                '--profile requires geometry, interface, analysis, or full'
            )

    declared = list(_dedupe(declared, key=lambda item: item.key))

    # ``-i`` and ``-c i`` are both interface requests.  An explicit request
    # wins over the implicit default but identical requests are de-duplicated.
    if all_interfaces:
        interface_requests.append(InterfaceRequest())
    interface_requests = list(_dedupe(interface_requests, key=lambda item: item.key))

    # Top-level ``-g`` always requests an independent group solve.  Group
    # selectors introduced by ``-i`` are interface definitions only and are
    # appended below with source ``interface``.  ``-c g`` covers declared
    # groups as one de-duplicated target.
    independent = [] if all_groups else list(declared)
    independent = list(_dedupe(independent, key=lambda item: item.key))

    # An interface selector can itself define groups.  Keep the declarations
    # in the plan so the command adapter can construct them before resolving
    # interface pairs.
    for interface in interface_requests:
        for operand in interface.operands:
            request = GroupRequest((operand,), source='interface')
            if request.key not in {group.key for group in declared}:
                declared.append(request)

    return CalculationPlan(
        input_source=input_source,
        independent_group_requests=tuple(independent),
        interface_requests=tuple(interface_requests),
        full_system_request=full_system,
        all_groups_request=all_groups,
        all_interfaces_request=all_interfaces,
        calculation_requests=_dedupe(calculations, key=lambda item: (item.target, item.selector)),
        declared_group_requests=tuple(declared),
        solver_settings=tuple(solver_settings),
        analysis_profile=analysis_profile,
        export_settings=tuple(exports),
        archive_settings=tuple(archive),
        blocks=blocks,
    )


__all__ = [
    'CalculationPlan', 'CalculationRequest', 'CommandBlock', 'GroupRequest',
    'InterfaceRequest', 'build_calculation_plan', 'parse_command_blocks',
]
