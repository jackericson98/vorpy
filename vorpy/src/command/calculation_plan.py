"""Parsing and planning for the calculation-oriented VorPy CLI.

The original command runner is intentionally retained as the execution
backend.  This module only turns argv into an ordered, validated description
of requested work.
"""

from dataclasses import dataclass
from itertools import combinations
from typing import Optional


_CONJUNCTIONS = frozenset({'&', 'and', 'nd', 'also', '+', '&&'})

# Export planning is deliberately kept in the command layer.  The output
# writers consume these canonical names but do not decide CLI precedence.
EXPORT_PRESETS = ('tiny', 'small', 'medium', 'large', 'all')
_PRESET_ALIASES = {
    'micro': 'tiny',
    'i': 'tiny',
    '0': 'tiny',
    'smallest': 'tiny',
    's': 'small',
    '1': 'small',
    'med': 'medium',
    '2': 'medium',
    'l': 'large',
    '3': 'large',
    'a': 'all',
    'everything': 'all',
}

# The values are cumulative and intentionally use the vocabulary frozen by
# NEXUSA, rather than writer keyword names such as ``surfs``.
_PRESET_COMPONENTS = {
    'system': {
        'tiny': ('info', 'edges', 'logs'),
        'small': ('info', 'edges', 'logs', 'surfaces'),
        'medium': (
            'info', 'edges', 'logs', 'surfaces', 'vertices',
            'source_structure', 'pymol_atoms',
        ),
        'large': (
            'info', 'edges', 'logs', 'surfaces', 'vertices',
            'source_structure', 'pymol_atoms', 'shell', 'shell_edges',
            'shell_vertices', 'atom_surfaces', 'atom_edges', 'atom_vertices',
        ),
        'all': (
            'info', 'edges', 'logs', 'surfaces', 'vertices',
            'source_structure', 'pymol_atoms', 'shell', 'shell_edges',
            'shell_vertices', 'atom_surfaces', 'atom_edges', 'atom_vertices',
            'separate_surfaces', 'separate_edges', 'separate_vertices',
            'surrounding_atoms', 'surrounding_residues', 'external_atoms',
        ),
    },
    'group': {
        'tiny': ('info', 'edges', 'logs'),
        'small': ('info', 'edges', 'logs', 'surfaces', 'group_atoms'),
        'medium': (
            'info', 'edges', 'logs', 'surfaces', 'group_atoms', 'vertices',
            'surrounding_atoms', 'surrounding_residues',
        ),
        'large': (
            'info', 'edges', 'logs', 'surfaces', 'group_atoms', 'vertices',
            'surrounding_atoms', 'surrounding_residues', 'shell',
            'shell_edges', 'shell_vertices', 'atom_surfaces', 'atom_edges',
            'atom_vertices',
        ),
        'all': (
            'info', 'edges', 'logs', 'surfaces', 'group_atoms', 'vertices',
            'surrounding_atoms', 'surrounding_residues', 'shell',
            'shell_edges', 'shell_vertices', 'atom_surfaces', 'atom_edges',
            'atom_vertices', 'separate_surfaces', 'separate_edges',
            'separate_vertices', 'external_atoms',
        ),
    },
    'interface': {
        'tiny': ('info', 'edges', 'interface_atoms'),
        'small': ('info', 'edges', 'interface_atoms', 'surfaces'),
        'medium': (
            'info', 'edges', 'interface_atoms', 'surfaces', 'vertices',
            'contact_tables', 'surrounding_atoms', 'surrounding_residues',
        ),
        'large': (
            'info', 'edges', 'interface_atoms', 'surfaces', 'vertices',
            'contact_tables', 'surrounding_atoms', 'surrounding_residues',
            'curvature', 'topology', 'shell', 'shell_edges', 'shell_vertices',
        ),
        'all': (
            'info', 'edges', 'interface_atoms', 'surfaces', 'vertices',
            'contact_tables', 'surrounding_atoms', 'surrounding_residues',
            'curvature', 'topology', 'shell', 'shell_edges', 'shell_vertices',
        ),
    },
}

_EXPORT_COMPONENT_ALIASES = {
    'group_info': 'info',
    'interface_info': 'info',
    'system_info': 'info',
    'edge': 'edges',
    'interface_edges': 'edges',
    'surfs': 'surfaces',
    'surface': 'surfaces',
    'interface_surfs': 'surfaces',
    'verts': 'vertices',
    'vertex': 'vertices',
    'interface_verts': 'vertices',
    'shell_surfs': 'shell',
    'shell_surfaces': 'shell',
    'shell_edge': 'shell_edges',
    'shell_verts': 'shell_vertices',
    'shell_vertex': 'shell_vertices',
    'atoms': 'group_atoms',
    'contact_atoms': 'interface_atoms',
    'surr_atoms': 'surrounding_atoms',
    'surr_resids': 'surrounding_residues',
    'surr_residues': 'surrounding_residues',
    'surrounding_resids': 'surrounding_residues',
    'ext_atoms': 'external_atoms',
    'sep_surfs': 'separate_surfaces',
    'sep_edges': 'separate_edges',
    'sep_verts': 'separate_vertices',
    'atom_surfs': 'atom_surfaces',
    'atom_verts': 'atom_vertices',
    'scheme_logs': 'logs',
    'classic_logs': 'logs',
    'lgs': 'logs',
    'canonical_results': 'results',
    'typed_results': 'results',
    'provenance_sidecar': 'provenance',
    'legacy_logs': 'legacy_results',
    'logs_csv': 'legacy_results',
    'contacts': 'contact_tables',
    'contact_table': 'contact_tables',
    'curv': 'curvature',
    'topo': 'topology',
    'solvent': 'water',
    'buried_water': 'water',
    'apollonius': 'dual',
    'alpha_dual': 'dual',
    'alpha_complex': 'alpha',
    'viz': 'visualization',
    'pymol': 'visualization',
    'pdb': 'source_structure',
    'set_atoms': 'pymol_atoms',
    'archive': 'checkpoint',
    'network_archive': 'checkpoint',
    'save_network': 'checkpoint',
    'input': 'source',
    'timing_report': 'timing',
    'profile': 'timing',
}

_EXPORT_FAMILIES = {
    'shell': ('shell', 'shell_edges', 'shell_vertices'),
    'atom_surfaces': ('atom_surfaces', 'atom_edges', 'atom_vertices'),
}
_LEGACY_EXPORT_COMMANDS = frozenset({
    'interface', 'default', 'nonpolar', 'dual', 'apollonius', 'visualize',
    'water', 'waters', 'mol', 'cif', 'xyz', 'txt',
})
_NO_EXPORT_COMMANDS = frozenset({'none', 'skip'})
_CANONICAL_EXPORT_COMPONENTS = frozenset(
    component
    for components in _PRESET_COMPONENTS.values()
    for preset_components in components.values()
    for component in preset_components
) | frozenset({
    'results', 'provenance', 'legacy_results', 'water', 'dual', 'alpha',
    'visualization', 'checkpoint', 'source', 'timing', 'mol', 'cif', 'xyz',
    'txt',
})


def canonical_export_component(token):
    """Normalize one CLI export component without applying precedence."""
    value = str(token).strip().lower().replace('-', '_')
    return _EXPORT_COMPONENT_ALIASES.get(value, value)


def _component_family(token):
    canonical = canonical_export_component(token)
    return _EXPORT_FAMILIES.get(canonical, (canonical,))


def _scope_components(scope, preset):
    if scope == 'mixed':
        components = []
        for scope_name in ('system', 'group', 'interface'):
            components.extend(_PRESET_COMPONENTS[scope_name][preset])
        return tuple(dict.fromkeys(components))
    try:
        return _PRESET_COMPONENTS[scope][preset]
    except KeyError as error:
        raise ValueError(
            f"Unknown export scope {scope!r}; expected system, group, interface, or mixed."
        ) from error


def _export_operation_tokens(exports):
    """Yield non-conjunction command tokens from raw or parsed ``-e`` blocks."""
    if exports is None:
        return
    if isinstance(exports, str):
        exports = (tuple(exports.split()),)
    elif isinstance(exports, (list, tuple)) and all(
        isinstance(value, str) for value in exports
    ):
        tokens = list(exports)
        if '-e' in tokens or '--export' in tokens:
            blocks = []
            current = []
            for token in tokens:
                if token in {'-e', '--export'}:
                    if current:
                        blocks.append(tuple(current))
                        current = []
                else:
                    current.append(token)
            if current:
                blocks.append(tuple(current))
            exports = tuple(blocks)
        else:
            exports = (tuple(tokens),)
    for block in exports:
        if isinstance(block, str):
            block = (block,)
        tokens = tuple(str(token) for token in block)
        current = []
        for token in tokens:
            if token.lower() in _CONJUNCTIONS:
                if current:
                    yield tuple(current)
                    current = []
            else:
                current.append(token)
        if current:
            yield tuple(current)


@dataclass(frozen=True)
class ExportPlan:
    """The normalized CLI export contract handed to output writers."""

    scope: str
    preset: str
    included_components: tuple[str, ...]
    excluded_components: tuple[str, ...]
    checkpoint_policy: str = 'shared'
    operations: tuple[tuple[str, ...], ...] = ()

    def as_dict(self):
        """Return a serialization-friendly view for summaries and adapters."""
        return {
            'scope': self.scope,
            'preset': self.preset,
            'included_components': self.included_components,
            'excluded_components': self.excluded_components,
            'checkpoint_policy': self.checkpoint_policy,
            'operations': self.operations,
        }


def build_export_plan(exports=(), *, scope='system'):
    """Normalize repeated ``-e`` blocks using the frozen CLI precedence rules."""
    scope = str(scope).strip().lower()
    if scope not in {'system', 'group', 'interface', 'mixed'}:
        raise ValueError(
            f"Unknown export scope {scope!r}; expected system, group, interface, or mixed."
        )

    preset = 'small'
    included = list(_scope_components(scope, preset))
    excluded = []
    operations = []
    checkpoint_policy = 'shared'
    directives = {'add', 'no', 'remove', 'only'}
    controls = {
        'checkpoint', 'no_checkpoint', 'no_archive', 'no_network',
        'no-archive', 'no-network', 'dir', 'directory', 'format',
        'file_format', 'file_type', 'mesh_format',
    }

    def add_values(values):
        for value in values:
            if value not in included:
                included.append(value)
            while value in excluded:
                excluded.remove(value)

    def exclude_values(values):
        for value in values:
            while value in included:
                included.remove(value)
            if value not in excluded:
                excluded.append(value)

    for command in _export_operation_tokens(exports):
        tokens = list(command)
        while tokens:
            head = tokens.pop(0).strip().lower().replace('-', '_')
            if head in _NO_EXPORT_COMMANDS or (head == 'no' and not tokens):
                # ``none`` suppresses optional files, not the mandatory
                # geometry checkpoint.
                included = []
                excluded = []
                operations.append(('none',))
                continue
            if head in controls:
                if head in {'no_archive', 'no_network', 'no_checkpoint'}:
                    # Checkpoints are durability records, not optional
                    # exports.  Keep the policy shared even for legacy
                    # suppression spellings.
                    operations.append((head,))
                elif head == 'checkpoint':
                    checkpoint_policy = 'shared'
                    operations.append(('checkpoint',))
                else:
                    operations.append(('control', head, *tokens))
                break

            if head in _PRESET_ALIASES or head in EXPORT_PRESETS:
                preset = _PRESET_ALIASES.get(head, head)
                preset_components = _scope_components(scope, preset)
                included = list(preset_components)
                # A later preset re-adds only components it explicitly
                # contains; unrelated prior removals remain explicit until
                # overridden.
                for value in preset_components:
                    while value in excluded:
                        excluded.remove(value)
                operations.append(('preset', preset))
                continue

            if head in _LEGACY_EXPORT_COMMANDS:
                operations.append(('legacy', head, *tokens))
                break

            if head in directives:
                operation = head
                values_tokens = []
                while tokens:
                    candidate = tokens[0].strip().lower().replace('-', '_')
                    if (
                        candidate in directives
                        or candidate in controls
                        or candidate in _PRESET_ALIASES
                        or candidate in EXPORT_PRESETS
                    ):
                        break
                    values_tokens.append(tokens.pop(0))
                if not values_tokens:
                    raise ValueError(
                        f"Export operation {head!r} requires a component."
                    )
                canonical = tuple(
                    canonical_export_component(token) for token in values_tokens
                )
                unknown = tuple(
                    token for token in canonical
                    if token not in _CANONICAL_EXPORT_COMPONENTS
                )
                if unknown:
                    raise ValueError(f"Unknown export component: {unknown[0]}")
                operations.append((operation, *canonical))
                values = []
                for token in canonical:
                    values.extend(_component_family(token))
                if operation == 'add':
                    add_values(values)
                elif operation in {'no', 'remove'}:
                    exclude_values(values)
                else:  # only
                    included = []
                    excluded = []
                    checkpoint_policy = 'explicit-only'
                    add_values(values)
                    if 'checkpoint' in values:
                        checkpoint_policy = 'shared'
                if 'checkpoint' in values and operation == 'add':
                    checkpoint_policy = 'shared'
                # A checkpoint cannot be removed from a successful solve.
                # ``-e no checkpoint`` remains parseable for compatibility,
                # but the mandatory shared checkpoint stays enabled.
                continue

            component = canonical_export_component(head)
            if component not in _CANONICAL_EXPORT_COMPONENTS:
                raise ValueError(f"Unknown export component: {head}")
            values = _component_family(component)
            operations.append(('add', component))
            add_values(values)
            if component == 'checkpoint':
                checkpoint_policy = 'shared'

    return ExportPlan(
        scope=scope,
        preset=preset,
        included_components=tuple(included),
        excluded_components=tuple(excluded),
        checkpoint_policy=checkpoint_policy,
        operations=tuple(operations),
    )


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
    '--timing-report': 'timing_report',
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

_ANALYSIS_PROFILE_ALIASES = {'simple': 'geometry'}


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


def expand_group_tokens(tokens):
    """Normalize one shell-quoted selector into legacy selector tokens.

    Windows shells commonly pass ``"c A and c B"`` as one argv value.  The
    legacy group interpreter consumes the equivalent unquoted token sequence.
    Group selectors do not support meaningful whitespace inside one token, so
    splitting only this command-block value is safe and keeps other CLI values
    unchanged.
    """
    expanded = []
    for token in tokens:
        expanded.extend(str(token).split())
    return tuple(expanded)


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

    @property
    def identity_key(self):
        """Return the unordered physical interface identity key."""
        return tuple(sorted(self.key))

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
    export_plan: Optional[ExportPlan] = None
    archive_settings: tuple[CommandBlock, ...] = ()
    blocks: tuple[CommandBlock, ...] = ()

    @property
    def interface_only(self):
        """Whether the plan requests interfaces without another network target."""
        return bool(self.interface_requests) and not (
            self.full_system_request
            or self.all_groups_request
            or self.independent_group_requests
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
            operands = request.key
            if len(operands) >= 3:
                for first, second in combinations(operands, 2):
                    pair = tuple(sorted((first, second)))
                    add(('interface', pair))
            else:
                add(('interface', request.identity_key))
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
            declared.append(
                GroupRequest(
                    _commands(expand_group_tokens(block.args)),
                    source='group',
                )
            )
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
        elif block.flag == 'timing_report':
            if len(block.args) != 1:
                raise ValueError('--timing-report requires one JSON output path.')
        elif block.flag == 'profile':
            if len(block.args) != 1:
                raise ValueError('--profile requires one profile name.')
            analysis_profile = _ANALYSIS_PROFILE_ALIASES.get(
                block.args[0].lower(), block.args[0].lower()
            )

    if analysis_profile is not None:
        analysis_profile = _ANALYSIS_PROFILE_ALIASES.get(
            str(analysis_profile).lower(), str(analysis_profile).lower()
        )
        if analysis_profile not in {'geometry', 'interface', 'analysis', 'full'}:
            raise ValueError(
                '--profile requires geometry, simple, interface, analysis, or full'
            )

    declared = list(_dedupe(declared, key=lambda item: item.key))

    # ``-i`` and ``-c i`` are both interface requests.  An explicit request
    # wins over the implicit default but identical requests are de-duplicated.
    if all_interfaces:
        interface_requests.append(InterfaceRequest())
    interface_requests = list(_dedupe(interface_requests, key=lambda item: item.key))

    # With ``-i``, top-level ``-g`` declarations define the ordered interface
    # sides.  They are selection containers, not independent group solves.
    # An explicit ``-c`` group target remains an independent calculation.
    independent = [] if all_groups else list(declared)
    if interface_requests:
        independent = [
            request for request in independent
            if request.source != 'group'
        ]
    independent = list(_dedupe(independent, key=lambda item: item.key))

    # An interface selector can itself define groups.  Keep the declarations
    # in the plan so the command adapter can construct them before resolving
    # interface pairs.
    for interface in interface_requests:
        for operand in interface.operands:
            request = GroupRequest((operand,), source='interface')
            if request.key not in {group.key for group in declared}:
                declared.append(request)

    # No target flag means the historical default: solve the complete input
    # system.  Keep that default explicit in the plan so diagnostics and
    # callers cannot confuse an ordinary full-system solve with "no work".
    if not (
        full_system
        or all_groups
        or all_interfaces
        or independent
        or interface_requests
    ):
        full_system = True

    if interface_requests and (
        full_system or all_groups or independent
    ):
        export_scope = 'mixed'
    elif interface_requests:
        export_scope = 'interface'
    elif all_groups or independent:
        export_scope = 'group'
    else:
        export_scope = 'system'

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
        export_plan=build_export_plan(exports, scope=export_scope),
        archive_settings=tuple(archive),
        blocks=blocks,
    )


__all__ = [
    'CalculationPlan', 'CalculationRequest', 'CommandBlock', 'ExportPlan',
    'EXPORT_PRESETS', 'GroupRequest', 'InterfaceRequest',
    'build_calculation_plan', 'build_export_plan', 'canonical_export_component',
    'expand_group_tokens', 'parse_command_blocks',
]
