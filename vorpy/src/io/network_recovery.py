"""Read-only recovery adapter for solved geometry checkpoints.

The adapter validates source identity before handing an already-loaded Network
to the existing NetworkView selection path.  It never constructs geometry.
"""

from dataclasses import dataclass
import hashlib
import json
from typing import Any, Mapping

from .network_archive import ArchiveError, validate_network_source


@dataclass(frozen=True)
class SourceGeneratorIdentity:
    """Exact source-table identity for one retained molecular generator."""

    system_index: int
    stable_id: str
    name: str
    residue_name: str
    residue_sequence: Any
    chain_id: str
    element: str
    coordinates: tuple[float, float, float]
    radius: float


@dataclass(frozen=True)
class SourceSystemIdentity:
    """Reconstructed source universe and checkpoint scientific status."""

    source_path: str
    source_sha256: str
    frame_index: Any
    generator_count: int
    generator_system_indices: tuple[int, ...]
    generator_stable_ids: tuple[str, ...]
    generators: tuple[SourceGeneratorIdentity, ...]
    scheme: str
    solvent_mode: str | None
    boundary_mode: str | None
    complete_cell_count: int | None
    complete_cell_evidence: bool
    geometry_status: str | None
    analysis_status: str | None
    solve_universe_id: str
    radii_digest: str
    generator_policy: Mapping[str, Any] | None
    source_network_id: str | None
    settings: Mapping[str, Any]


@dataclass(frozen=True)
class RecoveredSourceNetwork:
    """Validated solved Network plus its source identity.

    ``make_interface_view`` delegates only to ``Network.make_view``.  The
    existing view implementation snapshots local topology and records partial
    completeness; it does not call ``Network.build``.
    """

    network: Any
    source_identity: SourceSystemIdentity
    validation: Mapping[str, Any]

    def make_interface_view(
        self,
        group_a_system_ids,
        group_b_system_ids,
        *,
        source_network_id: str | None = None,
    ):
        """Create an ordered A -> B view from the reloaded physical network."""

        return self.network.make_view(
            iface_grps=(tuple(group_a_system_ids), tuple(group_b_system_ids)),
            source_network_id=source_network_id,
        )


def _scalar(value):
    if value is None:
        return None
    try:
        if value != value:  # NaN, without importing pandas into the adapter.
            return None
    except (TypeError, ValueError):
        pass
    return value.item() if hasattr(value, 'item') else value


def _normalized_setting(value):
    if hasattr(value, 'value'):
        value = value.value
    if isinstance(value, Mapping):
        return tuple(sorted(
            (str(key), _normalized_setting(item))
            for key, item in value.items()
        ))
    if isinstance(value, (list, tuple)):
        return tuple(_normalized_setting(item) for item in value)
    return value


def _radii_digest(generators):
    payload = [
        [generator.system_index, generator.radius]
        for generator in generators
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(',', ':'), allow_nan=False).encode('utf-8')
    ).hexdigest()


def recover_network_source(
    network,
    source_path=None,
    *,
    expected_scheme='pow',
    expected_solvent_mode='MOLECULE_ONLY',
    expected_settings=None,
    expected_solve_universe_id=None,
    expected_radii_digest=None,
    expected_source_network_id=None,
) -> RecoveredSourceNetwork:
    """Validate and wrap a loaded solved network for source-view recovery.

    ``source_path`` defaults to the checkpoint's recorded original structure.
    The default contract is the current POW molecular-only recovery path; the
    expected settings can be overridden for another explicitly selected mode.
    The source is parsed only for identity verification; no Network or Group
    build is performed.

    ``expected_solve_universe_id`` and ``expected_radii_digest`` are optional
    compatibility checks for callers that hold an independent producer record;
    a public recovery command should supply them rather than accepting a
    checkpoint-only identity.
    """

    metadata = dict(getattr(network, 'archive_metadata', {}) or {})
    source_path = source_path or metadata.get('original_structure_path')
    validation = validate_network_source(
        network,
        source_path,
        expected_scheme=expected_scheme,
        expected_solvent_mode=expected_solvent_mode,
    )

    if expected_settings:
        actual_settings = dict(validation['settings'])
        for key, expected in expected_settings.items():
            if (
                key not in actual_settings
                or _normalized_setting(actual_settings[key]) != _normalized_setting(expected)
            ):
                raise ArchiveError(
                    f'Recovery setting mismatch for {key!r}: '
                    f'expected={expected!r}, found={actual_settings.get(key)!r}'
                )

    try:
        solve_context = network.get_solve_context_key()
        solve_universe_id = str(solve_context.solve_universe_id)
    except (AttributeError, TypeError, ValueError) as error:
        raise ArchiveError(
            'Checkpoint lacks a valid solve-universe identity'
        ) from error
    archived_solve_universe_id = metadata.get('solve_universe_id')
    if (
        archived_solve_universe_id is not None
        and str(archived_solve_universe_id) != solve_universe_id
    ):
        raise ArchiveError(
            'Checkpoint solve-universe identity does not match its geometry: '
            f'recorded={archived_solve_universe_id!r}, found={solve_universe_id!r}'
        )
    if expected_solve_universe_id is not None and solve_universe_id != str(
        expected_solve_universe_id
    ):
        raise ArchiveError(
            'Recovery solve-universe identity does not match the producer record: '
            f'expected={expected_solve_universe_id!r}, found={solve_universe_id!r}'
        )
    source_network_id = getattr(network, 'source_network_id', None)
    if expected_source_network_id is not None and source_network_id != expected_source_network_id:
        raise ArchiveError(
            'Recovery source-network identity does not match the producer record: '
            f'expected={expected_source_network_id!r}, found={source_network_id!r}'
        )

    # Keep this local to the persistence boundary.  It is the same read-only
    # source parse used by validation and is never passed to a solver.
    from vorpy.src.system.system import System

    source_system = System(
        validation['source_path'], make_dir=False, print_actions=False,
        solvent_mode=expected_solvent_mode,
    )
    source_balls = source_system.balls
    generator_indices = tuple(network.generator_system_indices)
    stable_ids = tuple(str(value) for value in network.generator_stable_ids)
    generators = []
    for system_index, stable_id in zip(generator_indices, stable_ids):
        row = source_balls.iloc[system_index]
        coordinates = tuple(float(value) for value in row['loc'])
        generators.append(SourceGeneratorIdentity(
            system_index=system_index,
            stable_id=stable_id,
            name=str(_scalar(row.get('name', '')) or ''),
            residue_name=str(_scalar(row.get('res_name', '')) or ''),
            residue_sequence=_scalar(row.get('res_seq')),
            chain_id=str(_scalar(row.get('chain_id', '')) or ''),
            element=str(_scalar(row.get('element', '')) or ''),
            coordinates=coordinates,
            radius=float(row['rad']),
        ))
    radii_digest = _radii_digest(generators)
    if expected_radii_digest is not None and radii_digest != str(expected_radii_digest):
        raise ArchiveError(
            'Recovery generator-radii identity does not match the producer record: '
            f'expected={expected_radii_digest!r}, found={radii_digest!r}'
        )
    validation = dict(validation)
    validation.update({
        'solve_universe_id': solve_universe_id,
        'radii_digest': radii_digest,
        'source_network_id': source_network_id,
    })

    return RecoveredSourceNetwork(
        network=network,
        source_identity=SourceSystemIdentity(
            source_path=validation['source_path'],
            source_sha256=validation['source_sha256'],
            frame_index=getattr(getattr(network, 'sys', None), 'frame_index', None),
            generator_count=validation['generator_count'],
            generator_system_indices=generator_indices,
            generator_stable_ids=stable_ids,
            generators=tuple(generators),
            scheme=validation['scheme'],
            solvent_mode=validation['solvent_mode'],
            boundary_mode=validation['boundary_mode'],
            complete_cell_count=validation['complete_cell_count'],
            complete_cell_evidence=validation['complete_cell_evidence'],
            geometry_status=validation['geometry_status'],
            analysis_status=validation['analysis_status'],
            solve_universe_id=solve_universe_id,
            radii_digest=radii_digest,
            generator_policy=validation['generator_policy'],
            source_network_id=source_network_id,
            settings=dict(validation['settings']),
        ),
        validation=dict(validation),
    )


__all__ = [
    'RecoveredSourceNetwork',
    'SourceGeneratorIdentity',
    'SourceSystemIdentity',
    'recover_network_source',
]
