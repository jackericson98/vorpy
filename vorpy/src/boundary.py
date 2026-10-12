"""Boundary modes and geometry providers for finite AW network construction.

The virtual solvent shell is a numerical bounding aid only. Its sites are not
physical solvent and are never inserted into the molecular atom/residue tables.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from itertools import product
import hashlib
import json

import numpy as np


class BoundaryMode(str, Enum):
    NONE = "none"
    EXPLICIT = "explicit"
    SHELL = "shell"
    SES = "ses"
    SAS = "sas"
    PROBE = "probe"


class SolventMode(str, Enum):
    """Scientific generator-universe modes for molecular network solves."""

    MOLECULE_ONLY = "MOLECULE_ONLY"
    SOLVENT_CONTEXT = "SOLVENT_CONTEXT"
    FULL_SOLVENT = "FULL_SOLVENT"


def normalize_retained_system_indices(values, atom_count=None):
    """Return a deterministic explicit allowlist of retained System rows."""
    if values is None:
        return ()
    if isinstance(values, (str, bytes)):
        raise ValueError("retained_system_indices must be an iterable of integers")
    try:
        values = tuple(values)
    except TypeError as error:
        raise ValueError(
            "retained_system_indices must be an iterable of integers"
        ) from error

    normalized = set()
    for value in values:
        if isinstance(value, bool):
            raise ValueError("retained_system_indices must contain integers")
        try:
            index = int(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(
                "retained_system_indices must contain integers"
            ) from error
        if str(value).strip() not in {str(index), f"{index}.0"}:
            try:
                if float(value) != float(index):
                    raise ValueError
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError(
                    "retained_system_indices must contain integers"
                ) from error
        if index < 0 or (atom_count is not None and index >= int(atom_count)):
            if atom_count is None:
                bound = "non-negative"
            else:
                bound = f"0..{int(atom_count) - 1}"
            raise IndexError(
                f"retained System index {index} is outside {bound}"
            )
        normalized.add(index)
    return tuple(sorted(normalized))


def normalize_solvent_mode(value=None) -> SolventMode:
    """Normalize the explicit solvent context used by a geometry producer."""
    if value is None:
        return SolventMode.MOLECULE_ONLY
    if isinstance(value, SolventMode):
        return value
    normalized = str(value).strip().upper().replace("-", "_")
    try:
        return SolventMode(normalized)
    except ValueError as error:
        raise ValueError(
            "solvent_mode must be MOLECULE_ONLY, SOLVENT_CONTEXT, or FULL_SOLVENT"
        ) from error


# Core loader and Workbench solvent names are merged here to keep detection
# consistent. Ions count as explicit surrounding solvent/boundary sites too.
WATER_RESIDUES = frozenset({"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "TIP4", "TIP4P", "SPC", "SPCE"})
ION_RESIDUES = frozenset({
    "LI", "NA", "NA+", "SOD", "K", "K+", "POT", "RB", "CS", "MG", "MG2", "MG2+", "CA", "CA2", "CA2+", "SR", "BA", "ZN", "ZN2", "ZN2+", "CD",
    "FE", "FE2", "FE3", "MN", "CU", "CU1", "CU2", "CO", "NI", "AL",
    "F", "CL", "CL-", "BR", "I", "IOD", "SO4", "PO4", "NH4", "CLA", "ION", "OUT", "MN2",
})
SOLVENT_RESIDUES = WATER_RESIDUES | ION_RESIDUES

DEFAULT_PROBE_RADIUS = 1.4
# Slightly less than the 3 Å diameter of a 1.5 Å boundary site avoids exact
# tangencies between neighboring generators, while keeping the grid practical.
DEFAULT_SHELL_SPACING = 2.5
# A tiny outward perturbation avoids exact coplanar regular lattices, which can
# create degenerate AW junctions, while retaining a box-like enclosing shell.
SHELL_PERTURBATION_FRACTION = 0.02
# Two probe radii place the numerical boundary about one probe diameter beyond
# the van-der-Waals atomic union. This is not an SES/SAS construction.
DEFAULT_SHELL_OFFSET = 2.8
# Explicit HOH oxygen uses the core elemental O radius (1.50 Å).
DEFAULT_BOUNDARY_GENERATOR_RADIUS = 1.5


@dataclass(frozen=True)
class BoundaryConfig:
    mode: BoundaryMode | str | None = None
    probe_radius: float = DEFAULT_PROBE_RADIUS
    shell_spacing: float = DEFAULT_SHELL_SPACING
    shell_offset: float = DEFAULT_SHELL_OFFSET

    def __post_init__(self):
        if self.mode is not None:
            object.__setattr__(self, "mode", BoundaryMode(self.mode))
        for name in ("probe_radius", "shell_spacing", "shell_offset"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name.replace('_', ' ').capitalize()} must be positive")
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class BoundaryGenerator:
    position: tuple[float, float, float]
    radius: float = DEFAULT_BOUNDARY_GENERATOR_RADIUS
    is_boundary_generator: bool = True


@dataclass(frozen=True)
class NetworkGeometry:
    """Immutable generator inputs and membership evidence for one solve."""

    locations: tuple
    radii: tuple
    masses: tuple
    system_indices: tuple
    stable_ids: tuple
    solvent_local_indices: tuple
    solvent_generator_ids: tuple
    retained_system_indices: tuple
    retained_generator_ids: tuple
    boundary_indices: tuple
    boundary_generators: tuple
    boundary_config: BoundaryConfig | None
    solvent_mode: str
    solvent_cells_complete: bool

    @property
    def solvent_generator_digest(self):
        if not self.solvent_generator_ids:
            return None
        payload = json.dumps(
            list(self.solvent_generator_ids),
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    @property
    def participating_solvent_generator_digest(self):
        """Compatibility alias for the canonical contract field."""
        return self.solvent_generator_digest


def has_explicit_solvent(atoms) -> bool:
    """Return whether physical atom records contain a recognized solvent/ion residue."""
    if atoms is None:
        return False
    if hasattr(atoms, "columns"):
        if "res_name" not in atoms.columns:
            return False
        names = atoms["res_name"].fillna("")
    else:
        names = (getattr(atom, "residue_name", "") for atom in atoms)
    return any(str(name).strip().upper() in SOLVENT_RESIDUES for name in names)


def resolve_boundary(config: BoundaryConfig, solvent_present: bool) -> BoundaryConfig:
    """Resolve auto mode and reject unavailable or impossible requested modes."""
    mode = config.mode or BoundaryMode.NONE
    if mode in {BoundaryMode.SES, BoundaryMode.SAS, BoundaryMode.PROBE}:
        raise NotImplementedError(f"Boundary mode '{mode.value}' is planned but not implemented")
    if mode is BoundaryMode.EXPLICIT and not solvent_present:
        raise ValueError("Explicit solvent boundary requested, but no explicit solvent was detected.")
    return BoundaryConfig(mode, config.probe_radius, config.shell_spacing, config.shell_offset)


def virtual_shell(atoms, config: BoundaryConfig) -> tuple[BoundaryGenerator, ...]:
    """Generate a deterministic, deduplicated six-face box shell around atoms."""
    if atoms is None or len(atoms) == 0:
        raise ValueError("Cannot generate a boundary shell for an empty molecular structure")
    if hasattr(atoms, "columns"):
        coords = np.asarray(list(atoms["loc"]), dtype=float)
        radii = np.asarray(atoms["rad"], dtype=float)
    else:
        coords = np.asarray([atom.loc for atom in atoms], dtype=float)
        radii = np.asarray([atom.rad for atom in atoms], dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or not np.isfinite(coords).all():
        raise ValueError("Molecular coordinates must be finite 3D points")
    lower = np.min(coords - radii[:, None], axis=0) - config.shell_offset
    upper = np.max(coords + radii[:, None], axis=0) + config.shell_offset
    axes = [np.linspace(lo, hi, max(2, int(np.ceil((hi - lo) / config.shell_spacing)) + 1))
            for lo, hi in zip(lower, upper)]
    points = set()
    for axis in range(3):
        other = [i for i in range(3) if i != axis]
        for side in (lower[axis], upper[axis]):
            for a, b in product(axes[other[0]], axes[other[1]]):
                point = np.empty(3, dtype=float)
                point[axis] = side
                point[other[0]], point[other[1]] = a, b
                points.add(tuple(float(round(value, 10)) for value in point))
    center = 0.5 * (lower + upper)
    generators = []
    for raw_point in sorted(points):
        point = np.asarray(raw_point, dtype=float)
        direction = point - center
        direction /= np.linalg.norm(direction)
        phase = float(np.dot(point, np.array([12.9898, 78.233, 37.719])))
        fraction = 0.5 + 0.5 * np.sin(phase)
        point += direction * config.shell_spacing * SHELL_PERTURBATION_FRACTION * fraction
        generators.append(BoundaryGenerator(tuple(float(round(value, 10)) for value in point)))
    return tuple(generators)


def _row_value(row, name, default=None):
    if hasattr(row, "get"):
        return row.get(name, default)
    return getattr(row, name, default)


def _stable_row_id(row, system_index):
    for field in ("stable_id", "atom_stable_id", "generator_id"):
        value = _row_value(row, field)
        if value is not None and str(value) != "":
            return str(value)
    return f"system:{int(system_index)}"


def _solvent_row(row):
    return str(_row_value(row, "res_name", "")).strip().upper() in SOLVENT_RESIDUES


def _system_solvent_present(system, atoms):
    value = getattr(system, "has_explicit_solvent", None)
    if callable(value):
        return bool(value())
    if value is not None:
        return bool(value)
    return has_explicit_solvent(atoms)


def prepare_network_geometry(
    system,
    solvent_mode=None,
    retained_system_indices=None,
) -> NetworkGeometry:
    """Prepare one immutable physical-generator universe for a Network.

    ``System.balls`` is never filtered or rewritten.  The returned
    ``system_indices`` retain the source row for every physical generator;
    ``None`` identifies appended numerical boundary generators.
    """
    atoms = system.balls
    if atoms is None:
        raise ValueError("Cannot prepare network geometry without system balls")

    requested_mode = solvent_mode
    if requested_mode is None:
        requested_mode = getattr(system, "solvent_mode", None)
    # Direct legacy callers without a System mode retain the historical
    # all-physical-generator universe; System instances explicitly default to
    # MOLECULE_ONLY.
    if requested_mode is None:
        requested_mode = SolventMode.SOLVENT_CONTEXT
    mode = normalize_solvent_mode(requested_mode)

    rows = [atoms.iloc[index] if hasattr(atoms, "iloc") else atoms[index]
            for index in range(len(atoms))]
    requested_retained = (
        getattr(system, "retained_system_indices", ())
        if retained_system_indices is None
        else retained_system_indices
    )
    retained_system_indices = normalize_retained_system_indices(
        requested_retained,
        atom_count=len(rows),
    )
    retained_system_set = set(retained_system_indices)
    solvent_system_indices = tuple(
        index for index, row in enumerate(rows) if _solvent_row(row)
    )
    solvent_system_set = set(solvent_system_indices)
    if mode is SolventMode.MOLECULE_ONLY:
        physical_indices = tuple(
            index for index in range(len(rows))
            if index not in solvent_system_set or index in retained_system_set
        )
    else:
        physical_indices = tuple(range(len(rows)))

    physical_rows = [rows[index] for index in physical_indices]
    locs = [_row_value(row, "loc") for row in physical_rows]
    radii = [_row_value(row, "rad") for row in physical_rows]
    masses = [_row_value(row, "mass", 0.0) for row in physical_rows]
    stable_ids = [
        _stable_row_id(row, system_index)
        for row, system_index in zip(physical_rows, physical_indices)
    ]
    solvent_local_indices = tuple(
        local_index for local_index, system_index in enumerate(physical_indices)
        if system_index in solvent_system_set
        and system_index not in retained_system_set
    )
    solvent_generator_ids = tuple(
        stable_ids[local_index] for local_index in solvent_local_indices
    )
    retained_generator_ids = tuple(
        stable_ids[local_index]
        for local_index, system_index in enumerate(physical_indices)
        if system_index in retained_system_set
    )

    config = getattr(system, "boundary_config", None)
    resolved = None
    generators = ()
    if config is not None:
        resolved = resolve_boundary(config, _system_solvent_present(system, atoms))
        if resolved.mode is BoundaryMode.SHELL:
            if hasattr(atoms, "iloc"):
                selected_atoms = atoms.iloc[list(physical_indices)]
            else:
                selected_atoms = [rows[index] for index in physical_indices]
            generators = virtual_shell(selected_atoms, resolved)
        system.boundary_config = resolved
        system.boundary_mode = resolved.mode.value
        system.boundary_generators = generators

    boundary_start = len(physical_indices)
    boundary_indices = tuple(
        range(boundary_start, boundary_start + len(generators))
    )
    locs.extend(item.position for item in generators)
    radii.extend(item.radius for item in generators)
    masses.extend(0.0 for _ in generators)
    system_indices = physical_indices + (None,) * len(generators)
    stable_ids.extend(
        f"boundary:{index}" for index in range(len(generators))
    )

    return NetworkGeometry(
        locations=tuple(locs),
        radii=tuple(radii),
        masses=tuple(masses),
        system_indices=tuple(system_indices),
        stable_ids=tuple(stable_ids),
        solvent_local_indices=solvent_local_indices,
        solvent_generator_ids=solvent_generator_ids,
        retained_system_indices=retained_system_indices,
        retained_generator_ids=retained_generator_ids,
        boundary_indices=boundary_indices,
        boundary_generators=tuple(generators),
        boundary_config=resolved,
        solvent_mode=mode.value,
        solvent_cells_complete=(mode is SolventMode.FULL_SOLVENT),
    )


def network_geometry(system, solvent_mode=None):
    """Return geometry arrays augmented by the selected boundary provider.

    This compatibility wrapper keeps the historical five-value return shape;
    new producers should use :func:`prepare_network_geometry` for membership
    and solvent provenance.
    """
    geometry = prepare_network_geometry(system, solvent_mode=solvent_mode)
    return (
        list(geometry.locations),
        list(geometry.radii),
        list(geometry.masses),
        geometry.boundary_indices,
        geometry.boundary_config,
    )
