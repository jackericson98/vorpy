"""Boundary modes and geometry providers for finite AW network construction.

The virtual solvent shell is a numerical bounding aid only. Its sites are not
physical solvent and are never inserted into the molecular atom/residue tables.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from itertools import product

import numpy as np


class BoundaryMode(str, Enum):
    NONE = "none"
    EXPLICIT = "explicit"
    SHELL = "shell"
    SES = "ses"
    SAS = "sas"
    PROBE = "probe"


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


def network_geometry(system):
    """Return geometry arrays augmented by the selected boundary provider.

    Physical System tables are left untouched. The returned final indices are
    artificial generators and are recorded on the resulting Network.
    """
    atoms = system.balls
    config = getattr(system, "boundary_config", None)
    if config is None:
        return atoms["loc"], atoms["rad"], atoms["mass"], (), None
    resolved = resolve_boundary(config, system.has_explicit_solvent)
    generators = virtual_shell(atoms, resolved) if resolved.mode is BoundaryMode.SHELL else ()
    system.boundary_config = resolved
    system.boundary_mode = resolved.mode.value
    system.boundary_generators = generators
    locs = list(atoms["loc"]) + [item.position for item in generators]
    radii = list(atoms["rad"]) + [item.radius for item in generators]
    masses = list(atoms["mass"]) + [0.0] * len(generators)
    indices = tuple(range(len(atoms), len(atoms) + len(generators)))
    return locs, radii, masses, indices, resolved
