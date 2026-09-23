"""
Synthetic geometries for validating integrated surface curvatures.

This module generates deterministic XYZR point sets representing simple
geometries with analytically known integrated mean and Gaussian curvatures.

The generated systems are intended for:

1. Unit/regression testing of VorPy curvature calculations.
2. Numerical convergence studies.
3. Generation of validation data and figures for publications.

No VorPy network-building code should be imported here.  This module should
remain a lightweight geometry generator so that it can be imported by both
tests and analysis scripts without circular dependencies.

Curvature convention
--------------------
Mean curvature is defined as

    H = (k1 + k2) / 2

and the integrated mean curvature is

    C = integral(H dA).

Gaussian curvature is

    K = k1 * k2

and the integrated Gaussian curvature is

    X = integral(K dA).

For a closed orientable surface of genus g,

    integral(K dA) = 4*pi*(1 - g).

All coordinates and radii use the same arbitrary length unit.  In VorPy this
will normally be Angstroms.
"""

from __future__ import annotations


import argparse

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import sys

import numpy as np

# When invoked as ``python path/to/validation_shapes.py``, Python places the
# geometry directory first on sys.path and can otherwise import an unrelated
# installed VorPy package. Prefer this checkout for the lightweight PDB
# writer used by the generator.
if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from vorpy.src.output import make_pdb_line


# ============================================================================
# Data container
# ============================================================================


@dataclass
class ValidationShape:
    """
    Synthetic geometry with known analytic curvature values.

    Attributes
    ----------
    name
        Human-readable name of the geometry.

    xyzr
        Array with shape (N, 4). Columns are x, y, z, radius.

    expected_int_mean_curvature
        Analytic integrated mean curvature:

            integral(H dA)

        using H = (k1 + k2) / 2.

    expected_int_gaussian_curvature
        Analytic integrated Gaussian curvature:

            integral(K dA)

    parameters
        Parameters used to generate the system.

    notes
        Optional information describing the construction or interpretation.
    """

    name: str
    xyzr: np.ndarray

    expected_int_mean_curvature: float
    expected_int_gaussian_curvature: float

    parameters: dict[str, Any] = field(default_factory=dict)
    notes: str = ""

    @property
    def n_atoms(self) -> int:
        """Number of XYZR spheres in the system."""
        return int(self.xyzr.shape[0])

    def save_xyzr(
        self,
        filename: str | Path,
        precision: int = 8,
    ) -> Path:
        """
        Save the generated system as a whitespace-separated XYZR file.

        Parameters
        ----------
        filename
            Destination file.

        precision
            Number of digits after the decimal point.

        Returns
        -------
        pathlib.Path
            Path to the written file.
        """
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)

        fmt = f"%.{precision}f"

        np.savetxt(
            path,
            self.xyzr,
            fmt=fmt,
            header="x y z radius",
            comments="# ",
        )

        return path

    def save_pdb(
            self,
            filename: str | Path,
            residue_name: str = "VAL",
            chain: str = "A",
    ) -> Path:
        """
        Save the synthetic XYZR system as a PDB for visualization in PyMOL.

        The XYZR sphere radius is stored in the PDB B-factor field so that the
        original dummy-sphere radius remains available in the visualization file.

        Parameters
        ----------
        filename
            Destination PDB filename.

        residue_name
            Residue name assigned to the synthetic atoms.

        chain
            Chain identifier.

        Returns
        -------
        pathlib.Path
            Path to the written PDB file.
        """
        path = Path(filename)

        if path.suffix.lower() != ".pdb":
            path = path.with_suffix(".pdb")

        path.parent.mkdir(parents=True, exist_ok=True)

        with path.open("w") as pdb_file:
            pdb_file.write(
                f"HEADER    VORPY CURVATURE VALIDATION - "
                f"{self.name.upper()}\n"
            )

            pdb_file.write(
                f"REMARK    EXPECTED INT MEAN CURVATURE "
                f"{self.expected_int_mean_curvature:.12g}\n"
            )

            pdb_file.write(
                f"REMARK    EXPECTED INT GAUSSIAN CURVATURE "
                f"{self.expected_int_gaussian_curvature:.12g}\n"
            )

            for i, (x, y, z, radius) in enumerate(
                    self.xyzr,
                    start=1,
            ):
                interior_indices = self.parameters.get("interior_indices")
                interior_count = (
                    len(interior_indices)
                    if interior_indices is not None
                    else int(self.parameters.get("interior_count", 0))
                )
                is_interior = self.name in {
                    "spherocylinder_multicell", "torus_multicell",
                } and i <= interior_count
                is_center = (
                        np.isclose(x, 0.0)
                        and np.isclose(y, 0.0)
                        and np.isclose(z, 0.0)
                )

                if is_interior:
                    atom_name = "INT"
                    current_residue = "VAL"
                    current_res_seq = 1
                    current_chain = "I"
                else:
                    atom_name = "CTR" if is_center else "SUR"
                    current_residue = "CTR" if is_center else residue_name
                    current_res_seq = 1 if is_center else 2
                    current_chain = chain

                pdb_file.write(
                    make_pdb_line(
                        atom="HETATM",
                        ser_num=i,
                        name=atom_name,
                        res_name=current_residue,
                        chain=current_chain,
                        res_seq=current_res_seq,
                        x=float(x),
                        y=float(y),
                        z=float(z),
                        occ=1.0,
                        tfact=float(radius),
                        elem="C",
                    )
                )

            pdb_file.write("END\n")

        return path

    def save_pymol_script(
            self,
            filename: str | Path,
            pdb_filename: str | Path,
    ) -> Path:
        """
        Write a PyMOL script for viewing the generated validation system.

        Sphere radii are read from the PDB B-factor field.
        """
        path = Path(filename)

        if path.suffix.lower() != ".pml":
            path = path.with_suffix(".pml")

        path.parent.mkdir(parents=True, exist_ok=True)

        pdb_filename = Path(pdb_filename).name
        object_name = self.name

        script = f"""\
    # VorPy curvature validation visualization
    # Shape: {self.name}

    reinitialize

    load {pdb_filename}, {object_name}

    hide everything, {object_name}

    # The XYZR radius is stored in the B-factor.
    alter {object_name}, vdw=b
    rebuild

    show spheres, {object_name}

    # Surrounding geometry
    select surrounding, {object_name} and not resn CTR
    color gray70, surrounding

    # Central VorPy validation ball
    select center_ball, {object_name} and resn CTR
    color red, center_ball

    show spheres, center_ball
    show spheres, surrounding

    # Clean display
    bg_color white
    set ray_opaque_background, off
    set sphere_quality, 3

    orient {object_name}
    zoom {object_name}

    deselect
    """

        path.write_text(script)

        return path

    def save(
            self,
            directory: str | Path,
            basename: str | None = None,
    ) -> tuple[Path, Path, Path]:
        """
        Save a complete validation-shape package.

        Each generated system receives its own directory containing:

            <name>.xyzr
            <name>.pdb
            <name>.pml
        """
        directory = Path(directory)

        if basename is None:
            basename = self.name

        # Each validation system gets its own folder.
        shape_directory = directory / basename
        shape_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        xyzr_path = self.save_xyzr(
            shape_directory / f"{basename}.xyzr"
        )

        pdb_path = self.save_pdb(
            shape_directory / f"{basename}.pdb"
        )

        pymol_path = self.save_pymol_script(
            shape_directory / f"{basename}.pml",
            pdb_path,
        )

        return xyzr_path, pdb_path, pymol_path


# ============================================================================
# Internal helpers
# ============================================================================


def _validate_positive(value: float, name: str) -> float:
    value = float(value)

    if value <= 0.0:
        raise ValueError(f"{name} must be > 0; received {value}.")

    return value


def _validate_resolution(resolution: int, minimum: int = 1) -> int:
    resolution = int(resolution)

    if resolution < minimum:
        raise ValueError(
            f"resolution must be >= {minimum}; received {resolution}."
        )

    return resolution


def _make_xyzr(
        xyz: np.ndarray,
        sphere_radius: float,
    ) -> np.ndarray:
    """Append a common sphere radius to an (N, 3) coordinate array."""
    sphere_radius = _validate_positive(sphere_radius, "sphere_radius")

    xyz = np.asarray(xyz, dtype=float)

    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise ValueError(
            f"xyz must have shape (N, 3); received {xyz.shape}."
        )

    radii = np.full((xyz.shape[0], 1), sphere_radius, dtype=float)

    return np.hstack((xyz, radii))


def _unique_points(
        xyz: np.ndarray,
        decimals: int = 12,
    ) -> np.ndarray:
    """
    Remove duplicate Cartesian points.

    This is primarily useful for geometries constructed face-by-face, where
    edge and corner points may otherwise appear multiple times.
    """
    xyz = np.asarray(xyz, dtype=float)

    rounded = np.round(xyz, decimals=decimals)

    _, indices = np.unique(
        rounded,
        axis=0,
        return_index=True,
    )

    indices = np.sort(indices)

    return xyz[indices]


def _power_neighbor_distance(face_distance: float, center_radius: float,
                             neighbor_radius: float) -> float:
    """Return the generator spacing for a power bisector at ``face_distance``."""
    h = _validate_positive(face_distance, "face_distance")
    r0 = _validate_positive(center_radius, "center_radius")
    r1 = _validate_positive(neighbor_radius, "neighbor_radius")
    discriminant = h * h - r0 * r0 + r1 * r1
    if discriminant < 0.0:
        raise ValueError(
            "The requested power-cell face is incompatible with the supplied "
            "center and neighbor radii."
        )
    return h + float(np.sqrt(discriminant))


# ============================================================================
# Sphere
# ============================================================================


def sphere(
    radius: float,
    resolution: int,
    ball_radius: float = 1.0,
    center_radius: float = 1.0,
    include_center: bool = True,
) -> ValidationShape:
    """
    Generate a central validation ball surrounded by approximately uniformly
    distributed balls on a spherical shell.

    Parameters
    ----------
    radius
        Target radius of the central power cell. The surrounding generator
        centers are placed farther away as required by the power-bisector
        equation.

    resolution
        Number of surrounding balls.

    ball_radius
        Radius assigned to each surrounding ball.

    center_radius
        Radius assigned to the central validation ball.

    include_center
        If True, include the central ball at the origin.

    Returns
    -------
    ValidationShape
    """
    radius = _validate_positive(radius, "radius")
    ball_radius = _validate_positive(ball_radius, "ball_radius")
    center_radius = _validate_positive(center_radius, "center_radius")
    resolution = _validate_resolution(resolution, minimum=4)
    surrounding_distance = _power_neighbor_distance(
        radius, center_radius, ball_radius
    )

    indices = np.arange(resolution, dtype=float)

    golden_angle = np.pi * (3.0 - np.sqrt(5.0))

    z = 1.0 - 2.0 * (indices + 0.5) / resolution

    radial_xy = np.sqrt(
        np.maximum(0.0, 1.0 - z * z)
    )

    theta = golden_angle * indices

    x = surrounding_distance * radial_xy * np.cos(theta)
    y = surrounding_distance * radial_xy * np.sin(theta)
    z = surrounding_distance * z

    surrounding_xyz = np.column_stack((x, y, z))

    surrounding_xyzr = _make_xyzr(
        xyz=surrounding_xyz,
        sphere_radius=ball_radius,
    )

    if include_center:
        center_xyzr = np.array(
            [[0.0, 0.0, 0.0, center_radius]],
            dtype=float,
        )

        xyzr = np.vstack(
            (
                center_xyzr,
                surrounding_xyzr,
            )
        )
    else:
        xyzr = surrounding_xyzr

    return ValidationShape(
        name="sphere",
        xyzr=xyzr,
        expected_int_mean_curvature=4.0 * np.pi * radius,
        expected_int_gaussian_curvature=4.0 * np.pi,
        parameters={
            "radius": radius,
            "target_power_cell_radius": radius,
            "surrounding_generator_distance": surrounding_distance,
            "resolution": resolution,
            "ball_radius": ball_radius,
            "center_radius": center_radius,
            "include_center": include_center,
        },
        notes=(
            "One central validation ball surrounded by balls whose centers "
            "lie on a Fibonacci spherical shell. The shell radius is chosen "
            "so the power bisectors form the requested target sphere."
        ),
    )


# ============================================================================
# Rectangular box / cube
# ============================================================================


def box(
    length: float,
    width: float,
    height: float,
    resolution: int = 1,
    sphere_radius: float = 1.0,
    center_radius: float = 1.0,
    neighbor_radius: float | None = None,
    include_center: bool = True,
) -> ValidationShape:
    """
    Generate seven weighted generators for an exact central rectangular
    power cell. The six neighboring generators lie on the coordinate axes;
    their power bisectors with the central generator are the box faces.

    Parameters
    ----------
    length
        Full x dimension.

    width
        Full y dimension.

    height
        Full z dimension.

    resolution
        Retained for CLI/API compatibility. Exact boxes do not sample faces.

    sphere_radius
        Radius assigned to each neighboring generator.

    center_radius
        Radius assigned to the central generator.

    neighbor_radius
        Explicit neighboring-generator radius. If omitted, ``sphere_radius``
        is used.

    include_center
        If True, add one sphere at the origin.

    Returns
    -------
    ValidationShape

    Analytic values
    ---------------
    For a rectangular box with side lengths a, b, c:

        integral(H dA) = pi * (a + b + c)

        integral(K dA) = 4*pi

    using the convex-body convention and H = (k1 + k2)/2.
    """
    length = _validate_positive(length, "length")
    width = _validate_positive(width, "width")
    height = _validate_positive(height, "height")

    resolution = int(resolution)
    if resolution < 1:
        raise ValueError("resolution must be >= 1")
    center_radius = _validate_positive(center_radius, "center_radius")
    if neighbor_radius is None:
        neighbor_radius = sphere_radius
    neighbor_radius = _validate_positive(neighbor_radius, "neighbor_radius")

    hx = 0.5 * length
    hy = 0.5 * width
    hz = 0.5 * height

    dx = _power_neighbor_distance(hx, center_radius, neighbor_radius)
    dy = _power_neighbor_distance(hy, center_radius, neighbor_radius)
    dz = _power_neighbor_distance(hz, center_radius, neighbor_radius)
    xyzr = np.array([
        [0.0, 0.0, 0.0, center_radius],
        [dx, 0.0, 0.0, neighbor_radius],
        [-dx, 0.0, 0.0, neighbor_radius],
        [0.0, dy, 0.0, neighbor_radius],
        [0.0, -dy, 0.0, neighbor_radius],
        [0.0, 0.0, dz, neighbor_radius],
        [0.0, 0.0, -dz, neighbor_radius],
    ], dtype=float)
    if not include_center:
        xyzr = xyzr[1:]

    expected_mean = np.pi * (
        length + width + height
    )

    return ValidationShape(
        name="box",
        xyzr=xyzr,
        expected_int_mean_curvature=expected_mean,
        expected_int_gaussian_curvature=4.0 * np.pi,
        parameters={
            "length": length,
            "width": width,
            "height": height,
            "resolution": resolution,
            "sphere_radius": neighbor_radius,
            "center_radius": center_radius,
            "neighbor_radius": neighbor_radius,
            "include_center": include_center,
        },
        notes=(
            "Seven weighted generators define the exact central rectangular "
            "power cell; analytic values refer to that cell."
        ),
    )


def cube(
    side_length: float,
    resolution: int = 1,
    sphere_radius: float = 1.0,
    center_radius: float = 1.0,
    neighbor_radius: float | None = None,
    include_center: bool = True,
) -> ValidationShape:
    """
    Generate a cubic special case of :func:`box`.

    Analytic values
    ---------------

        integral(H dA) = 3*pi*L

        integral(K dA) = 4*pi
    """
    shape = box(
        length=side_length,
        width=side_length,
        height=side_length,
        resolution=resolution,
        sphere_radius=sphere_radius,
        center_radius=center_radius,
        neighbor_radius=neighbor_radius,
        include_center=include_center,
    )

    shape.name = "cube"

    shape.parameters = {
        "side_length": float(side_length),
        "resolution": int(resolution),
        "sphere_radius": float(sphere_radius),
        "center_radius": float(center_radius),
        "neighbor_radius": float(
            sphere_radius if neighbor_radius is None else neighbor_radius
        ),
        "include_center": bool(include_center),
    }

    return shape


# ============================================================================
# Spherocylinder
# ============================================================================


def spherocylinder(
    radius: float,
    cylinder_length: float,
    resolution: int,
    sphere_radius: float = 1.0,
    center_radius: float = 1.0,
    neighbor_radius: float | None = None,
    include_center: bool = True,
) -> ValidationShape:
    """
    Generate weighted generators supporting an exact central spherocylinder
    power cell.

    One generator is placed at the origin and ``resolution`` surrounding
    generators are placed along Fibonacci-distributed support normals. Their
    weighted bisectors with the center lie on the target spherocylinder.

    Parameters
    ----------
    radius
        Radius R.

    cylinder_length
        Length L of the cylindrical section only.

    resolution
        Approximate number of supporting planes and surrounding generators.

    sphere_radius
        Backward-compatible alias for ``neighbor_radius``.

    center_radius
        Radius of the central generator.

    neighbor_radius
        Radius of each surrounding generator. If omitted, ``sphere_radius``
        is used.

    include_center
        If True, include the central generator at index 0.

    Returns
    -------
    ValidationShape

    Analytic values
    ---------------
    For a spherocylinder:

        integral(H dA) = pi*L + 4*pi*R

        integral(K dA) = 4*pi
    """
    radius = _validate_positive(radius, "radius")
    cylinder_length = _validate_positive(
        cylinder_length,
        "cylinder_length",
    )

    resolution = _validate_resolution(resolution, minimum=4)
    center_radius = _validate_positive(center_radius, "center_radius")
    if neighbor_radius is None:
        neighbor_radius = sphere_radius
    neighbor_radius = _validate_positive(neighbor_radius, "neighbor_radius")

    indices = np.arange(resolution, dtype=float)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    nz = 1.0 - 2.0 * (indices + 0.5) / resolution
    radial_xy = np.sqrt(np.maximum(0.0, 1.0 - nz * nz))
    theta = golden_angle * indices
    normals = np.column_stack((
        radial_xy * np.cos(theta),
        radial_xy * np.sin(theta),
        nz,
    ))
    # Retain the Fibonacci distribution while including the three useful
    # validation directions exactly: both axial poles and one equatorial
    # normal. This makes the analytic extrema directly testable without
    # introducing duplicate pole samples.
    normals[0] = (0.0, 0.0, 1.0)
    normals[-1] = (0.0, 0.0, -1.0)
    equator_index = resolution // 2
    equator_angle = golden_angle * equator_index
    normals[equator_index] = (np.cos(equator_angle), np.sin(equator_angle), 0.0)
    support_distances = radius + 0.5 * cylinder_length * np.abs(normals[:, 2])
    discriminants = support_distances * support_distances - center_radius * center_radius + neighbor_radius * neighbor_radius
    if np.any(discriminants < 0.0):
        raise ValueError(
            "The requested spherocylinder support planes are incompatible "
            "with the supplied center and neighbor radii."
        )
    generator_distances = support_distances + np.sqrt(discriminants)
    surrounding_xyz = normals * generator_distances[:, None]
    surrounding_xyzr = _make_xyzr(surrounding_xyz, neighbor_radius)

    if include_center:
        xyzr = np.vstack((
            np.array([[0.0, 0.0, 0.0, center_radius]], dtype=float),
            surrounding_xyzr,
        ))
    else:
        xyzr = surrounding_xyzr

    expected_mean = (
        np.pi * cylinder_length
        + 4.0 * np.pi * radius
    )

    return ValidationShape(
        name="spherocylinder",
        xyzr=xyzr,
        expected_int_mean_curvature=expected_mean,
        expected_int_gaussian_curvature=4.0 * np.pi,
        parameters={
            "radius": radius,
            "cylinder_length": cylinder_length,
            "resolution": resolution,
            "sphere_radius": neighbor_radius,
            "center_radius": center_radius,
            "neighbor_radius": neighbor_radius,
            "include_center": include_center,
            "support_plane_count": resolution,
        },
        notes=(
            "One central weighted generator and Fibonacci-distributed "
            "supporting generators define the spherocylinder power cell."
        ),
    )


def spherocylinder_multicell(
    radius: float,
    cylinder_length: float,
    interior_count: int = 3,
    angular_resolution: int = 12,
    cap_resolution: int = 12,
    axial_ring_count: int = 3,
    sphere_radius: float = 1.0,
    center_radius: float = 1.0,
    neighbor_radius: float | None = None,
    include_center: bool = True,
) -> ValidationShape:
    """Generate a coarse multi-cell spherocylinder validation system.

    The interior generators form a selectable axial chain. Radial constraints
    are repeated along that chain to define the cylindrical side, while
    positive and negative cap normals are assigned to the corresponding end
    generator. Exterior generators are constraints only and are excluded from
    the interior group metadata.
    """
    radius = _validate_positive(radius, "radius")
    cylinder_length = _validate_positive(cylinder_length, "cylinder_length")
    interior_count = _validate_resolution(interior_count, minimum=2)
    angular_resolution = _validate_resolution(angular_resolution, minimum=8)
    cap_resolution = _validate_resolution(cap_resolution, minimum=4)
    axial_ring_count = _validate_resolution(axial_ring_count, minimum=1)
    center_radius = _validate_positive(center_radius, "center_radius")
    if neighbor_radius is None:
        neighbor_radius = sphere_radius
    neighbor_radius = _validate_positive(neighbor_radius, "neighbor_radius")

    z_positions = np.linspace(-0.5 * cylinder_length, 0.5 * cylinder_length, interior_count)
    interior_xyz = np.column_stack((
        np.zeros(interior_count),
        np.zeros(interior_count),
        z_positions,
    ))
    interior_radii = np.full((interior_count, 1), center_radius, dtype=float)
    interior_xyzr = np.hstack((interior_xyz, interior_radii))

    # Use a complete, non-coplanar Fibonacci shell around every interior cell.
    # This is intentionally redundant: each cell is independently bounded,
    # while the union of the axial cells forms a short, fat capsule.
    shell_count = angular_resolution * axial_ring_count + cap_resolution
    shell_indices = np.arange(shell_count, dtype=float)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    shell_z = 1.0 - 2.0 * (shell_indices + 0.5) / shell_count
    shell_xy = np.sqrt(np.maximum(0.0, 1.0 - shell_z * shell_z))
    shell_theta = golden_angle * shell_indices
    shell_normals = np.column_stack((
        shell_xy * np.cos(shell_theta),
        shell_xy * np.sin(shell_theta),
        shell_z,
    ))
    assignments: list[int] = []
    normals: list[np.ndarray] = []
    origins: list[np.ndarray] = []
    for cell_index, center in enumerate(interior_xyz):
        for normal in shell_normals:
            origins.append(interior_xyz[cell_index])
            normals.append(normal)
            assignments.append(cell_index)

    normals_array = np.asarray(normals, dtype=float)
    origins_array = np.asarray(origins, dtype=float)
    support_offsets = np.full(len(normals_array), radius, dtype=float)
    discriminants = support_offsets**2 - center_radius**2 + neighbor_radius**2
    if np.any(discriminants < 0.0):
        raise ValueError(
            "The requested multi-cell spherocylinder supports are incompatible "
            "with the supplied interior and exterior radii."
        )
    distances = support_offsets + np.sqrt(discriminants)
    exterior_xyz = origins_array + distances[:, None] * normals_array
    exterior_xyzr = _make_xyzr(exterior_xyz, neighbor_radius)

    if include_center:
        xyzr = np.vstack((interior_xyzr, exterior_xyzr))
        interior_indices = list(range(interior_count))
        exterior_start = interior_count
    else:
        xyzr = exterior_xyzr
        interior_indices = []
        exterior_start = 0

    return ValidationShape(
        name="spherocylinder_multicell",
        xyzr=xyzr,
        expected_int_mean_curvature=np.pi * cylinder_length + 4.0 * np.pi * radius,
        expected_int_gaussian_curvature=4.0 * np.pi,
        parameters={
            "radius": radius,
            "cylinder_length": cylinder_length,
            "interior_count": interior_count,
            "angular_resolution": angular_resolution,
            "cap_resolution": cap_resolution,
            "axial_ring_count": axial_ring_count,
            "sphere_radius": neighbor_radius,
            "center_radius": center_radius,
            "neighbor_radius": neighbor_radius,
            "include_center": include_center,
            "interior_indices": interior_indices,
            "interior_coordinates": interior_xyz,
            "exterior_start": exterior_start,
            "exterior_count": len(exterior_xyz),
            "constraint_assignments": assignments,
            "target_volume": np.pi * radius * radius * cylinder_length + (4.0 / 3.0) * np.pi * radius**3,
            "target_surface_area": 2.0 * np.pi * radius * cylinder_length + 4.0 * np.pi * radius**2,
        },
        notes=(
            "Interior axial generators form the selected group. Exterior "
            "generators are constraints only; internal group faces must be "
            "excluded from external surface measures."
        ),
    )


# ============================================================================
# Torus
# ============================================================================


def torus(
    major_radius: float,
    minor_radius: float,
    major_resolution: int,
    minor_resolution: int,
    sphere_radius: float = 1.0,
) -> ValidationShape:
    """
    Generate points on a standard ring torus.

    Parameters
    ----------
    major_radius
        Distance R from the center of the torus to the centerline of the tube.

    minor_radius
        Tube radius r.

    major_resolution
        Number of samples around the major circle.

    minor_resolution
        Number of samples around the tube.

    sphere_radius
        Radius assigned to every generated XYZR sphere.

    Returns
    -------
    ValidationShape

    Analytic values
    ---------------
    For a ring torus with R > r:

        integral(H dA) = 2*pi^2*R

        integral(K dA) = 0

    for the outward-oriented surface and H = (k1 + k2)/2.
    """
    major_radius = _validate_positive(
        major_radius,
        "major_radius",
    )

    minor_radius = _validate_positive(
        minor_radius,
        "minor_radius",
    )

    if major_radius <= minor_radius:
        raise ValueError(
            "A ring torus requires major_radius > minor_radius."
        )

    major_resolution = _validate_resolution(
        major_resolution,
        minimum=3,
    )

    minor_resolution = _validate_resolution(
        minor_resolution,
        minimum=3,
    )

    u_values = np.linspace(
        0.0,
        2.0 * np.pi,
        major_resolution,
        endpoint=False,
    )

    v_values = np.linspace(
        0.0,
        2.0 * np.pi,
        minor_resolution,
        endpoint=False,
    )

    points: list[list[float]] = []

    for u in u_values:
        cos_u = np.cos(u)
        sin_u = np.sin(u)

        for v in v_values:
            cos_v = np.cos(v)
            sin_v = np.sin(v)

            radial = (
                major_radius
                + minor_radius * cos_v
            )

            x = radial * cos_u
            y = radial * sin_u
            z = minor_radius * sin_v

            points.append([x, y, z])

    xyz = np.asarray(points, dtype=float)

    xyzr = _make_xyzr(
        xyz=xyz,
        sphere_radius=sphere_radius,
    )

    return ValidationShape(
        name="torus",
        xyzr=xyzr,
        expected_int_mean_curvature=(
            2.0 * np.pi**2 * major_radius
        ),
        expected_int_gaussian_curvature=0.0,
        parameters={
            "major_radius": major_radius,
            "minor_radius": minor_radius,
            "major_resolution": major_resolution,
            "minor_resolution": minor_resolution,
            "sphere_radius": sphere_radius,
        },
        notes=(
            "Standard ring torus. Gaussian curvature is positive on the "
            "outer region and negative on the inner region, but integrates "
            "to zero by Gauss-Bonnet because the surface has genus 1."
        ),
    )


def torus_multicell(
    major_radius: float,
    minor_radius: float,
    interior_count: int = 12,
    cross_section_resolution: int = 12,
    sphere_radius: float = 1.0,
    center_radius: float = 1.0,
    neighbor_radius: float | None = None,
    include_center: bool = True,
) -> ValidationShape:
    """Generate a multi-cell power-diagram validation ring.

    The first ``interior_count`` generators form a closed ring on the major
    circle and are tagged as the selectable interior group.  Each receives a
    modest set of weighted constraints in its local radial/vertical
    cross-section.  Their power bisectors are at ``minor_radius`` from the
    ring centerline, leaving the global origin (the torus hole) empty.

    This is deliberately a coarse solid-torus construction.  The selected
    cells, rather than the constraint shell, define the object to be tested;
    increasing ``interior_count`` and ``cross_section_resolution`` improves
    the polygonal approximation without changing the topology strategy.
    """
    major_radius = _validate_positive(major_radius, "major_radius")
    minor_radius = _validate_positive(minor_radius, "minor_radius")
    if major_radius <= minor_radius:
        raise ValueError("A ring torus requires major_radius > minor_radius.")
    interior_count = _validate_resolution(interior_count, minimum=3)
    cross_section_resolution = _validate_resolution(
        cross_section_resolution, minimum=4,
    )
    center_radius = _validate_positive(center_radius, "center_radius")
    if neighbor_radius is None:
        neighbor_radius = sphere_radius
    neighbor_radius = _validate_positive(neighbor_radius, "neighbor_radius")

    theta = 2.0 * np.pi * np.arange(cross_section_resolution) / cross_section_resolution
    phi = 2.0 * np.pi * np.arange(interior_count) / interior_count
    interior_xyz = np.column_stack((
        major_radius * np.cos(phi),
        major_radius * np.sin(phi),
        np.zeros(interior_count),
    ))

    # Every local support plane is h=minor_radius from its ring centerline.
    radicand = (
        minor_radius * minor_radius
        - center_radius * center_radius
        + neighbor_radius * neighbor_radius
    )
    if radicand < 0.0:
        raise ValueError(
            "The requested torus support planes are incompatible with the "
            "supplied center and neighbor radii."
        )
    constraint_distance = minor_radius + np.sqrt(radicand)
    exterior_xyz = []
    assignments = []
    for interior_index, angle in enumerate(phi):
        e_r = np.array([np.cos(angle), np.sin(angle), 0.0])
        e_z = np.array([0.0, 0.0, 1.0])
        # A small alternating phase avoids exact coincident three-plane
        # intersections between neighboring cells while retaining uniform
        # cross-sectional support at the requested tube radius.
        phase = (interior_index % 2) * (0.35 * np.pi / cross_section_resolution)
        local_theta = theta + phase
        normals = np.cos(local_theta)[:, None] * e_r + np.sin(local_theta)[:, None] * e_z
        points = interior_xyz[interior_index] + constraint_distance * normals
        exterior_xyz.extend(points.tolist())
        assignments.extend([interior_index] * len(points))

    interior_xyzr = _make_xyzr(interior_xyz, center_radius)
    exterior_xyzr = _make_xyzr(np.asarray(exterior_xyz), neighbor_radius)
    xyzr = np.vstack((interior_xyzr, exterior_xyzr)) if include_center else exterior_xyzr
    interior_indices = list(range(interior_count)) if include_center else []
    exterior_start = len(interior_xyzr) if include_center else 0

    return ValidationShape(
        name="torus_multicell",
        xyzr=xyzr,
        expected_int_mean_curvature=2.0 * np.pi**2 * major_radius,
        expected_int_gaussian_curvature=0.0,
        parameters={
            "major_radius": major_radius,
            "minor_radius": minor_radius,
            "interior_count": interior_count,
            "cross_section_resolution": cross_section_resolution,
            "sphere_radius": neighbor_radius,
            "center_radius": center_radius,
            "neighbor_radius": neighbor_radius,
            "include_center": include_center,
            "interior_indices": interior_indices,
            "interior_coordinates": interior_xyz,
            "exterior_start": exterior_start,
            "exterior_count": len(exterior_xyz),
            "constraint_distance": constraint_distance,
            "constraint_assignments": assignments,
            "target_volume": 2.0 * np.pi**2 * major_radius * minor_radius**2,
            "target_surface_area": 4.0 * np.pi**2 * major_radius * minor_radius,
        },
        notes=(
            "Interior ring generators form the selected group. Local radial "
            "and vertical weighted constraints bound their tube cells while "
            "preserving the central torus hole."
        ),
    )


# ============================================================================
# Convenience utilities
# ============================================================================


def analytic_summary(shape: ValidationShape) -> str:
    """Return a compact human-readable summary of a validation geometry."""
    return (
        f"{shape.name}\n"
        f"  XYZR spheres: {shape.n_atoms}\n"
        f"  Expected integrated mean curvature: "
        f"{shape.expected_int_mean_curvature:.12g}\n"
        f"  Expected integrated Gaussian curvature: "
        f"{shape.expected_int_gaussian_curvature:.12g}\n"
        f"  Parameters: {shape.parameters}"
    )



# ============================================================================
# Command-line interface
# ============================================================================


def _build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line parser for synthetic validation geometries.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate synthetic XYZR systems with analytically known "
            "integrated mean and Gaussian curvatures."
        )
    )

    # ------------------------------------------------------------------
    # Shape selection
    # ------------------------------------------------------------------

    parser.add_argument(
        "-s",
        "--shape",
        required=True,
        choices=[
            "sphere",
            "cube",
            "box",
            "spherocylinder",
            "spherocylinder_multicell",
            "torus",
            "torus_multicell",
        ],
        help="Synthetic geometry to generate.",
    )

    # ------------------------------------------------------------------
    # Resolution
    # ------------------------------------------------------------------

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=10,
        help=(
            "Primary geometry resolution. Interpretation depends on shape. "
            "Default: 10."
        ),
    )

    # ------------------------------------------------------------------
    # General sphere parameters
    # ------------------------------------------------------------------

    parser.add_argument(
        "--ball-radius",
        type=float,
        default=1.0,
        help=(
            "Radius assigned to each surrounding validation ball. "
            "Default: 1.0."
        ),
    )

    parser.add_argument(
        "--center-radius",
        type=float,
        default=1.0,
        help=(
            "Radius assigned to the central VorPy validation ball. "
            "Default: 1.0."
        ),
    )

    parser.add_argument(
        "--sphere-radius",
        type=float,
        default=1.0,
        help=(
            "Radius assigned to validation balls for cube, box, "
            "spherocylinder, and torus geometries. Default: 1.0."
        ),
    )

    parser.add_argument(
        "--neighbor-radius",
        type=float,
        default=None,
        help=(
            "Explicit radius for cube/box/spherocylinder neighboring "
            "generators. Defaults to the shape's sphere/ball radius."
        ),
    )

    parser.add_argument(
        "--no-center",
        action="store_true",
        help=(
            "Do not add the central validation ball. "
            "By default a central ball is included."
        ),
    )

    # ------------------------------------------------------------------
    # Sphere
    # ------------------------------------------------------------------

    parser.add_argument(
        "--radius",
        type=float,
        default=10.0,
        help=(
            "Target central power-cell radius for sphere; target radius "
            "for spherocylinder. "
            "Default: 10.0."
        ),
    )

    # ------------------------------------------------------------------
    # Cube / box
    # ------------------------------------------------------------------

    parser.add_argument(
        "--side-length",
        type=float,
        default=10.0,
        help="Cube side length. Default: 10.0.",
    )

    parser.add_argument(
        "--length",
        type=float,
        default=10.0,
        help="Box x dimension. Default: 10.0.",
    )

    parser.add_argument(
        "--width",
        type=float,
        default=10.0,
        help="Box y dimension. Default: 10.0.",
    )

    parser.add_argument(
        "--height",
        type=float,
        default=10.0,
        help="Box z dimension. Default: 10.0.",
    )

    # ------------------------------------------------------------------
    # Spherocylinder
    # ------------------------------------------------------------------

    parser.add_argument(
        "--cylinder-length",
        type=float,
        default=10.0,
        help=(
            "Length of the cylindrical section of a spherocylinder. "
            "Default: 10.0."
        ),
    )

    parser.add_argument(
        "--angular-resolution",
        type=int,
        default=None,
        help=(
            "Deprecated compatibility option; spherocylinder resolution "
            "now directly counts supporting generators."
        ),
    )

    parser.add_argument(
        "--interior-count",
        type=int,
        default=3,
        help="Number of axial interior cells for spherocylinder_multicell.",
    )

    parser.add_argument(
        "--cap-resolution",
        type=int,
        default=12,
        help="Number of constraints per cap for spherocylinder_multicell.",
    )

    # ------------------------------------------------------------------
    # Torus
    # ------------------------------------------------------------------

    parser.add_argument(
        "--major-radius",
        type=float,
        default=10.0,
        help="Torus major radius R. Default: 10.0.",
    )

    parser.add_argument(
        "--minor-radius",
        type=float,
        default=3.0,
        help="Torus minor radius r. Default: 3.0.",
    )

    parser.add_argument(
        "--minor-resolution",
        type=int,
        default=None,
        help=(
            "Torus minor-circle resolution. "
            "If omitted, --resolution is used."
        ),
    )

    parser.add_argument(
        "--cross-section-resolution",
        type=int,
        default=12,
        help=(
            "Number of local radial/vertical constraints per interior cell "
            "for torus_multicell. Default: 12."
        ),
    )

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    parser.add_argument(
        "-e",
        "--export",
        nargs="+",
        default=None,
        help=(
            "Export settings. Use '-e dir PATH' to specify the "
            "output directory."
        ),
    )

    parser.add_argument(
        "--name",
        type=str,
        default=None,
        help=(
            "Optional output basename. If omitted, a descriptive "
            "name is generated automatically."
        ),
    )

    return parser


def _get_output_directory(args) -> Path:
    """
    Interpret VorPy-style '-e dir PATH' output arguments.

    If no export directory is specified, validation systems are written to:

        ./data/validation_shapes
    """
    if args.export is None:
        return Path("data") / "validation_shapes"

    if len(args.export) >= 2 and args.export[0].lower() == "dir":
        return Path(" ".join(args.export[1:]))

    raise ValueError(
        "Invalid export syntax. Use:\n\n"
        "    -e dir PATH\n"
    )


def _generate_from_args(args) -> ValidationShape:
    """
    Generate the requested validation geometry from CLI arguments.
    """
    shape_name = args.shape.lower()

    if shape_name == "sphere":
        return sphere(
            radius=args.radius,
            resolution=args.resolution,
            ball_radius=args.ball_radius,
            center_radius=args.center_radius,
            include_center=not args.no_center,
        )

    if shape_name == "cube":
        return cube(
            side_length=args.side_length,
            resolution=args.resolution,
            sphere_radius=args.sphere_radius,
            center_radius=args.center_radius,
            neighbor_radius=args.neighbor_radius,
            include_center=not args.no_center,
        )

    if shape_name == "box":
        return box(
            length=args.length,
            width=args.width,
            height=args.height,
            resolution=args.resolution,
            sphere_radius=args.sphere_radius,
            center_radius=args.center_radius,
            neighbor_radius=args.neighbor_radius,
            include_center=not args.no_center,
        )

    if shape_name == "spherocylinder":
        return spherocylinder(
            radius=args.radius,
            cylinder_length=args.cylinder_length,
            resolution=args.resolution,
            sphere_radius=args.sphere_radius,
            center_radius=args.center_radius,
            neighbor_radius=args.neighbor_radius,
            include_center=not args.no_center,
        )

    if shape_name == "spherocylinder_multicell":
        angular_resolution = (
            args.angular_resolution
            if args.angular_resolution is not None
            else 12
        )
        return spherocylinder_multicell(
            radius=args.radius,
            cylinder_length=args.cylinder_length,
            interior_count=args.interior_count,
            angular_resolution=angular_resolution,
            cap_resolution=args.cap_resolution,
            sphere_radius=args.sphere_radius,
            center_radius=args.center_radius,
            neighbor_radius=args.neighbor_radius,
            include_center=not args.no_center,
        )

    if shape_name == "torus":
        minor_resolution = (
            args.minor_resolution
            if args.minor_resolution is not None
            else args.resolution
        )

        return torus(
            major_radius=args.major_radius,
            minor_radius=args.minor_radius,
            major_resolution=args.resolution,
            minor_resolution=minor_resolution,
            sphere_radius=args.sphere_radius,
        )

    if shape_name == "torus_multicell":
        return torus_multicell(
            major_radius=args.major_radius,
            minor_radius=args.minor_radius,
            interior_count=args.interior_count,
            cross_section_resolution=args.cross_section_resolution,
            sphere_radius=args.sphere_radius,
            center_radius=args.center_radius,
            neighbor_radius=args.neighbor_radius,
            include_center=not args.no_center,
        )

    raise ValueError(
        f"Unsupported validation shape: {shape_name}"
    )


def _default_basename(shape: ValidationShape) -> str:
    """
    Generate a descriptive output basename.
    """
    p = shape.parameters

    if shape.name == "sphere":
        return (
            f"sphere_R{p['radius']:g}"
            f"_br{p['ball_radius']:g}"
            f"_cr{p['center_radius']:g}"
            f"_n{p['resolution']}"
        )

    if shape.name == "cube":
        return (
            f"cube_L{p['side_length']:g}"
            f"_r{p['sphere_radius']:g}"
            f"_n{p['resolution']}"
        )

    if shape.name == "box":
        return (
            f"box_{p['length']:g}x"
            f"{p['width']:g}x"
            f"{p['height']:g}"
            f"_r{p['sphere_radius']:g}"
            f"_n{p['resolution']}"
        )

    if shape.name == "spherocylinder_multicell":
        return (
            f"spherocylinder_multicell_R{p['radius']:g}"
            f"_L{p['cylinder_length']:g}"
            f"_r{p['sphere_radius']:g}"
            f"_i{p['interior_count']}"
            f"_a{p['angular_resolution']}"
            f"_c{p['cap_resolution']}"
            f"_z{p['axial_ring_count']}"
        )

    if shape.name == "spherocylinder_multicell":
        return (
            f"spherocylinder_multicell_R{p['radius']:g}"
            f"_L{p['cylinder_length']:g}"
            f"_r{p['sphere_radius']:g}"
            f"_n{p['interior_count']}"
        )

    if shape.name == "torus":
        return (
            f"torus_R{p['major_radius']:g}"
            f"_r{p['minor_radius']:g}"
            f"_ball{p['sphere_radius']:g}"
            f"_n{p['major_resolution']}"
            f"x{p['minor_resolution']}"
        )

    if shape.name == "torus_multicell":
        return (
            f"torus_multicell_R{p['major_radius']:g}"
            f"_r{p['minor_radius']:g}"
            f"_i{p['interior_count']}"
            f"_x{p['cross_section_resolution']}"
        )

    return shape.name


def main() -> None:
    """
    Command-line entry point.
    """
    parser = _build_parser()
    args = parser.parse_args()

    try:
        shape = _generate_from_args(args)
        output_directory = _get_output_directory(args)

    except ValueError as exc:
        parser.error(str(exc))
        return

    basename = (
        args.name
        if args.name is not None
        else _default_basename(shape)
    )

    xyzr_path, pdb_path, pymol_path = shape.save(
        directory=output_directory,
        basename=basename,
    )

    print()
    print("=" * 72)
    print("VORPY CURVATURE VALIDATION SHAPE")
    print("=" * 72)

    print(analytic_summary(shape))

    print()
    print("Output")
    print(f"  Directory: {xyzr_path.parent}")
    print(f"  XYZR:      {xyzr_path.name}")
    print(f"  PDB:       {pdb_path.name}")
    print(f"  PyMOL:     {pymol_path.name}")

    print()
    print("PyMOL")
    print(f"  load {pdb_path}")
    print("  hide everything")
    print("  alter all, vdw=b")
    print("  rebuild")
    print("  show spheres")

    print("=" * 72)


if __name__ == "__main__":
    main()
