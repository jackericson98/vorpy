"""Common implicit-surface geometry interfaces for VorPy.

This module is intentionally additive.  Existing curvature construction still
uses :mod:`vorpy.src.calculations.curvature`, so introducing these objects does
not change any current network values.  Later boundary-edge calculations can
depend on this stable interface instead of branching on the network scheme.
"""

from abc import ABC, abstractmethod

import numpy as np


class SurfaceGeometry(ABC):
    """Interface exposed by a pairwise VorPy boundary surface."""

    @abstractmethod
    def implicit_value(self, point):
        """Return ``F(point)`` for the implicit surface ``F = 0``."""

    @abstractmethod
    def gradient(self, point):
        """Return the (not normalized) gradient of ``F`` at *point*."""

    def normal(self, point, orientation=1.0, tol=1e-12):
        """Return a unit normal, optionally reversed by ``orientation``.

        ``orientation`` must be positive for the native implicit-function
        orientation or negative for its reverse.  Group-relative outward
        orientation is deliberately selected by the boundary-classification
        layer, not inferred here.
        """
        gradient = np.asarray(self.gradient(point), dtype=float)
        magnitude = float(np.linalg.norm(gradient))
        if not np.isfinite(magnitude) or magnitude < tol:
            raise ValueError("Surface normal is undefined at this point.")
        sign = 1.0 if float(orientation) >= 0.0 else -1.0
        return sign * gradient / magnitude

    @abstractmethod
    def mean_curvature(self, point, tol=1e-12):
        """Return mean curvature at *point*."""

    @abstractmethod
    def gaussian_curvature(self, point, tol=1e-12):
        """Return Gaussian curvature at *point*."""


class QuadraticSurfaceGeometry(SurfaceGeometry):
    """Implicit quadratic surface using VorPy's 14-coefficient convention.

    The coefficient sequence is ``A, B, C, D, E, F, G, H, I, J, K,
    dx, dy, dz``.  The last three values are retained as construction metadata
    but are not part of the implicit polynomial, matching the existing
    curvature implementation.
    """

    def __init__(self, coefficients):
        values = np.asarray(coefficients, dtype=float)
        if values.shape != (14,):
            raise ValueError("Quadratic surfaces require exactly 14 coefficients.")
        self.coefficients = values

    def implicit_value(self, point):
        A, B, C, D, E, F, G, H, I, J, K, _, _, _ = self.coefficients
        x, y, z = np.asarray(point, dtype=float)
        return float(
            A*x*x + B*y*y + C*z*z + D*x*y + E*y*z + F*x*z
            + G*x + H*y + I*z + J
        )

    def gradient(self, point):
        A, B, C, D, E, F, G, H, I, _, _, _, _, _ = self.coefficients
        x, y, z = np.asarray(point, dtype=float)
        return np.array([
            2.0*A*x + D*y + F*z + G,
            2.0*B*y + D*x + E*z + H,
            2.0*C*z + F*x + E*y + I,
        ], dtype=float)

    def hessian(self):
        A, B, C, D, E, F = self.coefficients[:6]
        return np.array([
            [2.0*A, D, F],
            [D, 2.0*B, E],
            [F, E, 2.0*C],
        ], dtype=float)

    def mean_curvature(self, point, tol=1e-12):
        gradient = self.gradient(point)
        grad_mag = float(np.linalg.norm(gradient))
        if not np.isfinite(grad_mag) or grad_mag < tol:
            return 0.0
        hessian = self.hessian()
        numerator = (
            float(np.trace(hessian))*grad_mag**2
            - float(gradient @ hessian @ gradient)
        )
        denominator = 2.0*grad_mag**3
        value = numerator/denominator
        return float(value) if np.isfinite(value) else 0.0

    def gaussian_curvature(self, point, tol=1e-12):
        gradient = self.gradient(point)
        grad_mag = float(np.linalg.norm(gradient))
        if not np.isfinite(grad_mag) or grad_mag < tol:
            return 0.0
        hessian = self.hessian()
        # Cross-product construction is the adjugate for a 3x3 matrix and
        # remains valid when the Hessian is singular.
        adjugate = np.column_stack((
            np.cross(hessian[:, 1], hessian[:, 2]),
            np.cross(hessian[:, 2], hessian[:, 0]),
            np.cross(hessian[:, 0], hessian[:, 1]),
        ))
        value = float(gradient @ adjugate @ gradient)/(grad_mag**4)
        return value if np.isfinite(value) else 0.0


class PlaneSurfaceGeometry(SurfaceGeometry):
    """Implicit plane ``normal . (point - origin) = 0``."""

    def __init__(self, origin, normal):
        self.origin = np.asarray(origin, dtype=float)
        raw_normal = np.asarray(normal, dtype=float)
        magnitude = float(np.linalg.norm(raw_normal))
        if self.origin.shape != (3,) or raw_normal.shape != (3,):
            raise ValueError("Plane origin and normal must be three-vectors.")
        if not np.isfinite(magnitude) or magnitude < 1e-12:
            raise ValueError("Plane normal must be finite and nonzero.")
        self._normal = raw_normal/magnitude

    def implicit_value(self, point):
        return float(self._normal @ (np.asarray(point, dtype=float) - self.origin))

    def gradient(self, point):
        return self._normal.copy()

    def mean_curvature(self, point, tol=1e-12):
        return 0.0

    def gaussian_curvature(self, point, tol=1e-12):
        return 0.0


def surface_geometry_from_coefficients(coefficients):
    """Create the current AW quadratic geometry without changing callers."""
    return QuadraticSurfaceGeometry(coefficients)

def aw_clearance_difference(
        point,
        cell_location,
        cell_radius,
        neighbor_location,
        neighbor_radius):
    """Return additive-clearance difference d_cell - d_neighbor.

    Negative values lie on the cell side of the pairwise AW boundary,
    positive values lie on the neighboring-cell side, and zero defines
    the pairwise boundary surface.
    """
    point = np.asarray(point, dtype=float)
    cell_location = np.asarray(cell_location, dtype=float)
    neighbor_location = np.asarray(neighbor_location, dtype=float)

    if (
        point.shape != (3,)
        or cell_location.shape != (3,)
        or neighbor_location.shape != (3,)
    ):
        raise ValueError(
            "AW clearance inputs must be three-dimensional vectors."
        )

    if not (
        np.all(np.isfinite(point))
        and np.all(np.isfinite(cell_location))
        and np.all(np.isfinite(neighbor_location))
        and np.isfinite(float(cell_radius))
        and np.isfinite(float(neighbor_radius))
    ):
        raise ValueError("AW clearance inputs must be finite.")

    cell_clearance = (
        np.linalg.norm(point - cell_location)
        - float(cell_radius)
    )
    neighbor_clearance = (
        np.linalg.norm(point - neighbor_location)
        - float(neighbor_radius)
    )

    return float(cell_clearance - neighbor_clearance)


def aw_clearance_gradient(
        point,
        cell_location,
        neighbor_location,
        tol=1e-12):
    """Return grad(d_cell - d_neighbor).

    The gradient is oriented from the cell side toward the neighboring
    cell side of the pairwise additively weighted boundary.
    """
    point = np.asarray(point, dtype=float)
    cell_location = np.asarray(cell_location, dtype=float)
    neighbor_location = np.asarray(neighbor_location, dtype=float)

    if (
        point.shape != (3,)
        or cell_location.shape != (3,)
        or neighbor_location.shape != (3,)
    ):
        raise ValueError(
            "AW clearance-gradient inputs must be three-dimensional vectors."
        )

    if not (
        np.all(np.isfinite(point))
        and np.all(np.isfinite(cell_location))
        and np.all(np.isfinite(neighbor_location))
    ):
        raise ValueError("AW clearance-gradient inputs must be finite.")

    cell_vector = point - cell_location
    neighbor_vector = point - neighbor_location

    cell_distance = float(np.linalg.norm(cell_vector))
    neighbor_distance = float(np.linalg.norm(neighbor_vector))

    if cell_distance < tol or neighbor_distance < tol:
        raise ValueError(
            "AW clearance gradient is undefined at a generator center."
        )

    gradient = (
        cell_vector / cell_distance
        - neighbor_vector / neighbor_distance
    )

    if not np.all(np.isfinite(gradient)):
        raise ValueError("AW clearance gradient is non-finite.")

    return gradient


def aw_outward_normal(
        point,
        cell_location,
        neighbor_location,
        tol=1e-12):
    """Return the unit normal pointing outward from one AW cell.

    The returned direction points from the region satisfying

        d_cell <= d_neighbor

    toward the neighboring-cell side

        d_cell > d_neighbor.

    This orientation is independent of whatever generator ordering is
    required internally for surface construction.
    """
    gradient = aw_clearance_gradient(
        point,
        cell_location,
        neighbor_location,
        tol=tol,
    )

    magnitude = float(np.linalg.norm(gradient))

    if not np.isfinite(magnitude) or magnitude < tol:
        raise ValueError(
            "AW outward normal is undefined at this point."
        )

    return gradient / magnitude
