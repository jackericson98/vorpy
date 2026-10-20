"""
Curvature-dependent surface geometry and energy calculations.

The primary geometric descriptors are

    A = integral(dA)
    C = integral(H dA)
    Q = integral(H^2 dA)
    X = integral(K dA)

where H is mean curvature and K is Gaussian curvature.

When VorPy coordinates are in Angstroms:

    A : A^2
    H : A^-1
    K : A^-2
    C : A
    Q : dimensionless
    X : dimensionless

Energy conversion is kept separate from the geometric calculation so
that thermodynamic coefficients can be documented and parameterized
independently.
"""

import numpy as np

from vorpy.src.calculations.calcs import calc_tri
from vorpy.src.calculations.curvature import mean_curvature, gaussian_curvature


_CURVATURE_UNAVAILABLE = {
    "NOT_CALCULATED", "NOT_SUPPORTED", "UNRESOLVED", "STALE", "UNAVAILABLE",
}


def _status_value(value):
    value = getattr(value, "value", value)
    return str(value).upper() if value is not None else None


def _geometry_status(geometry):
    status = geometry.get("Curvature Status", geometry.get("curvature_status"))
    return _status_value(status) or "CERTIFIED"


def _finite_value(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _energy_status(geometry, required):
    status = _geometry_status(geometry)
    if status in _CURVATURE_UNAVAILABLE or status not in {"CERTIFIED", "CALCULATED"}:
        return status
    if any(_finite_value(geometry.get(name)) is None for name in required):
        return "NOT_CALCULATED"
    return "CERTIFIED"


def calc_surface_energy_geometry_from_curvatures(
        points,
        tris,
        mean_curvatures,
        gaussian_curvatures,
        area=None):
    """
    Integrate already-calculated triangle curvatures over one surface.

    This is the preferred function during VorPy network construction
    because build_surf() already calculates H and K at every triangle
    centroid. Reusing those values avoids recalculating curvature.

    Parameters
    ----------
    points : sequence of array-like
        Surface vertex coordinates.
    tris : sequence
        Triangle indices into points.
    mean_curvatures : sequence of float
        Mean curvature H for each triangle centroid.
    gaussian_curvatures : sequence of float
        Gaussian curvature K for each triangle centroid.
    area : float, optional
        Precalculated total surface area. If omitted, it is accumulated
        from the triangle areas.

    Returns
    -------
    dict
        Area, area-weighted average curvatures, and curvature integrals.
    """
    if len(tris) != len(mean_curvatures):
        raise ValueError("tris and mean_curvatures must have the same length.")

    if len(tris) != len(gaussian_curvatures):
        raise ValueError("tris and gaussian_curvatures must have the same length.")

    calculated_area = 0.0
    integrated_mean_curvature = 0.0
    integrated_mean_curvature_squared = 0.0
    integrated_gaussian_curvature = 0.0
    curvature_status = "CERTIFIED" if tris else "NOT_CALCULATED"

    for i, tri in enumerate(tris):
        p0 = np.asarray(points[tri[0]], dtype=float)
        p1 = np.asarray(points[tri[1]], dtype=float)
        p2 = np.asarray(points[tri[2]], dtype=float)

        tri_area = float(calc_tri(p0, p1, p2))

        if not np.isfinite(tri_area) or tri_area <= 0.0:
            curvature_status = "UNRESOLVED"
            continue

        H = _finite_value(mean_curvatures[i])
        K = _finite_value(gaussian_curvatures[i])
        if H is None or K is None:
            curvature_status = "UNRESOLVED"
            continue

        calculated_area += tri_area
        integrated_mean_curvature += H * tri_area
        integrated_mean_curvature_squared += H * H * tri_area
        integrated_gaussian_curvature += K * tri_area

    total_area = calculated_area if area is None else float(area)

    if curvature_status == "CERTIFIED" and total_area > 0.0:
        average_mean_curvature = integrated_mean_curvature / total_area
        average_gaussian_curvature = integrated_gaussian_curvature / total_area
    elif curvature_status == "CERTIFIED":
        average_mean_curvature = 0.0
        average_gaussian_curvature = 0.0
    else:
        integrated_mean_curvature = None
        integrated_mean_curvature_squared = None
        integrated_gaussian_curvature = None
        average_mean_curvature = None
        average_gaussian_curvature = None

    return {
        "Area": total_area,
        "Mean Curvature": average_mean_curvature,
        "Gaussian Curvature": average_gaussian_curvature,
        "Integrated Mean Curvature": integrated_mean_curvature,
        "Integrated Mean Curvature Squared": integrated_mean_curvature_squared,
        "Integrated Gaussian Curvature": integrated_gaussian_curvature,
        "Curvature Status": curvature_status,
    }


def calc_surface_energy_geometry(func, points, tris):
    """
    Calculate surface geometry when triangle curvatures are not available.

    Curvature is evaluated analytically at each triangle centroid, then
    passed through the same integration routine used during network builds.
    """
    mean_curvatures = []
    gaussian_curvatures = []

    for tri in tris:
        p0 = np.asarray(points[tri[0]], dtype=float)
        p1 = np.asarray(points[tri[1]], dtype=float)
        p2 = np.asarray(points[tri[2]], dtype=float)
        centroid = (p0 + p1 + p2) / 3.0

        mean_curvatures.append(mean_curvature(func, centroid))
        gaussian_curvatures.append(gaussian_curvature(func, centroid))

    return calc_surface_energy_geometry_from_curvatures(
        points=points,
        tris=tris,
        mean_curvatures=mean_curvatures,
        gaussian_curvatures=gaussian_curvatures
    )


def sum_surface_energy_geometry(surface_values):
    """
    Combine geometric descriptors from multiple surfaces.
    """
    area = 0.0
    mean_values = []
    mean_squared_values = []
    gaussian_values = []
    statuses = []

    for values in surface_values:
        area += float(values.get("Area", 0.0))
        statuses.append(_geometry_status(values))
        mean_values.append(_finite_value(values.get("Integrated Mean Curvature")))
        mean_squared_values.append(_finite_value(values.get("Integrated Mean Curvature Squared")))
        gaussian_values.append(_finite_value(values.get("Integrated Gaussian Curvature")))

    def complete_sum(values):
        return sum(values) if values and all(value is not None for value in values) else None

    integrated_mean_curvature = complete_sum(mean_values)
    integrated_mean_curvature_squared = complete_sum(mean_squared_values)
    integrated_gaussian_curvature = complete_sum(gaussian_values)

    status = "CERTIFIED"
    for candidate in statuses:
        if candidate in _CURVATURE_UNAVAILABLE or candidate not in {"CERTIFIED", "CALCULATED"}:
            status = candidate
            break
    if any(value is None for value in (
        integrated_mean_curvature,
        integrated_mean_curvature_squared,
        integrated_gaussian_curvature,
    )):
        status = status if status != "CERTIFIED" else "NOT_CALCULATED"

    if integrated_mean_curvature is not None and area > 0.0:
        average_mean_curvature = integrated_mean_curvature / area
    else:
        average_mean_curvature = None
    if integrated_gaussian_curvature is not None and area > 0.0:
        average_gaussian_curvature = integrated_gaussian_curvature / area
    else:
        average_gaussian_curvature = None

    return {
        "Area": area,
        "Mean Curvature": average_mean_curvature,
        "Gaussian Curvature": average_gaussian_curvature,
        "Integrated Mean Curvature": integrated_mean_curvature,
        "Integrated Mean Curvature Squared": integrated_mean_curvature_squared,
        "Integrated Gaussian Curvature": integrated_gaussian_curvature,
        "Curvature Status": status,
    }


def calc_morphometric_energy(
        geometry,
        surface_tension,
        mean_curvature_coefficient,
        gaussian_curvature_coefficient):
    """
    Calculate

        G = gamma*A + kappa*C + kappa_bar*X.
    """
    required = (
        "Area",
        "Integrated Mean Curvature",
        "Integrated Gaussian Curvature",
    )
    status = _energy_status(geometry, required)
    area = _finite_value(geometry.get("Area"))
    area_energy = surface_tension * area if area is not None else None
    if status != "CERTIFIED":
        return {
            "Area Energy": area_energy,
            "Mean Curvature Energy": None,
            "Gaussian Curvature Energy": None,
            "Total Energy": None,
            "Status": status,
        }
    mean_energy = (
        mean_curvature_coefficient * geometry["Integrated Mean Curvature"]
    )
    gaussian_energy = (
        gaussian_curvature_coefficient * geometry["Integrated Gaussian Curvature"]
    )

    return {
        "Area Energy": area_energy,
        "Mean Curvature Energy": mean_energy,
        "Gaussian Curvature Energy": gaussian_energy,
        "Total Energy": area_energy + mean_energy + gaussian_energy,
        "Status": "CERTIFIED",
    }


def calc_helfrich_energy(
        geometry,
        surface_tension,
        bending_modulus,
        gaussian_modulus,
        spontaneous_curvature=0.0):
    """
    Calculate the Helfrich-style surface energy

        G = gamma*A
            + (kappa_b/2) * integral((2H-C0)^2 dA)
            + kappa_G*integral(K dA).
    """
    required = (
        "Area",
        "Integrated Mean Curvature",
        "Integrated Mean Curvature Squared",
        "Integrated Gaussian Curvature",
    )
    status = _energy_status(geometry, required)
    area = _finite_value(geometry.get("Area"))
    int_mean = _finite_value(geometry.get("Integrated Mean Curvature"))
    int_mean_sq = _finite_value(geometry.get("Integrated Mean Curvature Squared"))
    int_gaussian = _finite_value(geometry.get("Integrated Gaussian Curvature"))

    if status != "CERTIFIED" or any(value is None for value in (area, int_mean, int_mean_sq, int_gaussian)):
        return {
            "Area Energy": surface_tension * area if area is not None else None,
            "Bending Integral": None,
            "Bending Energy": None,
            "Gaussian Curvature Energy": None,
            "Total Energy": None,
            "Status": status if status != "CERTIFIED" else "NOT_CALCULATED",
        }

    c0 = float(spontaneous_curvature)

    area_energy = surface_tension * area
    bending_integral = (
        4.0 * int_mean_sq
        - 4.0 * c0 * int_mean
        + c0 * c0 * area
    )
    bending_energy = 0.5 * bending_modulus * bending_integral
    gaussian_energy = gaussian_modulus * int_gaussian

    return {
        "Area Energy": area_energy,
        "Bending Integral": bending_integral,
        "Bending Energy": bending_energy,
        "Gaussian Curvature Energy": gaussian_energy,
        "Total Energy": area_energy + bending_energy + gaussian_energy,
        "Status": "CERTIFIED",
    }


def calc_representative_surface_energy(
        geometry,
        bending_modulus=1.0,
        spontaneous_curvature=0.0):
    """
    Calculate a representative curvature-dependent surface energy.

    This uses the mean-curvature bending term of the Helfrich model:

        E = (kappa_b / 2) * integral((2H - C0)^2 dA)

    By default:

        kappa_b = 1 kBT
        C0 = 0

    so:

        E / kBT = 2 * integral(H^2 dA)

    Parameters
    ----------
    geometry : dict
        Surface geometry containing:
            Integrated Mean Curvature
            Integrated Mean Curvature Squared
            Area

    bending_modulus : float, optional
        Bending modulus in units of kBT.
        Default = 1.0.

    spontaneous_curvature : float, optional
        Spontaneous curvature C0 in Angstrom^-1.
        Default = 0.0.

    Returns
    -------
    float
        Representative bending energy in units of kBT.
    """

    status = _energy_status(
        geometry,
        ('Area', 'Integrated Mean Curvature', 'Integrated Mean Curvature Squared'),
    )
    if status != "CERTIFIED":
        return None

    area = geometry['Area']
    int_mean_curv = geometry['Integrated Mean Curvature']
    int_mean_curv_sq = geometry['Integrated Mean Curvature Squared']

    c0 = float(spontaneous_curvature)

    bending_integral = (
        4.0 * int_mean_curv_sq
        - 4.0 * c0 * int_mean_curv
        + c0 ** 2 * area
    )

    energy = (
        0.5
        * bending_modulus
        * bending_integral
    )

    return energy
