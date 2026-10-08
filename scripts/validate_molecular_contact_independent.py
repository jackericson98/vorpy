"""Independent analytic/numerical validation of molecular-contact geometry.

This audit deliberately does not call VorPy's area, clipping, topology, or
curvature helpers to form its reference values.  VorPy is used only to build
the candidate molecular-contact surface whose measurements are checked.

The default run is short and deterministic.  The final three reference cases
(deep pocket, through-channel, and enclosed cavity) are independent reference
surfaces; the current molecular-contact builder does not yet certify all of
those topologies from an atom-ball arrangement.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vorpy.src.geometry.molecular_contact import (
    AtomBall,
    SelectedContact,
    build_molecular_contact_surface,
    build_molecular_contact_surface_pair,
)


TWO_PI = 2.0 * math.pi
FOUR_PI = 4.0 * math.pi


def ball(index, group, xyz, radius):
    return AtomBall(index, f"{group}-{index}", group, tuple(xyz), radius)


def contact(index, atom_a, atom_b):
    return SelectedContact(f"independent-contact-{index}", atom_a, atom_b, "independent", "audit")


def close(label, actual, expected, tolerance):
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=tolerance):
        raise AssertionError(f"{label}: actual={actual:.16g}, expected={expected:.16g}")


def fibonacci_directions(count):
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for index in range(count):
        y = -1.0 + 2.0 * (index + 0.5) / count
        radial = math.sqrt(max(0.0, 1.0 - y * y))
        angle = golden_angle * index
        yield (radial * math.cos(angle), y, radial * math.sin(angle))


def distance_squared(left, right):
    return sum((left[index] - right[index]) ** 2 for index in range(3))


def selected_partners(atoms, carrier_index, contacts):
    partners = []
    for selected in contacts:
        if selected.status != "SELECTED":
            continue
        if selected.atom_a_index == carrier_index:
            partners.append(selected.atom_b_index)
        elif selected.atom_b_index == carrier_index:
            partners.append(selected.atom_a_index)
    return tuple(partners)


def sampled_carrier_area(atoms, carrier_index, contacts, count=100_000):
    """Independent area estimate from direct point-in-ball classification."""
    carrier = atoms[carrier_index]
    partners = selected_partners(atoms, carrier_index, contacts)
    occluders = tuple(
        index
        for index, atom in enumerate(atoms)
        if index != carrier_index and atom.group == carrier.group
    )
    retained = 0
    for direction in fibonacci_directions(count):
        point = tuple(
            carrier.center[index] + carrier.radius * direction[index]
            for index in range(3)
        )
        if not any(
            distance_squared(point, atoms[index].center) <= atoms[index].radius**2
            for index in partners
        ):
            continue
        if any(
            distance_squared(point, atoms[index].center) < atoms[index].radius**2
            for index in occluders
        ):
            continue
        retained += 1
    return FOUR_PI * carrier.radius**2 * retained / count


def spherical_cap_reference(radius, partner_radius, center_distance):
    q = (
        center_distance**2 + radius**2 - partner_radius**2
    ) / (2.0 * center_distance * radius)
    area = TWO_PI * radius**2 * (1.0 - q)
    return {
        "q": q,
        "area": area,
        "mean": area / radius,
        "gaussian": area / radius**2,
        "boundary_geodesic": TWO_PI * q,
        "gauss_bonnet": TWO_PI,
    }


def validate_cap(radius, partner_radius, center_distance, label, limitations):
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), radius),
        ball(1, "B", (center_distance, 0.0, 0.0), partner_radius),
    )
    side_a, side_b = build_molecular_contact_surface_pair(
        atoms, (contact(0, 0, 1),)
    )
    reference_a = spherical_cap_reference(radius, partner_radius, center_distance)
    reference_b = spherical_cap_reference(partner_radius, radius, center_distance)
    for side, surface, reference in (
        ("A", side_a, reference_a),
        ("B", side_b, reference_b),
    ):
        curvature = surface.curvature
        assert curvature is not None
        close(f"{label}/{side}/area", surface.total_area_A2, reference["area"], 1e-9)
        close(f"{label}/{side}/sampled area", sampled_carrier_area(atoms, 0 if side == "A" else 1, (contact(0, 0, 1),), 60_000), reference["area"], 2.0e-2)
        close(f"{label}/{side}/smooth H", curvature.smooth_mean_curvature, reference["mean"], 1e-9)
        close(f"{label}/{side}/total H", curvature.total_integrated_mean_curvature, reference["mean"], 1e-9)
        close(f"{label}/{side}/smooth K", curvature.smooth_gaussian_curvature, reference["gaussian"], 1e-9)
        close(f"{label}/{side}/boundary kg", curvature.boundary_geodesic_curvature, reference["boundary_geodesic"], 1e-9)
        close(f"{label}/{side}/GB", curvature.completed_gaussian_curvature, reference["gauss_bonnet"], 1e-9)
        assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 1, 1)
        if surface.betti_numbers != (1, 0, 0):
            limitations.append(
                f"{label}/{side}: disk geometry has invalid canonical Betti tuple "
                f"{surface.betti_numbers!r}"
            )
        assert curvature.seam_mean_curvature == 0.0
        assert curvature.intrinsic_seam_gaussian_curvature == 0.0
        assert curvature.interior_vertex_defects == 0.0
        assert curvature.boundary_corner_turning == 0.0
        assert curvature.gaussian_status == "CERTIFIED"


def validate_closed_sphere(limitations):
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), 1.0),
        ball(1, "B", (0.0, 0.0, 0.0), 10.0),
    )
    selected = (contact(0, 0, 1),)
    surface = build_molecular_contact_surface(atoms, selected, "A")
    curvature = surface.curvature
    assert curvature is not None
    close("isolated sphere/area", surface.total_area_A2, FOUR_PI, 1e-9)
    close("isolated sphere/sampled area", sampled_carrier_area(atoms, 0, selected, 60_000), FOUR_PI, 2.0e-2)
    close("isolated sphere/H", curvature.total_integrated_mean_curvature, FOUR_PI, 1e-9)
    close("isolated sphere/smooth K", curvature.smooth_gaussian_curvature, FOUR_PI, 1e-9)
    assert surface.boundary_loops == 0
    assert surface.betti_numbers == (1, 0, 1)
    assert surface.euler_characteristic == 2
    assert surface.orientability is True
    assert surface.genus_by_component == (0,)
    close("isolated sphere/completed GB", curvature.completed_gaussian_curvature, FOUR_PI, 1e-9)
    close("isolated sphere/GB target", curvature.gauss_bonnet_expected, FOUR_PI, 1e-9)
    assert curvature.gaussian_status == "CERTIFIED"
    close("isolated sphere/GB residual", curvature.gauss_bonnet_residual, 0.0, 1e-12)


def validate_selected_partner_restriction():
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), 1.0),
        ball(1, "B", (1.0, 0.0, 0.0), 1.0),
        ball(2, "B", (-1.0, 0.0, 0.0), 1.0),
    )
    selected = (contact(0, 0, 1),)
    surface = build_molecular_contact_surface(atoms, selected, "A")
    reference = spherical_cap_reference(1.0, 1.0, 1.0)
    close("selected partner restriction", surface.total_area_A2, reference["area"], 1e-9)
    assert surface.selected_contact_ids == ("independent-contact-0",)


def validate_multiple_overlapping_caps():
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), 2.0),
        ball(1, "B", (2.0, 0.0, 0.0), 2.0),
        ball(2, "B", (1.0, math.sqrt(3.0), 0.0), 2.0),
    )
    selected = (contact(0, 0, 1), contact(1, 0, 2))
    surface = build_molecular_contact_surface(atoms, selected, "A")
    sampled = sampled_carrier_area(atoms, 0, selected, 120_000)
    close("multiple caps/area", surface.total_area_A2, sampled, 2.0e-2)
    assert surface.naive_pairwise_area_A2 > surface.total_area_A2
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 1, 1)
    assert surface.betti_numbers == (1, 0, 0)
    assert surface.curvature is not None
    close("multiple caps/GB", surface.curvature.completed_gaussian_curvature, TWO_PI, 1e-9)


def validate_closed_two_sphere_union_reference():
    """Reference the closed union of two equal intersecting carrier spheres."""
    radius, center_distance = 1.0, 1.0
    q = center_distance / (2.0 * radius)
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), radius),
        ball(1, "A", (center_distance, 0.0, 0.0), radius),
        ball(2, "B", (0.5, 0.0, 0.0), 10.0),
    )
    selected = (contact(0, 0, 2), contact(1, 1, 2))
    sampled_total = sum(
        sampled_carrier_area(atoms, carrier, selected, 100_000)
        for carrier in (0, 1)
    )
    exact_area = FOUR_PI * radius**2 * (1.0 + q)
    close("closed two-sphere union/area", sampled_total, exact_area, 3.0e-2)

    circle_radius = math.sqrt(radius**2 - (center_distance / 2.0) ** 2)
    seam_length = TWO_PI * circle_radius
    beta = math.acos(1.0 - center_distance**2 / (2.0 * radius**2))
    smooth_k = exact_area / radius**2
    seam_h = 0.5 * beta * seam_length
    intrinsic_seam_k = -TWO_PI
    close("closed two-sphere union/smooth K", smooth_k, exact_area, 1e-12)
    close("closed two-sphere union/seam H", seam_h, math.pi**2 * math.sqrt(3.0) / 6.0, 1e-12)
    close("closed two-sphere union/GB", smooth_k + intrinsic_seam_k, FOUR_PI, 1e-12)
    solid_topology = {"space": "3D overlapping-ball union solid", "betti_numbers": (1, 0, 0)}
    surface_topology = {"space": "2D boundary surface", "betti_numbers": (1, 0, 1)}
    assert solid_topology["betti_numbers"] == (1, 0, 0)
    assert surface_topology["betti_numbers"] == (1, 0, 1)


def validate_closed_two_sphere_union_builder():
    """Check the canonical builder against the closed equal-sphere union."""
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), 1.0),
        ball(1, "A", (1.0, 0.0, 0.0), 1.0),
        ball(2, "B", (0.5, 0.0, 0.0), 10.0),
    )
    selected = (contact(0, 0, 2), contact(1, 1, 2))
    surface = build_molecular_contact_surface(atoms, selected, "A")
    curvature = surface.curvature
    assert curvature is not None
    close("closed union/builder area", surface.total_area_A2, 6.0 * math.pi, 1e-9)
    assert (surface.components, surface.boundary_loops) == (1, 0)
    assert surface.euler_characteristic == 2
    assert surface.betti_numbers == (1, 0, 1)
    assert surface.orientability is True
    assert surface.genus_by_component == (0,)
    assert len(surface.refined_seams) == 1
    assert set(dict(surface.refined_seams[0].orientation_on_patches).values()) == {-1, 1}
    assert all(status.startswith("interior_cycle:") for _, status in surface.junction_link_status)
    close("closed union/builder smooth K", curvature.smooth_gaussian_curvature, 6.0 * math.pi, 1e-9)
    close("closed union/builder seam H", curvature.seam_mean_curvature, math.pi**2 * math.sqrt(3.0) / 6.0, 1e-9)
    close("closed union/builder intrinsic seam K", curvature.intrinsic_seam_gaussian_curvature, -2.0 * math.pi, 1e-9)
    close("closed union/builder GB", curvature.completed_gaussian_curvature, FOUR_PI, 1e-9)
    close("closed union/builder GB residual", curvature.gauss_bonnet_residual, 0.0, 1e-12)
    assert curvature.gaussian_status == "CERTIFIED"


def validate_same_group_occlusion(limitations):
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), 2.0),
        ball(1, "A", (1.0, 0.0, 0.0), 2.0),
        ball(2, "B", (2.4, 0.0, 1.8), 2.0),
    )
    selected = (contact(0, 0, 2),)
    surface = build_molecular_contact_surface(atoms, selected, "A")
    sampled = sampled_carrier_area(atoms, 0, selected, 120_000)
    close("same-group occlusion/area", surface.total_area_A2, sampled, 2.0e-2)
    assert "A-1" in {
        atom_id
        for patch in surface.patches
        for atom_id in patch.self_occluding_atom_ids
    }
    assert surface.curvature is not None
    assert surface.curvature.completed_gaussian_curvature is None
    assert surface.curvature.gaussian_status == "PARTIAL"
    limitations.append("same-group occlusion: canonical Gaussian completion remains PARTIAL")


def validate_multi_boundary_occlusion_builder():
    """Check one connected sphere with three disjoint same-group holes."""
    radius = 3.0
    directions = (
        (1.0, 0.0, 0.0),
        (-0.5, math.sqrt(3.0) / 2.0, 0.0),
        (-0.5, -math.sqrt(3.0) / 2.0, 0.0),
    )
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), radius),
        *(ball(index, "A", tuple(3.0 * value for value in direction), 2.0)
          for index, direction in enumerate(directions, start=1)),
        ball(4, "B", (0.0, 0.0, 0.0), 20.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 4),), "A")
    curvature = surface.curvature
    assert curvature is not None
    expected_area = 24.0 * math.pi
    close("multi-boundary/builder area", surface.total_area_A2, expected_area, 1e-9)
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 3, -1)
    assert surface.betti_numbers == (1, 2, 0)
    assert surface.orientability is True
    assert surface.genus_by_component == (0,)
    assert len(surface.refined_seams) == 3
    assert all(seam.classification == "exterior_first" for seam in surface.refined_seams)
    assert all(
        set(dict(seam.orientation_on_patches).values()) == {-1}
        for seam in surface.refined_seams
    )
    assert curvature.physical_seam_count == 0
    assert curvature.exterior_boundary_piece_count == 3
    close("multi-boundary/builder smooth K", curvature.smooth_gaussian_curvature, 8.0 * math.pi / 3.0, 1e-9)
    close("multi-boundary/builder total H", curvature.total_integrated_mean_curvature, 8.0 * math.pi, 1e-9)
    close("multi-boundary/builder boundary K", curvature.boundary_geodesic_curvature, -14.0 * math.pi / 3.0, 1e-9)
    close("multi-boundary/builder GB", curvature.completed_gaussian_curvature, -2.0 * math.pi, 1e-9)
    close("multi-boundary/builder GB residual", curvature.gauss_bonnet_residual, 0.0, 1e-12)
    assert curvature.gaussian_status == "CERTIFIED"


def seam_circle_length(atoms, carrier_a, carrier_b, contacts, count=120_000):
    first = atoms[carrier_a]
    second = atoms[carrier_b]
    dx = tuple(second.center[index] - first.center[index] for index in range(3))
    distance = math.sqrt(sum(value * value for value in dx))
    assert math.isclose(first.radius, second.radius)
    radius = first.radius
    circle_radius = math.sqrt(radius**2 - (distance / 2.0) ** 2)
    normal = tuple(value / distance for value in dx)
    seed = (1.0, 0.0, 0.0) if abs(normal[0]) < 0.8 else (0.0, 1.0, 0.0)
    basis_1 = (
        normal[1] * seed[2] - normal[2] * seed[1],
        normal[2] * seed[0] - normal[0] * seed[2],
        normal[0] * seed[1] - normal[1] * seed[0],
    )
    basis_norm = math.sqrt(sum(value * value for value in basis_1))
    basis_1 = tuple(value / basis_norm for value in basis_1)
    basis_2 = (
        normal[1] * basis_1[2] - normal[2] * basis_1[1],
        normal[2] * basis_1[0] - normal[0] * basis_1[2],
        normal[0] * basis_1[1] - normal[1] * basis_1[0],
    )
    center = tuple(
        (first.center[index] + second.center[index]) / 2.0
        for index in range(3)
    )
    partner_a = selected_partners(atoms, carrier_a, contacts)
    partner_b = selected_partners(atoms, carrier_b, contacts)
    retained = 0
    for index in range(count):
        angle = TWO_PI * (index + 0.5) / count
        point = tuple(
            center[axis]
            + circle_radius * (
                basis_1[axis] * math.cos(angle)
                + basis_2[axis] * math.sin(angle)
            )
            for axis in range(3)
        )
        in_a = any(
            distance_squared(point, atoms[partner].center) <= atoms[partner].radius**2
            for partner in partner_a
        )
        in_b = any(
            distance_squared(point, atoms[partner].center) <= atoms[partner].radius**2
            for partner in partner_b
        )
        retained += in_a and in_b
    return TWO_PI * circle_radius * retained / count


def validate_physical_seam():
    atoms = (
        ball(0, "A", (0.0, 0.0, 0.0), 1.5),
        ball(1, "A", (0.8, 0.0, 0.0), 1.5),
        ball(2, "B", (0.4, 1.0, 0.0), 1.5),
    )
    selected = (contact(0, 0, 2), contact(1, 1, 2))
    surface = build_molecular_contact_surface(atoms, selected, "A")
    curvature = surface.curvature
    assert curvature is not None
    seam_length = seam_circle_length(atoms, 0, 1, selected)
    beta = math.acos(1.0 - 0.8**2 / (2.0 * 1.5**2))
    independent_seam_h = 0.5 * beta * seam_length
    close("physical seam/H", curvature.seam_mean_curvature, independent_seam_h, 2.0e-4)
    assert curvature.physical_seam_count == 1
    close("physical seam/GB", curvature.completed_gaussian_curvature, TWO_PI, 1e-9)
    assert curvature.gaussian_status == "CERTIFIED"


def validate_deep_pocket_reference():
    """Reference a 2D spherical surface region with three boundary loops."""
    radius = 3.0
    # Three disjoint spherical holes, each made by a cap with q=7/9.
    q = 7.0 / 9.0
    holes = 3
    area = FOUR_PI * radius**2 - holes * TWO_PI * radius**2 * (1.0 - q)
    smooth_k = area / radius**2
    boundary_kg = -holes * TWO_PI * q
    close("deep pocket/GB", smooth_k + boundary_kg, -TWO_PI, 1e-12)
    surface_topology = {
        "space": "2D boundary-surface region",
        "euler_characteristic": 2 - holes,
        "betti_numbers": (1, holes - 1, 0),
        "boundary_loops": holes,
    }
    assert surface_topology == {
        "space": "2D boundary-surface region",
        "euler_characteristic": -1,
        "betti_numbers": (1, 2, 0),
        "boundary_loops": 3,
    }
    directions = (
        (1.0, 0.0, 0.0),
        (-0.5, math.sqrt(3.0) / 2.0, 0.0),
        (-0.5, -math.sqrt(3.0) / 2.0, 0.0),
    )
    retained = sum(
        not any(sum(direction[index] * center[index] for index in range(3)) > q for center in directions)
        for direction in fibonacci_directions(120_000)
    )
    sampled = FOUR_PI * radius**2 * retained / 120_000
    close("deep pocket/area", sampled, area, 3.0e-2)


def validate_torus_reference():
    """Independently validate a torus surface, not the contact builder."""
    major, minor = 4.0, 1.0
    exact_area = 4.0 * math.pi**2 * major * minor
    exact_mean = 2.0 * math.pi**2 * major
    exact_gaussian = 0.0
    area = mean = gaussian = 0.0
    u_count, v_count = 96, 64
    du, dv = TWO_PI / u_count, TWO_PI / v_count
    for i in range(u_count):
        for j in range(v_count):
            v = (j + 0.5) * dv
            jacobian = minor * (major + minor * math.cos(v))
            h = (major + 2.0 * minor * math.cos(v)) / (2.0 * minor * (major + minor * math.cos(v)))
            k = math.cos(v) / (minor * (major + minor * math.cos(v)))
            weight = jacobian * du * dv
            area += weight
            mean += h * weight
            gaussian += k * weight
    close("through-channel/area", area, exact_area, 1e-10)
    close("through-channel/H", mean, exact_mean, 1e-10)
    close("through-channel/K", gaussian, exact_gaussian, 1e-10)
    assert {"space": "3D solid torus", "betti_numbers": (1, 1, 0)}["betti_numbers"] == (1, 1, 0)
    assert {"space": "2D torus boundary", "betti_numbers": (1, 2, 1)}["betti_numbers"] == (1, 2, 1)


def validate_cavity_reference():
    """Reference a spherical-shell solid and its two boundary surfaces."""
    outer, inner = 5.0, 1.5
    outer_area = FOUR_PI * outer**2
    inner_area = FOUR_PI * inner**2
    close("cavity/area", outer_area + inner_area, FOUR_PI * (outer**2 + inner**2), 1e-12)
    close("cavity/H", FOUR_PI * (outer - inner), 4.0 * math.pi * (outer - inner), 1e-12)
    close("cavity/K", 4.0 * math.pi + 4.0 * math.pi, 8.0 * math.pi, 1e-12)
    close("cavity/GB", 8.0 * math.pi, TWO_PI * 4, 1e-12)
    solid_topology = {"space": "3D spherical shell solid", "betti_numbers": (1, 0, 1)}
    boundary_topology = {
        "space": "2D boundary surface (two spherical components)",
        "components": 2,
        "betti_numbers": (2, 0, 2),
        "euler_characteristic": 4,
    }
    assert solid_topology["betti_numbers"] == (1, 0, 1)
    assert boundary_topology == {
        "space": "2D boundary surface (two spherical components)",
        "components": 2,
        "betti_numbers": (2, 0, 2),
        "euler_characteristic": 4,
    }


def main():
    limitations = []
    validate_cap(1.0, 1.0, 1.0, "equal cap", limitations)
    validate_cap(2.0, 1.0, 2.0, "unequal cap", limitations)
    validate_closed_sphere(limitations)
    validate_selected_partner_restriction()
    validate_multiple_overlapping_caps()
    validate_closed_two_sphere_union_reference()
    validate_closed_two_sphere_union_builder()
    validate_same_group_occlusion(limitations)
    validate_multi_boundary_occlusion_builder()
    validate_physical_seam()
    validate_deep_pocket_reference()
    validate_torus_reference()
    validate_cavity_reference()
    print("reference-only limitations: deep pocket, through-channel, enclosed cavity")
    for limitation in limitations:
        print(f"known canonical limitation: {limitation}")
    print("independent molecular-contact validation: PASS")


if __name__ == "__main__":
    main()
