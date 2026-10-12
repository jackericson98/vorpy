import math

import numpy as np
import pytest

from vorpy.src.analyze.open_interface_curvature import (
    analyze_surface_topology,
    gauss_bonnet_accounting,
    integrated_mean_curvature,
    polygon_gaussian_terms,
    planar_crease_mean,
)
from vorpy.src.calculations.curvature import (
    integrate_quadratic_surface_gaussian_curvature,
)
from vorpy.src.network.fill import calc_surf_point_from_plane


def _subdivide_projected_quadratic_mesh(
        func, points, tris, plane_point, plane_normal, reference_point, levels):
    plane_point = np.asarray(plane_point, dtype=float)
    plane_normal = np.asarray(plane_normal, dtype=float)
    plane_normal /= np.linalg.norm(plane_normal)
    reference_point = np.asarray(reference_point, dtype=float)
    coefficients = np.asarray(func, dtype=float)
    leaves = []

    def split(first, second, third, depth):
        if depth == 0:
            leaves.append((first, second, third))
            return
        first_second = 0.5 * (first + second)
        second_third = 0.5 * (second + third)
        third_first = 0.5 * (third + first)
        split(first, first_second, third_first, depth - 1)
        split(first_second, second, second_third, depth - 1)
        split(third_first, second_third, third, depth - 1)
        split(first_second, second_third, third_first, depth - 1)

    for triangle in np.asarray(tris, dtype=int):
        actual = np.asarray(points, dtype=float)[triangle]
        projected = actual - np.outer(
            (actual - plane_point) @ plane_normal, plane_normal
        )
        split(*projected, int(levels))

    lifted = []
    subdivided_tris = []
    for triangle in leaves:
        indices = []
        for point in triangle:
            lifted_point = calc_surf_point_from_plane(
                point, plane_normal, coefficients, reference_point
            )
            if lifted_point is None:
                raise ValueError("Projected subdivision left the quadratic branch domain.")
            indices.append(len(lifted))
            lifted.append(lifted_point)
        subdivided_tris.append(tuple(indices))
    return np.asarray(lifted), np.asarray(subdivided_tris, dtype=int)


def test_planar_polygon_disk_gauss_bonnet_corner_turns():
    coords = {0: (0, 0, 0), 1: (2, 0, 0), 2: (2, 1, 0), 3: (0, 1, 0)}
    result = polygon_gaussian_terms([(0, 1, 2, 3)], coords)
    assert result["topology"]["euler_characteristic"] == 1
    assert result["smooth_K"] == 0
    assert result["boundary_geodesic_K"] == 0
    assert result["boundary_corner_K"] == pytest.approx(2 * math.pi)
    assert result["certified"]


def test_hemisphere_gauss_bonnet_exact_smooth_and_geodesic_terms():
    # Unit hemisphere: integral K=2pi; equator is a geodesic, so kg=0.
    result = gauss_bonnet_accounting(smooth_K=2 * math.pi,
        face_boundary_geodesic=0, vertex_corner_terms=0,
        euler_characteristic=1)
    assert result["total_GB_left"] == pytest.approx(2 * math.pi)
    assert result["certified"]


def test_full_sphere_and_cube_reference_decompose_gaussian_and_mean_terms():
    sphere = gauss_bonnet_accounting(
        smooth_K=4 * math.pi,
        face_boundary_geodesic=0,
        vertex_corner_terms=0,
        euler_characteristic=2,
    )
    sphere_mean = integrated_mean_curvature([4 * math.pi * 2.0], [])
    assert sphere["certified"]
    assert sphere["total_GB_left"] == pytest.approx(4 * math.pi)
    assert sphere["face_boundary_geodesic_K"] == 0
    assert sphere_mean["total_H"] == pytest.approx(8 * math.pi)

    coordinates = {
        0: (-1, -1, -1), 1: (1, -1, -1), 2: (1, 1, -1), 3: (-1, 1, -1),
        4: (-1, -1, 1), 5: (1, -1, 1), 6: (1, 1, 1), 7: (-1, 1, 1),
    }
    cube = polygon_gaussian_terms([
        (0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
        (3, 7, 6, 2), (0, 4, 7, 3), (1, 2, 6, 5),
    ], coordinates)
    assert cube["topology"]["vertices"] == 8
    assert cube["topology"]["edges"] == 12
    assert cube["topology"]["faces"] == 6
    assert cube["topology"]["euler_characteristic"] == 2
    assert cube["smooth_K"] == 0
    assert cube["boundary_geodesic_K"] == 0
    assert cube["interior_vertex_defect_K"] == pytest.approx(4 * math.pi)
    assert cube["intrinsic_total_K"] == pytest.approx(4 * math.pi)
    assert cube["certified"]

    # A side-2 convex cube has twelve edges of length 2 and exterior turn pi/2.
    cube_mean = integrated_mean_curvature([0.], [12 * 0.5 * 2.0 * (math.pi / 2)])
    assert cube_mean["total_H"] == pytest.approx(6 * math.pi)


def test_curved_quadratic_gaussian_integrator_spherical_cap_converges():
    radius = 2.0
    theta = math.pi / 3.0
    boundary_count = 128
    points = [(0.0, 0.0, radius)]
    for index in range(boundary_count):
        angle = 2.0 * math.pi * index / boundary_count
        points.append((
            radius * math.sin(theta) * math.cos(angle),
            radius * math.sin(theta) * math.sin(angle),
            radius * math.cos(theta),
        ))
    triangles = [
        (0, index + 1, (index + 1) % boundary_count + 1)
        for index in range(boundary_count)
    ]
    sphere = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -radius * radius, 0.0, 0.0, 0.0, 0.0]
    values = [integrate_quadratic_surface_gaussian_curvature(
        sphere, points, triangles, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
        (0.0, 0.0, radius), order=order,
    ) for order in (4, 8, 16)]
    expected = 2.0 * math.pi * (1.0 - math.cos(theta))
    assert values[2] == pytest.approx(expected, abs=2e-3)
    assert values[1] == pytest.approx(values[2], abs=1e-6)


def test_projected_domain_subdivision_reduces_fixed_face_quadrature_error():
    radius = 2.0
    theta = math.pi / 3.0
    boundary_count = 8
    points = [(0.0, 0.0, radius)]
    for index in range(boundary_count):
        angle = 2.0 * math.pi * index / boundary_count
        points.append((
            radius * math.sin(theta) * math.cos(angle),
            radius * math.sin(theta) * math.sin(angle),
            radius * math.cos(theta),
        ))
    triangles = [
        (0, index + 1, (index + 1) % boundary_count + 1)
        for index in range(boundary_count)
    ]
    sphere = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -radius * radius, 0.0, 0.0, 0.0, 0.0]
    fixed = {
        order: integrate_quadratic_surface_gaussian_curvature(
            sphere, points, triangles, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
            (0.0, 0.0, radius), order=order,
        )
        for order in (4, 8, 16)
    }
    subdivided_points, subdivided_tris = _subdivide_projected_quadratic_mesh(
        sphere, points, triangles, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
        (0.0, 0.0, radius), levels=1,
    )
    adaptive = {
        order: integrate_quadratic_surface_gaussian_curvature(
            sphere, subdivided_points, subdivided_tris,
            (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, radius), order=order,
        )
        for order in (4, 8)
    }
    print("G13 spherical-face comparison:", fixed, "adaptive=", adaptive)
    assert len(subdivided_tris) == 4 * len(triangles)
    assert abs(adaptive[4] - fixed[16]) < abs(fixed[4] - fixed[16])
    assert abs(fixed[8] - fixed[16]) < abs(adaptive[4] - fixed[16])
    assert abs(fixed[4] - fixed[8]) == pytest.approx(
        abs(fixed[4] - fixed[16]), rel=2e-3
    )


def test_curved_quadratic_gaussian_integrator_rejects_incomplete_surface_metadata():
    with pytest.raises(ValueError):
        integrate_quadratic_surface_gaussian_curvature(
            [1.0] * 11,
            [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
            [(0, 1, 2)],
            (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0),
        )


@pytest.mark.parametrize("theta", [math.pi / 6, math.pi / 3, math.pi / 2])
def test_spherical_cap_non_geodesic_boundary(theta):
    # Unit sphere cap, outward orientation; integral kg on its latitude circle
    # is 2pi*cos(theta), cancelling the cap's 2pi*(1-cos(theta)) K integral.
    smooth_k = 2 * math.pi * (1 - math.cos(theta))
    boundary_kg = 2 * math.pi * math.cos(theta)
    result = gauss_bonnet_accounting(smooth_K=smooth_k,
        face_boundary_geodesic=boundary_kg, vertex_corner_terms=0,
        euler_characteristic=1)
    assert result["certified"]
    assert result["total_GB_left"] == pytest.approx(2 * math.pi)


def test_planar_triangulation_does_not_change_gaussian_total():
    coords = {0: (0, 0, 0), 1: (1, 0, 0), 2: (1, 1, 0), 3: (0, 1, 0)}
    disk = polygon_gaussian_terms([(0, 1, 2, 3)], coords)
    triangulated = polygon_gaussian_terms([(0, 1, 2), (0, 2, 3)], coords)
    assert disk["topology"]["euler_characteristic"] == triangulated["topology"]["euler_characteristic"] == 1
    assert disk["boundary_corner_K"] == pytest.approx(triangulated["boundary_corner_K"])
    assert triangulated["topology"]["edges"] == 5  # diagonal is topological, not a boundary edge
    assert triangulated["certified"]


def test_folded_planar_patch_has_extrinsic_crease_mean_but_zero_intrinsic_gaussian():
    # Two triangles sharing a unit edge, folded by 60 degrees in R3.
    coords = {0: (0, 0, 0), 1: (1, 0, 0), 2: (0, 1, 0),
              3: (0, -0.5, math.sqrt(3) / 2)}
    cycles = [(0, 1, 2), (1, 0, 3)]
    gb = polygon_gaussian_terms(cycles, coords)
    crease = planar_crease_mean(*cycles, coords)
    mean = integrated_mean_curvature([0, 0], [crease])
    assert gb["certified"]
    assert gb["total_GB_left"] == pytest.approx(2 * math.pi)
    assert mean["internal_edge_H"] == pytest.approx(crease)
    assert mean["boundary_H"] == 0
    assert mean["total_H"] == pytest.approx(crease)
    assert gb["intrinsic_total_K"] == pytest.approx(0)
    assert planar_crease_mean(tuple(reversed(cycles[0])),
                             tuple(reversed(cycles[1])), coords) == pytest.approx(-crease)


def test_vertex_only_pinch_is_not_manifold_and_euler_is_cell_count():
    topo = analyze_surface_topology([(0, 1, 2), (0, 3, 4)])
    assert not topo["boundary_is_manifold"]
    assert topo["euler_characteristic"] == 5-6+2
    assert topo["genus_sum"] is None


def test_unresolved_or_nonfinite_mean_cannot_be_certified():
    assert integrated_mean_curvature([1], [2], unresolved_internal_edges=1)["total_H"] is None
    assert integrated_mean_curvature([float("nan")], [2])["total_H"] is None


def test_disk_fan_interior_vertex_has_no_intrinsic_curvature():
    coords = {0:(0,0,0), 1:(1,0,0), 2:(0,1,0), 3:(-1,0,0), 4:(0,-1,0)}
    result = polygon_gaussian_terms([(0,1,2),(0,2,3),(0,3,4),(0,4,1)], coords)
    assert result["interior_vertex_defect_K"] == pytest.approx(0)
    assert result["boundary_vertex_corner_K"] == pytest.approx(2*math.pi)
    assert result["certified"]


def test_nonmanifold_edge_is_rejected():
    assert not analyze_surface_topology([(0,1,2),(1,0,3),(0,1,4)])["boundary_is_manifold"]


@pytest.mark.parametrize("radius,theta", [(1,math.pi/2),(2,math.pi/3),(3,math.pi/6)])
def test_spherical_parameterization_integrals_and_boundary_mean(radius, theta):
    import numpy as np
    # X(u,v)=R(sin u cos v,sin u sin v,cos u). Dn=I/R,
    # |Xu cross Xv|=R^2 sin u. Integrate independently of GB.
    nodes, weights = np.polynomial.legendre.leggauss(32)
    u = (nodes+1)*theta/2
    dA = radius**2*np.sin(u)*theta/2*weights*2*math.pi
    smooth_h = float(np.sum(dA/radius))
    smooth_k = float(np.sum(dA/radius**2))
    # T=(-sin v,cos v,0); dT/ds=(-cos v,-sin v,0)/(R sin theta).
    # Positive boundary inward conormal n cross T gives kg=cos(theta)/(R sin(theta)).
    kg = math.cos(theta)/(radius*math.sin(theta))*2*math.pi*radius*math.sin(theta)
    mean = integrated_mean_curvature([smooth_h], [])
    assert mean["total_H"] == pytest.approx(2*math.pi*radius*(1-math.cos(theta)))
    assert mean["boundary_H"] == 0
    gb = gauss_bonnet_accounting(smooth_K=smooth_k, face_boundary_geodesic=kg,
        vertex_corner_terms=0, euler_characteristic=1)
    assert gb["certified"]


def test_smooth_planar_disk_has_boundary_geodesic_turn_but_no_H_or_K():
    radius = 2.0
    boundary_kg = (1/radius)*(2*math.pi*radius)
    result = gauss_bonnet_accounting(smooth_K=0, face_boundary_geodesic=boundary_kg,
        vertex_corner_terms=0, euler_characteristic=1)
    assert result["certified"]
    assert integrated_mean_curvature([0], [])["total_H"] == 0


def test_two_spherical_caps_glued_at_rim_have_intrinsic_seam_curvature():
    # The intrinsic metric of two identical latitude caps glued at their rims
    # is topologically a sphere. Their induced boundary kg terms add.
    theta = math.pi/3
    smooth_k = 4*math.pi*(1-math.cos(theta))
    seam_k = 4*math.pi*math.cos(theta)
    result = gauss_bonnet_accounting(smooth_K=smooth_k, face_boundary_geodesic=seam_k,
        vertex_corner_terms=0, euler_characteristic=2)
    assert seam_k != 0
    assert result["certified"]
    assert result["total_GB_left"] == pytest.approx(4*math.pi)


def test_free_boundary_term_is_not_mean_curvature():
    with pytest.raises(ValueError, match="not mean curvature"):
        integrated_mean_curvature([0], [], boundary_terms=[2 * math.pi])


def test_open_boundary_corner_is_gaussian_only_and_mean_total_can_certify():
    coords = {0: (0, 0, 0), 1: (2, 0, 0), 2: (2, 1, 0), 3: (0, 1, 0)}
    gaussian = polygon_gaussian_terms([(0, 1, 2, 3)], coords)
    mean = integrated_mean_curvature([0.], [])

    assert gaussian["boundary_corner_K"] == pytest.approx(2 * math.pi)
    assert mean["boundary_H"] == 0.
    assert mean["total_H"] == 0.
    assert mean["certified"] is True


def test_closed_regular_interface_has_no_boundary_or_corner_mean_term():
    tetrahedron = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
    topology = analyze_surface_topology(tetrahedron)
    mean = integrated_mean_curvature([0.], [])

    assert topology["boundary_components"] == 0
    assert topology["boundary_is_manifold"] and topology["orientable"]
    assert topology["genus_sum"] == 0
    assert mean["boundary_H"] == 0.
    assert mean["total_H"] == 0.
    assert mean["certified"] is True


def test_topology_components_boundary_loops_and_genus_for_annulus():
    cycles = [(0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    topo = analyze_surface_topology(cycles)
    assert topo["vertices"] == 8
    assert topo["edges"] == 12
    assert topo["faces"] == 4
    assert topo["euler_characteristic"] == 0
    assert topo["components"] == 1
    assert topo["boundary_components"] == 2
    assert topo["orientable"] and topo["boundary_is_manifold"]
    assert topo["genus_sum"] == 0


def test_group_reversal_flips_mean_but_not_gauss_bonnet_or_topology():
    forward = integrated_mean_curvature([1.25], [-2.0])
    reverse = integrated_mean_curvature([-1.25], [2.0])
    assert reverse["total_H"] == pytest.approx(-forward["total_H"])
    gb_a = gauss_bonnet_accounting(smooth_K=0.8, face_boundary_geodesic=1.2,
        vertex_corner_terms=2 * math.pi - 2.0, euler_characteristic=1)
    gb_b = gauss_bonnet_accounting(smooth_K=0.8, face_boundary_geodesic=1.2,
        vertex_corner_terms=2 * math.pi - 2.0, euler_characteristic=1)
    assert gb_a["total_GB_left"] == gb_b["total_GB_left"]
    assert gb_a["two_pi_chi"] == gb_b["two_pi_chi"]


def test_reversing_polygon_orientation_preserves_intrinsic_gaussian_and_topology():
    coords = {0: (0, 0, 0), 1: (1, 0, 0), 2: (1, 1, 0), 3: (0, 1, 0)}
    forward = polygon_gaussian_terms([(0, 1, 2, 3)], coords)
    reverse = polygon_gaussian_terms([(3, 2, 1, 0)], coords)
    assert forward["topology"]["euler_characteristic"] == reverse["topology"]["euler_characteristic"]
    assert forward["total_GB_left"] == pytest.approx(reverse["total_GB_left"])
    assert forward["boundary_corner_K"] == pytest.approx(reverse["boundary_corner_K"])
