import math

import pytest

from vorpy.src.analyze.open_interface_curvature import (
    analyze_surface_topology,
    gauss_bonnet_accounting,
    integrated_mean_curvature,
    polygon_gaussian_terms,
    planar_crease_mean,
)


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
