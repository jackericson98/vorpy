"""Synthetic regression tests for curved molecular contact surfaces."""

from __future__ import annotations

from dataclasses import replace
import math

import pytest

from vorpy.src.geometry.molecular_contact import (
    AtomBall,
    SelectedContact,
    audit_molecular_contact_surface_gauss_bonnet,
    build_molecular_contact_surface,
    build_molecular_contact_surface_pair,
    compute_molecular_contact_surface_curvature,
)
from vorpy.src.output.visualization import _iter_result_quantities
from vorpy.src.results import (
    Provenance,
    RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG,
    ResultIdentity,
    Quantity,
    Status,
)


def ball(index, group, xyz, radius):
    return AtomBall(index, f"{group}-{index}", group, tuple(xyz), radius)


def contact(index, a, b):
    return SelectedContact(f"contact-{index}", a, b, "network", "interface")


def one_pair(radius_a, radius_b, distance):
    atoms = (
        ball(0, "A", (0, 0, 0), radius_a),
        ball(1, "B", (distance, 0, 0), radius_b),
    )
    return build_molecular_contact_surface_pair(atoms, (contact(0, 0, 1),))


def test_equal_overlap_exact_cap_area():
    side_a, side_b = one_pair(1.0, 1.0, 1.0)
    assert side_a.status == "CERTIFIED"
    assert side_b.status == "CERTIFIED"
    assert side_a.total_area_A2 == pytest.approx(math.pi, rel=1e-10, abs=1e-10)
    assert side_b.total_area_A2 == pytest.approx(math.pi, rel=1e-10, abs=1e-10)


def test_unequal_overlap_exact_cap_area():
    side_a, side_b = one_pair(2.0, 1.0, 2.0)
    assert side_a.total_area_A2 == pytest.approx(math.pi, rel=1e-10, abs=1e-10)
    assert side_b.total_area_A2 == pytest.approx(1.5 * math.pi, rel=1e-10, abs=1e-10)


@pytest.mark.parametrize("distance", (3.0, 4.0))
def test_tangent_and_disjoint_have_no_positive_area(distance):
    side_a, side_b = one_pair(1.0, 2.0, distance)
    assert side_a.total_area_A2 == pytest.approx(0.0)
    assert side_b.total_area_A2 == pytest.approx(0.0)


def test_containment():
    side_a, side_b = one_pair(1.0, 3.0, 1.0)
    assert side_a.total_area_A2 == pytest.approx(4.0 * math.pi)
    assert side_b.total_area_A2 == pytest.approx(0.0)


def test_full_sphere_topology_and_gauss_bonnet_are_certified():
    surface = build_molecular_contact_surface(
        (ball(0, "A", (0, 0, 0), 1.0), ball(1, "B", (0, 0, 0), 10.0)),
        (contact(0, 0, 1),),
        "A",
    )
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 0, 2)
    assert surface.betti_numbers == (1, 0, 1)
    assert surface.orientability is True
    assert surface.genus_by_component == (0,)
    assert surface.curvature is not None
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(4.0 * math.pi)
    assert surface.curvature.gaussian_status == "CERTIFIED"


def test_disconnected_full_sphere_components_remain_partial_without_cell_poles():
    atoms = (
        ball(0, "A", (0, 0, 0), 1.0),
        ball(1, "A", (5, 0, 0), 1.0),
        ball(2, "B", (0, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(
        atoms, (contact(0, 0, 2), contact(1, 1, 2)), "A"
    )
    assert (surface.components, surface.boundary_loops) == (2, 0)
    assert surface.euler_characteristic is None
    assert surface.status == "PARTIAL"
    assert surface.curvature is not None
    assert surface.curvature.completed_gaussian_curvature is None
    assert surface.curvature.gaussian_status == "PARTIAL"


def test_closed_equal_sphere_union_topology_and_seam_are_certified():
    atoms = (
        ball(0, "A", (0, 0, 0), 1.0),
        ball(1, "A", (1, 0, 0), 1.0),
        ball(2, "B", (0.5, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(
        atoms, (contact(0, 0, 2), contact(1, 1, 2)), "A"
    )
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 0, 2)
    assert surface.betti_numbers == (1, 0, 1)
    assert len(surface.refined_seams) == 1
    assert set(dict(surface.refined_seams[0].orientation_on_patches).values()) == {-1, 1}
    assert surface.curvature is not None
    assert surface.curvature.intrinsic_seam_gaussian_curvature == pytest.approx(-2.0 * math.pi)
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(4.0 * math.pi)
    assert surface.curvature.gaussian_status == "CERTIFIED"


def test_three_same_group_occlusion_holes_have_certified_boundary_topology():
    directions = (
        (1.0, 0.0, 0.0),
        (-0.5, math.sqrt(3.0) / 2.0, 0.0),
        (-0.5, -math.sqrt(3.0) / 2.0, 0.0),
    )
    atoms = (
        ball(0, "A", (0, 0, 0), 3.0),
        *(ball(index, "A", tuple(3.0 * value for value in direction), 2.0)
          for index, direction in enumerate(directions, start=1)),
        ball(4, "B", (0, 0, 0), 20.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 4),), "A")
    assert surface.total_area_A2 == pytest.approx(24.0 * math.pi)
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 3, -1)
    assert surface.betti_numbers == (1, 2, 0)
    assert surface.curvature is not None
    assert surface.curvature.boundary_geodesic_curvature == pytest.approx(-14.0 * math.pi / 3.0)
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(-2.0 * math.pi)
    assert surface.curvature.gaussian_status == "CERTIFIED"


def test_topological_cuts_do_not_change_physical_curvature_terms():
    directions = (
        (1.0, 0.0, 0.0),
        (-0.5, math.sqrt(3.0) / 2.0, 0.0),
        (-0.5, -math.sqrt(3.0) / 2.0, 0.0),
    )
    atoms = (
        ball(0, "A", (0, 0, 0), 3.0),
        *(ball(index, "A", tuple(3.0 * value for value in direction), 2.0)
          for index, direction in enumerate(directions, start=1)),
        ball(4, "B", (0, 0, 0), 20.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 4),), "A")
    assert surface.topological_cuts
    without_cuts = replace(surface, topological_cuts=())
    cut_curvature = compute_molecular_contact_surface_curvature(surface)
    uncut_curvature = compute_molecular_contact_surface_curvature(without_cuts)
    for field_name in (
        "smooth_gaussian_curvature",
        "intrinsic_seam_gaussian_curvature",
        "boundary_geodesic_curvature",
        "interior_vertex_defects",
        "boundary_corner_turning",
        "completed_gaussian_curvature",
        "gauss_bonnet_residual",
    ):
        assert getattr(uncut_curvature, field_name) == pytest.approx(
            getattr(cut_curvature, field_name)
        )


def test_overlapping_partner_caps_are_not_pairwise_sum():
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "B", (2.0, 0, 0), 2.0),
        # The partner directions are 60 degrees apart, so the 60-degree
        # contact caps overlap with positive area rather than merely touching.
        ball(2, "B", (1.0, math.sqrt(3.0), 0), 2.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 1), contact(1, 0, 2)), "A")
    assert surface.total_area_A2 < surface.naive_pairwise_area_A2
    assert len(surface.patches) == 1


def test_same_group_occlusion_is_applied():
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "A", (1.0, 0, 0), 2.0),
        ball(2, "B", (0, 0, 3.0), 2.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 2),), "A")
    assert surface.total_area_A2 < surface.naive_pairwise_area_A2
    assert "A-1" in {atom for patch in surface.patches for atom in patch.self_occluding_atom_ids}


def test_intersecting_same_group_occlusion_boundaries_certify_with_junctions():
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "A", (1.0, 0, 0), 2.0),
        ball(2, "A", (0.5, math.sqrt(3.0) / 2.0, 0), 2.0),
        ball(3, "B", (0, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 3),), "A")
    assert surface.curvature is not None
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 1, 1)
    assert surface.betti_numbers == (1, 0, 0)
    assert surface.singular_vertices == 0
    assert surface.curvature.gaussian_status == "CERTIFIED"
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(2.0 * math.pi)
    assert surface.curvature.gauss_bonnet_residual == pytest.approx(0.0, abs=1e-12)


def test_multiple_overlapping_occlusion_caps_keep_component_topology_explicit():
    angle = math.radians(10.0)
    cluster_centers = (
        (2.0 * math.cos(angle), 0.0, 2.0 * math.sin(angle)),
        (2.0 * math.cos(angle), 0.0, -2.0 * math.sin(angle)),
        (-2.0 * math.cos(angle), 0.0, 2.0 * math.sin(angle)),
        (-2.0 * math.cos(angle), 0.0, -2.0 * math.sin(angle)),
    )
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        *(ball(index, "A", center, 1.0)
          for index, center in enumerate(cluster_centers, start=1)),
        ball(5, "B", (0, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 5),), "A")
    assert (surface.components, surface.boundary_loops) == (1, 2)
    assert surface.euler_characteristic == 0
    assert surface.betti_numbers == (1, 1, 0)
    assert surface.topological_cuts
    assert surface.junction_count > 0
    assert surface.singular_vertices == 0
    assert surface.curvature is not None
    assert surface.curvature.physical_seam_count == 0
    assert surface.curvature.gaussian_status == "CERTIFIED"
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(0.0, abs=1e-12)
    assert surface.curvature.gauss_bonnet_residual == pytest.approx(0.0, abs=1e-12)


def test_three_cap_belt_reports_two_disconnected_polar_components():
    # Each equal-radius occluder has carrier-sphere threshold 1/4, hence
    # cap half-angle acos(1/4).  Three centers 120 degrees apart close an
    # equatorial belt; the retained set has two polar components, not an annulus.
    directions = (
        (1.0, 0.0, 0.0),
        (-0.5, math.sqrt(3.0) / 2.0, 0.0),
        (-0.5, -math.sqrt(3.0) / 2.0, 0.0),
    )
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        *(ball(index, "A", direction, 2.0)
          for index, direction in enumerate(directions, start=1)),
        ball(4, "B", (0, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 4),), "A")
    assert surface.total_area_A2 == pytest.approx(2.686115311330023, abs=1e-12)
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (2, 2, 2)
    assert surface.betti_numbers == (2, 0, 0)
    assert surface.status == "CERTIFIED"
    assert surface.incidence_status == "CERTIFIED"
    assert surface.singular_vertices == 0
    assert len({patch.component_id for patch in surface.patches}) == 2
    assert all(len(patch.boundary_arc_ids) == 3 for patch in surface.patches)
    assert surface.topological_cuts == ()
    assert surface.curvature is not None
    assert surface.curvature.gaussian_status == "CERTIFIED"
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(4.0 * math.pi, abs=1e-12)
    assert surface.curvature.gauss_bonnet_residual == pytest.approx(0.0, abs=1e-12)
    patchwise = audit_molecular_contact_surface_gauss_bonnet(surface)
    assert patchwise["status"] == "CERTIFIED"
    assert patchwise["max_abs_residual"] == pytest.approx(0.0, abs=1e-12)


def test_two_disjoint_occlusion_components_have_certified_annular_winding():
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "A", (2.0, 0, 0), 1.0),
        ball(2, "A", (-2.0, 0, 0), 1.0),
        ball(3, "B", (0, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 3),), "A")
    assert (surface.components, surface.boundary_loops, surface.euler_characteristic) == (1, 2, 0)
    assert surface.curvature is not None
    assert surface.curvature.gaussian_status == "CERTIFIED"
    assert surface.curvature.completed_gaussian_curvature == pytest.approx(0.0, abs=1e-12)
    assert surface.curvature.gauss_bonnet_residual == pytest.approx(0.0, abs=1e-12)


def test_reversing_one_annular_loop_is_not_silently_certified():
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "A", (2.0, 0, 0), 1.0),
        ball(2, "A", (-2.0, 0, 0), 1.0),
        ball(3, "B", (0, 0, 0), 10.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 3),), "A")
    assert len(surface.refined_seams) == 2
    reversed_piece = replace(
        surface.refined_seams[0],
        orientation_on_patches=tuple(
            (patch_id, -direction)
            for patch_id, direction in surface.refined_seams[0].orientation_on_patches
        ),
    )
    reversed_surface = replace(
        surface,
        refined_seams=(reversed_piece,) + surface.refined_seams[1:],
    )
    curvature = compute_molecular_contact_surface_curvature(reversed_surface)
    assert curvature.gaussian_status == "PARTIAL"
    assert curvature.gauss_bonnet_residual is not None
    assert abs(curvature.gauss_bonnet_residual) > 1e-6


def test_power_exposed_patch_identity():
    atoms = (
        ball(0, "A", (0, 0, 0), 1.0),
        ball(1, "A", (1.0, 0, 0), 1.0),
        ball(2, "B", (0, 0, 0), 3.0),
    )
    surface = build_molecular_contact_surface(atoms, (contact(0, 0, 2),), "A")
    assert surface.total_area_A2 == pytest.approx(3.0 * math.pi, rel=1e-9, abs=1e-9)


def test_a_b_reversal_and_rigid_translation_rotation():
    atoms = (
        ball(0, "A", (0, 0, 0), 1.5),
        ball(1, "B", (1.4, 0.3, 0.2), 1.1),
    )
    first_a, first_b = build_molecular_contact_surface_pair(atoms, (contact(0, 0, 1),))
    transformed = (
        ball(0, "A", (5.0, -2.0, 3.0), 1.5),
        ball(1, "B", (5.0, -2.0, 3.0), 1.1),
    )
    # Rotate the original displacement with a fixed orthogonal map, then translate.
    transformed = (
        transformed[0],
        ball(1, "B", (5.0 + 0.2, -2.0 + 1.4, 3.0 + 0.3), 1.1),
    )
    second_a, second_b = build_molecular_contact_surface_pair(transformed, (contact(0, 0, 1),))
    distance = math.sqrt(1.4**2 + 0.3**2 + 0.2**2)

    def cap_area(radius, partner_radius):
        plane_distance = (distance**2 + radius**2 - partner_radius**2) / (2.0 * distance)
        return 2.0 * math.pi * radius * (radius - plane_distance)

    # Unequal carrier radii produce different spherical-cap areas; the
    # rigid-motion check below is the actual A/B reversal invariant.
    assert first_a.total_area_A2 == pytest.approx(cap_area(1.5, 1.1), rel=1e-10)
    assert first_b.total_area_A2 == pytest.approx(cap_area(1.1, 1.5), rel=1e-10)
    assert first_a.total_area_A2 != pytest.approx(first_b.total_area_A2, rel=1e-10)
    assert first_a.total_area_A2 == pytest.approx(second_a.total_area_A2, rel=1e-9)
    assert first_b.total_area_A2 == pytest.approx(second_b.total_area_A2, rel=1e-9)


def test_results_conversion_preserves_smooth_gaussian_and_unsupported_boundary_mean():
    atoms = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "B", (2.0, 0, 0), 2.0),
        ball(2, "B", (1.0, math.sqrt(3.0), 0), 2.0),
    )
    surface = build_molecular_contact_surface(
        atoms, (contact(0, 0, 1), contact(1, 0, 2)), "A",
        contact_selection_source="synthetic_contact_selection",
        contact_selection_id="synthetic:contacts",
    )
    assert surface.contact_selection_source == "synthetic_contact_selection"
    assert surface.summary()["contact_selection_source"] == "synthetic_contact_selection"
    identity = ResultIdentity(
        result_kind="molecular_contact_surface",
        system="synthetic",
        frame="static",
        network_id="synthetic:network",
        network_scope="interface",
        environment="dry",
        partition="power",
        representation="molecular_contact_surface",
        radius_configuration_id="synthetic:radii",
        system_id="synthetic",
        interface_id="synthetic:interface",
        group_A="A",
        group_B="B",
        side="A",
        orientation="group1 -> group2",
        partner_id="B",
        contact_selection_source=surface.contact_selection_source,
        contact_selection_id=surface.contact_selection_id,
        molecular_surface_convention="expanded_union_of_balls",
        definition_version="mcs-v1",
        source_network_id="synthetic:network",
        source_interface_id="synthetic:interface",
        interaction_id="synthetic:A:B",
        orientation_convention="outward_molecular_normal",
    )
    result = surface.to_contract(identity, Provenance("synthetic canonical contact"))

    assert result.metrics["area"].value == pytest.approx(surface.total_area_A2)
    assert result.metrics["area"].units == "Å²"
    assert result.metrics["area"].status is Status.CERTIFIED
    assert result.metrics["area"].method == surface.area_method
    assert result.metrics["area"].scope == "molecular_contact_surface"
    assert result.identity.representation == "molecular_contact_surface"
    assert result.identity.contact_selection_source == surface.contact_selection_source
    assert result.identity.contact_selection_id == surface.contact_selection_id
    assert result.metrics["area"].provenance is result.provenance
    assert [path for path, _ in _iter_result_quantities(result) if path == "metrics.area"] == [
        "metrics.area"
    ]

    partial = replace(surface, status="PARTIAL", area_status="PARTIAL")
    partial_result = partial.to_contract(identity, Provenance("partial molecular area"))
    assert partial_result.metrics["area"].value == pytest.approx(surface.total_area_A2)
    assert partial_result.metrics["area"].status is Status.PARTIAL

    legacy = replace(surface, status="PARTIAL", total_area_A2=None, area_status="PARTIAL")
    legacy_result = legacy.to_contract(identity, Provenance("legacy molecular area"))
    assert legacy_result.metrics["area"].value is None
    assert legacy_result.metrics["area"].status is Status.NOT_CALCULATED

    unresolved_result = replace(
        result,
        metrics={
            "area": Quantity(
                "area", None, "Å²", Status.UNRESOLVED,
                "unsupported legacy area", "molecular_contact_surface",
                provenance=result.provenance,
            )
        },
    )
    assert unresolved_result.metrics["area"].status is Status.UNRESOLVED

    for canonical_status in (Status.CERTIFIED, Status.PARTIAL, Status.STALE):
        status_surface = replace(surface, area_status=canonical_status.value)
        status_result = status_surface.to_contract(
            identity, Provenance(f"area status {canonical_status.value}")
        )
        assert status_result.metrics["area"].status is canonical_status

    for unavailable_status in (
        Status.UNRESOLVED,
        Status.NOT_CALCULATED,
        Status.NOT_SUPPORTED,
    ):
        unavailable_surface = replace(
            surface, total_area_A2=None, area_status=unavailable_status.value
        )
        unavailable_result = unavailable_surface.to_contract(
            identity, Provenance(f"area status {unavailable_status.value}")
        )
        assert unavailable_result.metrics["area"].value is None
        assert unavailable_result.metrics["area"].status is unavailable_status

    future_surface = replace(
        surface, total_area_A2=None, area_status="FUTURE_STATUS"
    )
    future_result = future_surface.to_contract(
        identity, Provenance("future area status")
    )
    assert future_result.metrics["area"].status is Status.UNRESOLVED

    for canonical_status in (Status.CERTIFIED, Status.PARTIAL, Status.STALE):
        result_status_surface = replace(surface, status=canonical_status.value)
        result_status = result_status_surface.to_contract(
            identity, Provenance(f"result status {canonical_status.value}")
        )
        assert result_status.status is canonical_status

    stale_partial_incidence = replace(
        surface, status=Status.STALE.value, incidence_status="PARTIAL"
    ).to_contract(identity, Provenance("stale partial incidence"))
    assert stale_partial_incidence.status is Status.STALE
    stale_quantities = dict(_iter_result_quantities(stale_partial_incidence))
    assert stale_quantities["topology.vertex_count"].status is Status.STALE

    assert result.mean_curvature.boundary_contribution.value is None
    assert result.mean_curvature.boundary_contribution.status is Status.NOT_SUPPORTED
    assert result.mean_curvature.edge_contribution.value == pytest.approx(
        surface.curvature.seam_mean_curvature
    )
    assert result.mean_curvature.complete_total.value == pytest.approx(
        surface.curvature.total_integrated_mean_curvature
    )
    intrinsic_total = (
        surface.curvature.smooth_gaussian_curvature
        + surface.curvature.intrinsic_seam_gaussian_curvature
        + surface.curvature.interior_vertex_defects
    )
    assert result.gaussian_curvature.integrated_gaussian_curvature.value == pytest.approx(
        intrinsic_total
    )
    assert result.gaussian_curvature.integrated_gaussian_curvature.scope == "intrinsic_surface"
    assert result.gaussian_curvature.smooth_integrated_gaussian_curvature.value == pytest.approx(
        surface.curvature.smooth_gaussian_curvature
    )
    assert result.gaussian_curvature.smooth_face_contribution.scope == "spherical_patches"
    assert result.gaussian_curvature.intrinsic_seam_contribution.value == pytest.approx(
        surface.curvature.intrinsic_seam_gaussian_curvature
    )
    assert result.gaussian_curvature.interior_vertex_defects.value == pytest.approx(
        surface.curvature.interior_vertex_defects
    )
    assert result.gaussian_curvature.boundary_geodesic_contribution.value == pytest.approx(
        surface.curvature.boundary_geodesic_curvature
    )
    assert result.gaussian_curvature.boundary_corner_contribution.value == pytest.approx(
        surface.curvature.boundary_corner_turning
    )
    assert result.gaussian_curvature.gauss_bonnet_completion.value == pytest.approx(
        surface.curvature.completed_gaussian_curvature
    )
    assert result.gaussian_curvature.integrated_gaussian_curvature.value != pytest.approx(
        result.gaussian_curvature.gauss_bonnet_completion.value
    )

    quantities = dict(_iter_result_quantities(result))
    logged = {
        RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG[path]: quantity
        for path, quantity in quantities.items()
        if path in RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG
    }
    assert logged['intrinsic_integrated_gaussian_curvature'].value == pytest.approx(
        intrinsic_total
    )
    assert logged['smooth_integrated_gaussian_curvature'].value == pytest.approx(
        surface.curvature.smooth_gaussian_curvature
    )
    assert logged['completed_gauss_bonnet'].value == pytest.approx(
        surface.curvature.completed_gaussian_curvature
    )
    assert 'intrinsic_total_gaussian_curvature' not in logged
    assert 'smooth_face_contribution' not in logged
