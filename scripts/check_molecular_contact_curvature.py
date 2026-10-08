"""Dependency-free curvature checks for canonical molecular contact geometry."""

from __future__ import annotations

import math
import sys
from dataclasses import asdict, replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vorpy.src.geometry.molecular_contact import (
    AtomBall,
    GeometryTolerance,
    RefinedSeamPiece,
    SelectedContact,
    TopologicalCut,
    build_molecular_contact_surface,
    build_molecular_contact_surface_pair,
    compute_molecular_contact_surface_curvature,
)


def ball(index, group, xyz, radius):
    return AtomBall(index, f"{group}-{index}", group, tuple(xyz), radius)


def contact(index, a, b):
    return SelectedContact(f"contact-{index}", a, b, "network", "interface")


def analytic_spherical_cycle_checks():
    """Check annular/multiple-boundary formulas without artificial geometry."""
    # On the unit sphere, a latitude circle u_z=a has integral kg ds=2*pi*a
    # for increasing longitude.  Reversing traversal changes only its sign.
    a_low, a_high = -0.2, 0.65
    smooth = 2.0 * math.pi * (a_high - a_low)
    boundary = -2.0 * math.pi * a_high + 2.0 * math.pi * a_low
    assert abs(boundary + smooth) < 1e-12  # annulus, chi=0
    assert abs((2.0 * math.pi * a_high) + (-2.0 * math.pi * a_high)) < 1e-12

    # Removing m disjoint spherical caps leaves chi=2-m.  Each cap boundary
    # is integrated from its own circle geometry, not from the target chi.
    caps = (0.2, 0.45, 0.7)
    remaining_smooth = 4.0 * math.pi - sum(2.0 * math.pi * (1.0 - a) for a in caps)
    hole_boundaries = sum(-2.0 * math.pi * a for a in caps)
    assert abs(remaining_smooth + hole_boundaries - 2.0 * math.pi * (2 - len(caps))) < 1e-12


def main():
    analytic_spherical_cycle_checks()

    # Two equal unit spheres at d=1: each retained cap has area pi and
    # geodesic boundary contribution pi, so C_H=pi and GB=2*pi.
    equal = (ball(0, "A", (0, 0, 0), 1.0), ball(1, "B", (1, 0, 0), 1.0))
    side_a, side_b = build_molecular_contact_surface_pair(equal, (contact(0, 0, 1),))
    for surface in (side_a, side_b):
        curvature = surface.curvature
        assert curvature is not None
        assert abs(curvature.smooth_mean_curvature - math.pi) < 1e-9
        assert curvature.seam_mean_curvature == 0.0
        assert abs(curvature.total_integrated_mean_curvature - math.pi) < 1e-9
        assert abs(curvature.smooth_gaussian_curvature - math.pi) < 1e-9
        assert abs(curvature.completed_gaussian_curvature - 2.0 * math.pi) < 1e-9
        assert curvature.gaussian_status == "CERTIFIED"

    # Tightening the geometric tolerance does not change the analytic cap
    # integrals or their Gauss--Bonnet closure.
    coarse = compute_molecular_contact_surface_curvature(
        side_a,
        GeometryTolerance(angular=1e-7, area=1e-8),
    )
    tight = compute_molecular_contact_surface_curvature(
        side_a,
        GeometryTolerance(angular=1e-11, area=1e-12),
    )
    assert abs(coarse.completed_gaussian_curvature - tight.completed_gaussian_curvature) < 1e-12

    # Unequal spheres: R=2, partner R=1, d=2 gives A-cap area pi and
    # smooth H=pi/2, smooth K=pi/4.
    unequal = (ball(0, "A", (0, 0, 0), 2.0), ball(1, "B", (2, 0, 0), 1.0))
    surface, _ = build_molecular_contact_surface_pair(unequal, (contact(0, 0, 1),))
    curvature = surface.curvature
    assert curvature is not None
    assert abs(surface.total_area_A2 - math.pi) < 1e-9
    assert abs(curvature.smooth_mean_curvature - math.pi / 2.0) < 1e-9
    assert abs(curvature.smooth_gaussian_curvature - math.pi / 4.0) < 1e-9

    # Two selected partner caps are unioned before integration.
    two_partners = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "B", (2.0, 0, 0), 2.0),
        ball(2, "B", (1.0, math.sqrt(3.0), 0), 2.0),
    )
    unioned = build_molecular_contact_surface(
        two_partners,
        (contact(0, 0, 1), contact(1, 0, 2)),
        "A",
    )
    assert unioned.curvature is not None
    assert unioned.total_area_A2 < unioned.naive_pairwise_area_A2

    # Same-group seam: mean seam curvature is present, while the intrinsic
    # seam/junction Gaussian closure remains explicitly partial until its
    # signed link accounting is certified.
    seam_case = (
        ball(0, "A", (0.0, 0.0, 0.0), 1.5),
        ball(1, "A", (0.8, 0.0, 0.0), 1.5),
        ball(2, "B", (0.4, 1.0, 0.0), 1.5),
    )
    seam_surface = build_molecular_contact_surface(
        seam_case,
        (contact(0, 0, 2), contact(1, 1, 2)),
        "A",
    )
    seam_curvature = seam_surface.curvature
    assert seam_curvature is not None
    assert seam_curvature.physical_seam_count == 1
    assert seam_curvature.seam_mean_curvature is not None
    assert seam_curvature.intrinsic_seam_gaussian_curvature is not None
    assert seam_curvature.gaussian_status == "CERTIFIED"
    assert abs(seam_curvature.completed_gaussian_curvature - 2.0 * math.pi) < 1e-9

    # Splitting one canonical seam interval changes only its cell subdivision.
    # The integrated physical seam term must be additive.
    seam_piece = next(piece for piece in seam_surface.refined_seams if piece.classification == "internal")
    begin, end = seam_piece.parameter_interval
    middle = 0.5 * (begin + end)
    midpoint = tuple((a + b) * 0.5 for a, b in zip(seam_piece.start_xyz, seam_piece.end_xyz))
    first_piece = replace(
        seam_piece,
        seam_id=seam_piece.seam_id + ":left",
        parameter_interval=(begin, middle),
        endpoint_junction_ids=(seam_piece.endpoint_junction_ids[0], "synthetic-mid"),
        end_xyz=midpoint,
    )
    second_piece = replace(
        seam_piece,
        seam_id=seam_piece.seam_id + ":right",
        parameter_interval=(middle, end),
        endpoint_junction_ids=("synthetic-mid", seam_piece.endpoint_junction_ids[1]),
        start_xyz=midpoint,
    )
    refined = replace(
        seam_surface,
        refined_seams=tuple(
            piece
            for piece in seam_surface.refined_seams
            if piece.seam_id != seam_piece.seam_id
        )
        + (first_piece, second_piece),
    )
    refined_curvature = compute_molecular_contact_surface_curvature(refined)
    assert abs(refined_curvature.total_integrated_mean_curvature - seam_curvature.total_integrated_mean_curvature) < 1e-12

    # Complete sphere: set the independently known topology χ=2 on the
    # geometry produced by the full-cover partner case.
    full_sphere = build_molecular_contact_surface(
        (ball(0, "A", (0, 0, 0), 1.0), ball(1, "B", (0, 0, 0), 10.0)),
        (contact(0, 0, 1),),
        "A",
    )
    full_sphere = replace(full_sphere, euler_characteristic=2)
    full_curvature = compute_molecular_contact_surface_curvature(full_sphere)
    assert abs(full_curvature.smooth_gaussian_curvature - 4.0 * math.pi) < 1e-9
    assert abs(full_curvature.completed_gaussian_curvature - 4.0 * math.pi) < 1e-9

    # Two intersecting carrier spheres forming a closed union boundary.  The
    # full seam has opposite induced boundary orientations on the two caps;
    # its intrinsic contribution is -2*pi, independently of topology.
    union_boundary = build_molecular_contact_surface(
        (
            ball(0, "A", (0, 0, 0), 1.0),
            ball(1, "A", (1, 0, 0), 1.0),
            ball(2, "B", (0.5, 0, 0), 10.0),
        ),
        (contact(0, 0, 2), contact(1, 1, 2)),
        "A",
    )
    first_patch, second_patch = union_boundary.patches
    seam = union_boundary.refined_seams[0]
    seam = replace(
        seam,
        orientation_on_patches=(
            (first_patch.patch_id, -1),
            (second_patch.patch_id, 1),
        ),
    )
    union_boundary = replace(
        union_boundary,
        refined_seams=(seam,),
        junctions=(),
        euler_characteristic=2,
    )
    union_curvature = compute_molecular_contact_surface_curvature(union_boundary)
    assert abs(union_curvature.intrinsic_seam_gaussian_curvature + 2.0 * math.pi) < 1e-9
    assert abs(union_curvature.completed_gaussian_curvature - 4.0 * math.pi) < 1e-9

    # Topology-only records are intentionally absent from curvature inputs.
    patch = seam_surface.patches[0]
    fake_cut = TopologicalCut(
        "synthetic-cut",
        patch.patch_id,
        patch.carrier_atom_id,
        (0, 1),
        ("cut-start", "cut-end"),
        (0.0, 0.0, patch.carrier_radius),
        (0.0, 0.0, -patch.carrier_radius),
    )
    with_cut = replace(seam_surface, topological_cuts=(fake_cut,))
    with_cut_curvature = compute_molecular_contact_surface_curvature(with_cut)
    assert asdict(with_cut_curvature) == asdict(seam_curvature)

    print("molecular contact curvature synthetic checks: PASS")


if __name__ == "__main__":
    main()
