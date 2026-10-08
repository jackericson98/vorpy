"""Small dependency-free regression check for molecular contact geometry."""

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


def ball(index, group, xyz, radius):
    return AtomBall(index, f"{group}-{index}", group, tuple(xyz), radius)


def contact(index, a, b):
    return SelectedContact(f"contact-{index}", a, b, "network", "interface")


def main():
    equal = (
        ball(0, "A", (0, 0, 0), 1.0),
        ball(1, "B", (1, 0, 0), 1.0),
    )
    side_a, side_b = build_molecular_contact_surface_pair(equal, (contact(0, 0, 1),))
    assert abs(side_a.total_area_A2 - math.pi) < 1e-9
    assert abs(side_b.total_area_A2 - math.pi) < 1e-9
    assert side_a.status == side_b.status == "CERTIFIED"

    unequal = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "B", (2, 0, 0), 1.0),
    )
    side_a, side_b = build_molecular_contact_surface_pair(unequal, (contact(0, 0, 1),))
    assert abs(side_a.total_area_A2 - math.pi) < 1e-9
    assert abs(side_b.total_area_A2 - 1.5 * math.pi) < 1e-9

    contained = (
        ball(0, "A", (0, 0, 0), 1.0),
        ball(1, "B", (1, 0, 0), 3.0),
    )
    side_a, side_b = build_molecular_contact_surface_pair(contained, (contact(0, 0, 1),))
    assert abs(side_a.total_area_A2 - 4.0 * math.pi) < 1e-9
    assert abs(side_b.total_area_A2) < 1e-12

    two_partners = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "B", (2.0, 0, 0), 2.0),
        ball(2, "B", (1.0, math.sqrt(3.0), 0), 2.0),
    )
    surface = build_molecular_contact_surface(two_partners, (contact(0, 0, 1), contact(1, 0, 2)), "A")
    assert surface.total_area_A2 < surface.naive_pairwise_area_A2
    assert len(surface.patches) == 1

    occluded = (
        ball(0, "A", (0, 0, 0), 2.0),
        ball(1, "A", (1.0, 0, 0), 2.0),
        ball(2, "B", (2.4, 0, 1.8), 2.0),
    )
    surface = build_molecular_contact_surface(occluded, (contact(0, 0, 2),), "A")
    assert surface.total_area_A2 < surface.naive_pairwise_area_A2
    assert "A-1" in {atom for patch in surface.patches for atom in patch.self_occluding_atom_ids}

    exposed = (
        ball(0, "A", (0, 0, 0), 1.0),
        ball(1, "A", (1.0, 0, 0), 1.0),
        ball(2, "B", (0, 0, 0), 3.0),
    )
    surface = build_molecular_contact_surface(exposed, (contact(0, 0, 2),), "A")
    assert abs(surface.total_area_A2 - 3.0 * math.pi) < 1e-9

    print("molecular contact synthetic checks: PASS")


if __name__ == "__main__":
    main()
