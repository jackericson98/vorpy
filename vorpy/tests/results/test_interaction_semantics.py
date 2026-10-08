import pytest

from results.model import (
    INTERACTION_COMPARISON_QUANTITIES,
    InteractionLinks,
    OrientationConvention,
)


def test_interaction_links_preserve_triplet_roles_and_provenance():
    links = InteractionLinks(
        "interaction-1",
        result_ids={
            "molecular_contact_surface:A": "result-a",
            "physical_partition_interface": "result-p",
            "molecular_contact_surface:B": "result-b",
        },
        contact_selection_source="power_alpha0_interface",
        contact_selection_id="contact-selection-1",
        orientation_conventions={
            "molecular_contact_surface:A": "outward_molecular_normal",
            "physical_partition_interface": "physical_interface_A_to_B",
            "molecular_contact_surface:B": "outward_molecular_normal",
        },
    )

    assert links.result_ids["physical_partition_interface"] == "result-p"
    assert links.contact_selection_id == "contact-selection-1"
    assert links.orientation_conventions["molecular_contact_surface:A"] == (
        OrientationConvention.MOLECULAR_OUTWARD_NORMAL.value
    )
    assert set(INTERACTION_COMPARISON_QUANTITIES) == {
        "area_A_contact",
        "area_partition_interface",
        "area_B_contact",
        "area_ratio_A_to_partition",
        "area_ratio_B_to_partition",
        "area_ratio_A_to_B",
        "integrated_mean_curvature_A",
        "integrated_mean_curvature_partition",
        "integrated_mean_curvature_B",
        "integrated_gaussian_curvature_A",
        "integrated_gaussian_curvature_partition",
        "integrated_gaussian_curvature_B",
    }


def test_interaction_links_reject_wrong_partition_orientation():
    with pytest.raises(ValueError):
        InteractionLinks(
            "interaction-1",
            orientation_conventions={
                "physical_partition_interface": "outward_molecular_normal"
            },
        )


def test_interaction_links_require_complete_selection_provenance():
    with pytest.raises(ValueError):
        InteractionLinks("interaction-1", contact_selection_id="selection-1")
