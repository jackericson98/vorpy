from types import SimpleNamespace

import pytest

from vorpy.src.network.views import NetworkView, NetworkViewError
from vorpy.src.output.export_scope import export_scope
from vorpy.tests.network.test_views import _network


@pytest.fixture
def wet_aw_reference():
    """Small full-source AW reference with water competition and an exterior defect."""
    source = _network(scheme="aw")
    source.solvent_mode = "FULL_SOLVENT"
    source.settings["solvent_mode"] = "FULL_SOLVENT"

    # Local generators 2 and 3 are competing water balls.  Generator 3 is an
    # unrelated incomplete exterior cell and must not contaminate A/B scope.
    source.balls.loc[2, "stable_id"] = "water-1"
    source.balls.loc[3, "stable_id"] = "water-2"
    source.balls.loc[3, "complete"] = False

    # Surfaces 0 and 2 are the known A/B facets, with water-1 in their
    # full-source incidence.  Surface 1 is unrelated exterior geometry.
    source.surfs.at[0, "balls"] = [0, 1, 2]
    source.surfs.at[0, "sa"] = 2.5
    selected = source.surfs.iloc[0].copy()
    selected["sa"] = 3.5
    source.surfs.loc[2] = selected
    source.surfs.at[1, "balls"] = [0, 3]

    return source


def test_full_source_wet_aw_interface_scope_is_closed_and_selective(
    wet_aw_reference,
):
    source = wet_aw_reference
    view = NetworkView.for_interface(source, [10], [11])
    owner = SimpleNamespace(source_network=source, source_view=view)
    scope = export_scope(owner)

    assert len(source.balls) == 5
    assert {10, 11, 12, 13}.issubset(set(source.group))
    assert view.source_ids("surface") == (0, 2)
    assert view.completeness.proven
    assert 2 in view.support_generator_local_indices
    assert 3 not in view.support_generator_local_indices
    assert scope.surfaces == (0, 2)
    assert scope.edges == (0, 1, 2, 3)
    assert scope.vertices == (0, 1, 2, 3)
    assert scope.support_generators == (0, 1, 2, 4)
    assert len(source.surfs) == 3


def test_full_source_wet_aw_interface_fails_closed_on_selected_incomplete_cell(
    wet_aw_reference,
):
    source = wet_aw_reference
    source.balls.loc[2, "complete"] = False

    with pytest.raises(NetworkViewError, match="completeness"):
        NetworkView.for_interface(source, [10], [11])
