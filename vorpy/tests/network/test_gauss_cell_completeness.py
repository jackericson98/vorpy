import pytest

from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.interfaces import build_full_interface
from vorpy.tests.geometry.test_interfaces import _tetra_network
from vorpy.tests.network.test_geometry_completeness import _closed_center_network


def test_solvent_and_boundary_adjacency_do_not_invalidate_a_closed_physical_cell():
    network = _closed_center_network()
    network.balls.loc[0, "complete"] = True
    network.balls["is_boundary_generator"] = [False, False, False, False, True]
    network.boundary_indices = (4,)
    network.solvent_mode = "SOLVENT_CONTEXT"
    network.solvent_local_indices = (1,)

    cell = network.diagnose_geometry_completeness()["cells"][0]

    assert cell["complete"] is True
    assert cell["classification"] == "COMPLETE"
    assert cell["bounded"] is True
    assert cell["missing_faces"] == ()
    assert cell["missing_edges"] == ()
    assert cell["closure_relations"]["open_edges"] == ()


def test_topology_defect_overrides_stale_producer_complete_flag():
    network = _closed_center_network()
    network.balls.loc[0, "complete"] = True
    network.balls.at[0, "edges"] = [0, 1, 2, 3, 99]

    cell = network.diagnose_geometry_completeness()["cells"][0]

    assert cell["complete"] is True
    assert cell["classification"] == "INCOMPLETE_CONSTRUCTION"
    assert (0, 99) in cell["missing_edges"]


@pytest.mark.parametrize("scheme", ["prm", "pow", "aw"])
def test_dual_interface_gate_is_scheme_specific_and_read_only(scheme):
    network = _tetra_network(scheme)
    network.balls["complete"] = [False] * len(network.balls)
    before = tuple(bool(value) for value in network.balls["complete"])

    dual = build_dual(network)
    interface = build_full_interface(network, dual, {0, 1}, {2, 3})

    assert tuple(bool(value) for value in network.balls["complete"]) == before
    if scheme == "aw":
        assert interface.selected_surfaces == []
        assert {row.status for row in interface.candidate_surfaces} == {
            "excluded_incomplete_cells"
        }
    else:
        assert len(interface.selected_surfaces) == 4

