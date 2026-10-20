from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.interface.export import summarize_interface_surface_areas
from vorpy.tests.geometry.test_alpha_interfaces import _network


@pytest.mark.parametrize("scheme", ["aw", "pow", "prm"])
def test_constructed_area_ledger_separates_direct_and_network_area(scheme):
    network = _network(scheme)
    network.surfs = pd.concat([
        network.surfs,
        pd.DataFrame(
            [{"balls": [10, 99], "sa": 6.0}],
            index=[23],
        ),
    ])
    interface = SimpleNamespace(
        net=network,
        group1_indices={10},
        group2_indices={13, 16},
    )

    ledger = summarize_interface_surface_areas(interface)

    assert ledger["direct_facet_count"] == 2
    assert ledger["direct_interface_area_A2"] == pytest.approx(7.0)
    assert ledger["internal_supporting_network_area_A2"] == pytest.approx(11.0)
    assert ledger["total_constructed_network_area_A2"] == pytest.approx(18.0)
    assert ledger["known_supported_partial_area_A2"] == pytest.approx(7.0)
    assert ledger["unresolved_direct_area_A2"] == pytest.approx(0.0)
    assert ledger["status"] == "VERIFIED"


def test_partial_direct_area_keeps_supported_and_unresolved_area_separate():
    network = _network("pow")
    network.surfs.loc[21, "complete"] = False
    interface = SimpleNamespace(
        net=network,
        group1_indices={10},
        group2_indices={13, 16},
    )

    ledger = summarize_interface_surface_areas(interface)

    assert ledger["direct_interface_area_A2"] == pytest.approx(7.0)
    assert ledger["known_supported_partial_area_A2"] == pytest.approx(3.0)
    assert ledger["unresolved_direct_area_A2"] == pytest.approx(4.0)
    assert ledger["status"] == "PARTIAL"
    assert ledger["supported_direct_facet_count"] == 1
    assert ledger["unresolved_direct_facet_count"] == 1
