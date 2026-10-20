from types import SimpleNamespace
from pathlib import Path

from vorpy.src.boundary import BoundaryConfig
from vorpy.src.group.layers import _build_index_maps
from vorpy.src.network.network import Network
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


EDTA = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"


def test_virtual_boundary_system_sentinel_is_not_a_system_mapping():
    network = Network(
        locs=[(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
        rads=[1.0, 1.0],
        generator_system_indices=(7, None),
        generator_stable_ids=("system:7", "boundary:0"),
        build_net=False,
    )
    network.boundary_config = BoundaryConfig(mode="shell")

    topology_to_row, topology_to_system, system_to_topology = _build_index_maps(
        SimpleNamespace(net=network)
    )

    assert topology_to_row == {0: 0, 1: 1}
    assert topology_to_system == {0: 7}
    assert system_to_topology == {7: 0}


def test_invalid_negative_system_reference_is_omitted_from_system_layer():
    network = Network(
        locs=[(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)],
        rads=[1.0, 1.0],
        balls=None,
        generator_system_indices=(0, None),
        generator_stable_ids=("system:0", "boundary:0"),
        build_net=False,
    )
    network.balls.loc[1, "system_num"] = -1

    _, topology_to_system, system_to_topology = _build_index_maps(
        SimpleNamespace(net=network)
    )

    assert topology_to_system == {0: 0}
    assert system_to_topology == {0: 0}


def test_shell_bounded_edta_keeps_boundary_only_in_topology_layers():
    result = VorPyBackend(
        VorPySolveSettings(network_type="aw", boundary_mode="shell")
    ).solve(
        EDTA,
        progress=lambda *args, **kwargs: None,
        is_cancelled=lambda: False,
    )
    group = result.export_group

    network = group.net
    boundary_topology_ids = {
        int(network.balls.iloc[index]["num"])
        for index in network.boundary_indices
    }
    layer_system_ids = {
        int(value)
        for layer in group.layer_atoms
        for value in layer
    }
    layer_topology_ids = {
        int(value)
        for layer in group.layer_net_atoms
        for value in layer
    }

    assert layer_system_ids
    assert min(layer_system_ids) >= 0
    assert -1 not in layer_system_ids
    assert boundary_topology_ids & layer_topology_ids
    assert boundary_topology_ids.isdisjoint(layer_system_ids)
