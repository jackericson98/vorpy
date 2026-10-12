from types import SimpleNamespace

import pandas as pd

from vorpy.src.boundary import BoundaryConfig, BoundaryGenerator
from vorpy.src.network.aw_geometry_reuse import (
    STATIC_FRAME_ID,
    _settings_for_identity,
    build_solve_context_key,
    build_solve_universe_key,
)
from vorpy.src.network.network import Network


def _network(
    *,
    scheme="aw",
    group=(10, 11),
    frame_index=4,
    locations=None,
    radii=None,
    boundary_config=None,
    boundary_indices=(),
    boundary_generators=(),
):
    locations = locations or [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    radii = radii or [1.0] * len(locations)
    boundary_indices = tuple(boundary_indices)
    balls = pd.DataFrame({
        "loc": locations,
        "rad": radii,
        "num": list(range(len(locations))),
        "system_num": [10, 11] + [None] * max(0, len(locations) - 2),
        "stable_id": ["atom-A", "atom-B"] + [None] * max(0, len(locations) - 2),
    })
    for index in boundary_indices:
        balls.loc[index, "is_boundary_generator"] = True
        balls.loc[index, "system_num"] = None
        balls.loc[index, "stable_id"] = None

    network = Network(
        locs=locations,
        rads=radii,
        group=list(group),
        settings={
            "net_type": scheme,
            "surf_res": 0.25,
            "max_vert": 5.0,
            "box_size": 1.5,
            "num_splits": 2,
            "build_type": "all",
        },
        balls=balls,
        system=SimpleNamespace(frame_index=frame_index),
    )
    network.boundary_config = boundary_config or BoundaryConfig(mode="none")
    network.boundary_indices = boundary_indices
    network.boundary_generators = tuple(boundary_generators)
    network.boundary_mode = network.boundary_config.mode.value
    return network


def test_identical_contexts_are_equal_and_preserve_system_mapping():
    first = _network()
    second = _network()

    first_key = build_solve_context_key(first)
    second_key = second.get_solve_context_key()

    assert first_key == second_key
    assert first_key.digest == second_key.digest
    assert first_key.scheme == "aw"
    assert first_key.frame_id == "4"
    assert first_key.physical_generator_ids == ("atom-A", "atom-B")
    assert first_key.local_to_system == ((0, 10), (1, 11))
    assert first_key.solve_universe_id == f"solve:{first_key.digest}"
    assert first.verts is None
    assert first.edges is None
    assert first.surfs is None


def test_settings_identity_normalizer_does_not_reference_network_state():
    settings = {
        "net_type": "aw",
        "max_vert": 5.0,
        "surf_col": "plasma",
        "print_metrics": True,
    }

    normalized = _settings_for_identity(settings)

    assert normalized == {
        "net_type": "aw",
        "max_vert": 5.0,
        "surf_col": "plasma",
    }
    assert build_solve_context_key(_network()).solve_universe_id.startswith("solve:")


def test_static_frame_sentinel_is_explicit_and_non_empty():
    key = build_solve_context_key(_network(), frame_id=STATIC_FRAME_ID)

    assert STATIC_FRAME_ID == "static"
    assert key.frame_id == STATIC_FRAME_ID
    assert key.frame_id


def test_context_changes_when_geometry_membership_scheme_or_frame_changes():
    baseline = build_solve_context_key(_network())

    changed = [
        _network(locations=[(0.1, 0.0, 0.0), (2.0, 0.0, 0.0)]),
        _network(radii=[1.1, 1.0]),
        _network(group=(10,)),
        _network(scheme="pow"),
        _network(frame_index=5),
        _network(),
    ]
    changed[-1].settings["surf_res"] = 0.5
    changed.append(_network())
    changed[-1].settings["max_vert"] = 9.0

    assert all(
        build_solve_context_key(network).digest != baseline.digest
        for network in changed
    )


def test_equal_value_environment_objects_do_not_use_object_identity():
    first = build_solve_context_key(
        _network(), environment=SimpleNamespace(mode="dry")
    )
    second = build_solve_context_key(
        _network(), environment=SimpleNamespace(mode="dry")
    )

    assert first == second


def test_explicit_build_options_are_part_of_context_identity():
    first = build_solve_context_key(
        _network(), build_options={"calculate_curvature": True}
    )
    second = build_solve_context_key(
        _network(), build_options={"calculate_curvature": False}
    )

    assert first.digest != second.digest
    assert first.build_options_digest != second.build_options_digest


def test_boundary_configuration_and_generators_are_identity_inputs():
    baseline = build_solve_context_key(_network())
    changed_config = build_solve_context_key(
        _network(boundary_config=BoundaryConfig(mode="shell"))
    )

    assert changed_config.digest != baseline.digest
    assert changed_config.boundary_config_digest != baseline.boundary_config_digest

    boundary = build_solve_context_key(
        _network(
            locations=[(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (4.0, 0.0, 0.0)],
            radii=[1.0, 1.0, 1.5],
            boundary_config=BoundaryConfig(mode="shell", shell_spacing=3.0),
            boundary_indices=(2,),
            boundary_generators=(
                BoundaryGenerator(position=(4.0, 0.0, 0.0), radius=1.5),
            ),
        )
    )

    assert boundary.digest != baseline.digest
    assert boundary.boundary_mode == "shell"
    assert boundary.boundary_generator_ids == ("boundary:2",)
    assert boundary.local_to_system == ((0, 10), (1, 11), (2, None))

    changed_generator = build_solve_context_key(
        _network(
            locations=[(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (4.0, 0.0, 0.0)],
            radii=[1.0, 1.0, 1.5],
            boundary_config=BoundaryConfig(mode="shell", shell_spacing=3.0),
            boundary_indices=(2,),
            boundary_generators=(
                BoundaryGenerator(position=(4.1, 0.0, 0.0), radius=1.5),
            ),
        )
    )

    assert changed_generator.digest != boundary.digest
    assert changed_generator.boundary_generator_digest != boundary.boundary_generator_digest


def test_generic_context_supports_power_and_primitive_but_legacy_aw_key_does_not():
    for scheme in ("pow", "prm"):
        key = build_solve_context_key(_network(scheme=scheme))
        assert key.scheme == scheme

        try:
            build_solve_universe_key(_network(scheme=scheme))
        except ValueError as error:
            assert "non-AW" in str(error)
        else:
            raise AssertionError("legacy AW key accepted a non-AW network")
