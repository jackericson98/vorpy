import os

import pandas as pd
import pytest

from vorpy.src.network.incremental import (
    NORMAL_BUILD_REQUIRED,
    VERIFIED_INCREMENTAL,
    IncrementalCompletionError,
    assess_incremental_group,
    complete_group_incrementally,
)
from vorpy.src.network.network import Network


LOCATIONS = [
    (0.0, 0.0, 0.0),
    (20.0, 0.0, 0.0),
    (-20.0, 0.0, 0.0),
    (0.0, 20.0, 0.0),
    (0.0, -20.0, 0.0),
    (0.0, 0.0, 20.0),
    (0.0, 0.0, -20.0),
]


def _network(group, scheme):
    settings = {
        "net_type": scheme,
        "surf_res": 0.25,
        "surf_scheme": "int_mean_curv",
        "max_vert": 40.0,
        "box_size": 1.5,
        "build_type": "all",
        "foam_box": None,
        "print_metrics": False,
        "sys_dir": os.getcwd(),
    }
    return Network(
        locs=LOCATIONS,
        rads=[1.0] * len(LOCATIONS),
        group=group,
        group_name="bounded-test",
        settings=settings,
        generator_system_indices=tuple(range(len(LOCATIONS))),
        generator_stable_ids=tuple(f"g{index}" for index in range(len(LOCATIONS))),
    )


@pytest.mark.parametrize("scheme", ("aw", "pow", "prm"))
def test_incremental_group_matches_clean_geometry_and_preserves_source(scheme):
    source = _network([0, 1], scheme)
    source.build()
    clean = _network([2], scheme)
    clean.build()
    source_vertices = source.verts.copy(deep=True)

    # Existing target-supported vertices are reused and topology is rebuilt.
    result = complete_group_incrementally(source, [2], clean_reference=clean)
    assert result.plan.decision == VERIFIED_INCREMENTAL
    assert "group_results" in result.plan.invalidated_caches
    assert result.comparison.equivalent
    assert len(result.network.verts) == len(clean.verts)
    pd.testing.assert_frame_equal(source.verts, source_vertices)

    # A damaged source copy with missing Group C vertices is still fail-closed
    # and exact when the bounded candidate search can recover them.
    source.verts = source.verts[
        ~source.verts["balls"].map(lambda balls: 2 in tuple(balls))
    ].reset_index(drop=True)
    recovered = complete_group_incrementally(source, [2], clean_reference=clean)
    assert recovered.plan.decision == VERIFIED_INCREMENTAL
    assert recovered.plan.missing_candidates
    assert recovered.comparison.equivalent


def test_incompatible_context_requires_explicit_normal_build():
    source = _network([0, 1], "aw")
    source.build()
    clean = _network([2], "aw")
    clean.build()
    clean.settings["surf_res"] = 0.5

    with pytest.raises(IncrementalCompletionError, match="context"):
        complete_group_incrementally(source, [2], clean_reference=clean)


def test_bounded_search_uses_explicit_normal_build_fallback():
    source = _network([0, 1], "aw")
    source.build()
    calls = []

    def normal_build(plan):
        calls.append(plan.decision)
        return "normal-build"

    result = complete_group_incrementally(
        source,
        [2],
        max_generators=6,
        normal_build=normal_build,
    )
    assert result.plan.decision == NORMAL_BUILD_REQUIRED
    assert result.network == "normal-build"
    assert calls == [NORMAL_BUILD_REQUIRED]


def test_unverified_candidates_fail_without_clean_reference():
    source = _network([0, 1], "prm")
    source.build()
    with pytest.raises(IncrementalCompletionError, match="clean_reference"):
        complete_group_incrementally(source, [2])
