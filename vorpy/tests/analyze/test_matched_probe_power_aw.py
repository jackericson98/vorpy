import numpy as np
import pytest

from vorpy.src.analyze.matched_probe_power_aw import (
    bicolor_pairs, pair_partition, validate_matched_inputs,
)


def test_matched_probe_input_identity_and_effective_radii():
    ids = [("A", 1, "", "GLY", "CA"), ("I", 2, "A", "LYS", "NZ")]
    xyz = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    radii = np.array([1.7, 1.9])
    assert validate_matched_inputs(ids, ids.copy(), xyz, xyz.copy(), radii, radii.copy())


@pytest.mark.parametrize("which", ["ids", "coordinates", "base_radii"])
def test_matched_probe_inputs_fail_before_geometry_on_mismatch(which):
    ids = [("A", 1, "", "GLY", "CA")]
    xyz = np.array([[0.0, 0.0, 0.0]])
    radii = np.array([1.8])
    values = [ids, ids.copy(), xyz, xyz.copy(), radii, radii.copy()]
    if which == "ids":
        values[1] = [("B", 1, "", "GLY", "CA")]
    elif which == "coordinates":
        values[3] = np.array([[0.0, 0.0, 0.1]])
    else:
        values[5] = np.array([1.81])
    with pytest.raises(AssertionError):
        validate_matched_inputs(*values)


def test_pair_partition_is_exact_and_disjoint():
    shared, power_only, aw_only = pair_partition({(1, 2), (2, 3)}, {(2, 1), (4, 5)})
    assert shared == {(1, 2)}
    assert power_only == {(2, 3)}
    assert aw_only == {(4, 5)}
    assert shared | power_only | aw_only == {(1, 2), (2, 3), (4, 5)}


def test_probe_is_not_double_counted():
    base = np.array([1.5, 1.8, 2.0])
    power_expanded = base + 1.4
    aw_input_radii = base.copy()
    aw_expanded_effective = aw_input_radii + 1.4
    assert np.array_equal(power_expanded, aw_expanded_effective)
    assert not np.array_equal(aw_input_radii + 2.8, power_expanded)


def test_group_reversal_preserves_bicolor_pair_membership():
    full = {(0, 1), (0, 2), (1, 3), (2, 3)}
    selected = {(0, 1), (0, 3)}
    first = bicolor_pairs(selected, {0, 2}, {1, 3})
    reversed_groups = bicolor_pairs(selected, {1, 3}, {0, 2})
    assert first == reversed_groups == {(0, 1), (0, 3)}
    assert bicolor_pairs(full, {0, 2}, {1, 3}) == bicolor_pairs(full, {1, 3}, {0, 2})


def test_rigid_transform_preserves_connection_distances():
    points = np.array([[0., 1., 2.], [-3., 4., 5.]])
    q, _ = np.linalg.qr(np.array([[1., 2., 3.], [4., 2., 1.], [2., 3., 5.]]))
    moved = points @ q + np.array([19., -8., 2.])
    assert np.isclose(np.linalg.norm(points[0] - points[1]),
                      np.linalg.norm(moved[0] - moved[1]), atol=1e-12)
