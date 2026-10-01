import numpy as np
import pytest

from vorpy.src.analyze.condition_b_geometry import (
    common_power_squared_radius,
    condition_b_decision,
)


def _pair_for_m_squared(q):
    # Equal unit balls centered two*sqrt(q+1) apart have the origin on their
    # bisector with common power q (when q >= -1).
    distance = 2.0 * np.sqrt(q + 1.0)
    return np.array([-distance/2, 0., 0.]), np.array([distance/2, 0., 0.])


@pytest.mark.parametrize("q", [2.25, 0.0, -0.25])
def test_two_ball_common_power_distinguishes_positive_zero_and_negative(q):
    pi, pj = _pair_for_m_squared(q)
    observed = common_power_squared_radius((0., 0., 0.), pi, 1., pj, 1.)
    assert observed == pytest.approx(q)
    decision = condition_b_decision([observed], 1., 5.)
    if q > 0:
        assert decision["m"] == pytest.approx(np.sqrt(q))
        assert decision["status"] == "positive_radius_candidate"
    elif q == 0:
        assert decision["m"] == 0.0
        assert decision["status"] == "zero_radius_candidate"
    else:
        assert decision["m"] is None
        assert decision["m_over_r"] is None
        assert decision["status"] == "no_real_orthogonal_ball"
        assert decision["accepted"] is True


@pytest.mark.parametrize("factor,accepted", [(1-1e-8, True), (1+1e-8, False)])
def test_condition_b_threshold_uses_strict_m_over_r_greater_than_M(factor, accepted):
    radius, threshold = 1.7, 5.0
    m = radius * threshold * factor
    decision = condition_b_decision([m*m], radius, threshold)
    assert decision["m_over_r"] == pytest.approx(threshold*factor)
    assert decision["accepted"] is accepted


def test_condition_b_is_invariant_to_translation_rotation_swapping_and_scaling():
    q, ri, rj = 2.25, 1.4, 1.8
    # Construct a radical-plane point for unequal radii on the x axis.
    pi, pj = np.array([-2., 0., 0.]), np.array([2., 0., 0.])
    x = np.array([(ri*ri-rj*rj)/(2*4.), 0., 0.])
    # Shift the point in y to obtain the requested common power where possible.
    current = common_power_squared_radius(x, pi, ri, pj, rj)
    x[1] = np.sqrt(q-current)
    q_observed = common_power_squared_radius(x, pi, ri, pj, rj)
    base = condition_b_decision([q_observed], min(ri, rj), 5.0)

    theta = 0.73
    rotation = np.array([[np.cos(theta), -np.sin(theta), 0.],
                         [np.sin(theta), np.cos(theta), 0.],
                         [0., 0., 1.]])
    shift = np.array([11., -7., 3.])
    transformed = common_power_squared_radius(rotation @ x + shift,
                                               rotation @ pi + shift, ri,
                                               rotation @ pj + shift, rj)
    swapped = common_power_squared_radius(x, pj, rj, pi, ri)
    assert transformed == pytest.approx(q_observed)
    assert swapped == pytest.approx(q_observed)

    scale = 3.25
    scaled_q = common_power_squared_radius(scale*x, scale*pi, scale*ri,
                                           scale*pj, scale*rj)
    scaled = condition_b_decision([scaled_q], scale*min(ri, rj), 5.0)
    assert scaled_q == pytest.approx(scale*scale*q_observed)
    assert scaled["m_over_r"] == pytest.approx(base["m_over_r"])
    assert scaled["accepted"] is base["accepted"]
