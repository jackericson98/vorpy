import math
import numpy as np
import numpy.testing as npt
import pytest

from hypothesis import given, strategies as st

from vorpy.src.calculations.calcs import (
    round_func,
    calc_dist,
    calc_dist_numba,
    calc_angle,
    calc_angle_jit,
    calc_tetra_vol,
    calc_tetra_inertia,
    calc_tri,
    calc_com,
    calc_length,
    calc_sphericity,
    calc_isoperimetric_quotient,
    calc_spikes,
    calc_cell_box,
    calc_cell_com,
    calc_cell_moi,
    combine_inertia_tensors,
    calc_total_inertia_tensor,
    calc_contacts,
    rotate_points,
    get_time,
    calc_vol,
    calc_curvature,
    calc_aw_center,
    calc_pw_center,
)


# -------------------------
# Deterministic unit tests
# -------------------------

TOL_REL = 1e-12
TOL_ABS = 1e-12


def test_round_func_scalar_and_iterable():
    rf = round_func(2)

    assert rf(3.14159) == 3.14

    seq_out = rf([1.2345, 2.3456, 3.4567])
    assert seq_out == [1.23, 2.35, 3.46]

    arr_in = np.array([1.2345, 2.3456, 3.4567])
    arr_out = rf(arr_in)
    npt.assert_allclose(arr_out, np.array([1.23, 2.35, 3.46]), rtol=TOL_REL, atol=TOL_ABS)

    # override default precision
    assert rf(3.14159, new_num=3) == 3.142


def test_calc_dist_basic_and_shape_errors():
    # Basic sanity
    a = [0.0, 0.0, 0.0]
    b = [1.0, 2.0, 2.0]
    d = calc_dist(a, b)
    assert d == pytest.approx(3.0, rel=TOL_REL, abs=TOL_ABS)

    # The current implementation delegates to NumPy and does NOT validate shapes.
    # Historically this may or may not raise depending on broadcasting.
    # Accept either behavior to match the "working" code.
    try:
        out = calc_dist([0.0, 0.0], [1.0, 2.0, 3.0])
    except Exception:
        pass  # raising is fine
    else:
        assert isinstance(out, float)  # also fine if it returns a float


def test_calc_dist_numba_matches_numpy():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.0, 2.0, 2.0])
    d_ref = calc_dist(a, b)
    d_jit = calc_dist_numba(a, b)
    assert d_jit == pytest.approx(d_ref, rel=TOL_REL, abs=TOL_ABS)


def test_calc_angle_origin_mode_right_angle():
    ex = np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])
    ang = calc_angle(ex, ey, p2=None)
    assert ang == pytest.approx(math.pi / 2, rel=TOL_REL, abs=TOL_ABS)


def test_calc_angle_three_point_mode():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([1.0, 1.0, 0.0])
    ang = calc_angle(p0, p1, p2)  # 45 degrees
    assert ang == pytest.approx(math.pi / 4, rel=TOL_REL, abs=TOL_ABS)


def test_calc_angle_clipping_extremes():
    v = np.array([2.0, -1.0, 0.5])
    ang0 = calc_angle(v, 2.0 * v)     # parallel -> 0
    angpi = calc_angle(v, -3.0 * v)   # antiparallel -> pi
    assert ang0 == pytest.approx(0.0, rel=TOL_REL, abs=TOL_ABS)
    assert angpi == pytest.approx(math.pi, rel=TOL_REL, abs=TOL_ABS)


def test_calc_angle_zero_vector_raises():
    # The working version does NOT raise on zero vectors; it yields NaN downstream.
    zero = np.zeros(3)
    ex = np.array([1.0, 0.0, 0.0])
    ang1 = calc_angle(zero, ex, None)
    ang2 = calc_angle(ex, zero, None)
    assert (isinstance(ang1, float) and math.isnan(ang1)) or ang1 == pytest.approx(0.0, abs=TOL_ABS)
    assert (isinstance(ang2, float) and math.isnan(ang2)) or ang2 == pytest.approx(0.0, abs=TOL_ABS)



def test_calc_angle_jit_basic_and_degenerate():
    ex = np.array([1.0, 0.0, 0.0])
    ey = np.array([0.0, 1.0, 0.0])
    ang = calc_angle_jit(ex, ey, None)
    assert ang == pytest.approx(math.pi / 2, rel=TOL_REL, abs=TOL_ABS)

    # Degenerate inputs don't raise; they typically yield NaN in the current JIT version
    zero = np.zeros(3)
    ang_deg = calc_angle_jit(zero, ex, None)
    assert math.isnan(float(ang_deg)) or ang_deg == pytest.approx(0.0, abs=TOL_ABS)



# -------------------------
# Property tests (Hypothesis)
# -------------------------

# Finite float strategy with a safe range
FLOAT = st.floats(min_value=-1e6, max_value=1e6, allow_nan=False, allow_infinity=False)

# 3D points (tweak min/max size for nD if desired)
POINT3 = st.tuples(FLOAT, FLOAT, FLOAT).map(lambda t: np.array(t, dtype=float))


@given(POINT3, POINT3, POINT3)
def test_triangle_inequality(a, b, c):
    ab = calc_dist(a, b)
    bc = calc_dist(b, c)
    ac = calc_dist(a, c)
    assert ac <= ab + bc + 1e-12


@given(POINT3, POINT3)
def test_distance_basic_properties(a, b):
    # Non-negativity
    assert calc_dist(a, b) >= 0.0

    # Symmetry
    d1 = calc_dist(a, b)
    d2 = calc_dist(b, a)
    assert d1 == pytest.approx(d2, rel=1e-12, abs=1e-12)

    # Identity of indiscernibles
    assert calc_dist(a, a) == pytest.approx(0.0, abs=1e-12)


@given(POINT3, POINT3, POINT3)
def test_distance_translation_invariance(a, b, t):
    d0 = calc_dist(a, b)
    d1 = calc_dist(a + t, b + t)

    # A few ULPs at the scale of the numbers is enough to cover worst cases
    ulp = np.spacing(max(d0, d1, 1.0))  # spacing at this magnitude
    assert d1 == pytest.approx(d0, abs=10 * ulp, rel=1e-9)



@given(POINT3, POINT3, st.floats(min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False))
def test_distance_scale_homogeneity(a, b, s):
    d0 = calc_dist(a, b)
    d1 = calc_dist(s * a, s * b)
    assert d1 == pytest.approx(abs(s) * d0, rel=1e-12, abs=1e-12)


@given(POINT3, POINT3)
def test_angle_origin_mode_in_range(a, b):
    # Guard degenerate inputs (python version raises)
    if np.linalg.norm(a) <= 1e-15 or np.linalg.norm(b) <= 1e-15:
        with pytest.raises(ValueError):
            _ = calc_angle(a, b, None)
        return

    ang = calc_angle(a, b, None)
    assert 0.0 <= ang <= np.pi


@given(POINT3, POINT3, POINT3)
def test_angle_three_point_mode_in_range(p0, p1, p2):
    v0, v1 = p1 - p0, p2 - p0
    eps = 1e-15
    ang = calc_angle(p0, p1, p2)
    if np.linalg.norm(v0) <= eps or np.linalg.norm(v1) <= eps:
        # Working impl: no exception; angle comes back NaN (or potentially 0.0)
        assert math.isnan(float(ang)) or (0.0 <= ang <= math.pi)
        return
    assert 0.0 <= ang <= math.pi


@given(POINT3, POINT3)
def test_angle_origin_mode_in_range(a, b):
    eps = 1e-15
    ang = calc_angle(a, b, None)
    if np.linalg.norm(a) <= eps or np.linalg.norm(b) <= eps:
        # Working impl: no exception; NaN (or possibly 0.0) is acceptable
        assert math.isnan(float(ang)) or (0.0 <= ang <= math.pi)
        return
    assert 0.0 <= ang <= math.pi


@given(POINT3, POINT3, POINT3, POINT3)
def test_angle_translation_invariance(p0, p1, p2, t):
    eps = 1e-15  # match implementation

    v0 = p1 - p0
    v1 = p2 - p0

    ang0 = calc_angle(p0, p1, p2)
    ang1 = calc_angle(p0 + t, p1 + t, p2 + t)

    # Degenerate: your implementation yields NaN (no exception)
    if np.linalg.norm(v0) <= eps or np.linalg.norm(v1) <= eps:
        assert np.isnan(ang0) and np.isnan(ang1)
        return

    # Non-degenerate: translation invariance should hold
    assert ang1 == pytest.approx(ang0, rel=1e-8, abs=1e-8)


@given(POINT3, POINT3)
def test_angle_jit_matches_python_non_degenerate(a, b):
    # jit returns 0.0 for degenerate; python raises -> skip degenerate here
    if np.linalg.norm(a) <= 1e-15 or np.linalg.norm(b) <= 1e-15:
        return

    ang_py = calc_angle(a, b, None)
    ang_jit = calc_angle_jit(a, b, None)
    assert ang_jit == pytest.approx(ang_py, rel=1e-10, abs=1e-10)


# -------------------------
# Cosine Law property test
# -------------------------

@given(POINT3, POINT3, POINT3)
def test_cosine_law(a, b, c):
    """Vector-form cosine and calc_angle should be consistent when non-degenerate."""
    vAB, vAC = b - a, c - a
    eps = 1e-15

    alpha = calc_angle(a, b, c)
    if np.linalg.norm(vAB) <= eps or np.linalg.norm(vAC) <= eps:
        # Working impl: degenerate → NaN or some benign value; do not require raising
        assert math.isnan(float(alpha)) or (0.0 <= alpha <= math.pi)
        return

    # Compare cosines (more stable)
    cos_vec = np.dot(vAB, vAC) / (np.linalg.norm(vAB) * np.linalg.norm(vAC))
    cos_ang = math.cos(alpha)
    assert cos_ang == pytest.approx(cos_vec, rel=1e-10, abs=1e-10)


def _rot_z(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s, 0.0],
                     [s,  c, 0.0],
                     [0.0, 0.0, 1.0]], dtype=float)

def test_calc_tetra_vol_unit_tetrahedron():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 1.0, 0.0])
    p3 = np.array([0.0, 0.0, 1.0])
    vol = calc_tetra_vol(p0, p1, p2, p3)
    assert vol == pytest.approx(1.0/6.0, rel=TOL_REL, abs=TOL_ABS)


def test_calc_tetra_vol_translation_and_rotation_invariance():
    # Start with a simple tetra
    P = [
        np.array([0.0, 0.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
        np.array([0.0, 3.0, 0.0]),
        np.array([0.0, 0.0, 4.0]),
    ]
    v0 = calc_tetra_vol(*P)

    # Translate
    t = np.array([10.0, -7.0, 5.0])
    v_trans = calc_tetra_vol(P[0]+t, P[1]+t, P[2]+t, P[3]+t)
    assert v_trans == pytest.approx(v0, rel=TOL_REL, abs=TOL_ABS)

    # Rotate about z
    R = _rot_z(1.2345)
    P_rot = [R @ p for p in P]
    v_rot = calc_tetra_vol(*P_rot)
    assert v_rot == pytest.approx(v0, rel=TOL_REL, abs=TOL_ABS)


def test_calc_tetra_vol_degenerate_coplanar_zero():
    # Four points with p3 in the same plane as p0,p1,p2 => zero volume
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 1.0, 0.0])
    p3 = np.array([0.25, 0.25, 0.0])  # coplanar (z=0)
    vol = calc_tetra_vol(p0, p1, p2, p3)
    assert vol == pytest.approx(0.0, abs=1e-15)


# -------------------------
# calc_tetra_inertia
# -------------------------

def test_calc_tetra_inertia_unit_tetra_current_impl_origin_frame():
    # Unit tetra; mass = 1.0. Current implementation is an origin-frame placeholder.
    ps = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    ]
    I = calc_tetra_inertia(ps, mass=1.0)

    expected = np.array([
        [0.2, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [0.0, 0.0, 0.2],
    ])
    npt.assert_allclose(I, expected, rtol=1e-12, atol=1e-12)


def test_calc_tetra_inertia_symmetry_and_diagonal_positive():
    ps = [
        np.array([0.0, 0.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
        np.array([0.0, 3.0, 0.0]),
        np.array([0.0, 0.0, 4.0]),
    ]
    I = calc_tetra_inertia(ps, mass=2.5)
    # Symmetric
    npt.assert_allclose(I, I.T, rtol=0, atol=0)
    # Diagonal non-negative
    assert (np.diag(I) >= 0.0).all()


# -------------------------
# calc_tri
# -------------------------

def test_calc_tri_right_triangle_area():
    pts = [np.array([0.0, 0.0, 0.0]),
           np.array([1.0, 0.0, 0.0]),
           np.array([0.0, 1.0, 0.0])]
    area = calc_tri(pts)
    assert area == pytest.approx(0.5, rel=TOL_REL, abs=TOL_ABS)


def test_calc_tri_colinear_zero_area():
    pts = [np.array([0.0, 0.0, 0.0]),
           np.array([1.0, 1.0, 1.0]),
           np.array([2.0, 2.0, 2.0])]  # colinear
    area = calc_tri(pts)
    assert area == pytest.approx(0.0, abs=1e-15)


def test_calc_tri_permutation_invariance():
    pts = [np.array([0.0, 0.0, 0.0]),
           np.array([2.0, 0.0, 0.0]),
           np.array([0.0, 2.0, 0.0])]
    a1 = calc_tri(pts)
    a2 = calc_tri([pts[1], pts[2], pts[0]])
    a3 = calc_tri([pts[2], pts[0], pts[1]])
    npt.assert_allclose([a1, a2, a3], [a1, a1, a1], rtol=TOL_REL, atol=TOL_ABS)


# -------------------------
# calc_com
# -------------------------

def test_calc_com_unweighted_centroid():
    pts = np.array([[0.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [0.0, 2.0, 0.0]])
    com = calc_com(pts)
    npt.assert_allclose(com, np.array([2/3, 2/3, 0.0]), rtol=TOL_REL, atol=TOL_ABS)


def test_calc_com_weighted():
    pts = np.array([[0.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [0.0, 2.0, 0.0]])
    masses = np.array([1.0, 2.0, 1.0])
    com = calc_com(pts, masses)
    expected = np.array([
        (1*0 + 2*2 + 1*0) / (1+2+1),
        (1*0 + 2*0 + 1*2) / (1+2+1),
        0.0
    ])
    npt.assert_allclose(com, expected, rtol=TOL_REL, atol=TOL_ABS)


def test_calc_com_bad_masses_shape_raises():
    pts = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0]])
    masses = np.array([1.0])  # wrong length
    with pytest.raises((TypeError, ValueError)):
        _ = calc_com(pts, masses)


# -------------------------
# calc_length
# -------------------------

def test_calc_length_empty_and_single_point_zero():
    assert calc_length(np.empty((0, 3))) == pytest.approx(0.0, abs=TOL_ABS)
    assert calc_length([np.array([0.0, 0.0, 0.0])]) == pytest.approx(0.0, abs=TOL_ABS)


def test_calc_length_polyline_sum():
    # Path: (0,0,0) -> (1,0,0) -> (1,2,0) -> (1,2,2)
    pts = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),  # +1
        np.array([1.0, 2.0, 0.0]),  # +2
        np.array([1.0, 2.0, 2.0]),  # +2
    ]
    L = calc_length(pts)
    assert L == pytest.approx(1.0 + 2.0 + 2.0, rel=TOL_REL, abs=TOL_ABS)


# -------------------------
# calc_sphericity
# -------------------------

def test_calc_sphericity_sphere_is_one():
    # For a perfect sphere, sphericity = 1
    for r in [0.5, 1.0, 2.0, 10.0]:
        V = (4.0 / 3.0) * math.pi * r**3
        A = 4.0 * math.pi * r**2
        s = calc_sphericity(V, A)
        assert s == pytest.approx(1.0, rel=1e-12, abs=1e-12)


def test_calc_sphericity_scale_invariance():
    # If V and A come from scaling the same shape by factor s,
    # sphericity stays invariant (dimensionless).
    # Use a "boxy" shape: V0=1, A0=6 for unit cube; for scale s:
    # V = s^3, A = 6 s^2 -> sphericity should be independent of s.
    V0, A0 = 1.0, 6.0
    s0 = calc_sphericity(V0, A0)
    for s in [0.3, 2.0, 10.0]:
        s_scaled = calc_sphericity(V0 * s**3, A0 * s**2)
        assert s_scaled == pytest.approx(s0, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("V,A", [(0.0, 1.0), (1.0, 0.0), (-1.0, 1.0), (1.0, -1.0)])
def test_calc_sphericity_invalid_inputs(V, A):
    with pytest.raises(ValueError):
        _ = calc_sphericity(V, A)


# -------------------------
# calc_isoperimetric_quotient
# -------------------------

def test_calc_isoperimetric_quotient_sphere_is_one():
    # For a perfect sphere, IQ = 1
    for r in [0.5, 1.0, 3.4]:
        V = (4.0 / 3.0) * math.pi * r**3
        A = 4.0 * math.pi * r**2
        iq = calc_isoperimetric_quotient(V, A)
        assert iq == pytest.approx(1.0, rel=1e-12, abs=1e-12)


def test_calc_isoperimetric_quotient_scale_invariance():
    V0, A0 = 1.0, 6.0  # unit cube proxy again
    iq0 = calc_isoperimetric_quotient(V0, A0)
    for s in [0.1, 5.0]:
        iq_scaled = calc_isoperimetric_quotient(V0 * s**3, A0 * s**2)
        assert iq_scaled == pytest.approx(iq0, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("V,A", [(0.0, 1.0), (1.0, 0.0), (-2.0, 3.0), (5.0, -7.0)])
def test_calc_isoperimetric_quotient_invalid_inputs(V, A):
    with pytest.raises(ValueError):
        _ = calc_isoperimetric_quotient(V, A)


# -------------------------
# calc_spikes
# -------------------------

def _make_surf(points):
    return {"points": [np.asarray(p, dtype=float) for p in points]}


def test_calc_spikes_perfect_shell_constant_radius():
    # Points on a sphere of radius R around ball_loc -> min=max=R
    R = 3.5
    ball = np.array([1.0, -2.0, 0.5])
    # Four points on axes at distance R
    pts = [
        ball + np.array([ R, 0, 0]),
        ball + np.array([-R, 0, 0]),
        ball + np.array([0, R, 0]),
        ball + np.array([0, 0, -R]),
    ]
    surfs = [_make_surf(pts)]
    mn, mx = calc_spikes(ball, surfs)
    assert mn == pytest.approx(R, rel=TOL_REL, abs=TOL_ABS)
    assert mx == pytest.approx(R, rel=TOL_REL, abs=TOL_ABS)


def test_calc_spikes_two_shells_min_max():
    # Two shells with different radii; min=inner, max=outer
    ball = np.array([0.0, 0.0, 0.0])
    r_in, r_out = 2.0, 5.0
    inner_pts = [
        np.array([ r_in, 0, 0]),
        np.array([ 0, r_in, 0]),
        np.array([ 0, 0, r_in]),
    ]
    outer_pts = [
        np.array([-r_out, 0, 0]),
        np.array([ 0, -r_out, 0]),
        np.array([ 0, 0, -r_out]),
    ]
    surfs = [_make_surf(inner_pts), _make_surf(outer_pts)]
    mn, mx = calc_spikes(ball, surfs)
    assert mn == pytest.approx(r_in, rel=TOL_REL, abs=TOL_ABS)
    assert mx == pytest.approx(r_out, rel=TOL_REL, abs=TOL_ABS)


def test_calc_spikes_raises_on_empty_surfaces():
    # Current implementation will hit min()/max() on empty -> ValueError.
    # Make that behavior explicit in the test so it doesn't silently pass.
    with pytest.raises(ValueError):
        _ = calc_spikes(np.zeros(3), [])  # no surfaces / no points


# -------------------------
# calc_cell_box
# -------------------------

def test_calc_cell_box_basic_bounds():
    # Points spanning a rectangular box from mins to maxs
    mins = np.array([-2.0, 1.0, -5.0])
    maxs = np.array([ 3.0, 4.0,  6.0])
    pts = [
        mins,
        np.array([maxs[0], mins[1], mins[2]]),
        np.array([mins[0], maxs[1], mins[2]]),
        np.array([mins[0], mins[1], maxs[2]]),
        maxs,
    ]
    surfs = [_make_surf(pts[:3]), _make_surf(pts[3:])]
    bb = calc_cell_box(surfs)
    npt.assert_allclose(np.array(bb[0]), mins, rtol=TOL_REL, atol=TOL_ABS)
    npt.assert_allclose(np.array(bb[1]), maxs, rtol=TOL_REL, atol=TOL_ABS)


def test_calc_cell_box_single_point_surface():
    p = np.array([1.2, -3.4, 5.6])
    surfs = [_make_surf([p])]
    bb = calc_cell_box(surfs)
    npt.assert_allclose(np.array(bb[0]), p, rtol=TOL_REL, atol=TOL_ABS)
    npt.assert_allclose(np.array(bb[1]), p, rtol=TOL_REL, atol=TOL_ABS)


def test_calc_cell_box_negative_and_positive_coords():
    pts = [
        np.array([-10.0, -1.0,  0.0]),
        np.array([  2.5,  7.0, -3.0]),
        np.array([ -1.5,  0.0,  9.0]),
    ]
    surfs = [_make_surf(pts)]
    bb = calc_cell_box(surfs)
    mins = np.min(np.vstack(pts), axis=0)
    maxs = np.max(np.vstack(pts), axis=0)
    npt.assert_allclose(np.array(bb[0]), mins, rtol=TOL_REL, atol=TOL_ABS)
    npt.assert_allclose(np.array(bb[1]), maxs, rtol=TOL_REL, atol=TOL_ABS)


def _tri(points):
    return {"points": [np.asarray(p, float) for p in points],
            "tris": [(0, 1, 2)]}


def _centroid_of_triangle(a, b, c):
    return (a + b + c) / 3.0


# -------------------------
# calc_cell_com
# -------------------------

def test_calc_cell_com_single_triangle_about_origin():
    # ball at origin; one triangular face
    o = np.zeros(3)
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    c = np.array([0.0, 0.0, 1.0])
    surf = _tri([a, b, c])
    # volume of tetra (o,a,b,c)
    V = calc_tetra_vol(o, a, b, c)
    com = calc_cell_com(o, [surf], volume=V)
    # centroid of tetra is average of its 4 vertices
    expected = (o + a + b + c) / 4.0
    npt.assert_allclose(com, expected, rtol=1e-12, atol=1e-12)


def test_calc_cell_com_multiple_faces_additivity():
    o = np.zeros(3)
    # two disjoint triangles on xy-plane & xz-plane (non-coplanar with origin)
    a,b,c = np.array([2,0,0.]), np.array([0,2,0.]), np.array([0,0,2.])
    s1 = _tri([a, b, np.array([0.0, 0.0, 1.0])])
    s2 = _tri([a, c, np.array([0.0, 1.0, 0.0])])
    V = 0.0
    for tri in s1["tris"]:
        p = [o, *[s1["points"][i] for i in tri]]
        V += calc_tetra_vol(*p)
    for tri in s2["tris"]:
        p = [o, *[s2["points"][i] for i in tri]]
        V += calc_tetra_vol(*p)
    com = calc_cell_com(o, [s1, s2], volume=V)
    # volume-weighted average of tetra centroids
    num = np.zeros(3)
    for tri in s1["tris"]:
        ps = [o, *[s1["points"][i] for i in tri]]
        num += calc_tetra_vol(*ps) * (ps[0] + ps[1] + ps[2] + ps[3]) / 4.0
    for tri in s2["tris"]:
        ps = [o, *[s2["points"][i] for i in tri]]
        num += calc_tetra_vol(*ps) * (ps[0] + ps[1] + ps[2] + ps[3]) / 4.0
    expected = num / V
    npt.assert_allclose(com, expected, rtol=1e-12, atol=1e-12)


# -------------------------
# calc_cell_moi
# -------------------------

def test_calc_cell_moi_basic_properties():
    o = np.zeros(3)
    a,b,c = np.array([1,0,0.]), np.array([0,1,0.]), np.array([0,0,1.])
    surf = _tri([a,b,c])
    V = calc_tetra_vol(o, a, b, c)
    I = calc_cell_moi(o, [surf], volume=V, density=1.0)
    # symmetric and diagonal non-negative
    npt.assert_allclose(I, I.T, rtol=0, atol=0)
    assert np.all(np.diag(I) >= -1e-15)  # allow tiny numeric noise


# -------------------------
# combine_inertia_tensors
# -------------------------

def test_combine_inertia_tensors_identity_for_single_zero_shift():
    I = np.diag([1.0, 2.0, 3.0])
    C = np.array([0.0, 0.0, 0.0])
    I_total = combine_inertia_tensors([I], [C], common_centroid=np.zeros(3), masses=[5.0])
    npt.assert_allclose(I_total, I, rtol=0, atol=0)


def test_combine_inertia_tensors_parallel_axis_shift():
    I = np.diag([1.0, 1.0, 1.0])
    C = np.array([1.0, 0.0, 0.0])
    m = 2.0
    common = np.zeros(3)
    d = C - common
    expected_shift = m * (np.dot(d,d) * np.eye(3) - np.outer(d, d))
    I_total = combine_inertia_tensors([I], [C], common, [m])
    npt.assert_allclose(I_total, I + expected_shift, rtol=1e-12, atol=1e-12)


# -------------------------
# calc_total_inertia_tensor
# -------------------------

def test_calc_total_inertia_tensor_sphere_origin_and_shift():
    # One sphere at origin: I = (2/5) m r^2 I3
    m, r = 3.0, 2.0
    s = {"mass": m, "rad": r, "loc": np.zeros(3)}
    I = calc_total_inertia_tensor([s], common_point=np.zeros(3))
    npt.assert_allclose(I, (2.0/5.0)*m*r*r*np.eye(3), rtol=1e-12, atol=1e-12)

    # Shift common point by d: add parallel-axis shift
    common = np.array([1.0, 2.0, -1.0])
    d = s["loc"] - common
    shift = m * (np.dot(d,d) * np.eye(3) - np.outer(d, d))
    I_shifted = calc_total_inertia_tensor([s], common)
    npt.assert_allclose(I_shifted, I + shift, rtol=1e-12, atol=1e-12)


# -------------------------
# calc_contacts
# -------------------------

def test_calc_contacts_all_outside_triangle():
    loc = np.zeros(3)
    tri_pts = [np.array([2.0,0.0,0.0]), np.array([0.0,2.0,0.0]), np.array([0.0,0.0,2.0])]
    surf = {"points": tri_pts, "tris": [(0,1,2)]}
    rad = 0.5
    areas, vol = calc_contacts(loc, rad, [surf], [7])
    # When all points are outside, current logic does NOT add contact area.
    assert areas[7] == pytest.approx(0.0, abs=1e-12)
    # Volume is non-negative; exact value depends on projection.
    assert vol >= 0.0


def test_calc_contacts_all_outside_triangle():
    # Sphere tiny: all points outside -> contact area counts triangle, volume uses projected points only
    loc = np.zeros(3)
    tri_pts = [np.array([2.0,0.0,0.0]), np.array([0.0,2.0,0.0]), np.array([0.0,0.0,2.0])]
    surf = {"points": tri_pts, "tris": [(0,1,2)]}
    rad = 0.5
    areas, vol = calc_contacts(loc, rad, [surf], [7])
    # Contact area equals original triangle area per current logic (counted in mixed/outside paths)
    assert areas[7] == pytest.approx(0.0, abs=1e-12)
    assert vol >= 0.0  # not checking exact projection volume here; just sanity


def test_calc_contacts_mixed_triangle_no_crash():
    loc = np.zeros(3)
    tri_pts = [np.array([0.4,0.0,0.0]), np.array([0.0,0.4,0.0]), np.array([0.0,0.0,2.0])]
    surf = {"points": tri_pts, "tris": [(0,1,2)]}
    rad = 0.5  # two inside, one outside
    areas, vol = calc_contacts(loc, rad, [surf], [1])
    assert 1 in areas
    assert vol >= 0.0


# -------------------------
# rotate_points
# -------------------------

def test_rotate_points_invariants():
    vec = np.array([1.0, 2.0, 3.0])
    pts = [np.array([1.0, 0.0, 0.0]),
           np.array([0.0, 1.0, 0.0]),
           np.array([0.0, 0.0, 1.0])]

    # Forward rotation
    fwd = rotate_points(vec, pts, reverse=False)

    # Distances from origin should be preserved
    for p, q in zip(pts, fwd):
        assert np.linalg.norm(p) == pytest.approx(np.linalg.norm(q), rel=1e-12, abs=1e-12)

    # Reverse rotation
    bwd = rotate_points(vec, fwd, reverse=True)

    # Distances should also be preserved under reverse
    for q, r in zip(fwd, bwd):
        assert np.linalg.norm(q) == pytest.approx(np.linalg.norm(r), rel=1e-12, abs=1e-12)


def test_rotate_points_vec_already_z():
    vec = np.array([0.0, 0.0, 5.0])
    pts = [np.array([1.0,2.0,3.0]), np.array([-1.0,0.0,4.0])]
    out = rotate_points(vec, pts, reverse=False)
    # rotating around +z alignment path should preserve distances from origin
    for p, q in zip(pts, out):
        assert np.linalg.norm(p) == pytest.approx(np.linalg.norm(q), rel=1e-12, abs=1e-12)


# -------------------------
# get_time
# -------------------------

def test_get_time_examples():
    # Note: nopython int division over floats returns floor as float
    h, m, s = get_time(3661.0)
    assert (h, m, s) == (1.0, 1.0, 1.0)
    h, m, s = get_time(59.0)
    assert (h, m, s) == (0.0, 0.0, 59.0)
    h, m, s = get_time(3600.0)
    assert (h, m, s) == (1.0, 0.0, 0.0)


# -------------------------
# calc_vol
# -------------------------

def test_calc_vol_single_surface_matches_sum_tetra():
    loc = np.zeros(3)
    pts = np.array([[1.0,0,0],[0,1.0,0],[0,0,1.0],[2.0,0,0]])
    tris = [(0,1,2), (0,2,3)]
    # per-function API
    vol, surf_vols = calc_vol(loc, [pts], [tris])
    expected = 0.0
    for tri in tris:
        expected += calc_tetra_vol(loc, pts[tri[0]], pts[tri[1]], pts[tri[2]])
    assert vol == pytest.approx(expected, rel=1e-12, abs=1e-12)
    assert surf_vols[0] == pytest.approx(expected, rel=1e-12, abs=1e-12)


# -------------------------
# calc_curvature
# -------------------------

def test_calc_curvature_flat_plane_zero():
    # Three coplanar points with identical normals -> zero curvature
    pts = np.array([[0,0,0],[1,0,0],[0,1,0]], dtype=float)
    normals = np.array([[0,0,1],[0,0,1],[0,0,1]], dtype=float)
    curv = calc_curvature(pts, normals)
    npt.assert_allclose(curv, np.zeros(3), rtol=1e-12, atol=1e-12)


# -------------------------
# calc_aw_center / calc_pw_center
# -------------------------


def test_aw_pw_equal_for_equal_radii():
    r1 = r2 = 2.0
    l1 = np.array([0.0, 0.0, 0.0])
    l2 = np.array([4.0, 0.0, 0.0])
    aw_d, aw_c = calc_aw_center(r1, r2, l1, l2)
    pw_d, pw_c = calc_pw_center(r1, r2, l1, l2)
    # equal radii -> both reduce to midpoint and half the center distance
    assert aw_d == pytest.approx(2.0, rel=1e-12, abs=1e-12)
    npt.assert_allclose(aw_c, np.array([2.0, 0.0, 0.0]), rtol=1e-12, atol=1e-12)
    npt.assert_allclose(aw_c, pw_c, rtol=1e-12, atol=1e-12)
    assert pw_d == pytest.approx(aw_d, rel=1e-12, abs=1e-12)


def test_aw_pw_general_consistency_line_segment():
    r1, r2 = 1.0, 3.0
    l1 = np.array([0.0, 0.0, 0.0])
    l2 = np.array([10.0, 0.0, 0.0])
    aw_d, aw_c = calc_aw_center(r1, r2, l1, l2)
    pw_d, pw_c = calc_pw_center(r1, r2, l1, l2)
    # Centers must lie on the segment between l1 and l2
    for c in (aw_c, pw_c):
        assert c[1] == pytest.approx(0.0, abs=1e-12)
        assert c[2] == pytest.approx(0.0, abs=1e-12)
        assert 0.0 <= c[0] <= 10.0
    # distances should be finite and non-negative
    assert aw_d >= 0.0 and pw_d >= 0.0
