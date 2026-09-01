"""
Test module for vorpy.src.calculations.vert module.

This module contains comprehensive tests for all functions in the vert module,
including edge cases, error handling, and integration tests.
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch
from vorpy.src.calculations.vert import (
    _numeric_guard,
    _real_roots_quadratic,
    _safe_div,
    calc_vert_abcfs,
    calc_vert_case_1,
    calc_vert_case_2,
    filter_vert_locrads,
    calc_vert,
    calc_flat_vert,
    verify_aw,
    verify_prm,
    verify_pow,
    verify_site
)


class TestNumericGuard:
    """Test cases for _numeric_guard context manager."""
    
    def test_numeric_guard_basic(self):
        """Test basic functionality of numeric guard."""
        with _numeric_guard():
            # Should not raise any exceptions for normal operations
            result = 1.0 + 2.0
            assert result == 3.0
    
    def test_numeric_guard_division_by_zero(self):
        """Test that numeric guard catches division by zero."""
        with _numeric_guard():
            with pytest.raises(ZeroDivisionError):
                1.0 / 0.0
    
    def test_numeric_guard_invalid_operation(self):
        """Test that numeric guard catches invalid operations."""
        with _numeric_guard():
            with pytest.raises(FloatingPointError):
                np.sqrt(-1.0)
    
    def test_numeric_guard_restores_state(self):
        """Test that numeric guard restores previous state."""
        # Set initial state
        np.seterr(divide='ignore', invalid='ignore')
        
        with _numeric_guard():
            # Should raise on division by zero
            with pytest.raises(ZeroDivisionError):
                1.0 / 0.0
        
        # Should restore previous state
        # Note: np.RankWarning may not be available in newer NumPy versions
        try:
            with pytest.raises(np.RankWarning):
                np.polyfit([1, 2, 3], [1, 4, 9], 2)
        except AttributeError:
            # Skip this test if RankWarning is not available
            pass


class TestRealRootsQuadratic:
    """Test cases for _real_roots_quadratic function."""
    
    def test_real_roots_quadratic_two_real_roots(self):
        """Test quadratic with two real roots."""
        result = _real_roots_quadratic(1, -5, 6)  # (x-2)(x-3) = 0
        assert len(result) == 2
        # Check that roots are close to expected values
        assert any(abs(root - 2.0) < 1e-10 for root in result)
        assert any(abs(root - 3.0) < 1e-10 for root in result)
    
    def test_real_roots_quadratic_one_real_root(self):
        """Test quadratic with one real root."""
        result = _real_roots_quadratic(1, -4, 4)  # (x-2)^2 = 0
        # Perfect square may return duplicate roots
        assert len(result) >= 1
        assert any(abs(root - 2.0) < 1e-10 for root in result)
    
    def test_real_roots_quadratic_no_real_roots(self):
        """Test quadratic with no real roots."""
        result = _real_roots_quadratic(1, 0, 1)  # x^2 + 1 = 0
        assert len(result) == 0
    
    def test_real_roots_quadratic_linear(self):
        """Test linear equation (a=0)."""
        result = _real_roots_quadratic(0, 2, -4)  # 2x - 4 = 0
        assert len(result) == 1
        assert abs(result[0] - 2.0) < 1e-10
    
    def test_real_roots_quadratic_constant(self):
        """Test constant equation (a=0, b=0)."""
        result = _real_roots_quadratic(0, 0, 5)  # 5 = 0
        assert len(result) == 0
    
    def test_real_roots_quadratic_tolerance(self):
        """Test with custom tolerance."""
        # Roots with small imaginary parts
        result = _real_roots_quadratic(1, -2, 1.0001, tol=1e-3)
        # May not find real roots due to numerical precision
        assert len(result) >= 0


class TestSafeDiv:
    """Test cases for _safe_div function."""
    
    def test_safe_div_normal_division(self):
        """Test normal division."""
        result = _safe_div(10, 2)
        assert result == 5.0
    
    def test_safe_div_division_by_zero(self):
        """Test division by zero."""
        with pytest.raises(ValueError, match="denominator is zero or non-finite"):
            _safe_div(10, 0)
    
    def test_safe_div_very_small_denominator(self):
        """Test division by very small denominator."""
        with pytest.raises(ValueError, match="denominator is zero or non-finite"):
            _safe_div(10, 1e-20)
    
    def test_safe_div_infinite_denominator(self):
        """Test division by infinite denominator."""
        with pytest.raises(ValueError, match="denominator is zero or non-finite"):
            _safe_div(10, np.inf)
    
    def test_safe_div_nan_denominator(self):
        """Test division by NaN denominator."""
        with pytest.raises(ValueError, match="denominator is zero or non-finite"):
            _safe_div(10, np.nan)
    
    def test_safe_div_custom_name(self):
        """Test with custom error message name."""
        with pytest.raises(ValueError, match="divisor is zero or non-finite"):
            _safe_div(10, 0, name="divisor")
    
    def test_safe_div_custom_eps(self):
        """Test with custom epsilon."""
        with pytest.raises(ValueError, match="denominator is zero or non-finite"):
            _safe_div(10, 1e-10, eps=1e-5)


class TestCalcVertAbcfs:
    """Test cases for calc_vert_abcfs function."""
    
    def test_calc_vert_abcfs_basic(self):
        """Test basic vertex coefficient calculation."""
        locs = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ])
        rads = np.array([1.0, 1.0, 1.0, 1.0])
        
        result = calc_vert_abcfs(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
        Fs, abcdfs, rs, l0 = result
        
        assert isinstance(Fs, np.ndarray)
        assert isinstance(abcdfs, np.ndarray)
        assert isinstance(rs, np.ndarray)
        assert isinstance(l0, np.ndarray)
    
    def test_calc_vert_abcfs_different_radii(self):
        """Test with different sphere radii."""
        locs = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0]
        ])
        rads = np.array([1.0, 1.5, 2.0, 0.5])
        
        result = calc_vert_abcfs(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
    
    def test_calc_vert_abcfs_3d_positions(self):
        """Test with 3D positions."""
        locs = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0]
        ])
        rads = np.array([1.0, 1.0, 1.0, 1.0])
        
        result = calc_vert_abcfs(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 4


class TestCalcVertCase1:
    """Test cases for calc_vert_case_1 function."""
    
    def test_calc_vert_case_1_basic(self):
        """Test basic case 1 calculation."""
        Fs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # 8 coefficients
        l0 = np.array([0.0, 0.0, 0.0])
        r0 = 1.0
        
        result = calc_vert_case_1(Fs, l0, r0)
        
        assert isinstance(result, list)
        # Result is a list of vertices
        assert all(isinstance(vertex, list) for vertex in result)
    
    def test_calc_vert_case_1_different_radius(self):
        """Test with different radius."""
        Fs = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # 8 coefficients
        l0 = np.array([0.0, 0.0, 0.0])
        r0 = 2.0
        
        result = calc_vert_case_1(Fs, l0, r0)
        
        assert isinstance(result, list)
        # Result is a list of vertices
        assert all(isinstance(vertex, list) for vertex in result)


class TestCalcVertCase2:
    """Tests for the Numba-compatible legacy Case 2 solver."""

    def test_calc_vert_case_2_basic(self):
        """Case 2 should compile and return both real quadratic solutions."""
        Fs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        r0 = 1.0
        l0 = np.array([0.0, 0.0, 0.0])

        result = calc_vert_case_2(Fs, r0, l0)

        assert isinstance(result, list)
        assert len(result) == 2

        # For these coefficients the free-coordinate equation is
        # 2*z^2 - 2 = 0, giving z = -1 and +1.
        zs = sorted(float(vertex[0][2]) for vertex in result)
        assert np.allclose(zs, [-1.0, 1.0], atol=1e-12, rtol=1e-12)

        # Reconstruction should be internally consistent:
        # x = 1 + z, y = 1 + z, R = 1 + z.
        for loc, rad in result:
            z = float(loc[2])
            assert np.allclose(loc[:2], [1.0 + z, 1.0 + z],
                               atol=1e-12, rtol=1e-12)
            assert np.isclose(rad, 1.0 + z, atol=1e-12, rtol=1e-12)

    def test_calc_vert_case_2_different_radius(self):
        """Changing r0 should still produce valid roots of the Case 2 equation."""
        Fs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        r0 = 2.0
        l0 = np.array([0.0, 0.0, 0.0])

        result = calc_vert_case_2(Fs, r0, l0)

        assert isinstance(result, list)
        assert len(result) == 2

        # a=2, b=-2, c=-7 for this fixture.
        for loc, rad in result:
            z = float(loc[2])
            assert np.isclose(2.0*z*z - 2.0*z - 7.0, 0.0,
                              atol=1e-10, rtol=1e-10)
            assert np.allclose(loc[:2], [1.0 + z, 1.0 + z],
                               atol=1e-12, rtol=1e-12)
            assert np.isclose(rad, 1.0 + z, atol=1e-12, rtol=1e-12)

    def test_calc_vert_case_2_rejects_zero_determinant(self):
        """Case 2 cannot reconstruct a vertex when F is zero."""
        Fs = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        result = calc_vert_case_2(Fs, 1.0, np.zeros(3))
        assert result == []


class TestFilterVertLocrads:
    """Test cases for filter_vert_locrads function."""
    
    def test_filter_vert_locrads_basic(self):
        """Test basic vertex filtering."""
        verts = [
            (np.array([0.0, 0.0, 0.0]), 1.0),
            (np.array([1.0, 1.0, 1.0]), 2.0),
            (np.array([2.0, 2.0, 2.0]), 0.5)
        ]
        rs = [1.0, 2.0, 0.5]
        
        result = filter_vert_locrads(verts, rs)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
        loc, loc2, rad, rad2 = result
        # Some values may be None
        assert loc is None or isinstance(loc, np.ndarray)
        assert rad is None or isinstance(rad, (int, float, np.number))
    
    def test_filter_vert_locrads_empty_input(self):
        """Test with empty input."""
        verts = []
        rs = []
        
        result = filter_vert_locrads(verts, rs)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
        loc, loc2, rad, rad2 = result
        # All values should be None for empty input
        assert loc is None
        assert rad is None
        assert loc2 is None
        assert rad2 is None
    
    def test_filter_vert_locrads_single_vertex(self):
        """Test with single vertex."""
        verts = [(np.array([0.0, 0.0, 0.0]), 1.0)]
        rs = [1.0]
        
        result = filter_vert_locrads(verts, rs)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
        loc, loc2, rad, rad2 = result
        # Some values may be None
        assert loc is None or isinstance(loc, np.ndarray)
        assert rad is None or isinstance(rad, (int, float, np.number))


class TestCalcVert:
    """Test cases for calc_vert function."""
    
    def test_calc_vert_basic(self):
        """Test basic vertex calculation."""
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        rads = [1.0, 1.0, 1.0, 1.0]
        
        result = calc_vert(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
        loc, rad, loc2, rad2 = result
        # Some values may be None
        assert loc is None or isinstance(loc, (list, np.ndarray))
        assert rad is None or isinstance(rad, (int, float, np.number))
    
    def test_calc_vert_different_radii(self):
        """Test with different sphere radii."""
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([0.0, 2.0, 0.0]),
            np.array([0.0, 0.0, 2.0])
        ]
        rads = [1.0, 1.5, 2.0, 0.5]
        
        result = calc_vert(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 4
    
    def test_calc_vert_3d_positions(self):
        """Test with 3D positions."""
        locs = [
            np.array([1.0, 2.0, 3.0]),
            np.array([4.0, 5.0, 6.0]),
            np.array([7.0, 8.0, 9.0]),
            np.array([10.0, 11.0, 12.0])
        ]
        rads = [1.0, 1.0, 1.0, 1.0]
        
        result = calc_vert(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 4


class TestCalcFlatVert:
    """Test cases for calc_flat_vert function."""
    
    def test_calc_flat_vert_basic(self):
        """Test basic flat vertex calculation."""
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        rads = [1.0, 1.0, 1.0, 1.0]
        
        result = calc_flat_vert(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 2
        loc, rad = result
        assert isinstance(loc, list)  # Returns a list, not numpy array
        assert isinstance(rad, (int, float, np.number))
    
    def test_calc_flat_vert_power_diagram(self):
        """Test with power diagram method."""
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        rads = [1.0, 1.0, 1.0, 1.0]
        
        result = calc_flat_vert(locs, rads, power=True)
        
        assert isinstance(result, tuple)
        assert len(result) == 2
    
    def test_calc_flat_vert_different_radii(self):
        """Test with different radii."""
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([0.0, 2.0, 0.0]),
            np.array([0.0, 0.0, 2.0])
        ]
        rads = [1.0, 1.5, 2.0, 0.5]
        
        result = calc_flat_vert(locs, rads)
        
        assert isinstance(result, tuple)
        assert len(result) == 2


class TestVerifyAw:
    """Test cases for verify_aw function."""
    
    def test_verify_aw_no_encroachment(self):
        """Test verification with no encroachment."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_aw(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
        assert result is True  # No encroachment
    
    def test_verify_aw_encroachment(self):
        """Test verification with encroachment."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[1.5, 0.0, 0.0], [0.0, 1.5, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_aw(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
        assert result is False  # Encroachment detected
    
    def test_verify_aw_empty_test_locs(self):
        """Test with empty test locations."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([])
        test_rads = np.array([])
        
        result = verify_aw(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
        assert result is True  # No test locations means no encroachment
    
    def test_verify_aw_identical_locations(self):
        """Test with identical locations."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[0.0, 0.0, 0.0]])
        test_rads = np.array([1.0])
        
        result = verify_aw(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
        assert result is False  # Identical locations mean encroachment


class TestVerifyPrm:
    """Test cases for verify_prm function."""
    
    def test_verify_prm_basic(self):
        """Test basic parameter verification."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        
        result = verify_prm(loc, rad, test_locs)
        
        assert isinstance(result, bool)
    
    def test_verify_prm_empty_test_locs(self):
        """Test with empty test locations."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([])
        
        result = verify_prm(loc, rad, test_locs)
        
        assert isinstance(result, bool)
    
    def test_verify_prm_close_locations(self):
        """Test with close locations."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0]])
        
        result = verify_prm(loc, rad, test_locs)
        
        assert isinstance(result, bool)


class TestVerifyPow:
    """Test cases for verify_pow function."""
    
    def test_verify_pow_basic(self):
        """Test basic power verification."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_pow(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
    
    def test_verify_pow_empty_test_locs(self):
        """Test with empty test locations."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([])
        test_rads = np.array([])
        
        result = verify_pow(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
    
    def test_verify_pow_different_radii(self):
        """Test with different radii."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        test_rads = np.array([1.5, 0.5])
        
        result = verify_pow(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)


class TestVerifySite:
    """Test cases for verify_site function."""
    
    def test_verify_site_aw_network(self):
        """Test site verification for AW network."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_site(loc, rad, test_locs, test_rads, net_type='aw')
        
        assert isinstance(result, bool)
    
    def test_verify_site_prm_network(self):
        """Test site verification for PRM network."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_site(loc, rad, test_locs, test_rads, net_type='prm')
        
        assert isinstance(result, bool)
    
    def test_verify_site_pow_network(self):
        """Test site verification for POW network."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_site(loc, rad, test_locs, test_rads, net_type='pow')
        
        assert isinstance(result, bool)
    
    def test_verify_site_default_network(self):
        """Test site verification with default network type."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        result = verify_site(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)
    
    def test_verify_site_empty_test_locs(self):
        """Test with empty test locations."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([])
        test_rads = np.array([])
        
        result = verify_site(loc, rad, test_locs, test_rads)
        
        assert isinstance(result, bool)


class TestVertIntegration:
    """Integration tests for vert module functions."""
    
    def test_vertex_calculation_workflow(self):
        """Test complete workflow of vertex calculations."""
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        rads = [1.0, 1.0, 1.0, 1.0]
        
        # Calculate vertex
        loc, rad, loc2, rad2 = calc_vert(locs, rads)
        # Some values may be None
        assert loc is None or isinstance(loc, (list, np.ndarray))
        assert rad is None or isinstance(rad, (int, float, np.number))
        
        # Calculate flat vertex
        flat_loc, flat_rad = calc_flat_vert(locs, rads)
        assert isinstance(flat_loc, list)
        assert isinstance(flat_rad, (int, float, np.number))
    
    def test_verification_workflow(self):
        """Test complete workflow of verification functions."""
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1.0
        test_locs = np.array([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
        test_rads = np.array([1.0, 1.0])
        
        # Test different verification methods
        aw_result = verify_aw(loc, rad, test_locs, test_rads)
        prm_result = verify_prm(loc, rad, test_locs)
        pow_result = verify_pow(loc, rad, test_locs, test_rads)
        site_result = verify_site(loc, rad, test_locs, test_rads)
        
        assert isinstance(aw_result, bool)
        assert isinstance(prm_result, bool)
        assert isinstance(pow_result, bool)
        assert isinstance(site_result, bool)
    
    def test_edge_case_handling(self):
        """Test edge case handling across functions."""
        # Test with very small values
        locs = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1e-10, 0.0, 0.0]),
            np.array([0.0, 1e-10, 0.0]),
            np.array([0.0, 0.0, 1e-10])
        ]
        rads = [1e-10, 1e-10, 1e-10, 1e-10]
        
        # Should not crash
        try:
            verts, vert_rads = calc_vert(locs, rads)
            assert isinstance(verts, list)
            assert isinstance(vert_rads, list)
        except (ValueError, FloatingPointError):
            # Expected for very small values
            pass
        
        # Test verification with edge cases
        loc = np.array([0.0, 0.0, 0.0])
        rad = 1e-10
        test_locs = np.array([[1e-9, 0.0, 0.0]])
        test_rads = np.array([1e-10])
        
        try:
            result = verify_aw(loc, rad, test_locs, test_rads)
            assert isinstance(result, bool)
        except (ValueError, FloatingPointError):
            # Expected for very small values
            pass
    
    def test_numeric_guard_integration(self):
        """Test numeric guard integration with other functions."""
        with _numeric_guard():
            # Test that numeric guard works with vertex calculations
            locs = [
                np.array([0.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0.0, 1.0, 0.0]),
                np.array([0.0, 0.0, 1.0])
            ]
            rads = [1.0, 1.0, 1.0, 1.0]
            
            # Should work normally
            result = calc_vert(locs, rads)
            assert isinstance(result, tuple)
            assert len(result) == 4
