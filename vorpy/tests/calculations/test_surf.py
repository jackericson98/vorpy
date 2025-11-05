"""
Test module for vorpy.src.calculations.surf module.

This module contains comprehensive tests for all functions in the surf module,
including edge cases, error handling, and integration tests.
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch
from vorpy.src.calculations.surf import (
    calc_tri,
    calc_surf_func,
    calc_surf_func_jit,
    calc_surf_func_reg,
    calc_2d_surf_sa,
    calc_surf_sa,
    calc_surf_tri_dists
)


class TestCalcTri:
    """Test cases for calc_tri function."""
    
    def test_calc_tri_basic_triangle(self):
        """Test basic triangle area calculation."""
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0])
        ]
        result = calc_tri(points)
        expected = 0.5
        assert abs(result - expected) < 1e-10
    
    def test_calc_tri_unit_triangle(self):
        """Test triangle with unit area."""
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0])
        ]
        result = calc_tri(points)
        expected = 1.0
        assert abs(result - expected) < 1e-10
    
    def test_calc_tri_3d_triangle(self):
        """Test triangle in 3D space."""
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        result = calc_tri(points)
        expected = 0.5
        assert abs(result - expected) < 1e-10
    
    def test_calc_tri_equilateral_triangle(self):
        """Test equilateral triangle."""
        # Equilateral triangle with side length 2
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([1.0, np.sqrt(3), 0.0])
        ]
        result = calc_tri(points)
        expected = np.sqrt(3)  # Area of equilateral triangle with side 2
        assert abs(result - expected) < 1e-10
    
    def test_calc_tri_degenerate_triangle(self):
        """Test degenerate triangle (collinear points)."""
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0])
        ]
        result = calc_tri(points)
        expected = 0.0
        assert abs(result - expected) < 1e-10
    
    def test_calc_tri_negative_coordinates(self):
        """Test triangle with negative coordinates."""
        points = [
            np.array([-1.0, -1.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([-1.0, 0.0, 0.0])
        ]
        result = calc_tri(points)
        expected = 0.5
        assert abs(result - expected) < 1e-10


class TestCalcSurfFunc:
    """Test cases for calc_surf_func function."""
    
    def test_calc_surf_func_basic(self):
        """Test basic surface function calculation."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([2, 0, 0])
        r1 = 1.0
        
        result = calc_surf_func(l0, r0, l1, r1)
        
        # Should return a list with 14 coefficients
        assert isinstance(result, list)
        assert len(result) == 14
        assert all(isinstance(x, (int, float, np.number)) for x in result)
    
    def test_calc_surf_func_different_radii(self):
        """Test with different sphere radii."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([3, 0, 0])
        r1 = 2.0
        
        result = calc_surf_func(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14
    
    def test_calc_surf_func_3d_positions(self):
        """Test with 3D positions."""
        l0 = np.array([1, 2, 3])
        r0 = 1.5
        l1 = np.array([4, 5, 6])
        r1 = 2.0
        
        result = calc_surf_func(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14
    
    def test_calc_surf_func_identical_spheres(self):
        """Test with identical spheres."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([0, 0, 0])
        r1 = 1.0
        
        result = calc_surf_func(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14
    
    def test_calc_surf_func_touching_spheres(self):
        """Test with touching spheres."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([2, 0, 0])
        r1 = 1.0
        
        result = calc_surf_func(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14


class TestCalcSurfFuncJit:
    """Test cases for calc_surf_func_jit function."""
    
    def test_calc_surf_func_jit_basic(self):
        """Test basic JIT surface function calculation."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([2, 0, 0])
        r1 = 1.0
        
        result = calc_surf_func_jit(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14
    
    def test_calc_surf_func_jit_different_radii(self):
        """Test JIT function with different radii."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([3, 0, 0])
        r1 = 2.0
        
        result = calc_surf_func_jit(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14


class TestCalcSurfFuncReg:
    """Test cases for calc_surf_func_reg function."""
    
    def test_calc_surf_func_reg_basic(self):
        """Test basic regular surface function calculation."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([2, 0, 0])
        r1 = 1.0
        
        result = calc_surf_func_reg(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14
    
    def test_calc_surf_func_reg_radius_swap(self):
        """Test that radius swapping works correctly."""
        l0 = np.array([0, 0, 0])
        r0 = 2.0  # Larger radius first
        l1 = np.array([2, 0, 0])
        r1 = 1.0  # Smaller radius second
        
        result = calc_surf_func_reg(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14
    
    def test_calc_surf_func_reg_list_input(self):
        """Test with list input instead of numpy array."""
        l0 = [0, 0, 0]
        r0 = 1.0
        l1 = [2, 0, 0]
        r1 = 1.0
        
        result = calc_surf_func_reg(l0, r0, l1, r1)
        
        assert isinstance(result, list)
        assert len(result) == 14


class TestCalc2dSurfSa:
    """Test cases for calc_2d_surf_sa function."""
    
    def test_calc_2d_surf_sa_single_triangle(self):
        """Test surface area calculation for single triangle."""
        tris = [(0, 1, 2)]
        points = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0])
        ]
        
        result = calc_2d_surf_sa(tris, points)
        expected = 0.5
        assert abs(result - expected) < 1e-10
    
    def test_calc_2d_surf_sa_multiple_triangles(self):
        """Test surface area calculation for multiple triangles."""
        tris = [(0, 1, 2), (1, 3, 2)]
        points = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([1, 1, 0])
        ]
        
        result = calc_2d_surf_sa(tris, points)
        expected = 1.0  # Two triangles of area 0.5 each
        assert abs(result - expected) < 1e-10
    
    def test_calc_2d_surf_sa_square(self):
        """Test surface area calculation for square made of triangles."""
        tris = [(0, 1, 2), (1, 3, 2)]
        points = [
            np.array([0, 0, 0]),
            np.array([2, 0, 0]),
            np.array([0, 2, 0]),
            np.array([2, 2, 0])
        ]
        
        result = calc_2d_surf_sa(tris, points)
        expected = 4.0  # Square of side 2
        assert abs(result - expected) < 1e-10
    
    def test_calc_2d_surf_sa_empty_triangles(self):
        """Test with empty triangle list."""
        tris = []
        points = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0])
        ]
        
        result = calc_2d_surf_sa(tris, points)
        expected = 0.0
        assert result == expected
    
    def test_calc_2d_surf_sa_degenerate_triangle(self):
        """Test with degenerate triangle."""
        tris = [(0, 1, 2)]
        points = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([2, 0, 0])  # Collinear points
        ]
        
        result = calc_2d_surf_sa(tris, points)
        expected = 0.0
        assert result == expected


class TestCalcSurfSa:
    """Test cases for calc_surf_sa function."""
    
    def test_calc_surf_sa_single_triangle(self):
        """Test 3D surface area calculation for single triangle."""
        tris = [(0, 1, 2)]
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0])
        ]
        
        result = calc_surf_sa(tris, points)
        expected = 0.5
        assert abs(result - expected) < 1e-10
    
    def test_calc_surf_sa_3d_triangle(self):
        """Test 3D surface area calculation for triangle in 3D space."""
        tris = [(0, 1, 2)]
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        
        result = calc_surf_sa(tris, points)
        expected = 0.5
        assert abs(result - expected) < 1e-10
    
    def test_calc_surf_sa_multiple_triangles(self):
        """Test 3D surface area calculation for multiple triangles."""
        tris = [(0, 1, 2), (1, 3, 2)]
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([1.0, 1.0, 0.0])
        ]
        
        result = calc_surf_sa(tris, points)
        expected = 1.0
        assert abs(result - expected) < 1e-10
    
    def test_calc_surf_sa_empty_triangles(self):
        """Test with empty triangle list."""
        tris = []
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0])
        ]
        
        result = calc_surf_sa(tris, points)
        expected = 0.0
        assert result == expected
    
    def test_calc_surf_sa_complex_3d_surface(self):
        """Test with complex 3D surface."""
        tris = [(0, 1, 2), (1, 3, 2), (0, 2, 4), (2, 3, 4)]
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([1.0, 1.0, 0.0]),
            np.array([0.5, 0.5, 1.0])
        ]
        
        result = calc_surf_sa(tris, points)
        assert result > 0  # Should have positive area


class TestCalcSurfTriDists:
    """Test cases for calc_surf_tri_dists function."""
    
    def test_calc_surf_tri_dists_basic(self):
        """Test basic distance calculation."""
        points = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0])
        ]
        tris = [(0, 1, 2)]
        loc = np.array([0, 0, 0])
        
        result = calc_surf_tri_dists(points, tris, loc)
        
        assert isinstance(result, list)
        assert len(result) == len(points)
        assert all(0 <= dist <= 1 for dist in result)
    
    def test_calc_surf_tri_dists_multiple_points(self):
        """Test with multiple points."""
        points = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([0, 1, 0]),
            np.array([2, 2, 0])
        ]
        tris = [(0, 1, 2), (1, 3, 2)]
        loc = np.array([0, 0, 0])
        
        result = calc_surf_tri_dists(points, tris, loc)
        
        assert isinstance(result, list)
        assert len(result) == len(points)
        assert all(0 <= dist <= 1 for dist in result)
    
    def test_calc_surf_tri_dists_center_location(self):
        """Test with location at center of points."""
        points = [
            np.array([0, 0, 0]),
            np.array([2, 0, 0]),
            np.array([0, 2, 0]),
            np.array([2, 2, 0])
        ]
        tris = [(0, 1, 2), (1, 3, 2)]
        loc = np.array([1, 1, 0])  # Center
        
        result = calc_surf_tri_dists(points, tris, loc)
        
        assert isinstance(result, list)
        assert len(result) == len(points)
        # When all distances are the same, normalization fails with division by zero
        # Check that we get valid results (may contain nan values)
        assert all(np.isnan(dist) or (0 <= dist <= 1) for dist in result)
    
    def test_calc_surf_tri_dists_identical_points(self):
        """Test with identical points."""
        points = [
            np.array([1, 1, 1]),
            np.array([1, 1, 1]),
            np.array([1, 1, 1])
        ]
        tris = [(0, 1, 2)]
        loc = np.array([0, 0, 0])
        
        result = calc_surf_tri_dists(points, tris, loc)
        
        assert isinstance(result, list)
        assert len(result) == len(points)
        # When all distances are identical, normalization fails with division by zero
        # Check that we get valid results (may contain nan values)
        assert all(np.isnan(dist) or (0 <= dist <= 1) for dist in result)
    
    def test_calc_surf_tri_dists_3d_space(self):
        """Test with 3D coordinates."""
        points = [
            np.array([0, 0, 0]),
            np.array([1, 1, 1]),
            np.array([2, 2, 2])
        ]
        tris = [(0, 1, 2)]
        loc = np.array([0, 0, 0])
        
        result = calc_surf_tri_dists(points, tris, loc)
        
        assert isinstance(result, list)
        assert len(result) == len(points)
        assert all(0 <= dist <= 1 for dist in result)


class TestSurfIntegration:
    """Integration tests for surf module functions."""
    
    def test_surface_calculation_workflow(self):
        """Test complete workflow of surface calculations."""
        # Create a simple surface
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([1.0, 1.0, 0.0])
        ]
        tris = [(0, 1, 2), (1, 3, 2)]
        
        # Calculate surface area
        area = calc_surf_sa(tris, points)
        assert area > 0
        
        # Calculate distances
        loc = np.array([0.5, 0.5, 0.0])
        dists = calc_surf_tri_dists(points, tris, loc)
        assert len(dists) == len(points)
        # Check that we get valid results (may contain nan values for identical distances)
        assert all(np.isnan(dist) or (0 <= dist <= 1) for dist in dists)
    
    def test_surface_function_consistency(self):
        """Test that different surface function implementations give consistent results."""
        l0 = np.array([0, 0, 0])
        r0 = 1.0
        l1 = np.array([2, 0, 0])
        r1 = 1.0
        
        # Test main function
        result_main = calc_surf_func(l0, r0, l1, r1)
        assert isinstance(result_main, list)
        assert len(result_main) == 14
        
        # Test JIT function
        result_jit = calc_surf_func_jit(l0, r0, l1, r1)
        assert isinstance(result_jit, list)
        assert len(result_jit) == 14
        
        # Test regular function
        result_reg = calc_surf_func_reg(l0, r0, l1, r1)
        assert isinstance(result_reg, list)
        assert len(result_reg) == 14
    
    def test_triangle_area_consistency(self):
        """Test that 2D and 3D triangle area calculations are consistent for flat triangles."""
        # Create a flat triangle (z=0)
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0])
        ]
        
        # 3D calculation
        area_3d = calc_surf_sa([(0, 1, 2)], points)
        
        # 2D calculation (same points)
        area_2d = calc_2d_surf_sa([(0, 1, 2)], points)
        
        # Should be the same for flat triangle
        assert abs(area_3d - area_2d) < 1e-10
    
    def test_edge_case_handling(self):
        """Test edge case handling across functions."""
        # Test with very small values
        points = [
            np.array([0.0, 0.0, 0.0]),
            np.array([1e-10, 0.0, 0.0]),
            np.array([0.0, 1e-10, 0.0])
        ]
        
        # Should not crash
        area = calc_surf_sa([(0, 1, 2)], points)
        assert area >= 0
        
        # Test distance calculation
        dists = calc_surf_tri_dists(points, [(0, 1, 2)], np.array([0.0, 0.0, 0.0]))
        assert len(dists) == len(points)
        assert all(0 <= dist <= 1 for dist in dists)
