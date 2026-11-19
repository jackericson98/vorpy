"""Nondestructive smoke tests for canonical F6 multiscale modules."""
import importlib

import pytest


BASE = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F6_Proteins_Interfaces_Multiscale'


@pytest.mark.parametrize('name', [
    'F3_common', 'F3A_multiscale_summary', 'F3B_atomic_groups',
    'F3C_residue_groups', 'F3D_geometry_examples', 'F3E_solvent_exposure',
    'F3_build_figure',
])
def test_canonical_f6_modules_import_without_execution(name):
    importlib.import_module(f'{BASE}.{name}')


def test_multiscale_difference_definitions_are_unchanged():
    common = importlib.import_module(f'{BASE}.F3_common')
    assert common.signed_percent_difference(12, 10) == pytest.approx(20)
    assert common.absolute_percent_difference(8, 10) == pytest.approx(20)
    assert common.mean_absolute_percent_difference([12, 8], [10, 10]) == pytest.approx(20)
