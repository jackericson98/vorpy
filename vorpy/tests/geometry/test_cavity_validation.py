"""Regression checks for enclosed components, independently of handles."""
import numpy as np
import pytest

from vorpy.src.geometry.validation_shapes import sphere_with_cavities
from vorpy.src.geometry.topology_diagnostic import (
    inspect_group, boundary_component_reports, cavity_acceptance_errors,
    export_planar_boundary,
)


def test_cavity_generator_selects_surrounding_cells_only(tmp_path):
    from vorpy.src.geometry.validation_shapes import _build_parser, _generate_from_args

    shape = _generate_from_args(_build_parser().parse_args([
        '--shape', 'sphere_with_cavities', '--n-cavities', '1', '-r', '32',
    ]))
    assert shape.parameters['interior_indices'] == list(range(32))
    assert shape.parameters['cavity_indices'] == [32]
    assert shape.xyzr[32, :3] == pytest.approx([0, 0, 0])
    assert np.linalg.norm(shape.xyzr[:32, :3], axis=1) == pytest.approx(np.full(32, 6.))
    assert np.linalg.norm(shape.xyzr[33:, :3], axis=1) == pytest.approx(np.full(65, 14.))
    lines = shape.save_pdb(tmp_path / 'cavity.pdb').read_text().splitlines()
    atoms = [line for line in lines if line.startswith('HETATM')]
    assert all(line[12:16].strip() == 'INT' and line[21] == 'I' for line in atoms[:32])
    assert all(line[12:16].strip() != 'INT' and line[21] != 'I' for line in atoms[32:])
    with pytest.raises(ValueError, match='at least'):
        sphere_with_cavities(outer_radius=8)


@pytest.mark.parametrize("n_cavities", [0, 1, 2, 3])
@pytest.mark.timeout(180)
def test_cavity_measured_boundary_and_orientation(tmp_path, n_cavities):
    from vorpy.src.command.set import sett
    from vorpy.src.group.group import Group
    from vorpy.src.system.system import System

    shape = sphere_with_cavities(n_cavities=n_cavities)
    path = shape.save_pdb(tmp_path / 'cavity.pdb')
    system = System(file=str(path), make_dir=False, print_actions=False)
    selected = [int(row['num']) for _, row in system.balls.iterrows()
                if str(row['name']).strip() == 'INT']
    group = Group(system, name='cavity_test', atoms=selected,
                  settings=sett('mv', ['30']), make_net=True)
    group.build()
    group.get_info()
    report = inspect_group(group)
    report['component_reports'] = boundary_component_reports(group)
    assert cavity_acceptance_errors(report, n_cavities) == []
    assert report['complete_selected_cells'] == report['selected_cells'] == shape.parameters['interior_count']
    assert report['selected_cell_components'] == 1
    assert report['boundary_components'] == n_cavities + 1
    assert report['euler_characteristic'] == 2 * (n_cavities + 1)
    assert report['integrated_gaussian_curvature'] == pytest.approx(4 * np.pi * (n_cavities + 1), abs=1e-8, rel=0)
    for part in report['component_reports']:
        assert part['chi'] == 2
        assert part['integrated_gaussian_curvature'] == pytest.approx(4 * np.pi, abs=1e-8, rel=0)
    # Inspect exported assembled surfaces, not merely the generator spheres.
    export_planar_boundary(group, report['component_reports'], tmp_path)
    obj = (tmp_path / 'boundary.obj').read_text().splitlines()
    assert sum(line.startswith('o ') for line in obj) == n_cavities + 1
    assert sum(line.startswith('f ') for line in obj) == report['boundary_faces']
    compile((tmp_path / 'boundary.py').read_text(), 'boundary.py', 'exec')


@pytest.mark.parametrize("count", [0, 1, 2, 3, 4, 5, 12, 25])
def test_cavity_placement_clearance_and_reproducibility(count):
    shape = sphere_with_cavities(count)
    assert np.array_equal(shape.xyzr, sphere_with_cavities(count).xyzr)
    p = shape.parameters
    centers = p['cavity_centers']
    assert len(centers) == count
    assert set(p['cavity_indices']).isdisjoint(p['interior_indices'])
    assert np.all(np.linalg.norm(centers, axis=1) + 3 * p['cavity_radius'] <= p['outer_radius'] + 1e-10)
    if count > 1:
        distances = np.linalg.norm(centers[:, None] - centers[None, :], axis=2)
        np.fill_diagonal(distances, np.inf)
        assert distances.min() == pytest.approx(6 * p['cavity_radius'])
        with pytest.raises(ValueError, match='at least'):
            sphere_with_cavities(count, outer_radius=p['minimum_outer_radius'] - 0.1)


@pytest.mark.parametrize("count", [-1, 1.5, True])
def test_invalid_cavity_counts(count):
    with pytest.raises(ValueError, match='nonnegative integer'):
        sphere_with_cavities(count)


def test_cavity_cli_defaults_and_unique_output_names():
    from vorpy.src.geometry.validation_shapes import _build_parser, _generate_from_args, _default_basename
    shapes = [_generate_from_args(_build_parser().parse_args([
        '--shape', 'sphere_with_cavities', '--n-cavities', str(n),
    ])) for n in (1, 5, 10)]
    assert all(s.parameters['resolution'] == 32 for s in shapes)
    assert len({_default_basename(s) for s in shapes}) == 3
    assert shapes[-1].parameters['outer_radius'] > shapes[0].parameters['outer_radius']


def test_zero_cavities_accepts_small_body_and_radius_cli_alias():
    from vorpy.src.geometry.validation_shapes import _build_parser, _generate_from_args
    shape = _generate_from_args(_build_parser().parse_args([
        '--shape', 'sphere_with_cavities', '--n-cavities', '0', '--radius', '2',
    ]))
    assert shape.parameters['outer_radius'] == 2
    assert shape.parameters['cavity_indices'] == []
    assert np.linalg.norm(shape.xyzr[1:, :3], axis=1) == pytest.approx(np.full(65, 4.))
