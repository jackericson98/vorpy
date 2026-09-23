import csv

import pandas as pd
import pytest

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.output.logs import EDGE_CURVATURE_HEADERS, edge_curvature_log_values


def write_section(writer, name, values):
    writer.writerow([name])
    # Column order must not affect interpretation.
    writer.writerow(list(values)[::-1])
    writer.writerow(list(values.values())[::-1])


@pytest.mark.parametrize('return_dict', [True, False])
def test_current_log_curvature_and_boundary_fields(tmp_path, return_dict):
    path = tmp_path / 'current.csv'
    breakdown = {
        'Integrated Mean Curvature (Face)': -1.5,
        'Integrated Mean Curvature (Edge)': 2.5,
        'Integrated Mean Curvature (Total)': 1.,
        'Integrated Gaussian Curvature (Face)': -3.,
        'Integrated Gaussian Curvature (Edge)': 4.,
        'Integrated Gaussian Curvature (Vertex)': 5.,
        'Integrated Gaussian Curvature (Total)': 6.,
    }
    group = {
        'Name': 'sample', 'Integrated Mean Curvature': 1.,
        'Integrated Mean Curvature Squared': 2., 'Integrated Gaussian Curvature': 6.,
        **breakdown,
        'Euler Characteristic': 2, 'Gauss-Bonnet Expected': 12.56,
        'Gauss-Bonnet Error': -0.01, 'Gauss-Bonnet Relative Error': -0.001,
        'Boundary Complete': True, 'Boundary Closed': False,
        'Boundary Manifold': True, 'Boundary Orientable': False,
        'Boundary Components': 1, 'Missing Faces': 2, 'Missing Edges': 3,
        'Missing Vertices': 4, 'Nonmanifold Edges': 5, 'Nonmanifold Vertices': 6,
    }
    edge_values = edge_curvature_log_values(
        [10, 20, 30], {10: -1.25, 20: 0., 30: 2.},
        {10: 3., 30: 4.}, {(10, 20): -0.5, (30, 10): 0.}, lambda value: value,
    )
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        write_section(writer, 'build informaiton', {'Name': 'sample', 'vorPy version': '3.5'})
        write_section(writer, 'group information', group)
        write_section(writer, 'Atoms', {'Index': 10, 'Name': 'CA', **breakdown})
        write_section(writer, 'Surfaces', {
            'Index': 1, 'Ball 1': 10, 'Ball 2': 20, 'Surface Area': 7.,
            'Integrated Gaussian Curvature (Face)': -3.,
        })
        write_section(writer, 'Edges', {
            'Index': 2, 'Ball 1': 10, 'Ball 2': 20, 'Ball 3': 30, 'Length': 8.,
            **dict(zip(EDGE_CURVATURE_HEADERS, edge_values)),
        })
        write_section(writer, 'Vertices', {
            'Index': 3, 'Ball 1': 10, 'Ball 2': 20, 'Ball 3': 30, 'Ball 4': 40,
            'x': 1., 'y': 2., 'z': 3., 'r': 4.,
            'Integrated Mean Curvature': 0., 'Integrated Gaussian Curvature': -2.,
        })

    result = read_logs2(path, return_dict=return_dict)
    assert result['group data'] == group
    assert result['group data']['Boundary Closed'] is False
    assert isinstance(result['group data']['Missing Edges'], int)

    def row(section):
        table = result[section]
        return table[0] if return_dict else table.iloc[0].to_dict()

    assert {key: row('atoms')[key] for key in breakdown} == breakdown
    assert row('surfs')['Integrated Gaussian Curvature (Face)'] == -3.
    assert row('verts')['Integrated Gaussian Curvature'] == -2.
    assert row('verts')['Integrated Mean Curvature'] == 0.
    assert row('verts')['loc'] == [1., 2., 3.]
    edge = row('edges')
    assert edge['Balls'] == [10, 20, 30]
    for header, value in zip(EDGE_CURVATURE_HEADERS, edge_values):
        if value == '':
            assert header not in edge
        else:
            assert edge[header] == value


def test_legacy_logs_selection_and_missing_vertex_values(tmp_path):
    path = tmp_path / 'legacy.csv'
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        write_section(writer, 'Atoms', {
            'Index': 1, 'Name': 'CA', 'Maximum Curvature': 2.,
            'Average Surface Curvature': 3., 'Neighbors': '[2, 3]',
        })
        writer.writerow([2, 'OW', 0., 0., '[]'][::-1])
        write_section(writer, 'Edges', {
            'Index': 4, 'Ball 1': 1, 'Ball 2': 2, 'Ball 3': 3, 'Length': 2.,
        })
        writer.writerow(['Vertices'])
        writer.writerow(['Index', 'Ball 1', 'Ball 2', 'Ball 3', 'Ball 4', 'x', 'y', 'z', 'r',
                         'Integrated Mean Curvature', 'Integrated Gaussian Curvature'])
        # Interface logs currently leave these trailing curvature values absent.
        writer.writerow([0, 1, 2, 3, 4, 0., 0., 0., 1.])

    result = read_logs2([path, path], return_dict=True, no_sol=True)
    assert list(result) == ['legacy', 'legacy0']
    atom = result['legacy']['atoms'][0]
    assert len(result['legacy']['atoms']) == 1
    assert atom['Maximum Mean Curvature'] == 2.
    assert atom['Average Mean Surface Curvature'] == 3.
    assert atom['Neighbors'] == [2, 3]
    assert result['legacy']['edges'] == [{'Index': 4, 'Balls': [1, 2, 3], 'Length': 2.}]
    assert 'Integrated Gaussian Curvature' not in result['legacy']['verts'][0]
    selected = read_logs2(path, all_=False, edges=True)
    assert isinstance(selected['edges'], pd.DataFrame)
    assert len(selected['edges']) == 1
    assert selected['atoms'].empty and selected['verts'].empty
