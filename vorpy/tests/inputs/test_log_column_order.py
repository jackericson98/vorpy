import csv

import pytest

from vorpy.src.inputs.logs import read_logs
from vorpy.src.analyze.tools.compare.read_logs import read_logs as read_analysis_logs
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.log_columns import LogWriter


SECTIONS = {
    'build informaiton': {
        'Name': 'sample', 'Location': 'sample.pdb', 'Completion Date': '2026-10-06',
        'Network Type': 'aw', 'Surface Resolution': .2, 'Box Size': 1.25,
        'Maximum Allowable Vertex': 5., 'Maximum Found Vertex': 2.,
        'Total Time': 10., 'Vertex Time': 1., 'Connect Time': 2.,
        'Surface Building Time': 3., 'Analysis time': 4., 'vorPy version': '3.5.2',
    },
    'group information': {
        'Name': 'sample', 'Volume': 20., 'Surface Area': 30., 'Mass': 12.,
        'Density': .6, 'Center of Mass': [1., 2., 3.], 'VDW Volume': 10.,
        'VDW Center of Mass': [4., 5., 6.], 'Moment of Inertia': [[1., 0.], [0., 2.]],
        'Spatial Moment of Inertia': [[3., 0.], [0., 4.]],
        'Integrated Mean Curvature': -1., 'Integrated Gaussian Curvature': 12.,
        'Boundary Components': 2, 'Boundary Closed': False, 'Missing Edges': 3,
    },
    'Atoms': {
        'Index': 0, 'Name': 'CA', 'Residue': 'ALA', 'Residue Sequence': 7,
        'Chain': 'B', 'Mass': 12., 'X': 1., 'Y': 2., 'Z': 3., 'Radius': 1.8,
        'Volume': 5., 'Surface Area': 6., 'Van Der Waals Volume': 4.,
        'Complete Cell?': True, 'Maximum Mean Curvature': -.2,
        'Integrated Mean Curvature': -2., 'Integrated Mean Curvature (Edge)': 0.,
        'Integrated Gaussian Curvature (Vertex)': 12., 'Non-Overlap Volume': 4.,
        'neighbors': [1, 2], 'Number of Neighbors': 2, 'Closest Neighbor': 1,
        'Moment of Inertia Tensor': [[1., 2.], [3., 4.]],
    },
    'Surfaces': {
        'Index': 0, 'Ball 1': 0, 'Ball 2': 1, 'Surface Area': 7.,
        'Mean Curvature': -.1, 'Gaussian Curvature': -.2,
        'Integrated Gaussian Curvature (Face)': -.3,
        'Ball 1 Volume Contribution': 2., 'Ball 2 Volume Contribution': 3.,
        'Contact Area': 4., 'Overlap': 5.,
    },
    'Edges': {
        'Index': 0, 'Ball 1': 0, 'Ball 2': 1, 'Ball 3': 2, 'Length': 8.,
        'Integrated Mean Curvature (Ball 1)': -.4,
        'Integrated Gaussian Curvature (Ball 1, Face with Ball 2)': 0.,
        'Integrated Gaussian Curvature (Ball 2)': '',
    },
    'Vertices': {
        'Index': 0, 'Ball 1': 0, 'Ball 2': 1, 'Ball 3': 2, 'Ball 4': 3,
        'x': 1., 'y': 2., 'z': 3., 'r': 4., 'Integrated Gaussian Curvature': 12.,
    },
}


def write_log(path, writer_type, reverse=False):
    with path.open('w', newline='') as stream:
        writer = writer_type(stream)
        for section, values in SECTIONS.items():
            pairs = list(values.items())
            if reverse:
                pairs.reverse()
            writer.writerow([section])
            writer.writerow([h for h, _ in pairs])
            writer.writerow([v for _, v in pairs])


@pytest.mark.parametrize('reader', [read_logs, read_logs2, read_analysis_logs])
@pytest.mark.parametrize('reverse', [False, True])
def test_reordering_preserves_import_and_analysis_values(tmp_path, reader, reverse):
    original = tmp_path / 'original.csv'
    reordered = tmp_path / 'reordered.csv'
    write_log(original, csv.writer)
    write_log(reordered, LogWriter, reverse=reverse)
    assert reader(str(original), return_dict=True) == reader(str(reordered), return_dict=True)
    result = reader(str(reordered), return_dict=True)
    assert result['data']['max_vertex'] == 2.
    if reader is not read_analysis_logs:
        assert result['group data']['Boundary Components'] == 2
        assert result['surfs'][0]['Ball Volumes'] == [2., 3.]
        assert result['atoms'][0]['Integrated Mean Curvature (Edge)'] == 0.
        assert 'Integrated Gaussian Curvature (Ball 2)' not in result['edges'][0]


def test_writer_rejects_misaligned_row(tmp_path):
    with (tmp_path / 'bad.csv').open('w', newline='') as stream:
        writer = LogWriter(stream)
        writer.writerow(['Atoms'])
        writer.writerow(['Index', 'Name'])
        with pytest.raises(ValueError, match='widths differ'):
            writer.writerow([0])
