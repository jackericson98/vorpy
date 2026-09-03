import xml.etree.ElementTree as ET

import numpy as np

from mesh import MeshData, combine_mesh_parts, write_mesh


def sample_mesh():
    return combine_mesh_parts(
        [np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
         np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1]])],
        [np.array([[0, 1, 2]]), np.array([[0, 1, 2]])],
        [[[1, 0, 0]], [[0, 1, 0]]],
        {'surface_index': [[4], [7]], 'surf_energy': [[2.5], [3.5]]},
    )


def test_combination_offsets_triangle_indices():
    mesh = sample_mesh()
    assert mesh.points.shape == (6, 3)
    assert mesh.triangles.tolist() == [[0, 1, 2], [3, 4, 5]]
    assert mesh.face_data['surface_index'].tolist() == [4, 7]


def test_writers_preserve_counts_colors_and_vtp_metadata(tmp_path):
    mesh = sample_mesh()
    off_path = write_mesh(mesh, 'sample', 'off', tmp_path)
    ply_path = write_mesh(mesh, 'sample', 'ply', tmp_path)
    vtp_path = write_mesh(mesh, 'sample', 'vtp', tmp_path)

    assert off_path.read_text().splitlines()[1] == '6 2 0'
    ply = ply_path.read_text()
    assert 'element vertex 6' in ply
    assert 'element face 2' in ply
    assert 'property uchar red' in ply

    root = ET.parse(vtp_path).getroot()
    piece = root.find('./PolyData/Piece')
    assert piece.attrib['NumberOfPoints'] == '6'
    assert piece.attrib['NumberOfPolys'] == '2'
    names = {item.attrib.get('Name') for item in root.findall('./PolyData/Piece/CellData/DataArray')}
    assert {'RGB', 'surface_index', 'surf_energy'} <= names


def test_mesh_rejects_wrong_color_count():
    try:
        MeshData(np.zeros((3, 3)), np.zeros((1, 3), dtype=int), np.zeros((2, 3)))
    except ValueError:
        pass
    else:
        raise AssertionError('invalid face-color count was accepted')
