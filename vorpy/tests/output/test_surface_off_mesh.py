from types import SimpleNamespace

import numpy as np
import pandas as pd

from vorpy.src.output.mesh import MeshData, write_mesh
from vorpy.src.output.surfs import prepare_surfs


def test_surface_off_contains_faces_and_valid_indices(tmp_path):
    net = SimpleNamespace(
        settings={"net_type": "pow", "surf_scheme": "solid"},
        surfs=pd.DataFrame({
            "points": [np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])],
            "tris": [np.array([[0, 1, 2]])],
        }),
    )
    mesh = prepare_surfs(net, [0])
    assert len(mesh.points) == 3
    assert len(mesh.triangles) == 1
    assert np.all(mesh.triangles >= 0)
    assert np.all(mesh.triangles < len(mesh.points))
    path = write_mesh(mesh, "surfs", directory=tmp_path)
    header = path.read_text().splitlines()[1].split()
    assert int(header[0]) > 0
    assert int(header[1]) > 0


def test_mesh_rejects_out_of_range_face_indices():
    try:
        MeshData(np.zeros((3, 3)), np.array([[0, 1, 3]]))
    except ValueError as error:
        assert "triangle index" in str(error)
    else:
        raise AssertionError("invalid mesh indices were accepted")
