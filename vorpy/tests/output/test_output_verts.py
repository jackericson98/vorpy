import numpy as np
import pandas as pd

from vorpy.src.output.verts import write_off_verts


class DummyNet:
    def __init__(self):
        self.verts = pd.DataFrame(
            {
                "loc": [
                    np.array([0.0, 0.0, 0.0]),
                    np.array([1.0, 0.0, 0.0]),
                ]
            }
        )


def test_write_off_verts_header_matches_written_geometry(tmp_path):
    net = DummyNet()

    write_off_verts(
        net,
        [0, 1],
        "shell_verts",
        directory=str(tmp_path),
        color="red",
    )

    lines = (tmp_path / "shell_verts.off").read_text().splitlines()

    assert lines[0] == "OFF"
    n_vertices, n_faces, n_edges = map(int, lines[1].split())

    body = [line for line in lines[2:] if line.strip()]
    vertex_lines = [line for line in body if len(line.split()) == 3]
    face_lines = [line for line in body if line.split()[0] == "3" and len(line.split()) >= 4]

    # Each displayed Voronoi vertex is an octahedron:
    # 6 coordinate vertices and 8 triangular faces.
    assert n_vertices == 12
    assert n_faces == 16
    assert n_edges == 0

    assert len(vertex_lines) == n_vertices
    assert len(face_lines) == n_faces
