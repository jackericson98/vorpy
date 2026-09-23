from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from vorpy.src.output.output import _export_nonpolar_geometry


def test_large_export_geometry_helper_writes_analysis_subdirectory(tmp_path):
    balls = pd.DataFrame([
        {"num": 0, "name": "C", "res_name": "ALA", "res_seq": 1,
         "chain_name": "A", "element": "C", "rad": 1.7,
         "vol": 1., "sa": 2., "complete": True},
        {"num": 1, "name": "C", "res_name": "LEU", "res_seq": 2,
         "chain_name": "B", "element": "C", "rad": 1.7,
         "vol": 1., "sa": 2., "complete": True},
    ])
    network = SimpleNamespace(
        balls=balls,
        surfs=pd.DataFrame([{"balls": [0, 1], "sa": 2., "mean_curv": .1,
                             "gauss_curv": .2}]),
        iface_grps=({0}, {1}),
        group_name="iface",
        sys=SimpleNamespace(balls=balls),
    )
    system = SimpleNamespace(name="sample", frame_index=0)
    _export_nonpolar_geometry(system, network, str(tmp_path), "iface")
    output = Path(tmp_path) / "nonpolar_interface_geometry"
    assert (output / "nonpolar_interface_geometry_surfaces.csv").is_file()
    assert (output / "nonpolar_interface_geometry_atoms.csv").is_file()
    assert (output / "nonpolar_interface_geometry_metadata.json").is_file()
    assert "interface_label" in (output / "nonpolar_interface_geometry_surfaces.csv").read_text()
