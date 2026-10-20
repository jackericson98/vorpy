from pathlib import Path

from vorpy.src.geometry.validation_shapes import spherocylinder_multicell
from vorpy.src.group.export import export_info
from vorpy.src.group.group import Group
from vorpy.src.system.system import System


def _synthetic_group(tmp_path):
    pdb = tmp_path / "synthetic.pdb"
    spherocylinder_multicell(
        10.0, 10.0, interior_count=2, angular_resolution=8,
        cap_resolution=8, axial_ring_count=2,
    ).save_pdb(pdb)
    system = System(file=str(pdb), make_dir=False, print_actions=False)
    group = Group(
        sys=system, name="synthetic_group", chains=[system.chains[0]],
        make_net=True, print_metrics=False,
    )
    group.build()
    group.get_info()
    group.dir = str(tmp_path / "group")
    Path(group.dir).mkdir()
    return group


def test_group_info_preserves_missing_status_and_hides_energy_by_default(tmp_path, monkeypatch):
    group = _synthetic_group(tmp_path)
    monkeypatch.setattr(group, "get_info", lambda: None)
    group.int_mean_curv_face = None
    group.int_mean_curv_edge = None
    group.int_mean_curv = None
    group.int_gauss_curv_face = None
    group.int_gauss_curv_edge = None
    group.int_gauss_curv_vertex = None
    group.int_gauss_curv = None
    group.int_mean_curv_sq = 4.0
    group.boundary_is_complete = False
    monkeypatch.chdir(tmp_path)

    export_info(group)
    report = next(Path(group.dir).glob("*network_info.txt"), None)
    if report is None:
        report = Path(group.dir) / "group_info.txt"
    text = report.read_text(encoding="utf-8")

    assert "SURFACE ENERGY ESTIMATE" not in text
    assert "UNRESOLVED" in text
    assert "NOT_APPLICABLE" in text
    assert "Supported boundary (not certified complete)" in text

    export_info(group, include_surface_energy=True)
    text = report.read_text(encoding="utf-8")
    assert "SURFACE ENERGY ESTIMATE" in text
    assert "Representative Surface Energy:       8.00000 kBT" in text
