import pytest
from types import SimpleNamespace

from vorpy.src.command.interpret import get_set, get_val
from vorpy.src.command.set import sett
from vorpy.src.group import export as group_export


def test_edge_width_and_vertex_size_are_cli_settings():
    settings = sett("ew", ["0.12"])
    settings = sett("vs", ["0.24"], settings)

    assert settings["edge_width"] == pytest.approx(0.12)
    assert settings["vertex_size"] == pytest.approx(0.24)
    assert settings["max_vert"] == 40
    assert get_set("ew") == "ew"
    assert get_set("vs") == "vs"
    assert get_val("ew", ["0.3"]) == pytest.approx(0.3)
    assert get_val("vs", ["0.4"]) == pytest.approx(0.4)


def test_vertex_and_edge_color_cli_settings_are_applied():
    settings = sett("vc", ["blue"])
    settings = sett("ec", ["orange"], settings)

    assert settings["vert_col"] == "blue"
    assert settings["edge_col"] == "orange"


@pytest.mark.parametrize("setting", ["ew", "vs"])
@pytest.mark.parametrize("value", ["0", "-0.1", "nan", "inf", "bad"])
def test_mesh_width_settings_reject_nonpositive_or_nonfinite_values(setting, value):
    settings = sett(setting, [value])
    key, default = ("edge_width", 0.05) if setting == "ew" else ("vertex_size", 0.08)
    assert settings[key] == pytest.approx(default)


def test_group_exports_apply_configured_edge_and_vertex_sizes(tmp_path, monkeypatch):
    system_dir = tmp_path / "system"
    group_dir = system_dir / "group"
    group_dir.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    settings = sett("vc", ["blue"])
    settings = sett("ec", ["orange"], settings)
    settings.update({
            "file_type": "off",
            "surf_col": "plasma",
            "surf_scheme": None,
            "edge_width": 0.17,
            "vertex_size": 0.29,
        })
    group = SimpleNamespace(
        settings=settings,
        net=SimpleNamespace(
            settings={"surf_scheme": None}, balls=None, surfs=[0],
            edges=[], verts=[],
        ),
        sys=SimpleNamespace(files={"dir": str(system_dir)}),
        dir=str(group_dir),
        ball_ndxs=[],
        layer_edges=None,
        layer_verts=None,
    )
    observed = {}
    monkeypatch.setattr(
        group_export, "write_edges",
        lambda *args, **kwargs: observed.update(
            edge_radius=kwargs["radius"], edge_color=kwargs["color"],
            edge_scheme=kwargs["color_scheme"],
        ),
    )
    monkeypatch.setattr(
        group_export, "write_off_verts",
        lambda *args, **kwargs: observed.update(
            vertex_radius=kwargs["vert_rad"], vertex_color=kwargs["color"],
            vertex_scheme=kwargs["color_scheme"],
        ),
    )

    group_export.group_exports(group, edges=True, verts=True)

    assert observed == {
        "edge_radius": 0.17,
        "vertex_radius": 0.29,
        "edge_color": "orange",
        "vertex_color": "blue",
        "edge_scheme": "fixed",
        "vertex_scheme": "fixed",
    }
