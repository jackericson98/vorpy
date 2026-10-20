from pathlib import Path
from types import SimpleNamespace

import pytest

from vorpy.src.command.calculation_plan import build_calculation_plan
from vorpy.src.command.group import ggroup
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.interface import Interface
from vorpy.src.network.network import Network
from vorpy.src.system.system import System

from vorpy.tests.network.test_views import _network


DATA = Path(__file__).resolve().parents[2] / "data"


class _Group:
    def __init__(self, name, system_ids):
        self.name = name
        self.group_id = name
        self.ball_ndxs = list(system_ids)
        self.settings = {"net_type": "aw"}
        self.net = None
        self.interfaces = {}
        self.interface_metadata = {}

    def register_interface(self, interface_id, interface=None, other_group=None,
                           side=None, interface_name=None):
        self.interfaces[interface_id] = interface
        self.interface_metadata[interface_id] = {
            "interface_id": interface_id,
            "interface_name": interface_name,
            "group_name": self.name,
            "other_group_id": getattr(other_group, "group_id", None),
            "other_group_name": getattr(other_group, "name", None),
            "side": side,
            "status": "defined",
            "built": False,
        }


def _command_with_source(source, groups):
    system = SimpleNamespace(
        name="synthetic",
        frame_index=1,
        groups=list(groups),
        ifaces=[],
        files={"dir": None},
    )
    command = Command(sys=system)
    command.base_file = Path("synthetic.pdb")
    command.calculation_plan = build_calculation_plan(
        [
            "synthetic.pdb",
            "-g", "c", "0",
            "-g", "c", "1",
            "-i",
            "-c", "system",
        ],
        input_source="synthetic.pdb",
    )
    command.groups = {0: [["c", "0"]], 1: [["c", "1"]]}
    command._full_system_group = SimpleNamespace(net=source)
    return command


@pytest.mark.parametrize("scheme", ["aw", "prm"])
def test_full_system_groups_and_interfaces_reuse_one_source_network(monkeypatch, scheme):
    source = _network(scheme=scheme)
    source.balls["vol"] = [1.0] * len(source.balls)
    source.settings["solvent_mode"] = "SOLVENT_CONTEXT"
    source.solvent_mode = "SOLVENT_CONTEXT"
    source.refresh_solvent_provenance()

    group_a = _Group("A", [10])
    group_b = _Group("B", [11])
    command = _command_with_source(source, [group_a, group_b])

    def fail_build(*args, **kwargs):
        raise AssertionError("source-derived full-system views rebuilt geometry")

    monkeypatch.setattr(Network, "build", fail_build)
    source_tables = (source.verts, source.edges, source.surfs)

    command._attach_full_system_group_views()
    command._attach_full_system_interface_views([(group_a, group_b)])

    assert all(group.net is None for group in (group_a, group_b))
    assert all(group.source_view.source_network_id == source.source_network_id
               for group in (group_a, group_b))
    assert tuple(group_a.ball_ndxs) == (10,)
    assert tuple(group_b.ball_ndxs) == (11,)
    assert source.generator_system_indices == (10, 11, 12, 13, None)
    assert source.verts is source_tables[0]
    assert source.edges is source_tables[1]
    assert source.surfs is source_tables[2]

    interface = command.sys.ifaces[0]
    assert interface.net is source
    assert interface.source_view.selection_kind == "interface"
    assert interface.source_view.source_ids("surface") == (0,)
    assert interface.geometry_analysis.metadata["source_geometry"] == (
        "full_system_network_view"
    )
    assert interface.geometry_analysis.metadata["geometry_status"] == (
        "GEOMETRY_COMPLETE"
    )
    assert interface.geometry_analysis.metadata["analysis_status"] == (
        "NOT_CALCULATED"
    )
    assert interface.geometry_analysis.metadata["alpha_status"] == "NOT_CALCULATED"
    assert interface.geometry_analysis.alpha_selection["status"] == "not_calculated"
    assert interface.geometry_analysis.metadata["additional_network_solves"] == 0
    assert interface.geometry_analysis.voronoi_side_1["source_surface_ids"] == (0,)
    assert interface.geometry_analysis.voronoi_side_1["source_generator_pairs"] == (
        (10, 11),
    )
    assert interface.geometry_analysis.voronoi_side_1["topology_faces"] == 1
    assert group_a.interface_metadata[interface.interface_id]["network_created"] is False
    assert group_a.interface_metadata[interface.interface_id]["source_view_id"] == (
        interface.source_view.selection_id
    )


def test_interface_only_does_not_materialize_unrelated_group_views():
    source = _network(scheme="aw")
    group_a = _Group("A", [10])
    group_b = _Group("B", [11])
    command = _command_with_source(source, [group_a, group_b])
    command.calculation_plan = build_calculation_plan(
        ["synthetic.pdb", "-g", "c", "0", "-g", "c", "1", "-i"],
        input_source="synthetic.pdb",
    )

    command._attach_full_system_group_views()

    assert not hasattr(group_a, "source_view")
    assert not hasattr(group_b, "source_view")


def test_failed_source_interface_is_reported_without_discarding_checkpoint():
    source = _network(scheme="aw")
    source.balls.loc[0, "complete"] = False
    group_a = _Group("A", [10])
    group_b = _Group("B", [11])
    command = _command_with_source(source, [group_a, group_b])

    command._attach_full_system_interface_views([(group_a, group_b)])

    assert command.sys.ifaces == []
    assert len(command.interface_result_errors) == 1
    assert "generator cells 0" in command.interface_result_errors[0]


def test_partial_source_interface_preserves_unresolved_surface_area():
    source = _network(scheme="aw")
    source.surfs["sa"] = 0.0
    extra = source.surfs.iloc[1].copy()
    extra["balls"] = [0, 1, 2]
    extra["sa"] = 1.5
    source.surfs.loc[2] = extra
    source.balls.loc[2, "complete"] = False
    group_a = _Group("A", [10])
    group_b = _Group("B", [11])
    command = _command_with_source(source, [group_a, group_b])

    command._attach_full_system_interface_views([(group_a, group_b)])

    interface = command.sys.ifaces[0]
    view = interface.source_view
    analysis = interface.geometry_analysis
    assert view.source_ids("surface") == (0,)
    assert view.completeness.proven is False
    assert view.completeness.unresolved_surface_ids == (2,)
    assert analysis.voronoi_side_1["area_selection_complete"] is False
    assert analysis.voronoi_side_1["unresolved_surface_ids"] == (2,)
    assert analysis.voronoi_side_1["unresolved_surface_area"] == pytest.approx(1.5)
    assert analysis.metadata["result_status"] == "PARTIAL"


@pytest.mark.parametrize("scheme", ["aw", "pow"])
def test_interface_only_pipeline_builds_no_full_system_or_group_networks(
    monkeypatch, scheme, capsys, tmp_path
):
    source = _network(scheme=scheme)
    source.balls["vol"] = [1.0] * len(source.balls)
    source.settings["solvent_mode"] = "SOLVENT_CONTEXT"
    source.solvent_mode = "SOLVENT_CONTEXT"
    source.refresh_solvent_provenance()
    group_a = _Group("A", [10])
    group_b = _Group("B", [11])
    system = SimpleNamespace(
        name="synthetic",
        frame_index=1,
        balls=source.balls,
        groups=[group_a, group_b],
        ifaces=[],
        files={"dir": tmp_path},
    )
    command = Command(sys=system)
    command.base_file = Path("synthetic.input")
    command.groups = {0: [["c", "0"]], 1: [["c", "1"]]}
    command.settings_dict = {"net_type": scheme}
    command.interface_mode = True
    command._build_group_networks = False
    command.calculation_plan = build_calculation_plan(
        ["synthetic.input", "-g", "c", "0", "-g", "c", "1", "-i"],
        input_source="synthetic.input",
    )
    builds = []
    interface_calls = []

    class FullGroup:
        def __init__(self, **kwargs):
            self.net = source

        def build(self, **kwargs):
            builds.append(kwargs)

    monkeypatch.setattr("vorpy.src.group.Group", FullGroup)
    monkeypatch.setattr(Network, "build", lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("a derived view called Network.build")
    ))
    monkeypatch.setattr(Interface, "build", lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("a derived interface called Interface.build")
    ))
    command.load_files = lambda: None
    command.apply_settings = lambda: None
    command.create_groups = lambda: None
    system.make_interfaces = lambda pairs, *, analyze_waters: interface_calls.append(
        (pairs, analyze_waters)
    )
    command.run_exports = lambda: None
    command._save_archive = lambda: None
    command._make_cli_interface_results = lambda pairs: ()
    command._attach_full_system_group_views = lambda: pytest.fail(
        "interface-only source reuse materialized a group view"
    )
    monkeypatch.setattr(command, "_persist_geometry_checkpoint", lambda *args, **kwargs: None)

    command._process_system()

    output = capsys.readouterr().out

    assert builds == []
    assert len(interface_calls) == 1
    assert all(group.net is None for group in (group_a, group_b))
    assert system.ifaces == []
    assert f"Solving {scheme.upper()} interface network..." in output
    assert "Extracting " not in output
    assert '"A_B_interface" built - ' not in output
    assert "Exporting results..." in output
    assert "Exporting interface results..." not in output
    assert "Interface Result status:" not in output
    assert "<vorpy." not in output
    assert "0x" not in output


def test_canonical_1aay_groups_keep_membership_and_zinc_without_geometry_build():
    settings = {"net_type": "aw"}
    baseline = System(file=str(DATA / "1AAY_no_sol.pdb"), make_dir=False)
    ggroup(
        baseline,
        {0: [["dna"]], 1: [["protein"]]},
        settings=settings,
        make_net=False,
    )
    expected = [tuple(group.ball_ndxs) for group in baseline.groups]

    system = System(file=str(DATA / "1AAY_no_sol.pdb"), make_dir=False)
    command = Command(sys=system)
    command.groups = {0: [["dna"]], 1: [["protein"]]}
    command.settings_dict = settings.copy()
    command.interface_mode = True
    command._build_group_networks = False
    command.create_groups()

    assert [tuple(group.ball_ndxs) for group in system.groups] == expected
    assert sum(
        str(value).upper() == "ZN"
        for value in system.balls["res_name"]
    ) == 3
    assert all(group.net is None for group in system.groups)


def test_full_system_non_interface_groups_do_not_create_subgroup_networks(monkeypatch):
    source = _network(scheme="aw")
    command = _command_with_source(
        source,
        [_Group("A", [10]), _Group("B", [11])],
    )
    command.interface_mode = False
    calls = []
    monkeypatch.setattr(
        "vorpy.src.command.vpy_cmnd.ggroup",
        lambda *args, make_net, **kwargs: calls.append(make_net),
    )

    command.create_groups()

    assert calls == [False]
