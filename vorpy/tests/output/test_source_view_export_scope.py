import csv
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.group.export import group_exports
from vorpy.src.interface import Interface
from vorpy.src.interface.export import interface_exports
from vorpy.src.network.views import NetworkView
from vorpy.src.output.edges import write_edges
from vorpy.src.output.export_scope import ExportScope, export_scope
from vorpy.src.output.mesh import write_prepared_mesh
from vorpy.tests.network.test_views import _network


def _register(group, interface_id, interface=None, other_group=None, side=None,
              interface_name=None):
    group.interfaces[interface_id] = interface
    group.interface_metadata[interface_id] = {
        "interface_id": interface_id,
        "interface_name": interface_name,
        "group_name": group.name,
        "other_group_name": getattr(other_group, "name", None),
        "side": side,
    }


def _case(tmp_path, scheme, solvent_mode="SOLVENT_CONTEXT"):
    source = _network(scheme=scheme)
    source.balls["vol"] = [1.0] * len(source.balls)
    source.settings["solvent_mode"] = solvent_mode
    source.solvent_mode = solvent_mode
    source.refresh_solvent_provenance()

    # Add ligand-water, water-water, and DNA-water faces. They must remain
    # available to the source archive but never enter the A/B export.
    water_surface = source.surfs.iloc[1].copy()
    water_surface["balls"] = [2, 3]
    dna_water_surface = source.surfs.iloc[1].copy()
    dna_water_surface["balls"] = [1, 2]
    source.surfs = pd.concat(
        [source.surfs,
         pd.DataFrame([water_surface, dna_water_surface], index=[2, 3])]
    )

    group_a = SimpleNamespace(
        name="ligand", group_id="ligand", ball_ndxs=[10],
        settings=dict(source.settings), net=None, interfaces={},
        interface_metadata={},
    )
    group_b = SimpleNamespace(
        name="dna", group_id="dna", ball_ndxs=[11],
        settings=dict(source.settings), net=None, interfaces={},
        interface_metadata={},
    )
    group_a.register_interface = lambda **kwargs: _register(group_a, **kwargs)
    group_b.register_interface = lambda **kwargs: _register(group_b, **kwargs)
    system = SimpleNamespace(
        name="synthetic", frame_index=1, groups=[group_a, group_b], ifaces=[],
        balls=source.balls, round_to=3, files={"dir": str(tmp_path), "base_file": None},
    )
    group_a.sys = system
    group_b.sys = system
    interface = Interface(system, group_a, group_b)
    interface.net = source
    interface.source_network = source
    interface.source_view = source.make_view(iface_grps=([10], [11]))
    interface.group1_indices = {10}
    interface.group2_indices = {11}
    interface.group2_name = "dna"
    interface.dir = str(tmp_path / "interface")
    system.ifaces.append(interface)
    return source, group_a, interface


@pytest.mark.parametrize("scheme", ["aw", "pow"])
@pytest.mark.parametrize("solvent_mode", ["SOLVENT_CONTEXT", "MOLECULE_ONLY"])
def test_materialized_interface_exports_only_selected_source_geometry(tmp_path, scheme, solvent_mode):
    source, _, interface = _case(tmp_path, scheme, solvent_mode)

    assert interface.source_view.source_ids("surface") == (0,)
    assert len(source.surfs) == 4

    scope = export_scope(interface)
    assert scope.source_ids("surface") == (0,)
    assert scope.source_ids("edge") == tuple(
        source.edges.index[index] for index in source.surfs.iloc[0]["edges"]
    )
    assert scope.source_ids("vertex") == tuple(
        source.verts.index[index] for index in source.surfs.iloc[0]["verts"]
    )

    interface_exports(interface, surfs=True, edges=True, verts=True, info=True, logs=True)

    surface_header = (tmp_path / "interface" / "interface_surfs.off").read_text().splitlines()[1]
    assert surface_header == "4 1 0"
    assert (tmp_path / "interface" / "interface_edges.off").is_file()
    assert (tmp_path / "interface" / "interface_verts.off").is_file()
    for name in ("interface_edges.off", "interface_verts.off"):
        counts = (tmp_path / "interface" / name).read_text(encoding="ascii").splitlines()[1].split()
        assert len(counts) == 3 and all(int(value) >= 0 for value in counts)
    assert not (tmp_path / "interface" / "surfs.off").exists()
    assert not (tmp_path / "interface" / "visualization_manifest.json").exists()

    info = (tmp_path / "interface" / "interface_info.txt").read_text(encoding="utf-8")
    assert "Physical export selection: selected_source_view" in info
    assert "Direct Group 1–Group 2 surfaces: 1" in info
    assert "Supporting/unclassified surfaces: 0" in info

    with (tmp_path / "interface" / "interface_surfaces.csv").open(newline="", encoding="utf-8") as stream:
        surface_rows = list(csv.DictReader(stream))
    assert len(surface_rows) == 1
    assert surface_rows[0]["Surface ID"] == "0"
    assert surface_rows[0]["Group A Generator ID"] == "atom-A"
    assert surface_rows[0]["Group B Generator ID"] == "atom-B"
    assert surface_rows[0]["Surface Area (Å²)"] == ""

    with (tmp_path / "interface" / "interface_edges.csv").open(newline="", encoding="utf-8") as stream:
        edge_rows = list(csv.DictReader(stream))
    with (tmp_path / "interface" / "interface_vertices.csv").open(newline="", encoding="utf-8") as stream:
        vertex_rows = list(csv.DictReader(stream))
    assert {row["Edge ID"] for row in edge_rows} == {"0", "1", "2", "3"}
    assert {row["Vertex ID"] for row in vertex_rows} == {"0", "1", "2", "3"}
    assert all("2" not in row["Incident Selected Interface Surfaces"] for row in edge_rows)

    logs = (tmp_path / "interface" / f"{scheme}_logs.csv").read_text(encoding="utf-8")
    assert "0,0,1," in logs
    assert "1,0,2," not in logs
    assert "2,2,3," not in logs
    assert "3,1,2," not in logs
    assert interface.net is source

    from vorpy.src.output.visualization import _compact_physical_interface
    compact = _compact_physical_interface(interface, tmp_path / "compact")
    assert compact["surface_count"] == 1
    assert compact["selection_mode"] == "selected_source_view"
    assert (tmp_path / "compact" / "interface_surfs.off").is_file()
    assert (tmp_path / "compact" / "interface_edges.off").is_file()
    assert (tmp_path / "compact" / "interface_verts.off").is_file()
    vertex_lines = (tmp_path / "compact" / "interface_verts.off").read_text().splitlines()
    assert vertex_lines[0] == "OFF"
    vertex_count, face_count, _ = (int(value) for value in vertex_lines[1].split()[:3])
    assert vertex_count > 0 and face_count > 0


def test_empty_interface_exports_valid_empty_layers_and_unresolved_info(tmp_path):
    source, _, interface = _case(tmp_path, "aw", "MOLECULE_ONLY")
    interface.group1_indices = {13}
    interface.group2_indices = {10}
    interface.source_view = None
    interface.source_network = None
    interface.dir = str(tmp_path / "empty_interface")

    interface_exports(
        interface,
        surfs=True,
        edges=True,
        verts=True,
        logs=True,
        info=True,
    )

    for name in ("interface_surfs.off", "interface_edges.off", "interface_verts.off"):
        lines = (tmp_path / "empty_interface" / name).read_text(encoding="ascii").splitlines()
        assert lines[:2] == ["OFF", "0 0 0"]

    info = (tmp_path / "empty_interface" / "interface_info.txt").read_text(encoding="utf-8")
    assert "Completeness: NOT_CALCULATED" in info
    assert "No direct A/B physical facets selected" in info
    assert "Interface area (Å²): Not available" in info
    assert not ("Completeness: COMPLETE" in info)

    for name in (
        "interface_surfaces.csv", "interface_edges.csv", "interface_vertices.csv",
    ):
        assert (tmp_path / "empty_interface" / name).is_file()


def test_classic_log_keeps_incomplete_cells_and_blank_unavailable_curvature(tmp_path):
    source, _, interface = _case(tmp_path, "aw")
    source.balls["sa"] = [1.0] * len(source.balls)
    source.balls.loc[1, "complete"] = False

    interface_exports(interface, logs=True)

    with (tmp_path / "interface" / "aw_logs.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.reader(stream))
    header = rows[rows.index(["Atoms"]) + 1]
    atom_rows = [dict(zip(header, row)) for row in rows[rows.index(header) + 1:]
                 if row and len(row) == len(header)]
    incomplete = next(row for row in atom_rows if row["Index"] == "1")
    assert incomplete["Complete Cell?"] == "False"
    assert incomplete["Maximum Mean Curvature"] == ""
    assert incomplete["Integrated Gaussian Curvature"] == ""


def test_interface_rechecks_broad_materialized_view_before_export(tmp_path):
    source, _, interface = _case(tmp_path, "pow", "SOLVENT_CONTEXT")
    # A group view contains the A-internal surface as well as the A/B facet.
    interface.source_view = source.make_view(group_ids=[10])

    interface_exports(interface, surfs=True, edges=True, verts=True, info=True)

    assert (tmp_path / "interface" / "interface_surfs.off").read_text(encoding="ascii").splitlines()[1] == "4 1 0"
    info = (tmp_path / "interface" / "interface_info.txt").read_text(encoding="utf-8")
    assert "Direct Group 1–Group 2 surfaces: 1" in info


@pytest.mark.parametrize("scheme", ["aw", "pow"])
def test_partial_materialized_interface_exports_supported_subset(tmp_path, scheme):
    source, _, interface = _case(tmp_path, scheme)
    source.surfs["sa"] = [2.5, 8.0, 9.0, 10.0]
    unresolved = source.surfs.iloc[1].copy()
    unresolved["balls"] = [0, 1, 3]
    unresolved["sa"] = 1.5
    source.surfs.loc[4] = unresolved
    source.balls.loc[3, "complete"] = False
    interface.source_view = NetworkView.for_interface(source, [10], [11])

    scope = export_scope(interface)
    assert scope.source_ids("surface") == (0,)
    assert interface.source_view.completeness.proven is False

    interface_exports(interface, surfs=True, edges=True, verts=True, info=True, logs=True)

    directory = tmp_path / "interface"
    assert (directory / "interface_surfs.off").read_text(encoding="ascii").splitlines()[1] == "4 1 0"
    with (directory / "interface_surfaces.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert [row["Surface ID"] for row in rows] == ["0"]
    assert rows[0]["Surface Area (Å²)"] == "2.5"

    info = (directory / "interface_info.txt").read_text(encoding="utf-8")
    assert "Completeness: PARTIAL" in info
    assert "Exported physical subset:" in info
    assert "Missing physical surfaces: 1" in info
    assert "Missing surface IDs: 4" in info
    assert "Missing physical area (Å²): 1.500000" in info
    assert "Interface area (Å²): 2.500000" in info


def test_partial_materialized_interface_rejects_incomplete_retained_feature(tmp_path):
    source, _, interface = _case(tmp_path, "pow")
    unresolved = source.surfs.iloc[1].copy()
    unresolved["balls"] = [0, 1, 3]
    source.surfs.loc[4] = unresolved
    source.balls.loc[3, "complete"] = False
    interface.source_view = NetworkView.for_interface(source, [10], [11])
    source.surfs.at[0, "complete"] = False

    with pytest.raises(ValueError, match="locally incomplete materialized surface"):
        export_scope(interface)


def test_interface_preserves_area_and_leaves_partial_curvature_blank(tmp_path):
    source, _, interface = _case(tmp_path, "aw", "SOLVENT_CONTEXT")
    source.surfs["sa"] = [2.5, 8.0, 9.0, 10.0]
    source.surfs["mean_curv"] = [0.25, 0.5, 0.5, 0.5]
    source.surfs["gauss_curv"] = [float("nan")] * len(source.surfs)

    interface_exports(interface, surfs=True, info=True, logs=True)

    with (tmp_path / "interface" / "interface_surfaces.csv").open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream))
    assert row["Surface Area (Å²)"] == "2.5"
    assert row["Mean Curvature (Å⁻¹)"] == "0.25"
    assert row["Gaussian Curvature (Å⁻²)"] == ""
    info = (tmp_path / "interface" / "interface_info.txt").read_text(encoding="utf-8")
    assert "Interface area (Å²): 2.500000" in info


def test_optional_interface_profile_reports_scope_counts_and_stage_timings(tmp_path, monkeypatch):
    _, _, interface = _case(tmp_path, "pow", "SOLVENT_CONTEXT")
    monkeypatch.setenv("VORPY_INTERFACE_EXPORT_PROFILE", "1")

    profile = interface_exports(
        interface,
        surfs=True,
        edges=True,
        verts=True,
        info=True,
        logs=True,
    )

    assert profile is interface.interface_export_profile
    assert profile["counts"]["source_generators"] == 5
    assert profile["counts"]["source_surfaces"] == 4
    assert profile["counts"]["selected_surfaces"] == 1
    assert profile["counts"]["selected_surface_points"] == 4
    assert profile["counts"]["selected_surface_triangles"] == 1
    assert profile["counts"]["selected_edges"] == 4
    assert profile["counts"]["selected_vertices"] == 4
    assert profile["measurements"]["selected_surface_area_A2"] is None
    assert profile["timings"]["scope_selection"] >= 0.0
    assert profile["timings"]["surface_mesh"] >= 0.0
    assert profile["timings"]["total"] >= profile["timings"]["scope_selection"]
    assert not (tmp_path / "interface" / "interface_export_profile.json").exists()


def test_mesh_only_interface_export_does_not_materialize_info_proxy(tmp_path, monkeypatch):
    _, _, interface = _case(tmp_path, "aw")

    monkeypatch.setattr(
        ExportScope,
        "scoped_network",
        lambda self: pytest.fail("mesh-only export materialized an info proxy"),
    )
    interface_exports(interface, surfs=True)


def test_interface_info_survives_missing_optional_dual_records(tmp_path, monkeypatch):
    _, _, interface = _case(tmp_path, "aw")
    interface.source_view = None
    interface.source_network = None
    interface._geometry_representation = object()

    import importlib
    dual_module = importlib.import_module("vorpy.src.geometry.visualization.interface")
    def missing_dual(*args, **kwargs):
        raise KeyError("missing simplex")

    monkeypatch.setattr(dual_module, "write_interface_dual_summary", missing_dual)
    interface_exports(interface, info=True)

    info = (tmp_path / "interface" / "interface_info.txt").read_text(encoding="utf-8")
    assert "Apollonius dual visualization: unavailable (missing simplex)" in info


@pytest.mark.parametrize("scheme", ["aw", "pow"])
def test_materialized_group_export_uses_group_view_surface_ids(tmp_path, scheme):
    source, group_a, _ = _case(tmp_path, scheme)
    group_a.source_network = source
    group_a.source_view = source.make_view(group_ids=[10])
    group_a.source_view_analysis = SimpleNamespace(area=2.0, volume=1.0)
    group_a.dir = str(tmp_path / "group")

    scope = export_scope(group_a)
    assert scope.source_ids("surface") == (0, 1)
    group_exports(group_a, surfs=True, edges=True, verts=True, info=True)

    surface_header = (tmp_path / "group" / "surfs.off").read_text().splitlines()[1]
    assert surface_header == "8 2 0"
    info = (tmp_path / "group" / "network_info.txt").read_text(encoding="utf-8")
    assert "Source view ID:" in info
    assert "  0: generators=[atom-A, atom-B, boundary:4]" in info
    assert "  2:" not in info
    assert "  3:" not in info


@pytest.mark.parametrize("scheme", ["aw", "pow"])
def test_dedicated_interface_exports_only_direct_physical_facets(tmp_path, scheme):
    source, _, interface = _case(tmp_path, scheme)
    interface.source_network = None
    interface.source_view = None
    interface.dir = str(tmp_path / "dedicated_interface")

    interface_exports(interface, surfs=True, edges=True, verts=True, info=True, logs=True)

    surface_header = (tmp_path / "dedicated_interface" / "interface_surfs.off").read_text().splitlines()[1]
    assert surface_header == "4 1 0"
    with (tmp_path / "dedicated_interface" / "interface_surfaces.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert [row["Surface ID"] for row in rows] == ["0"]
    assert "2,2,3," not in (tmp_path / "dedicated_interface" / f"{scheme}_logs.csv").read_text(encoding="utf-8")


def test_oversized_mesh_is_skipped_before_combined_allocation(tmp_path):
    source = _network(scheme="aw")
    source.export_memory_budget_bytes = 1

    assert write_edges(
        source, [0], "edges", directory=tmp_path,
    ) is None
    assert not (tmp_path / "edges.off").exists()
    assert source._export_mesh_skips[0]["component"] == "edges"


def test_optional_mesh_skip_preserves_source_archive_and_summary(tmp_path):
    source, _, interface = _case(tmp_path, "aw")
    source.export_memory_budget_bytes = 1
    archive = tmp_path / "synthetic.vpy"
    archive.write_bytes(b"already-saved-network")

    interface_exports(interface, surfs=True, edges=True, verts=True, info=True, logs=True)

    assert archive.read_bytes() == b"already-saved-network"
    assert (tmp_path / "interface" / "network_info.txt").is_file()
    assert (tmp_path / "interface" / "aw_logs.csv").is_file()


def test_allocation_memory_error_is_reported_as_optional_mesh_skip(tmp_path):
    source = _network(scheme="aw")

    def fail_prepare():
        raise MemoryError("synthetic allocation failure")

    assert write_prepared_mesh(
        source, fail_prepare, "surfs", "off", tmp_path, 100,
    ) is None
    assert source._export_mesh_skips[-1] == {
        "component": "surfs", "reason": "MemoryError",
    }


def test_failed_optional_compact_mesh_preserves_saved_archive(tmp_path, monkeypatch):
    source, _, interface = _case(tmp_path, "aw")
    archive = tmp_path / "synthetic.vpy"
    archive.write_bytes(b"already-saved-network")

    from vorpy.src.output import visualization
    def fail_mesh(*args, **kwargs):
        raise MemoryError("synthetic")

    monkeypatch.setattr(visualization, "_surface_mesh", fail_mesh)
    result = visualization._compact_physical_interface(interface, tmp_path / "compact")

    assert result["state"] == "UNRESOLVED"
    assert archive.read_bytes() == b"already-saved-network"


@pytest.mark.parametrize("scheme", ["aw", "pow"])
@pytest.mark.parametrize("solvent_mode", ["SOLVENT_CONTEXT", "MOLECULE_ONLY"])
def test_interface_atom_exports_use_selected_contact_generators(tmp_path, monkeypatch, scheme, solvent_mode):
    source, _, interface = _case(tmp_path, scheme, solvent_mode)
    calls = []

    from vorpy.src.interface import export as interface_export
    monkeypatch.setattr(
        interface_export,
        "write_pdb",
        lambda atoms, file_name, sys, directory=None: calls.append((file_name, list(atoms))),
    )

    interface_exports(interface, atoms=True)

    selected = dict(calls)
    assert selected["group_a_atoms"] == [10]
    assert selected["group_b_atoms"] == [11]
    assert selected["surr_atoms"] == []
    assert selected["surr_residues"] == []
    assert "interface_atoms" not in selected


def test_compact_visualization_preflights_selected_physical_edges(tmp_path):
    source, _, interface = _case(tmp_path, "aw")
    source.export_memory_budget_bytes = 1

    from vorpy.src.output.visualization import _compact_physical_interface
    result = _compact_physical_interface(interface, tmp_path / "compact")

    assert result["surface_count"] == 1
    assert not (tmp_path / "compact" / "interface_edges.off").exists()
