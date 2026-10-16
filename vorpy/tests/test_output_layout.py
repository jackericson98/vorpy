import csv
from pathlib import Path
from types import SimpleNamespace

from vorpy.src.interface.interface import Interface
from vorpy.src.output.layout import interface_folder_name
from vorpy.src.output.output import _set_group_directory
from vorpy.src.output.set_pymol_atoms import set_pymol_atoms
from vorpy.src.results import (
    INTERFACE_LOG_COLUMNS,
    Provenance,
    RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG,
    RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG,
    Status,
)


def test_full_system_network_uses_root_and_selected_network_gets_one_folder(tmp_path):
    system = SimpleNamespace(files={"dir": str(tmp_path)}, balls=[0, 1, 2], groups=[])
    full = SimpleNamespace(name="2KAI", ball_ndxs=[0, 1, 2], dir=None)
    selected = SimpleNamespace(name="group_1", ball_ndxs=[0], dir=None)
    system.groups[:] = [full, selected]

    _set_group_directory(system, full)
    _set_group_directory(system, selected)

    assert full.dir == str(tmp_path)
    assert selected.dir == str(tmp_path / "group_1")
    assert not (tmp_path / "system").exists()


def test_interface_directory_is_named_and_collision_safe(tmp_path):
    first_group = SimpleNamespace(name="group_1")
    second_group = SimpleNamespace(name="group_2")
    system = SimpleNamespace(files={"dir": str(tmp_path)}, ifaces=[])

    first = Interface.__new__(Interface)
    first.sys, first.group1, first.group2, first.dir = system, first_group, second_group, None
    Interface.set_dir(first)
    second = Interface.__new__(Interface)
    second.sys, second.group1, second.group2, second.dir = system, first_group, second_group, None
    system.ifaces.append(first)
    Interface.set_dir(second)

    assert first.dir.endswith(interface_folder_name(first_group, second_group))
    assert second.dir.endswith("interface_group_1_group_2_2")
    assert (tmp_path / "interface_group_1_group_2").is_dir()
    assert (tmp_path / "interface_group_1_group_2_2").is_dir()


def test_system_pymol_script_loads_root_pdb(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    system = SimpleNamespace(
        type="mol", name="2KAI", residues=[], element_radii={},
        balls=[],
    )

    set_pymol_atoms(system)

    assert (tmp_path / "set_atoms.pml").read_text().splitlines()[0] == "load 2KAI.pdb, 2KAI"


def test_system_exports_stay_at_root(tmp_path, monkeypatch):
    from vorpy.src.output import sys as output_sys

    system = SimpleNamespace(
        files={"dir": str(tmp_path)}, name="2KAI", type="balls", balls=[0],
        residues=[], element_radii={},
    )

    monkeypatch.setattr(
        output_sys,
        "write_pdb",
        lambda atoms, file_name, sys: Path(sys.files["dir"], f"{file_name}.pdb").write_text(
            "ATOM\n", encoding="utf-8"
        ),
    )
    monkeypatch.setattr(
        output_sys,
        "export_sys_info",
        lambda sys: Path("system_info.txt").write_text("system\n", encoding="utf-8"),
    )
    output_sys.export_sys(system, pdb=True, set_atoms=True, info=True)

    assert (tmp_path / "2KAI.pdb").is_file()
    assert (tmp_path / "system_info.txt").is_file()
    assert (tmp_path / "set_atoms.pml").read_text(encoding="utf-8").startswith(
        "load 2KAI.pdb, 2KAI\n"
    )


def test_canonical_interface_log_includes_provenance():
    assert INTERFACE_LOG_COLUMNS[-1] == "provenance"
    assert len(INTERFACE_LOG_COLUMNS) == 8


def test_compact_log_sorts_canonical_rows_and_rewrites_without_duplicates(tmp_path, monkeypatch):
    from vorpy.src.output import visualization

    def row(interface_id, representation, side, quantity):
        return {
            "interface_id": interface_id,
            "representation": representation,
            "side": side,
            "quantity": quantity,
            "value": 1,
            "units": "1",
            "status": "CERTIFIED",
            "provenance": "{}",
        }

    unordered = [
        row("B", "physical_partition_interface", "shared", "zeta"),
        row("A", "molecular_contact_surface", "B", "selected_ab_contact_count"),
        row("A", "molecular_contact_surface", "A", "selected_ab_contact_count"),
    ]
    monkeypatch.setattr(visualization, "_compact_result_rows", lambda iface: unordered)
    path = tmp_path / "logs.csv"

    visualization._compact_log(None, path)
    visualization._compact_log(None, path)

    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    assert reader.fieldnames == list(INTERFACE_LOG_COLUMNS)
    keys = [tuple(item[column] for column in ("interface_id", "representation", "side", "quantity"))
            for item in rows]
    assert list(keys) == sorted(keys)
    assert len(rows) == len(unordered)


def test_compact_rows_use_canonical_results_mappings(monkeypatch):
    from vorpy.src.output import visualization

    provenance = Provenance("mapping regression")
    identity = SimpleNamespace(interface_id="I", result_id="I-result", alpha_value=None)
    surface_identity = SimpleNamespace(interface_id="I", result_id="I-surface", alpha_value=None)
    interface_result = SimpleNamespace(identity=identity, provenance=provenance)
    surface_result = SimpleNamespace(identity=surface_identity, provenance=provenance)

    def quantity(name, value, status=Status.CERTIFIED):
        return SimpleNamespace(
            name=name,
            value=value,
            units="Å²" if name == "area" else "Å" if name == "complete_total" else "1",
            status=status,
            provenance=provenance,
        )

    def quantities(result):
        if result is interface_result:
            return iter((
                ("metrics.area", quantity("area", 2.5)),
                ("metrics.n_contacts", quantity("n_contacts", 2)),
                ("metrics.n_residue_contacts", quantity("n_residue_contacts", 1)),
                ("mean_curvature.complete_total", quantity("complete_total", 4.0)),
                ("gaussian_curvature.gauss_bonnet_completion",
                 quantity("gauss_bonnet_completion", 6.0)),
                ("topology.vertices", quantity("vertices", 7)),
                ("topology.boundary_loops", quantity("boundary_loops", 2)),
                ("selection.selected_features", quantity(
                    "selected_features", None, Status.UNRESOLVED,
                )),
            ))
        return iter((
            ("metrics.area", quantity("area", 3.5)),
            ("molecular_selection.selected_contacts", quantity("selected_contacts", 2)),
            ("molecular_selection.unresolved_contacts", quantity("unresolved_contacts", 1)),
        ))

    monkeypatch.setattr(visualization, "_cached_interface_result", lambda iface: interface_result)
    monkeypatch.setattr(visualization, "_cached_alpha_selection_result", lambda iface: None)
    monkeypatch.setattr(visualization, "_cached_surface_result", lambda iface, surface: surface_result)
    monkeypatch.setattr(visualization, "_iter_result_quantities", quantities)
    iface = SimpleNamespace(
        geometry_analysis=None,
        water_topology=None,
        representation_caches={
            "molecular_contact_surface:A": object(),
            "molecular_contact_surface:B": object(),
        },
    )

    rows = visualization._compact_result_rows(iface)
    logged = {(row["representation"], row["side"], row["quantity"]): row for row in rows}
    physical = "physical_partition_interface", "shared"
    assert logged[physical + (RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG["metrics.area"],)]["value"] == "2.5"
    assert logged[physical + (RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG["metrics.n_contacts"],)]["value"] == "2"
    assert logged[physical + (RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG["metrics.n_residue_contacts"],)]["value"] == "1"
    assert logged[physical + (RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG["mean_curvature.complete_total"],)]["value"] == "4.0"
    assert logged[physical + (RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG["gaussian_curvature.gauss_bonnet_completion"],)]["value"] == "6.0"
    assert logged[("physical_partition_interface", "shared", "vertex_count")]["value"] == "7"
    assert logged[("physical_partition_interface", "shared", "boundary_component_count")]["value"] == "2"
    assert logged[("physical_partition_interface", "shared", "alpha_selected_count")]["status"] == "UNRESOLVED"
    assert logged[("physical_partition_interface", "shared", "alpha_selected_count")]["value"] == ""
    for side in ("A", "B"):
        assert logged[("molecular_contact_surface", side, "area")]["value"] == "3.5"
        assert logged[("molecular_contact_surface", side, "selected_ab_contact_count")]["value"] == "2"
        assert logged[("molecular_contact_surface", side, "unresolved_contact_count")]["value"] == "1"
    assert not {row["quantity"] for row in rows} & {"complete_total", "selected_contacts", "unresolved_contacts"}
    assert all(row["provenance"] for row in rows)

    missing_surface_result = SimpleNamespace(identity=surface_identity, provenance=provenance)
    monkeypatch.setattr(visualization, "_cached_result_identity", lambda iface: identity)
    monkeypatch.setattr(visualization, "unavailable_interface", lambda *args, **kwargs: missing_surface_result)
    monkeypatch.setattr(visualization, "_iter_result_quantities",
                        lambda result: iter((("metrics.area", SimpleNamespace(
                            name="area", value=None, units="Å²", status=Status.NOT_CALCULATED,
                            provenance=provenance,
                        )),)) if result is missing_surface_result else iter(()))
    missing_rows = visualization._compact_result_rows(SimpleNamespace(
        geometry_analysis=None, water_topology=None, representation_caches={}
    ))
    missing_area = [row for row in missing_rows if row["quantity"] == "area"]
    assert missing_area and all(row["status"] == "NOT_CALCULATED" and row["value"] == ""
                                for row in missing_area)



def test_compact_interface_bundle_has_no_generic_wrappers(tmp_path, monkeypatch):
    from vorpy.src.output import visualization
    from vorpy.src.output.visualization import export_compact_visualization_bundle
    from vorpy.tests.interface.test_dual_visualization import cached_interface

    iface = cached_interface(tmp_path)
    iface.group1.ball_ndxs = [10]
    iface.group2.ball_ndxs = [13, 16]
    iface.dir = None
    system = SimpleNamespace(
        name="fixture", files={"dir": str(tmp_path)}, groups=[], ifaces=[iface],
        balls=None, update_progress=lambda **kwargs: None, verbose=False,
    )
    iface.sys = system
    monkeypatch.setattr(
        visualization,
        "write_pdb",
        lambda atoms, file_name, sys, directory=None: Path(directory, f"{file_name}.pdb").write_text(
            "fixture\n", encoding="utf-8"
        ),
    )

    export_compact_visualization_bundle(system)

    root = tmp_path / "interface_A_B"
    assert (root / "interface_surfs.off").is_file()
    assert (root / "interface_edges.off").is_file()
    assert (root / "interface_verts.off").is_file()
    assert (root / "aw_dual_edges.off").is_file()
    assert not (root / "interface.pml").exists()
    assert (root / "interface_info.txt").is_file()
    info = (root / "interface_info.txt").read_text(encoding="utf-8")
    assert "Interface network topology:" in info
    assert "Surface classification:" in info
    assert "Physical export selection: direct_inter_group" in info
    with (root / "logs.csv").open(newline="", encoding="utf-8") as stream:
        assert csv.DictReader(stream).fieldnames == list(INTERFACE_LOG_COLUMNS)
    assert not (tmp_path / "groups").exists()
    assert not (tmp_path / "interfaces").exists()


def test_compact_physical_mesh_preserves_selected_cached_facets(tmp_path):
    from vorpy.src.output import visualization
    from vorpy.tests.interface.test_dual_visualization import cached_interface

    iface = cached_interface(tmp_path)
    iface.group1.ball_ndxs = [10]
    iface.group2.ball_ndxs = [13, 16]
    result = visualization._compact_physical_interface(iface, tmp_path / "bundle")
    off = tmp_path / "bundle" / "interface_surfs.off"
    lines = [line for line in off.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    vertex_count, facet_count, _ = (int(value) for value in lines[1].split()[:3])
    direct = iface.net.surfs.iloc[[0, 1]]
    expected_facets = sum(len(tris) for tris in direct["tris"])
    faces = lines[2 + vertex_count:2 + vertex_count + facet_count]
    assert result["surface_count"] == len(direct)
    assert result["selection_mode"] == "direct_inter_group"
    assert result["area_A2"] == float(direct["sa"].sum())
    assert facet_count == expected_facets
    assert len(faces) == facet_count
    assert all(
        fields[0] == "3" and all(0 <= int(index) < vertex_count for index in fields[1:4])
        for fields in (face.split() for face in faces)
    )


def test_compact_physical_mesh_scopes_legacy_interface_to_direct_surfaces(tmp_path):
    from vorpy.src.output import visualization
    from vorpy.tests.interface.test_dual_visualization import cached_interface

    iface = cached_interface(tmp_path)
    iface.group1.ball_ndxs = [10]
    iface.group2.ball_ndxs = [13, 16]
    result = visualization._compact_physical_interface(iface, tmp_path / "bundle")

    off = tmp_path / "bundle" / "interface_surfs.off"
    lines = [line for line in off.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    vertex_count, facet_count, _ = (int(value) for value in lines[1].split()[:3])
    direct = iface.net.surfs.iloc[[0, 1]]
    assert result["surface_count"] == 2
    assert result["selection_mode"] == "direct_inter_group"
    assert result["area_A2"] == float(direct["sa"].sum())
    assert facet_count == sum(len(tris) for tris in direct["tris"])
    assert vertex_count == 8


def test_compact_interface_bundle_writes_both_input_group_pdbs(tmp_path, monkeypatch):
    from vorpy.src.output import visualization
    from vorpy.tests.interface.test_dual_visualization import cached_interface

    iface = cached_interface(tmp_path)
    iface.dir = None
    iface.group1 = SimpleNamespace(name="A", ball_ndxs=[10])
    iface.group2 = SimpleNamespace(name="B", ball_ndxs=[13, 16])
    system = SimpleNamespace(
        name="fixture", files={"dir": str(tmp_path)}, groups=[], ifaces=[iface],
        balls=None, update_progress=lambda **kwargs: None, verbose=False,
    )
    iface.sys = system

    def write_pdb(atoms, file_name, sys, directory=None):
        Path(directory, f"{file_name}.pdb").write_text(
            ",".join(str(atom) for atom in atoms), encoding="utf-8"
        )

    monkeypatch.setattr(visualization, "write_pdb", write_pdb)
    visualization.export_compact_visualization_bundle(system)

    root = tmp_path / "interface_A_B"
    assert (root / "group_1_atoms.pdb").read_text(encoding="utf-8") == "10"
    assert (root / "group_2_atoms.pdb").read_text(encoding="utf-8") == "13,16"
