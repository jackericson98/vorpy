import csv
from types import SimpleNamespace

import numpy as np
import pandas as pd

from vorpy.src.analyze.apollonius import PrimalFeatureRef
from vorpy.src.command.vpy_cmnd import visualization_cli_options
from vorpy.src.geometry.duals.model import DualComplex, DualSimplex
from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex
from vorpy.src.geometry.interfaces.model import (
    InterfaceRepresentation,
    InterfaceSurface,
)
from vorpy.src.geometry.visualization import (
    export_alpha_complex_visualization,
    export_dual_visualization,
)


def _case(scheme="aw"):
    ids = (10, 11, 12, 13)
    coordinates = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                   (0.0, 2.0, 0.0), (0.0, 0.0, 2.0))
    network = SimpleNamespace(
        settings={"net_type": scheme},
        balls=pd.DataFrame({
            "loc": [np.asarray(point) for point in coordinates],
            "rad": [1.0] * 4,
            "name": ["C"] * 4,
            "res": ["GLY"] * 4,
            "chn": ["A"] * 4,
        }, index=ids),
    )
    simplices = []
    for index, generator_id in enumerate(ids):
        simplices.append(DualSimplex(0, (generator_id,), (), True, True, True,
                                     "supported", (coordinates[index],), (1.0,)))
    for pair in ((10, 11), (10, 12), (11, 12), (10, 13), (11, 13), (12, 13)):
        points = tuple(coordinates[ids.index(generator_id)] for generator_id in pair)
        simplices.append(DualSimplex(1, pair, (), True, True, True, "supported", points, (1.0, 1.0)))
    simplices.append(DualSimplex(2, (10, 11, 12), (), True, True, True, "supported",
                                 (coordinates[0], coordinates[1], coordinates[2]), (1.0,) * 3))
    simplices.append(DualSimplex(3, ids, (), True, True, True, "supported", coordinates, (1.0,) * 4))
    feature_map = {("cell", generator_id): [f"0:{generator_id}"] for generator_id in ids}
    dual = DualComplex(network, scheme, simplices, feature_map)
    records = {dimension: {} for dimension in range(4)}
    for generator_id in ids:
        item = AlphaSimplex(0, (generator_id,), 0.0, 0.0, "synthetic", True)
        records[0][item.generator_tuple] = item
    for pair in ((10, 11), (10, 12)):
        item = AlphaSimplex(1, pair, 0.0, 0.0, "synthetic", True)
        records[1][item.generator_tuple] = item
    unresolved = AlphaSimplex(2, (10, 11, 12), None, None, "unresolved", False,
                              notes="synthetic unresolved AW triangle")
    records[2][unresolved.generator_tuple] = unresolved
    filtration = AlphaFiltration(
        network, scheme, records,
        {"native_alpha_units": "A" if scheme == "aw" else "A^2",
         "filtration_kind": "synthetic test filtration", "experimental": scheme == "aw"},
        blocked=[unresolved],
    )
    return network, dual, filtration


def test_dual_export_writes_center_coordinates_and_simplices(tmp_path):
    network, dual, _ = _case("pow")
    before = network.balls.copy(deep=True)
    export_dual_visualization(network, dual, tmp_path)
    target = tmp_path / "pow" / "dual"
    with (target / "dual_full_vertices.csv").open(newline="", encoding="utf-8") as stream:
        vertices = list(csv.DictReader(stream))
    assert len(vertices) == 4
    assert [float(vertices[0][key]) for key in ("x0_A", "y0_A", "z0_A")] == [0.0, 0.0, 0.0]
    with (target / "dual_full_triangles.csv").open(newline="", encoding="utf-8") as stream:
        triangles = list(csv.DictReader(stream))
    assert len(triangles) == 1
    assert tuple(float(triangles[0][key]) for key in ("x0_A", "y0_A", "z0_A")) == (0.0, 0.0, 0.0)
    with (target / "dual_full_tetrahedra.csv").open(newline="", encoding="utf-8") as stream:
        tetrahedra = list(csv.DictReader(stream))
    assert len(tetrahedra) == 1
    assert tuple(float(tetrahedra[0][key]) for key in ("x3_A", "y3_A", "z3_A")) == (0.0, 0.0, 2.0)
    assert (target / "dual_full_edges.off").exists()
    assert (target / "dual_full_triangles.off").exists()
    assert (target / "dual_full_tetrahedra_edges.off").exists()
    assert "run \"" in (tmp_path / "pow" / "visualize.pml").read_text(encoding="utf-8")
    assert (tmp_path / "pow" / "visualize.py").exists()
    pd.testing.assert_frame_equal(network.balls, before)


def test_alpha_export_uses_query_and_does_not_render_unresolved_aw_faces(tmp_path):
    network, dual, filtration = _case("aw")
    before = {dimension: dict(rows) for dimension, rows in filtration.records.items()}
    result = export_alpha_complex_visualization(network, filtration, 0.0, tmp_path, dual=dual)
    target = tmp_path / "aw" / "alpha_0"
    assert result["alpha_counts"] == {"0": 4, "1": 2, "2": 0, "3": 0}
    assert result["alpha_blocked_records_not_rendered"] == 1
    assert (target / "alpha_edges.off").exists()
    assert not (target / "alpha_triangles.off").exists()
    with (target / "alpha_edges.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert {row["simplex_id"] for row in rows} == {"1:10,11", "1:10,12"}
    assert filtration.records == before
    assert result["aw_warning"].startswith("AW dual simplices are visualization only")


def test_visualization_cli_options_are_distinct_from_existing_alpha_flags():
    remaining, requested, alpha = visualization_cli_options([
        "2KAI.pdb", "--export-dual-visualization", "--export-alpha-visualization", "0",
        "--alpha-complex", "1.4"
    ])
    assert remaining == ["2KAI.pdb", "--alpha-complex", "1.4"]
    assert requested is True
    assert alpha == 0.0


def test_alpha_interface_export_maps_to_solved_surface_without_mutation(tmp_path):
    from vorpy.src.geometry.visualization import export_dual_voronoi_mapping

    network, base_dual, _ = _case("pow")
    network.settings["surf_scheme"] = "none"
    surface_points = np.asarray(((0.5, 0.0, 0.0), (0.5, 1.0, 0.0), (0.5, 0.0, 1.0)))
    network.surfs = pd.DataFrame({
        "balls": [(10, 11)], "points": [surface_points], "tris": [[(0, 1, 2)]],
        "sa": [0.5],
    }, index=[501])
    network.edges = pd.DataFrame({"points": []})
    network.verts = pd.DataFrame()
    features = (PrimalFeatureRef("surf", 501, (10, 11), True, True),)
    pair = DualSimplex(1, (10, 11), features, True, True, True, "supported",
                       ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0)), (1.0, 1.0))
    all_simplices = [s for dim in range(4) for s in base_dual.simplices[dim].values()]
    all_simplices = [s for s in all_simplices if not (s.dimension == 1 and s.generator_ids == (10, 11))]
    dual = DualComplex(network, "pow", [*all_simplices, pair],
                       {("surf", 501): [pair.simplex_id]})
    surface = InterfaceSurface(501, (10, 11), pair.simplex_id, 0.5, True, True,
                               True, True, True, "selected")
    interface = InterfaceRepresentation(
        "pow", "alpha", 0.0, frozenset({10}), frozenset({11}),
        surfaces=[surface], selection_mappings=[{
            "generator_i": 10, "generator_j": 11, "surface_id": 501,
            "mapped": True, "selected": True,
        }],
    )
    before_surface = network.surfs.copy(deep=True)
    before_selected = list(interface.selected_surfaces)
    exported = export_dual_voronoi_mapping(network, dual, interface, tmp_path)
    assert exported["interface_pairs"] == 1
    assert exported["interface_area_A2"] == 0.5
    assert (tmp_path / "interface_alpha_surfaces.off").exists()
    assert (tmp_path / "alpha_interface_surface_mapping.csv").exists()
    pd.testing.assert_frame_equal(network.surfs, before_surface)
    assert interface.selected_surfaces == before_selected
