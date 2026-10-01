from types import SimpleNamespace

import numpy as np
import pandas as pd

from vorpy.src.analyze.aw_pipeline_audit import (
    audit_interface_surfaces, compare_surface_tables, stable_atom_identity,
)


def _tiny_filtered_network():
    balls = pd.DataFrame({
        "loc": [(0., 0., 0.), (2., 0., 0.), (0., 2., 0.)],
        "rad": [1., 1., 1.], "num": [10, 11, 12], "complete": [False, False, False],
        "chain_name": ["A", "B", "A"], "res_name": ["GLY"] * 3,
        "res_seq": [1, 2, 3], "pdb_ins_code": [""] * 3,
        "name": ["CA", "CA", "CB"], "element": ["C"] * 3,
    })
    surfs = pd.DataFrame({
        "balls": [[0, 1]], "edges": [[0, 1, 2]], "verts": [[0, 1, 2]],
        "points": [np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])],
        "tris": [np.array([[0, 1, 2]])], "sa": [0.5], "func": [object()],
    })
    edges = pd.DataFrame({
        "verts": [[0, 1], [1, 2], [2, 0]],
        "surfs": [[0], [0], [0]],
    })
    verts = pd.DataFrame({"surfs": [[0], [0], [0]]})
    return SimpleNamespace(balls=balls, surfs=surfs, edges=edges, verts=verts)


def test_feature_boundary_can_be_complete_with_incomplete_generator_cells(monkeypatch):
    from vorpy.src.network import edge_geometry_diagnostics

    monkeypatch.setattr(
        edge_geometry_diagnostics, "resolve_aw_network_edge",
        lambda network, edge_id: SimpleNamespace(
            geometry=SimpleNamespace(t_min=0., t_max=1.), length=1.
        ),
    )
    net = _tiny_filtered_network()
    rows = audit_interface_surfaces(net, {0}, {1})
    assert len(rows) == 1
    assert rows[0]["both_cells_complete"] is False
    assert rows[0]["boundary_closed_cycles"] is True
    assert rows[0]["all_boundary_edges_supported"] is True
    assert rows[0]["feature_complete"] is True


def test_open_surface_boundary_is_not_feature_complete(monkeypatch):
    from vorpy.src.network import edge_geometry_diagnostics

    monkeypatch.setattr(
        edge_geometry_diagnostics, "resolve_aw_network_edge",
        lambda network, edge_id: SimpleNamespace(
            geometry=SimpleNamespace(t_min=0., t_max=1.), length=1.
        ),
    )
    net = _tiny_filtered_network()
    net.surfs.at[0, "edges"] = [0, 1]
    rows = audit_interface_surfaces(net, {0}, {1})
    assert rows[0]["boundary_closed_cycles"] is False
    assert rows[0]["feature_complete"] is False
    assert "boundary_not_disjoint_closed_cycles" in rows[0]["feature_errors"]


def test_missing_analytic_boundary_edge_rejects_feature_completeness(monkeypatch):
    from vorpy.src.network import edge_geometry_diagnostics

    monkeypatch.setattr(
        edge_geometry_diagnostics, "resolve_aw_network_edge",
        lambda network, edge_id: (_ for _ in ()).throw(ValueError("unresolved edge")),
    )
    rows = audit_interface_surfaces(_tiny_filtered_network(), {0}, {1})
    assert rows[0]["feature_complete"] is False
    assert "edge_geometry:0:ValueError" in rows[0]["feature_errors"]


def test_stable_generator_identity_uses_structure_identity_not_row_id():
    net = _tiny_filtered_network()
    first = stable_atom_identity(net, 0)
    assert first.startswith("A|GLY|1||CA|C|")
    # The canonical identity is based on chain/residue/atom plus position.
    reordered = net.balls.rename(index={0: 20, 1: 21, 2: 22})
    net.balls = reordered
    assert stable_atom_identity(net, 20) == first


def test_surface_matching_uses_stable_atoms_across_network_index_changes():
    full_balls = pd.DataFrame({
        "loc": [(0., 0., 0.), (2., 0., 0.)], "rad": [1., 1.],
        "num": [100, 200], "complete": [True, True],
        "chain_name": ["A", "B"], "res_name": ["GLY", "SER"],
        "res_seq": [1, 2], "pdb_ins_code": ["", ""],
        "name": ["CA", "CA"], "element": ["C", "C"],
    }, index=[20, 21])
    full = SimpleNamespace(
        balls=full_balls,
        sys=SimpleNamespace(balls=full_balls),
        surfs=pd.DataFrame({"balls": [[20, 21]], "sa": [3.5]}, index=[8]),
    )
    stable = tuple(sorted(stable_atom_identity(full, i) for i in (20, 21)))
    filtered_rows = [{
        "stable_generator_ids": " || ".join(stable), "surface_id": 0,
        "surface_area": 3.5, "both_cells_complete": False,
        "feature_complete": True,
    }]
    rows = compare_surface_tables(filtered_rows, full, {20}, {21})
    assert len(rows) == 1
    assert rows[0]["full_geometry_present"] is True
    assert rows[0]["full_surface_id"] == 8
