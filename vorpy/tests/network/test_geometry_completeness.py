import io
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd

from vorpy.src.boundary import SolventMode, prepare_network_geometry
from vorpy.src.network.network import Network
from vorpy.src.system.system import System


ROOT = Path(__file__).resolve().parents[2]
EDTA = ROOT / "data" / "EDTA.pdb"


def test_weighted_reference_distinguishes_hull_and_interior_cells():
    # The four tetrahedron corners have unbounded weighted cells; the center
    # has a bounded cell even before any primal topology is supplied.
    locations = [
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (0.0, 0.0, 0.0),
    ]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [0.0] * 5,
        "num": list(range(5)),
        "complete": [False] * 5,
    })
    empty_verts = pd.DataFrame(columns=["balls", "edges", "surfs", "loc"])
    empty_edges = pd.DataFrame(columns=["balls", "verts", "surfs"])
    empty_surfs = pd.DataFrame(columns=["balls", "verts", "edges"])
    network = Network(
        locs=locations,
        rads=[0.0] * 5,
        group=list(range(5)),
        settings={"net_type": "aw", "build_type": "all"},
        balls=balls,
        verts=empty_verts,
        edges=empty_edges,
        surfs=empty_surfs,
    )

    audit = network.diagnose_geometry_completeness()
    cells = {cell["local_id"]: cell for cell in audit["cells"]}

    assert cells[0]["bounded"] is False
    assert cells[3]["classification"] == "UNBOUNDED"
    assert cells[4]["bounded"] is True
    assert cells[4]["classification"] == "INCOMPLETE_CONSTRUCTION"
    assert audit["counts"]["classifications"] == {
        "INCOMPLETE_CONSTRUCTION": 1,
        "UNBOUNDED": 4,
    }

    # Generator incidence is keyed by the explicit ``num`` field, not by a
    # caller's incidental DataFrame index.
    network.balls.index = range(10, 15)
    reindexed = network.diagnose_geometry_completeness()
    assert reindexed["counts"] == audit["counts"]
    assert [cell["classification"] for cell in reindexed["cells"]] == [
        cell["classification"] for cell in audit["cells"]
    ]


def _closed_center_network():
    locations = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    ]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [0.0] * 5,
        "num": list(range(5)),
        "verts": [list(range(4)), [], [], [], []],
        "edges": [list(range(6)), [], [], [], []],
        "surfs": [list(range(4)), [], [], [], []],
        "complete": [False] * 5,
    })
    vertices = [
        ([0, 1, 2, 3], [0, 1, 3], [0, 1, 2]),
        ([0, 1, 2, 4], [0, 2, 4], [0, 1, 3]),
        ([0, 1, 3, 4], [1, 2, 5], [0, 2, 3]),
        ([0, 2, 3, 4], [3, 4, 5], [1, 2, 3]),
    ]
    edges = [
        ([0, 1, 2], [0, 1]),
        ([0, 1, 3], [0, 2]),
        ([0, 1, 4], [1, 2]),
        ([0, 2, 3], [0, 3]),
        ([0, 2, 4], [1, 3]),
        ([0, 3, 4], [2, 3]),
    ]
    surfaces = [
        ([0, 1], [0, 1, 2], [0, 1, 2]),
        ([0, 2], [0, 3, 4], [0, 1, 3]),
        ([0, 3], [1, 3, 5], [0, 2, 3]),
        ([0, 4], [2, 4, 5], [1, 2, 3]),
    ]
    verts = pd.DataFrame({
        "balls": [row[0] for row in vertices],
        "edges": [row[1] for row in vertices],
        "surfs": [row[2] for row in vertices],
        "loc": [(0.0, 0.0, 0.0)] * 4,
    })
    edge_table = pd.DataFrame({
        "balls": [row[0] for row in edges],
        "verts": [row[1] for row in edges],
        "surfs": [[]] * 6,
    })
    surface_table = pd.DataFrame({
        "balls": [row[0] for row in surfaces],
        "edges": [row[1] for row in surfaces],
        "verts": [row[2] for row in surfaces],
    })
    return Network(
        locs=locations,
        rads=[0.0] * 5,
        group=[0],
        settings={"net_type": "aw", "build_type": "all"},
        balls=balls,
        verts=verts,
        edges=edge_table,
        surfs=surface_table,
    )


def test_completeness_audit_separates_analytic_and_conservative_statuses():
    analytic = _closed_center_network()
    analytic._aw_unresolved_vertex_curvature = {
        (0, 0, "unsupported analytic vertex branch"),
    }
    analytic_cell = analytic.diagnose_geometry_completeness()["cells"][0]
    assert analytic_cell["classification"] == "UNSUPPORTED_ANALYTIC_GEOMETRY"
    assert analytic_cell["missing_faces"] == ()
    assert analytic_cell["missing_edges"] == ()
    assert analytic_cell["closure_relations"]["open_edges"] == ()

    conservative = _closed_center_network()
    conservative_cell = conservative.diagnose_geometry_completeness()["cells"][0]
    assert conservative_cell["classification"] == "OVERLY_CONSERVATIVE"


def test_edta_molecule_only_audit_separates_unbounded_and_incomplete_cells():
    system = System(file=str(EDTA), make_dir=False)
    geometry = prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)
    assert len(geometry.system_indices) == 32
    assert geometry.system_indices == tuple(range(32))

    from vorpy.workbench.services.vorpy_backend import (
        VorPyBackend,
        VorPySolveSettings,
    )

    with redirect_stdout(io.StringIO()):
        result = VorPyBackend(
            VorPySolveSettings(
                boundary_mode="none",
                max_vertices=5,
                surface_resolution=0.2,
            )
        ).solve(EDTA, lambda *args: None, lambda: False)
    network = result.export_group.net
    before = (
        network.balls.copy(deep=True),
        network.verts.copy(deep=True),
        network.edges.copy(deep=True),
        network.surfs.copy(deep=True),
    )

    audit = network.diagnose_geometry_completeness()
    cells = {cell["system_id"]: cell for cell in audit["cells"]}

    assert network.solvent_mode == SolventMode.MOLECULE_ONLY.value
    assert audit["counts"] == {
        "cells": 32,
        "complete": 2,
        "classifications": {
            "COMPLETE": 2,
            "INCOMPLETE_CONSTRUCTION": 10,
            "UNBOUNDED": 20,
        },
    }
    assert {cell_id for cell_id, cell in cells.items() if cell["bounded"]} == {
        0, 1, 3, 4, 6, 8, 10, 12, 14, 16, 18, 19,
    }
    assert cells[0]["classification"] == "COMPLETE"
    assert cells[0]["missing_faces"] == ()
    assert cells[0]["missing_edges"] == ()
    assert cells[1]["classification"] == "INCOMPLETE_CONSTRUCTION"
    assert cells[1]["bounded"] is True
    assert cells[1]["missing_faces"]
    assert cells[1]["missing_edges"]
    assert cells[2]["classification"] == "UNBOUNDED"
    assert cells[2]["bounded"] is False
    assert cells[2]["closure_relations"]["open_edges"]

    pd.testing.assert_frame_equal(network.balls, before[0])
    pd.testing.assert_frame_equal(network.verts, before[1])
    pd.testing.assert_frame_equal(network.edges, before[2])
    pd.testing.assert_frame_equal(network.surfs, before[3])


def test_edta_mg_is_explicitly_in_the_current_solvent_classification():
    system = System(file=str(EDTA), make_dir=False)
    assert str(system.balls.iloc[32]["res_name"]).upper() == "MG"

    molecule_only = prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)
    context = prepare_network_geometry(system, SolventMode.SOLVENT_CONTEXT)

    assert molecule_only.system_indices == tuple(range(32))
    assert 32 not in molecule_only.system_indices
    assert 32 in context.system_indices
    assert context.solvent_local_indices[0] == 32
    assert len(context.solvent_local_indices) == 580
