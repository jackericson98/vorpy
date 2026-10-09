from types import SimpleNamespace

import pandas as pd

from vorpy.src.network.aw_geometry_envelope import extract_aw_geometry_envelope
from vorpy.src.results.model import Status


def _network(*, invalid_edge=False):
    """A small table-backed AW system with an exact equal-radius edge."""
    balls = pd.DataFrame({
        "loc": [
            (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
            (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
        ],
        "rad": [1.0, 1.0, 1.0, 1.0],
        "complete": [True, True, True, True],
        "stable_id": ["a", "b", "c", "d"],
    })
    verts = pd.DataFrame({
        "balls": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "loc": [(0.5, 0.5, 0.5), (0.6, 0.5, 0.5)],
        "rad": [-0.1339745962, -0.0630145813],
    })
    edges = pd.DataFrame({
        "balls": [[0, 1, 2]],
        "verts": [[0] if invalid_edge else [0, 1]],
        "surfs": [[0]],
        "points": [[(0.5, 0.5, 0.5), (0.6, 0.5, 0.5)]],
        "length": [0.1],
    })
    surfs = pd.DataFrame({
        "balls": [[0, 1]],
        "verts": [[0, 1]],
        "edges": [[0]],
        # Omit ``func`` intentionally: extraction must perform the existing
        # exact generator reconstruction rather than fit this mesh.
        "points": [[(1.0, 0.0, 0.0)]],
        "tris": [[(0, 0, 0)]],
    })
    return SimpleNamespace(
        settings={"net_type": "aw", "surf_res": 0.25},
        balls=balls,
        verts=verts,
        edges=edges,
        surfs=surfs,
        group=[0, 1, 2, 3],
        iface_grps=None,
        frame="synthetic",
    )


def test_envelope_connects_atoms_cells_surfaces_edges_and_vertices_deterministically():
    first = extract_aw_geometry_envelope(_network())
    second = extract_aw_geometry_envelope(_network())

    assert first.coordinate_units == "Å"
    assert [atom.cell_id for atom in first.atoms] == [cell.cell_id for cell in first.cells]
    assert [record.surface_id for record in first.surfaces] == [record.surface_id for record in second.surfaces]
    assert [record.edge_id for record in first.edges] == [record.edge_id for record in second.edges]
    assert [record.vertex_id for record in first.vertices] == [record.vertex_id for record in second.vertices]

    surface = first.surfaces[0]
    edge = first.edges[0]
    assert surface.status is Status.CERTIFIED
    assert surface.coefficient_source == "exact_generator_reconstruction"
    assert surface.geometry_kind == "quadratic"
    assert len(surface.coefficients) == 14
    assert surface.spherical_patch_status is Status.NOT_SUPPORTED
    assert len(surface.owners) == 2
    assert {owner.cell_id for owner in surface.owners} == set(surface.cell_ids)
    assert all(owner.normal_convention.startswith("aw_clearance_gradient") for owner in surface.owners)
    assert edge.status is Status.CERTIFIED
    assert edge.analytic["kind"] == "line"
    assert edge.component_id == surface.component_id
    assert first.components[0].component_id == surface.component_id

    # Table-backed synthetic vertices have no solver provenance, so the
    # envelope exposes the usable analytic subset without claiming a complete
    # geometry solve or a downstream curvature certification.
    assert first.completeness.status is Status.PARTIAL
    assert "solver_completeness_not_proven" in first.completeness.reasons


def test_unresolved_edge_is_explicit_and_never_replaced_by_a_zero_or_polyline():
    envelope = extract_aw_geometry_envelope(_network(invalid_edge=True))
    edge = envelope.edges[0]

    assert edge.status is Status.UNRESOLVED
    assert edge.analytic is None
    assert "unsupported_topology" in edge.reason
    assert envelope.completeness.counts["unresolved_primitives"] >= 1
