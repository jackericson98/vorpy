from itertools import combinations
from types import SimpleNamespace
from pathlib import Path
import math
import os
import time

import pandas as pd
import pytest

from vorpy.src.analyze.apollonius import ApolloniusComplex


@pytest.fixture(scope="module")
def edta_network():
    from vorpy.workbench.services.vorpy_backend import VorPyBackend

    source = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"
    old_cwd = Path.cwd()
    try:
        result = VorPyBackend().solve(source, lambda *args: None, lambda: False)
        yield result.export_group.net
    finally:
        os.chdir(old_cwd)


def _network(n_generators, surfaces=(), edges=(), vertices=(), complete=True):
    locs = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
    ][:n_generators]
    balls = pd.DataFrame({
        "loc": locs,
        "rad": [1.0 + 0.1 * i for i in range(n_generators)],
        "num": list(range(10, 10 + n_generators)),
        "complete": [complete] * n_generators,
    })
    verts = pd.DataFrame(vertices, columns=["balls", "loc"])
    edges = pd.DataFrame(edges, columns=["balls", "verts"])
    surfaces = pd.DataFrame(surfaces, columns=["balls", "edges", "verts"])
    return SimpleNamespace(
        balls=balls, surfs=surfaces, edges=edges, verts=verts, group_name="synthetic"
    )


def _closed_triangle_network():
    vertices = [([0, 1, 2], [0.0, 0.0, 0.0]), ([0, 1, 2], [0.0, 0.0, 1.0])]
    edges = [([0, 1, 2], [0, 1])]
    surfaces = [
        ([0, 1], [0], [0, 1]),
        ([0, 2], [0], [0, 1]),
        ([1, 2], [0], [0, 1]),
    ]
    return _network(3, surfaces, edges, vertices)


def _closed_tetra_network(duplicate_edge=False):
    # Four weighted spheres can have two distinct common tangent-sphere
    # centers. Keeping both rows also exercises multiplicity on the dual tet.
    vertices = [
        ([0, 1, 2, 3], [0.1, 0.1, 0.1]),
        ([0, 1, 2, 3], [0.2, 0.2, 0.2]),
    ]
    edge_rows = [(list(triple), [0, 1]) for triple in combinations(range(4), 3)]
    if duplicate_edge:
        edge_rows.append(([0, 1, 2], [0]))
    surfaces = [
        (list(pair), list(range(4)), [0]) for pair in combinations(range(4), 2)
    ]
    return _network(4, surfaces, edge_rows, vertices)


def test_two_generators_define_one_straight_dual_edge():
    network = _network(2, surfaces=[([0, 1], [0], [0, 1])])
    complex_ = ApolloniusComplex(network)
    edge = complex_.simplices[1][(0, 1)]
    assert edge.geometry == ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert edge.generator_radii == (1.0, 1.1)
    assert edge.primal_feature_ids == ("surf:0",)
    assert complex_.validate() == []


def test_three_generator_edge_defines_triangle_and_closure():
    complex_ = ApolloniusComplex(_closed_triangle_network())
    assert list(complex_.simplices[2]) == [(0, 1, 2)]
    assert len(complex_.simplices[1]) == 3
    assert complex_.validate() == []


def test_four_generator_vertex_defines_tetrahedron_and_closure():
    complex_ = ApolloniusComplex(_closed_tetra_network())
    assert list(complex_.simplices[3]) == [(0, 1, 2, 3)]
    assert len(complex_.simplices[2]) == 4
    assert len(complex_.simplices[1]) == 6
    assert complex_.assert_valid()
    assert "Euler characteristic:" in complex_.summary()


def test_repeated_primal_components_are_preserved_on_canonical_simplex():
    complex_ = ApolloniusComplex(_closed_tetra_network(duplicate_edge=True))
    triangle = complex_.simplices[2][(0, 1, 2)]
    assert triangle.num_primal_components == 2
    assert triangle.primal_feature_ids == ("edge:0", "edge:4")
    assert len(complex_.simplices[2]) == 4
    assert {simplex.dimension for simplex in complex_.multiple_primal_simplices} == {2, 3}


def test_disconnected_surface_mesh_components_remain_distinct_primal_references():
    network = _network(2, surfaces=[(
        [0, 1], [0, 1], [0, 1, 2, 3, 4, 5],
    )])
    network.surfs["points"] = [[
        (0., 0., 0.), (1., 0., 0.), (0., 1., 0.),
        (3., 0., 0.), (4., 0., 0.), (3., 1., 0.),
    ]]
    network.surfs["tris"] = [[[0, 1, 2], [3, 4, 5]]]
    complex_ = ApolloniusComplex(network)
    edge = complex_.simplices[1][(0, 1)]
    assert edge.num_primal_components == 2
    assert edge.primal_feature_ids == (
        "surf:0#component-0", "surf:0#component-1"
    )


def test_higher_order_vertex_record_is_reported_not_split():
    network = _network(5, vertices=[([0, 1, 2, 3, 4], [0.0, 0.0, 0.0])])
    complex_ = ApolloniusComplex(network)
    assert not complex_.simplices[3]
    assert complex_.unsupported_features[0].generator_ids == (0, 1, 2, 3, 4)
    assert "Degenerate/unsupported:     1" in complex_.summary()


def test_incomplete_boundary_is_exposed_without_inventing_closure():
    network = _network(3, surfaces=[([0, 1], [0], [0, 1])], complete=False)
    complex_ = ApolloniusComplex(network)
    assert not complex_.simplices[1][(0, 1)].complete
    assert len(complex_.incomplete_features) == 4  # 3 cells + the surface
    assert complex_.validate() == []


def test_aw_edge_birth_uses_exact_clearance_parameter_not_vertex_value():
    from vorpy.tests.network.test_build_net import _analytic_diagnostic_network

    network = _analytic_diagnostic_network()
    complex_ = ApolloniusComplex(network)
    complex_.calculate_alpha_births()
    triangle = complex_.simplices[2][(0, 1, 2)]
    assert triangle.alpha_status == "exact_aw_edge_minimum"
    # The exact branch parameter is the common AW clearance; its minimum is
    # the lower endpoint rho=2, while the stored endpoint vertex has rho=4.
    assert triangle.alpha_birth == 2.0
    assert min(triangle.diagnostics["alpha"]["components"][0]["candidate_clearances"]) == 2.0


def test_equal_radius_aw_edge_uses_straight_segment_clearance():
    network = SimpleNamespace(
        settings={"net_type": "aw"},
        balls=pd.DataFrame({
            "loc": [(0., 0., 0.), (2., 0., 0.), (0., 2., 0.)],
            "rad": [1., 1., 1.],
        }),
        verts=pd.DataFrame({
            "loc": [(1., 1., -1.), (1., 1., 1.)],
        }),
        edges=pd.DataFrame({
            "balls": [[0, 1, 2]], "verts": [[0, 1]],
            "points": [[(1., 1., -1.), (1., 1., 0.), (1., 1., 1.)]],
        }),
        surfs=pd.DataFrame(), group_name="equal-radius",
    )
    complex_ = ApolloniusComplex(network)
    complex_.calculate_alpha_births()
    triangle = complex_.simplices[2][(0, 1, 2)]
    assert triangle.alpha_status == "exact_aw_edge_minimum"
    assert abs(triangle.alpha_birth - (2. ** 0.5 - 1.)) < 1e-12


def test_aw_tetrahedron_birth_checks_all_four_generator_clearances():
    network = _closed_tetra_network()
    network.balls["loc"] = [
        (1., 1., 1.), (1., -1., -1.), (-1., 1., -1.), (-1., -1., 1.)
    ]
    network.balls["rad"] = [1., 1., 1., 1.]
    network.verts["loc"] = [(0., 0., 0.), (0., 0., 0.)]
    network.verts["rad"] = [3. ** 0.5 - 1., 3. ** 0.5 - 1.]
    complex_ = ApolloniusComplex(network)
    complex_.calculate_alpha_births()
    tetrahedron = complex_.simplices[3][(0, 1, 2, 3)]
    assert tetrahedron.alpha_status == "exact_aw_vertex_clearance"
    assert abs(tetrahedron.alpha_birth - (3. ** 0.5 - 1.)) < 1e-12
    for component in tetrahedron.diagnostics["alpha"]["components"]:
        assert component["spread"] < 1e-12
        assert max(component["clearances"]) - min(component["clearances"]) < 1e-12


def test_cell_birth_requires_generator_center_to_belong_to_its_aw_cell():
    network = _network(2)
    network.balls.loc[0, "rad"] = 1.
    network.balls.loc[1, "rad"] = 3.
    network.settings = {"net_type": "aw"}
    complex_ = ApolloniusComplex(network)
    complex_.calculate_alpha_births()
    hidden_center = complex_.simplices[0][(0,)]
    assert hidden_center.alpha_birth is None
    assert hidden_center.alpha_status == "cell_minimum_unresolved"
    assert complex_.simplices[0][(1,)].alpha_birth == -3.


def test_alpha_complex_is_filtered_and_closed_over_explicit_birth_values():
    complex_ = ApolloniusComplex(_closed_tetra_network())
    for dimension, table in complex_.simplices.items():
        for simplex in table.values():
            simplex.alpha_birth = dimension * 0.5
            simplex.alpha_status = "test_exact"
    complex_.alpha_tolerance = 1e-6
    alpha = complex_.alpha_complex(1.0)
    assert alpha.counts == {0: 4, 1: 6, 2: 4, 3: 0}
    assert alpha.validate() == []
    assert "Alpha (Å): 1" in alpha.summary()


def test_filtration_reports_substantive_monotonicity_violations():
    complex_ = ApolloniusComplex(_closed_triangle_network())
    for simplex in complex_.simplices[0].values():
        simplex.alpha_birth = 0.0
    complex_.simplices[1][(0, 1)].alpha_birth = 2.0
    complex_.simplices[1][(0, 2)].alpha_birth = 0.0
    complex_.simplices[1][(1, 2)].alpha_birth = 0.0
    complex_.simplices[2][(0, 1, 2)].alpha_birth = 1.0
    issues = complex_.validate_filtration()
    assert any(issue["kind"] == "monotonicity_violation" for issue in issues)


def test_edta_real_network_keeps_only_incidence_closed_supported_simplices(edta_network):
    complex_ = ApolloniusComplex(edta_network, "EDTA")
    diagnostics = complex_.incidence_diagnostics()

    # Raw primal tuples remain discoverable, but only records with actual
    # lower-dimensional boundary incidence enter the supported subcomplex.
    assert tuple(len(complex_.simplices[d]) for d in range(4)) == (612, 424, 890, 499)
    assert complex_.supported_counts == {0: 612, 1: 424, 2: 422, 3: 138}
    assert complex_.euler_characteristic == 472
    assert complex_.raw_euler_characteristic == 579
    assert len(diagnostics["missing_face_relations"]) == 1474
    assert len(diagnostics["missing_face_causes"]["incomplete_network_geometry"]) == 1474
    assert not diagnostics["missing_face_causes"]["missing_primal_feature_despite_valid_incidence"]
    assert not diagnostics["missing_face_causes"]["unsupported_higher_order_incidence"]
    assert not diagnostics["causes"]["canonical_id_mismatch"]
    assert len(complex_.multiple_primal_simplices) == 4
    assert complex_.validate() == []

    # This concrete triple-edge has two linked pair surfaces. Its third pair
    # is absent from the network and the edge is unsupported due incomplete
    # generator-cell geometry; the dual extractor must not invent it.
    triangle = complex_.simplices[2][(1, 199, 381)]
    assert triangle.primal_feature_ids == ("edge:73",)
    assert edta_network.edges.loc[73, "surfs"] == [25, 26]
    assert (199, 381) not in complex_.simplices[1]
    assert not complex_._supported(triangle)


def test_edta_interface_local_faces_are_distinct_from_global_closure_gaps(edta_network):
    complex_ = ApolloniusComplex(edta_network, "EDTA")
    diagnostics = complex_.incidence_diagnostics()
    # F10's default molecule-only producer no longer guarantees the old
    # archive row numbers (and may exclude an explicitly listed ion).  Pick a
    # real generator interface from the current network instead of coupling
    # this local-completeness diagnostic to legacy row IDs.
    interface_generator = next(
        int(index)
        for index, row in edta_network.surfs.iterrows()
        if len(row["balls"]) == 2 and 0 in row["balls"]
    )
    interface_surfaces = {
        int(index)
        for index, row in edta_network.surfs.iterrows()
        if len(row["balls"]) == 2 and interface_generator in row["balls"]
    }
    interface_pairs = {
        tuple(sorted(map(int, edta_network.surfs.loc[sid, "balls"])))
        for sid in interface_surfaces
    }
    candidate_edges = []
    for edge_id, row in edta_network.edges.iterrows():
        triangle = tuple(sorted(map(int, row["balls"])))
        incident = tuple(sorted(
            int(sid) for sid in row["surfs"] if int(sid) in interface_surfaces
        ))
        if interface_generator in triangle and incident:
            candidate_edges.append((triangle, incident))

    missing = diagnostics["missing_face_relations"]
    missing_pairs = {
        tuple(relation["missing_face"])
        for relation in missing
        if relation["face_dimension"] == 1
    }
    assert not interface_pairs & missing_pairs

    seam_triangles = {
        triangle for triangle, incident in candidate_edges if len(incident) == 2
    }
    boundary_triangles = {
        triangle for triangle, incident in candidate_edges if len(incident) == 1
    }
    missing_boundary_relations = [
        relation for relation in missing
        if relation["simplex"] in {
            "2:" + ",".join(map(str, triangle))
            for triangle in boundary_triangles
        }
    ]
    assert seam_triangles
    # The selected pair faces are present, while unrelated cofaces may still
    # be absent from the globally incomplete generator complex.
    assert not boundary_triangles
    assert not missing_boundary_relations
    assert diagnostics["missing_face_relations"]
    assert missing_pairs - interface_pairs


def test_edta_molecule_only_interface_is_local_complete_but_global_unresolved():
    from vorpy.src.analyze.aw_interface_curvature import _integrate_edge_beta
    from vorpy.src.analyze.open_interface_curvature import (
        _aw_supported_gauss_bonnet,
        _network_face_cycle,
        analyze_surface_topology,
    )
    from vorpy.src.network.aw_geometry_reuse import (
        build_aw_interface_view,
        certify_aw_interface_completeness,
        extract_aw_geometry,
    )
    from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge
    from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings

    source = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"
    result = VorPyBackend(
        VorPySolveSettings(boundary_mode="none", max_vertices=5, surface_resolution=0.2)
    ).solve(source, lambda *args: None, lambda: False)
    network = result.export_group.net
    complex_ = ApolloniusComplex(network, "EDTA-molecule-only")
    diagnostics = complex_.incidence_diagnostics()

    assert network.solvent_mode == "MOLECULE_ONLY"
    assert len(network.balls) < 42
    assert not network.boundary_indices
    assert int(network.balls["complete"].sum()) < len(network.balls)
    interface_generator = 0
    interface_surface_ids = {
        int(index) for index, row in network.surfs.iterrows()
        if len(row["balls"]) == 2 and interface_generator in row["balls"]
    }
    assert interface_surface_ids
    pair_complete = [
        bool(network.balls.loc[list(map(int, network.surfs.loc[sid, "balls"])), "complete"].all())
        for sid in interface_surface_ids
    ]
    assert not all(pair_complete)

    interface_pairs = {
        tuple(sorted(map(int, network.surfs.loc[sid, "balls"])))
        for sid in interface_surface_ids
    }
    candidate_edges = []
    for edge_id, row in network.edges.iterrows():
        triangle = tuple(sorted(map(int, row["balls"])))
        incident = tuple(sorted(
            int(sid) for sid in row["surfs"] if int(sid) in interface_surface_ids
        ))
        if interface_generator in triangle and incident:
            candidate_edges.append((int(edge_id), triangle, incident))
    assert candidate_edges
    assert {len(incident) for _, _, incident in candidate_edges} == {2}

    topology = analyze_surface_topology([
        _network_face_cycle(network, network.surfs.loc[sid])
        for sid in sorted(interface_surface_ids)
    ])
    assert topology["faces"] == len(interface_surface_ids)
    assert topology["euler_characteristic"] == 2
    assert topology["boundary_edges"] == topology["boundary_components"] == 0
    assert topology["boundary_is_manifold"] and topology["orientable"]

    # The selected patch is topologically closed, but the conservative local
    # certificate follows its shared dual support component.  On this real
    # molecule that component reaches every retained generator, so incomplete
    # unrelated cells must still block certification.
    store = extract_aw_geometry(network, allow_incomplete=True)
    view = build_aw_interface_view(
        store, {interface_generator}, set(range(1, len(network.balls)))
    )
    certificate = certify_aw_interface_completeness(view)
    assert len(view.surface_ids) == len(interface_surface_ids)
    assert not certificate.proven
    assert not certificate.globally_complete
    assert "incomplete_relevant_competitor" in certificate.reasons
    assert certificate.counts["relevant_generators"] == len(network.balls)

    missing = diagnostics["missing_face_relations"]
    assert not interface_pairs & {
        tuple(relation["missing_face"])
        for relation in missing
        if relation["face_dimension"] == 1
    }
    simplex_ids = {
        "2:" + ",".join(map(str, triangle))
        for _, triangle, _ in candidate_edges
    }
    seam_missing = [relation for relation in missing if relation["simplex"] in simplex_ids]
    assert seam_missing
    assert all(
        interface_generator not in tuple(relation["missing_face"])
        for relation in seam_missing
    )
    assert 0 < sum(
        complex_.is_supported(2, triangle)
        for _, triangle, _ in candidate_edges
    ) < len(candidate_edges)

    resolved_edges = []
    signed = []
    for edge_id, _, incident in candidate_edges:
        resolved = resolve_aw_network_edge(network, edge_id)
        assert resolved.status in {"matched_nonsingular", "matched_curve"}
        resolved_edges.append(edge_id)
        forward = _integrate_edge_beta(
            network, resolved, incident, set(range(1, len(network.balls))), {interface_generator},
            order=8, max_order=64, rtol=1e-9, atol=1e-11,
        )
        reverse = _integrate_edge_beta(
            network, resolved, incident, {interface_generator}, set(range(1, len(network.balls))),
            order=8, max_order=64, rtol=1e-9, atol=1e-11,
        )
        assert forward["status"] == reverse["status"] == "converged"
        assert forward["length"] > 0.0
        assert forward["signed"] == pytest.approx(-reverse["signed"])
        signed.append(forward["signed"])
    assert resolved_edges
    assert len(resolved_edges) == len(candidate_edges)
    assert sum(signed) != 0.0

    face_records = [
        (
            sid,
            tuple(sorted(map(int, network.surfs.loc[sid, "balls"]))),
            _network_face_cycle(network, network.surfs.loc[sid]),
        )
        for sid in sorted(interface_surface_ids)
    ]
    gaussian_16 = _aw_supported_gauss_bonnet(
        network, face_records, interface_surface_ids, set(range(1, len(network.balls))), order=16
    )
    gaussian_32 = _aw_supported_gauss_bonnet(
        network, face_records, interface_surface_ids, set(range(1, len(network.balls))), order=32
    )
    gaussian_baseline = _aw_supported_gauss_bonnet(
        network, face_records, interface_surface_ids, set(range(1, len(network.balls))),
        order=16, smooth_face_method="triangulated",
    )
    gaussian_reversed = _aw_supported_gauss_bonnet(
        network, face_records, interface_surface_ids, set(range(1, len(network.balls))), order=16
    )
    assert gaussian_16["unresolved_edges"] == 0
    assert gaussian_16["unresolved_vertices"] == 0
    assert gaussian_16["unresolved_faces"] == 0
    assert not gaussian_16["certified"]
    assert not gaussian_32["certified"]
    assert not gaussian_16["global_complete"]
    assert gaussian_16["smooth_face_method"] == "curved"
    assert abs(gaussian_16["residual"]) < abs(gaussian_baseline["residual"])
    for key in (
        "smooth_K", "internal_seam_face_geodesic_K", "free_boundary_face_geodesic_K",
        "interior_vertex_defect_K", "boundary_vertex_corner_K", "total_GB_left", "residual",
    ):
        assert gaussian_reversed[key] == pytest.approx(gaussian_16[key], abs=1e-12)
    assert gaussian_16["free_boundary_face_geodesic_K"] == pytest.approx(0.0)
    assert gaussian_16["boundary_vertex_corner_K"] == pytest.approx(0.0)
    assert len(gaussian_16["face_diagnostics"]) == len(interface_surface_ids)
    assert all(row["complete"] for row in gaussian_16["face_diagnostics"])
    assert sum(row["face_GB_residual"] for row in gaussian_16["face_diagnostics"]) == pytest.approx(
        gaussian_16["residual"], abs=1e-12
    )
    assert len(gaussian_16["edge_diagnostics"]) == len(candidate_edges)
    assert all(
        row["complete"] and len(row["surface_ids"]) == 2
        and len(row["incidence_terms"]) == 2 and row["sum"] is not None
        for row in gaussian_16["edge_diagnostics"]
    )
    assert len(gaussian_16["vertex_terms"]) == topology["vertices"]
    assert gaussian_16["total_GB_left"] == pytest.approx(
        gaussian_16["smooth_K"]
        + gaussian_16["internal_seam_face_geodesic_K"]
        + gaussian_16["free_boundary_face_geodesic_K"]
        + gaussian_16["interior_vertex_defect_K"]
        + gaussian_16["boundary_vertex_corner_K"],
        abs=1e-12,
    )
    assert abs(gaussian_16["residual"]) > 1e-6
    # Validated conic rematching removes the clearance-branch conditioning
    # issue; the persistent nonzero residual is an accounting blocker.
    assert gaussian_32["residual"] == pytest.approx(
        gaussian_16["residual"], abs=1e-12
    )


def test_edta_curved_gaussian_surface_resolution_and_quadrature_converge():
    from vorpy.src.analyze.open_interface_curvature import (
        _aw_curved_surface_gaussian,
        _aw_supported_gauss_bonnet,
        _network_face_cycle,
    )
    from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings

    source = Path(__file__).resolve().parents[2] / "data" / "EDTA.pdb"
    study = {}
    for resolution in (0.2, 0.1):
        network = VorPyBackend(
            VorPySolveSettings(boundary_mode="none", max_vertices=5,
                               surface_resolution=resolution)
        ).solve(source, lambda *args: None, lambda: False).export_group.net
        interface_generator = 0
        surface_ids = {
            int(index) for index, row in network.surfs.iterrows()
            if len(row["balls"]) == 2 and interface_generator in row["balls"]
        }
        face_records = [
            (sid, tuple(sorted(map(int, network.surfs.loc[sid, "balls"]))),
             _network_face_cycle(network, network.surfs.loc[sid]))
            for sid in sorted(surface_ids)
        ]
        group_a = set(range(1, len(network.balls)))
        baseline_started = time.perf_counter()
        fixed_terms = _aw_supported_gauss_bonnet(
            network, face_records, surface_ids, group_a,
            order=16, smooth_face_method="triangulated",
        )
        baseline_seconds = time.perf_counter() - baseline_started
        non_smooth = fixed_terms["total_GB_left"] - fixed_terms["smooth_K"]
        values = {}
        curved_seconds = {}
        for smooth_order in (2, 4):
            curved_started = time.perf_counter()
            values[smooth_order] = non_smooth + sum(
                _aw_curved_surface_gaussian(network, sid, order=smooth_order)
                for sid in surface_ids
            ) - fixed_terms["two_pi_chi"]
            curved_seconds[smooth_order] = time.perf_counter() - curved_started
        study[resolution] = {
            "values": values,
            "baseline_residual": fixed_terms["residual"],
            "baseline_seconds": baseline_seconds,
            "curved_seconds": curved_seconds,
        }

    print("G12 EDTA residual study:", study)
    for resolution, result in study.items():
        values = result["values"]
        assert abs(values[4]) < abs(result["baseline_residual"])
        assert result["baseline_seconds"] >= 0.0
        assert all(seconds >= 0.0 for seconds in result["curved_seconds"].values())
        if resolution == 0.2:
            assert abs(values[2] - values[4]) > 1e-6
        else:
            assert abs(values[2] - values[4]) < 1e-7
    assert abs(study[0.1]["values"][4]) < abs(study[0.2]["values"][4])
    observed_order = math.log(
        abs(study[0.2]["values"][4] / study[0.1]["values"][4]), 2.0
    )
    assert observed_order > 1.0


def test_edta_alpha_stays_closed_and_reports_unresolved_births(edta_network):
    complex_ = ApolloniusComplex(edta_network, "EDTA")
    complex_.calculate_alpha_births()
    alpha = complex_.alpha_complex(1.4)

    assert len(complex_.alpha_unresolved) == 785
    assert alpha.counts == {0: 612, 1: 0, 2: 0, 3: 0}
    assert alpha.validate() == []
    filtration = complex_.validate_filtration()
    assert len(filtration) == 1266
    assert {issue["kind"] for issue in filtration} == {"unresolved_face"}


def test_aw_interface_curvature_round_trips_on_existing_network_archive(edta_network, tmp_path, monkeypatch):
    from vorpy.src.analyze.aw_interface_curvature import analyze_aw_interface_curvature
    from vorpy.src.io import load_network, save_network
    import numpy as np

    # The supported project environment exposes this NumPy 1.x private helper;
    # NumPy 2.5 removed it. Restore its tiny version dispatch for this test so
    # this remains a serialization regression without changing archive code.
    if not hasattr(np.lib.format, "_read_array_header"):
        def read_array_header(stream, version):
            if version == (1, 0):
                return np.lib.format.read_array_header_1_0(stream)
            if version == (2, 0):
                return np.lib.format.read_array_header_2_0(stream)
            raise ValueError(f"Unsupported NumPy array header version: {version}")
        monkeypatch.setattr(np.lib.format, "_read_array_header", read_array_header, raising=False)

    supported = ApolloniusComplex(edta_network)
    pair = next(key for key in supported.simplices[1] if supported.is_supported(1, key))
    before = analyze_aw_interface_curvature(edta_network, {pair[0]}, {pair[1]})
    archive = save_network(edta_network, tmp_path / "edta.vpy")
    restored = load_network(archive)
    after = analyze_aw_interface_curvature(restored, {pair[0]}, {pair[1]})

    assert [record.surface_id for record in before.selected_surfaces] == [
        record.surface_id for record in after.selected_surfaces
    ]
    assert after.interface_area == pytest.approx(before.interface_area)
    assert after.surface_curvature == pytest.approx(before.surface_curvature)
    assert after.edge_curvature == pytest.approx(before.edge_curvature)


def test_apollonius_cli_options_are_opt_in_and_alpha_implies_complex():
    from vorpy.src.command.vpy_cmnd import apollonius_cli_options

    remaining, requested, alpha = apollonius_cli_options(
        ["molecule.pdb", "-s", "mv", "5", "--alpha-complex", "1.4"]
    )
    assert remaining == ["molecule.pdb", "-s", "mv", "5"]
    assert requested and alpha == 1.4
    remaining, requested, alpha = apollonius_cli_options(["molecule.pdb", "--apollonius"])
    assert remaining == ["molecule.pdb"] and requested and alpha is None

    from vorpy.src.command.vpy_cmnd import Command
    command = Command()
    command._arguments = ["molecule.pdb", "--alpha-complex", "1.4"]
    command.parse_commands()
    assert command.apollonius_requested and command.alpha_value == 1.4


def test_dual_export_contains_primal_multiplicity_and_straight_coordinates(tmp_path):
    import csv
    import json
    from vorpy.src.analyze.apollonius_export import export_apollonius

    complex_ = ApolloniusComplex(_closed_tetra_network(duplicate_edge=True))
    paths = export_apollonius(complex_, tmp_path)
    with paths[1].open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    row = next(row for row in rows if row["generator_ids"] == "[0,1]")
    assert json.loads(row["generator_coordinates"]) == [[0., 0., 0.], [1., 0., 0.]]
    assert int(row["num_primal_components"]) == 1
    assert paths["summary"].exists()
    assert "CYLINDER" in paths["visualization"].read_text(encoding="utf-8")


def test_boundary_generators_remain_labeled_in_dual_metadata():
    network = _network(2, surfaces=[([0, 1], [0], [0, 1])])
    network.balls["is_boundary_generator"] = [False, True]
    complex_ = ApolloniusComplex(network)
    assert complex_.simplices[1][(0, 1)].diagnostics["boundary_generator_flags"] == (False, True)


def test_cli_analysis_writes_full_dual_and_alpha_outputs(tmp_path):
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.tests.network.test_build_net import _analytic_diagnostic_network

    network = _analytic_diagnostic_network()
    system = SimpleNamespace(
        groups=[SimpleNamespace(name="synthetic", net=network)],
        interfaces=[], files={"dir": str(tmp_path)},
    )
    command = Command(sys=system)
    command.apollonius_requested = True
    command.alpha_value = 2.1
    command.run_apollonius()
    output = tmp_path / "apollonius" / "synthetic"
    assert (output / "apollonius_edges.csv").exists()
    assert (output / "apollonius_summary.txt").exists()
    assert (output / "alpha_complex_2.1_summary.txt").exists()
