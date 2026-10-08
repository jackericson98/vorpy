from results.model import (
    MolecularContactSurfaceGeometry,
    MolecularGeometryCounts,
    Quantity,
    Status,
    TopologicalCellComplex,
)


def q(name, value, units="1", status=Status.CERTIFIED):
    return Quantity(name, value, units, status, "frozen Results fixture", "test fixture")


def test_geometry_and_topology_counts_are_independent():
    surface = MolecularContactSurfaceGeometry(
        quantities={"area": q("area", 10.0, "Å²")},
        status=Status.CERTIFIED,
    )
    geometry = MolecularGeometryCounts(
        quantities={"spherical_patch_count": q("spherical_patch_count", 111)}
    )
    topology = TopologicalCellComplex(
        quantities={"face_count": q("face_count", 70)},
        status=Status.PARTIAL,
    )

    assert surface.quantities["area"].status is Status.CERTIFIED
    assert geometry.quantities["spherical_patch_count"].value == 111
    assert topology.quantities["face_count"].value == 70


def test_unknown_orientability_is_not_nonorientable():
    unknown = TopologicalCellComplex(
        quantities={"orientable": q("orientable", None, status=Status.NOT_CALCULATED)}
    )
    nonorientable = TopologicalCellComplex(
        quantities={"orientable": q("orientable", False)}
    )

    assert unknown.quantities["orientable"].value is None
    assert unknown.quantities["orientable"].status is Status.NOT_CALCULATED
    assert nonorientable.quantities["orientable"].value is False


def test_unavailable_euler_and_genus_remain_null_and_status_bearing():
    topology = TopologicalCellComplex(
        quantities={
            "euler_characteristic": q(
                "euler_characteristic", None, status=Status.UNRESOLVED
            ),
            "genus": q("genus", None, status=Status.NOT_CALCULATED),
        },
        status=Status.PARTIAL,
    )

    assert topology.quantities["euler_characteristic"].value is None
    assert topology.quantities["genus"].value is None


def test_boundary_arc_count_is_not_boundary_component_count():
    geometry = MolecularGeometryCounts(
        quantities={"boundary_arc_count": q("boundary_arc_count", 12)}
    )
    topology = TopologicalCellComplex(
        quantities={"boundary_component_count": q("boundary_component_count", 5)}
    )

    assert geometry.quantities["boundary_arc_count"].value == 12
    assert topology.quantities["boundary_component_count"].value == 5


def test_legacy_patch_count_is_not_promoted_to_topological_face_count():
    legacy = MolecularContactSurfaceGeometry(
        quantities={"patch_count": q("patch_count", 111)},
        status=Status.PARTIAL,
    )

    assert legacy.quantities["patch_count"].value == 111


def test_certified_2kai_power_topology_regressions_keep_edges_and_cuts_distinct():
    a_geometry = MolecularGeometryCounts(
        quantities={"geometric_edge_count": q("geometric_edge_count", 784)}
    )
    a = TopologicalCellComplex(
        quantities={
            "topological_edge_count": q("topological_edge_count", 786),
            "topological_cut_count": q("topological_cut_count", 2),
            "patch_attachment_record_count": q("patch_attachment_record_count", 111),
            "vertex_count": q("vertex_count", 670),
            "face_count": q("face_count", 111),
            "euler_characteristic": q("euler_characteristic", -5),
            "beta0": q("beta0", 1),
            "beta1": q("beta1", 6),
            "beta2": q("beta2", 0),
            "connected_component_count": q("connected_component_count", 1),
            "boundary_component_count": q("boundary_component_count", 5),
            "orientable": q("orientable", True),
            "genus": q("genus", 1),
            "nonorientable_genus": q(
                "nonorientable_genus", None, status=Status.NOT_CALCULATED
            ),
            "manifold_vertex_count": q("manifold_vertex_count", 670),
            "singular_vertex_count": q("singular_vertex_count", 0),
            "chain_complex_valid": q("chain_complex_valid", True),
            "boundary_of_boundary_zero": q("boundary_of_boundary_zero", True),
        },
        status=Status.CERTIFIED,
    )
    b_geometry = MolecularGeometryCounts(
        quantities={"geometric_edge_count": q("geometric_edge_count", 754)}
    )
    b = TopologicalCellComplex(
        quantities={
            "topological_edge_count": q("topological_edge_count", 757),
            "topological_cut_count": q("topological_cut_count", 3),
            "patch_attachment_record_count": q("patch_attachment_record_count", 70),
            "vertex_count": q("vertex_count", 684),
            "face_count": q("face_count", 70),
            "euler_characteristic": q("euler_characteristic", -3),
            "beta0": q("beta0", 1),
            "beta1": q("beta1", 4),
            "beta2": q("beta2", 0),
            "connected_component_count": q("connected_component_count", 1),
            "boundary_component_count": q("boundary_component_count", 5),
            "orientable": q("orientable", True),
            "genus": q("genus", 0),
            "nonorientable_genus": q(
                "nonorientable_genus", None, status=Status.NOT_CALCULATED
            ),
            "manifold_vertex_count": q("manifold_vertex_count", 684),
            "singular_vertex_count": q("singular_vertex_count", 0),
            "chain_complex_valid": q("chain_complex_valid", True),
            "boundary_of_boundary_zero": q("boundary_of_boundary_zero", True),
        },
        status=Status.CERTIFIED,
    )

    assert a_geometry.quantities["geometric_edge_count"].value == 784
    assert a.quantities["topological_edge_count"].value == 786
    assert a.quantities["topological_cut_count"].value == 2
    assert a.quantities["euler_characteristic"].value == -5
    assert a.quantities["genus"].value == 1
    assert b_geometry.quantities["geometric_edge_count"].value == 754
    assert b.quantities["topological_edge_count"].value == 757
    assert b.quantities["topological_cut_count"].value == 3
    assert b.quantities["euler_characteristic"].value == -3
    assert b.quantities["genus"].value == 0
