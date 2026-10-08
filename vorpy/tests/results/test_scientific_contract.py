from dataclasses import fields

import pytest

from vorpy.src.results import (
    GaussianCurvature,
    INTERFACE_LOG_BASE_COLUMNS,
    INTERFACE_LOG_COLUMNS,
    INTERFACE_LOG_NATIVE_UNIT_QUANTITIES,
    INTERFACE_LOG_QUANTITY_UNITS,
    INTERFACE_LOG_REPRESENTATIONS,
    INTERFACE_LOG_STATUS_VALUES,
    InterfaceLogRow,
    MeanCurvature,
    Provenance,
    Quantity,
    RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG,
    RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG,
    RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG,
    Status,
    canonical_interface_log_representation,
)


FROZEN_2KAI = {
    "A": {
        "area": 638.9827835376809,
        "smooth_k": 67.582108181661,
        "intrinsic_discrete": 109.29089804245606,
        "intrinsic_seam": -207.98185596219884,
        "boundary_geodesic": 50.837993210148376,
        "boundary_corner": -51.14507000796407,
        "completed_gb": -31.41592653589747,
        "expected_gb": -31.41592653589793,
        "gb_residual": 4.618527782440651e-13,
        "smooth_h": 207.3596759415419,
        "seam_h": 232.63605596302926,
        "total_h": 439.99573190457113,
    },
    "B": {
        "area": 774.8019193578036,
        "smooth_k": 80.6368577799271,
        "intrinsic_discrete": 42.325534352671525,
        "intrinsic_seam": -117.17541458536404,
        "boundary_geodesic": 51.12287325648051,
        "boundary_corner": -75.75940672525387,
        "completed_gb": -18.849555921538773,
        "expected_gb": -18.84955592153876,
        "gb_residual": -1.4210854715202004e-14,
        "smooth_h": 249.41721941700987,
        "seam_h": 158.29373032675855,
        "total_h": 407.7109497437684,
    },
}


def _quantity(name, value, units, status, provenance, scope):
    return Quantity(name, value, units, status, "frozen 2KAI curvature audit", scope, provenance=provenance)


@pytest.mark.parametrize("side", ["A", "B"])
def test_frozen_2kai_curvature_fields_are_separate_and_status_bearing(side):
    values = FROZEN_2KAI[side]
    provenance = Provenance(
        f"frozen 2KAI Power molecular contact surface {side}",
        source_revision="2kai-power-contact-v1",
    )
    mean = MeanCurvature(
        _quantity("integrated_mean_curvature", values["total_h"], "Å", Status.CERTIFIED, provenance, "canonical_surface"),
        _quantity("complete_mean_curvature", values["total_h"], "Å", Status.CERTIFIED, provenance, "canonical_surface"),
        _quantity("smooth_face_contribution", values["smooth_h"], "Å", Status.CERTIFIED, provenance, "spherical_patches"),
        _quantity("edge_contribution", values["seam_h"], "Å", Status.CERTIFIED, provenance, "physical_refined_seams"),
        _quantity("boundary_contribution", None, "Å", Status.NOT_SUPPORTED, provenance, "free_boundary_convention"),
        "H=(k1+k2)/2; smooth faces + 1/2 signed physical seam dihedral",
    )
    intrinsic_total = values["smooth_k"] + values["intrinsic_discrete"] + values["intrinsic_seam"]
    gaussian = GaussianCurvature(
        _quantity("integrated_gaussian_curvature", intrinsic_total, "1", Status.CERTIFIED, provenance, "intrinsic_surface"),
        _quantity("smooth_face_contribution", values["smooth_k"], "1", Status.CERTIFIED, provenance, "spherical_patches"),
        _quantity("intrinsic_discrete_contribution", values["intrinsic_discrete"], "1", Status.CERTIFIED, provenance, "interior_junctions"),
        _quantity("intrinsic_seam_contribution", values["intrinsic_seam"], "1", Status.CERTIFIED, provenance, "physical_refined_seams"),
        _quantity("interior_vertex_defects", values["intrinsic_discrete"], "1", Status.CERTIFIED, provenance, "interior_junctions"),
        _quantity("boundary_geodesic_contribution", values["boundary_geodesic"], "1", Status.CERTIFIED, provenance, "physical_boundary_arcs"),
        _quantity("boundary_corner_contribution", values["boundary_corner"], "1", Status.CERTIFIED, provenance, "physical_boundary_junctions"),
        _quantity("gauss_bonnet_completion", values["completed_gb"], "1", Status.CERTIFIED, provenance, "gauss_bonnet_completed_surface"),
        _quantity("gauss_bonnet_expected", values["expected_gb"], "1", Status.CERTIFIED, provenance, "certified_topology"),
        _quantity("gauss_bonnet_residual", values["gb_residual"], "1", Status.CERTIFIED, provenance, "gauss_bonnet_audit"),
        "intrinsic surface total; exterior boundary terms complete Gauss-Bonnet",
    )
    area = InterfaceLogRow(
        "AB_I",
        "molecular_contact_surface",
        side,
        "area",
        values["area"],
        "Å²",
        Status.CERTIFIED,
        provenance,
    )

    assert area.value == pytest.approx(values["area"])
    assert '"source_revision":"2kai-power-contact-v1"' in area.as_csv_row()["provenance"]
    assert mean.integrated_mean_curvature.value == pytest.approx(values["total_h"])
    assert mean.complete_total.value == pytest.approx(values["total_h"])
    assert mean.smooth_face_contribution.value == pytest.approx(values["smooth_h"])
    assert mean.edge_contribution.value == pytest.approx(values["seam_h"])
    assert mean.boundary_contribution.value is None
    assert mean.boundary_contribution.status is Status.NOT_SUPPORTED
    assert gaussian.smooth_face_contribution.value == pytest.approx(values["smooth_k"])
    assert gaussian.smooth_integrated_gaussian_curvature.value == pytest.approx(values["smooth_k"])
    assert gaussian.intrinsic_total_gaussian_curvature.value == pytest.approx(intrinsic_total)
    assert gaussian.integrated_gaussian_curvature.value != pytest.approx(values["completed_gb"])
    assert gaussian.gauss_bonnet_completion.value == pytest.approx(values["completed_gb"])
    assert gaussian.gauss_bonnet_residual.value == pytest.approx(values["gb_residual"])
    assert all(quantity.units == "Å" for quantity in vars(mean).values() if isinstance(quantity, Quantity))
    assert all(quantity.units == "1" for quantity in vars(gaussian).values() if isinstance(quantity, Quantity))
    assert all(quantity.provenance is provenance for quantity in (*vars(mean).values(), *vars(gaussian).values()) if isinstance(quantity, Quantity))


def test_interface_log_contract_is_canonical_and_deterministic():
    assert INTERFACE_LOG_COLUMNS[:7] == INTERFACE_LOG_BASE_COLUMNS
    assert "provenance" in INTERFACE_LOG_COLUMNS
    assert len(INTERFACE_LOG_COLUMNS) == 8
    assert len(set(INTERFACE_LOG_COLUMNS)) == len(INTERFACE_LOG_COLUMNS)
    assert "AVAILABLE" not in INTERFACE_LOG_STATUS_VALUES
    assert set(INTERFACE_LOG_REPRESENTATIONS) == {
        "physical_partition_interface",
        "molecular_contact_surface",
        "dual_contact_complex",
        "dual_generator_complex",
        "dual_separator",
        "alpha_selection",
        "water_analysis",
        "timing",
    }
    assert canonical_interface_log_representation("voronoi") == "physical_partition_interface"
    assert INTERFACE_LOG_QUANTITY_UNITS["area"] == "Å²"
    assert INTERFACE_LOG_NATIVE_UNIT_QUANTITIES == frozenset({"alpha_value"})
    assert RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG[
        "gaussian_curvature.gauss_bonnet_completion"
    ] == "completed_gauss_bonnet"
    row = InterfaceLogRow(
        "AB_I",
        "molecular_contact_surface",
        "A",
        "alpha_selected_generator_pairs",
        {"z": [2, 3], "a": [0, 1]},
        "1",
        Status.CERTIFIED,
        {"partition": "power", "source": "frozen test"},
    )
    csv_row = row.as_csv_row()
    assert tuple(csv_row) == INTERFACE_LOG_COLUMNS
    assert csv_row["value"] == '{"a":[0,1],"z":[2,3]}'
    assert csv_row["status"] == "CERTIFIED"
    assert csv_row["provenance"] == '{"partition":"power","source":"frozen test"}'


def test_curvature_fields_have_unique_canonical_log_names():
    paths = {
        f"mean_curvature.{field.name}"
        for field in fields(MeanCurvature)
        if field.name != "convention"
    } | {
        f"gaussian_curvature.{field.name}"
        for field in fields(GaussianCurvature)
        if field.name != "convention"
    }
    assert paths == set(RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG)
    log_names = tuple(RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG.values())
    assert len(log_names) == len(set(log_names))
    assert "smooth_integrated_gaussian_curvature" not in {
        field.name for field in fields(GaussianCurvature)
    }
    assert RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG == {
        "molecular_selection.selected_contacts": "selected_ab_contact_count",
        "molecular_selection.unresolved_contacts": "unresolved_contact_count",
    }
    assert RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG == {
        "metrics.area": "area",
        "metrics.n_contacts": "contact_count",
        "metrics.n_residue_contacts": "residue_contact_count",
    }
    assert len(set(RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG.values())) == len(
        RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG
    )
    assert INTERFACE_LOG_QUANTITY_UNITS["mean_curvature_per_area"] == "Å^-1"
    assert INTERFACE_LOG_QUANTITY_UNITS["gaussian_curvature_per_area"] == "Å^-2"
    molecular_log_names = tuple(RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG.values())
    assert len(molecular_log_names) == len(set(molecular_log_names))
    assert canonical_interface_log_representation("dual_a") == "dual_generator_complex"
    assert canonical_interface_log_representation("dual_b") == "dual_generator_complex"


@pytest.mark.parametrize(
    ("status", "value"),
    [
        (Status.CERTIFIED, 1.0),
        (Status.PARTIAL, 1.0),
        (Status.STALE, 1.0),
        (Status.STALE, None),
        (Status.UNRESOLVED, None),
        (Status.NOT_CALCULATED, None),
        (Status.NOT_SUPPORTED, None),
    ],
)
def test_interface_log_status_value_pairs_are_explicit(status, value):
    row = InterfaceLogRow(
        "AB_I", "physical_partition_interface", "shared", "area", value, "Å²", status
    )
    assert row.as_csv_row()["status"] == status.value
    assert row.as_csv_row()["value"] == ("" if value is None else str(value))


@pytest.mark.parametrize(
    ("representation", "quantity", "units", "value"),
    [
        ("physical_partition_interface", "mean_curvature_per_area", "Å^-1", -0.25),
        ("physical_partition_interface", "gaussian_curvature_per_area", "Å^-2", 0.01),
        ("alpha_selection", "alpha_candidate_count", "1", 2),
        ("alpha_selection", "alpha_value", "Å²", 0.0),
        ("water_analysis", "water_component_count", "1", 1),
        ("timing", "timing_total", "s", 0.25),
    ],
)
def test_non_curvature_log_families_have_explicit_units(
    representation, quantity, units, value
):
    row = InterfaceLogRow("AB_I", representation, "shared", quantity, value, units, Status.PARTIAL)
    assert row.units == units

    if quantity != "alpha_value":
        wrong_units = "s" if units != "s" else "1"
        with pytest.raises(ValueError):
            InterfaceLogRow(
                "AB_I", representation, "shared", quantity, value, wrong_units, Status.PARTIAL
            )


def test_interface_log_missing_is_not_zero_and_stale_may_retain_value():
    unavailable = InterfaceLogRow(
        "AB_I", "physical_partition_interface", "shared", "area", None, "Å²", Status.NOT_CALCULATED
    )
    assert unavailable.as_csv_row()["value"] == ""
    zero = InterfaceLogRow(
        "AB_I", "physical_partition_interface", "shared", "area", 0.0, "Å²", Status.CERTIFIED
    )
    assert zero.as_csv_row()["value"] == "0.0"
    stale = InterfaceLogRow(
        "AB_I", "physical_partition_interface", "shared", "area", 3.0, "Å²", Status.STALE
    )
    assert stale.as_csv_row()["status"] == "STALE"
    with pytest.raises(ValueError):
        InterfaceLogRow(
            "AB_I", "physical_partition_interface", "shared", "area", 0.0, "Å²", Status.NOT_SUPPORTED
        )
