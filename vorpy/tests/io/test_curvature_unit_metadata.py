"""Regression tests for pointwise versus integrated curvature units."""

from vorpy.src.results.model import (
    INTERFACE_LOG_QUANTITY_UNITS,
    InterfaceLogRow,
    Status,
)


def test_curvature_units_distinguish_pointwise_and_integrated_quantities():
    assert INTERFACE_LOG_QUANTITY_UNITS["mean_curvature_per_area"] == "Å^-1"
    assert INTERFACE_LOG_QUANTITY_UNITS["gaussian_curvature_per_area"] == "Å^-2"
    assert INTERFACE_LOG_QUANTITY_UNITS["total_integrated_mean_curvature"] == "Å"
    assert INTERFACE_LOG_QUANTITY_UNITS["intrinsic_integrated_gaussian_curvature"] == "1"


def test_integrated_curvature_units_survive_log_serialization():
    mean = InterfaceLogRow(
        "iface", "physical_partition_interface", "shared",
        "total_integrated_mean_curvature", 1.25, "Å", Status.CERTIFIED,
    ).as_csv_row()
    gaussian = InterfaceLogRow(
        "iface", "physical_partition_interface", "shared",
        "intrinsic_integrated_gaussian_curvature", 2.0, "1", Status.CERTIFIED,
    ).as_csv_row()

    assert mean["units"] == "Å"
    assert gaussian["units"] == "1"
