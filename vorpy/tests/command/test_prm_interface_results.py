"""Bounded PRM coverage for the typed CLI interface Results seam."""

import pandas as pd

from vorpy.src.boundary import SolventMode
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.io.network_archive import _restore_interface_result
from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.results import (
    InterfaceResult,
    Provenance,
    ResultIdentity,
    Status,
    adapt_interface_analysis,
)
from vorpy.src.results.model import dictionary
from vorpy.src.system.system import System
from vorpy.workbench.services.vorpy_backend import (
    InterfaceSolveSelection,
    VorPyBackend,
    VorPySolveSettings,
)


def _prm_interface_system():
    locations = [
        (-1.8, -0.9, -0.4), (-1.1, 0.7, 0.8),
        (-0.5, -0.8, 1.1), (-0.2, 0.9, -1.2),
        (0.5, -0.2, -0.8), (1.0, 1.2, 0.3),
        (1.7, -1.0, 0.9), (1.9, 0.4, -1.1),
    ]
    system = System(
        spheres=pd.DataFrame({
            "loc": locations,
            "rad": [0.7] * len(locations),
            "mass": [12.0] * len(locations),
            "num": list(range(len(locations))),
            "res_name": ["ALA"] * len(locations),
            "name": ["CA"] * len(locations),
        }),
        name="synthetic-prm-interface",
        make_dir=False,
        solvent_mode=SolventMode.MOLECULE_ONLY,
        retained_system_indices=(),
    )
    selection = InterfaceSolveSelection(
        "synthetic-A__synthetic-B", "synthetic-A", (0, 1, 2, 3),
        "synthetic-B", (4, 5, 6, 7),
    )
    interface = VorPyBackend(VorPySolveSettings(
        network_type="prm",
        max_vertices=5,
        surface_resolution=0.2,
        source_network_id="producer:synthetic-prm-interface",
    )).make_interface_network(system, selection)
    return system, interface


def test_cli_prm_interface_results_preserve_scheme_and_missing_curvature(tmp_path):
    system, interface = _prm_interface_system()
    interface.build(analyze_waters=False, calculate_curvature=False)
    system.ifaces.append(interface)

    command = Command(sys=system)
    command.base_file = tmp_path / "synthetic-prm.pdb"
    handoffs = command._make_cli_interface_results(
        [(interface.group1, interface.group2)]
    )

    assert len(handoffs) == 1
    result = handoffs[0].result
    assert isinstance(result, InterfaceResult)
    assert result.identity.partition == "prm"
    assert result.identity.orientation == "group1 -> group2"
    assert result.identity.alpha_method == "squared_radius"
    assert result.identity.alpha_units == "A^2"
    assert result.mean_curvature.integrated_mean_curvature.value is None
    assert result.mean_curvature.integrated_mean_curvature.status is Status.NOT_CALCULATED
    assert result.gaussian_curvature.integrated_gaussian_curvature.value is None
    assert result.gaussian_curvature.integrated_gaussian_curvature.status is Status.NOT_CALCULATED


def test_prm_typed_result_restore_keeps_partition_and_unresolved_values():
    identity = ResultIdentity(
        result_kind="interface",
        system="synthetic-prm",
        frame="static",
        network_id="archive:synthetic-prm",
        network_scope="interface",
        environment="dry",
        partition="prm",
        representation="voronoi",
        radius_configuration_id="radii:synthetic-prm",
        interface_id="A__B",
        group_A="A",
        group_B="B",
        side="shared",
        orientation="group1 -> group2",
    )
    analysis = InterfaceGeometryAnalysis(
        metadata={
            "scheme": "prm",
            "curvature_requested": False,
            "result_status": "PARTIAL",
        },
        voronoi_side_1={},
    )
    result = adapt_interface_analysis(
        analysis,
        identity,
        provenance=Provenance("synthetic PRM interface"),
    )
    restored = _restore_interface_result({
        "format": "vorpy.interface_result_snapshot",
        "schema_version": 1,
        "result_schema_version": dictionary(result)["schema_version"],
        "payload": dictionary(result),
    })

    assert restored.identity.partition == "prm"
    assert restored.metrics["area"].value is None
    assert restored.metrics["area"].status is Status.NOT_CALCULATED
    assert restored.gaussian_curvature.integrated_gaussian_curvature.value is None
    assert restored.gaussian_curvature.integrated_gaussian_curvature.status is Status.NOT_CALCULATED
