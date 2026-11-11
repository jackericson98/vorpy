from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.boundary import (
    BoundaryConfig,
    BoundaryGenerator,
    SolventMode,
    prepare_network_geometry,
)
from vorpy.src.group.shared_view import (
    analyze_group_view,
    attach_group_view_results,
    make_group_view,
    shared_network_identity,
    system_ids_for_group,
)
from vorpy.src.io.network_archive import load_network, save_network
from vorpy.src.network.network import Network
from vorpy.src.network.views import NetworkView, NetworkViewError
from vorpy.src.results.model import Status
from vorpy.src.system.system import System
from vorpy.workbench.domain import AnalysisResult
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings


CAMBRIN_ARCHIVE = (
    Path(__file__).resolve().parents[3] / "output" / "cambrin_33" / "cambrin.vpy"
)
EDTA_PDB = Path(__file__).resolve().parents[3] / "vorpy" / "data" / "EDTA.pdb"


def _solved_network():
    locations = [
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 2.0),
        (4.0, 4.0, 4.0),
    ]
    balls = pd.DataFrame({
        "loc": locations,
        "rad": [1.0] * 5,
        "num": list(range(5)),
        "system_num": [10, 11, 12, 13, None],
        "stable_id": ["atom-A", "atom-B", "atom-C", "atom-D", None],
        "is_boundary_generator": [False, False, False, False, True],
        "complete": [True] * 5,
        "vol": [2.0, 3.0, 5.0, 7.0, 0.0],
    })
    verts = pd.DataFrame({
        "balls": [[0, 1, 4]] * 4,
        "loc": [(0.2, 0.2, 0.2)] * 4,
        "rad": [0.2] * 4,
        "complete": [True] * 4,
    })
    edges = pd.DataFrame({
        "balls": [[0, 1, 4]] * 4,
        "verts": [[0, 1], [1, 2], [2, 3], [3, 0]],
        "points": [[(0.0, 0.0, 0.0), (0.1, 0.1, 0.1)]] * 4,
        "surfs": [[0, 1]] * 4,
        "complete": [True] * 4,
    })
    surfaces = pd.DataFrame({
        "balls": [[0, 1, 4], [0, 2]],
        "verts": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "edges": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "points": [[(0.0, 0.0, 0.0)] * 4] * 2,
        "tris": [[(0, 1, 2)]] * 2,
        "sa": [2.5, 1.5],
        "complete": [True, True],
    })
    network = Network(
        locs=locations,
        rads=[1.0] * 5,
        group=[10, 11, 12, 13],
        settings={
            "net_type": "prm",
            "surf_res": 0.25,
            "max_vert": 5.0,
            "box_size": 1.5,
            "num_splits": 2,
            "build_type": "all",
        },
        balls=balls,
        verts=verts,
        edges=edges,
        surfs=surfaces,
        system=SimpleNamespace(frame_index=7),
    )
    network.boundary_config = BoundaryConfig(mode="shell")
    network.boundary_mode = "shell"
    network.boundary_indices = (4,)
    network.boundary_generators = (
        BoundaryGenerator(position=locations[4], radius=1.0),
    )
    return network


def _molecule_only_network():
    """A solved universe with input-global IDs that leave solvent IDs absent."""
    network = _solved_network()
    system = System(
        spheres=pd.DataFrame({
            "loc": [
                (0.0, 0.0, 0.0),
                (0.5, 0.5, 0.5),
                (2.0, 0.0, 0.0),
                (2.5, 0.5, 0.5),
                (0.0, 2.0, 0.0),
                (0.5, 2.5, 0.5),
                (0.0, 0.0, 2.0),
                (0.5, 0.5, 2.5),
                (4.0, 4.0, 4.0),
            ],
            "rad": [1.0] * 9,
            "mass": [1.0] * 9,
            "num": list(range(100, 109)),
            "res_name": ["ALA", "HOH", "GLY", "HOH", "SER", "HOH", "THR", "HOH", "LEU"],
            "name": ["CA", "O", "CA", "O", "CA", "O", "CA", "O", "CA"],
            "res": [0] * 9,
            "chn": [0] * 9,
        }),
        name="molecule-only",
        make_dir=False,
        boundary_config=BoundaryConfig("none"),
    )
    geometry = prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)
    assert geometry.system_indices == (0, 2, 4, 6, 8)
    network.balls["system_num"] = list(geometry.system_indices)
    network.balls["stable_id"] = list(geometry.stable_ids)
    network.balls["is_boundary_generator"] = False
    network.group = list(range(len(geometry.system_indices)))
    network.sys = system
    network.boundary_indices = geometry.boundary_indices
    network.boundary_config = geometry.boundary_config
    network.boundary_mode = system.boundary_mode
    network.boundary_generators = geometry.boundary_generators
    network.source_network_id = "molecule-only-test-source"
    network.generator_system_indices = geometry.system_indices
    network.generator_stable_ids = geometry.stable_ids
    network.solvent_mode = geometry.solvent_mode
    network.solvent_local_indices = geometry.solvent_local_indices
    network.solvent_generator_ids = geometry.solvent_generator_ids
    network.refresh_solvent_provenance()
    network.system = SimpleNamespace(
        frame_index=7,
        balls=pd.DataFrame({
            "num": [2, 3, 4, 5, 6, 7, 8, 9, 10],
            "res_name": ["C", "HOH", "C", "HOH", "C", "HOH", "C", "HOH", "C"],
        }),
    )
    return network


def test_two_groups_consume_one_solved_network_without_rebuild(monkeypatch):
    network = _solved_network()
    balls_before = network.balls.copy(deep=True)
    surfaces_before = network.surfs.copy(deep=True)

    def fail(*args, **kwargs):
        raise AssertionError("shared-view analysis rebuilt the source network")

    monkeypatch.setattr(Network, "build", fail)

    group_a = analyze_group_view(NetworkView.for_group(network, [10]))
    group_b = analyze_group_view(NetworkView.for_group(network, [11]))

    assert group_a.source_network_id == group_b.source_network_id
    assert group_a.solve_universe_id == group_b.solve_universe_id
    assert group_a.selection_id != group_b.selection_id
    assert group_a.system_ids == (10,)
    assert group_b.system_ids == (11,)
    assert group_a.volume == pytest.approx(2.0)
    assert group_b.volume == pytest.approx(3.0)
    assert group_a.area == pytest.approx(4.0)
    assert group_b.area == pytest.approx(2.5)
    assert group_a.competitor_system_ids == (11, 12)
    assert group_b.competitor_system_ids == (10,)
    assert group_a.boundary_surface_ids == (0, 1)
    assert group_b.boundary_surface_ids == (0,)
    assert group_a.complete and group_b.complete
    pd.testing.assert_frame_equal(network.balls, balls_before)
    pd.testing.assert_frame_equal(network.surfs, surfaces_before)


def test_molecule_only_views_preserve_input_global_ids_and_exclude_solvent(monkeypatch):
    network = _molecule_only_network()
    balls_before = network.balls.copy(deep=True)
    surfaces_before = network.surfs.copy(deep=True)

    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("molecule-only view rebuilt Network"),
    )

    context = network.get_solve_context_key()
    assert context.local_to_system == (
        (0, 0),
        (1, 2),
        (2, 4),
        (3, 6),
        (4, 8),
    )
    assert context.physical_generator_ids == (
        "system:0",
        "system:2",
        "system:4",
        "system:6",
        "system:8",
    )
    assert network.solvent_provenance == {
        "solvent_mode": "MOLECULE_ONLY",
        "solvent_generator_digest": None,
        "solvent_cells_complete": False,
        "solvent_local_indices": (),
        "solvent_generator_ids": (),
        "generator_membership": (
            (0, 0, "system:0", False, False),
            (1, 2, "system:2", False, False),
            (2, 4, "system:4", False, False),
            (3, 6, "system:6", False, False),
            (4, 8, "system:8", False, False),
        ),
    }

    view_a = make_group_view(SimpleNamespace(ball_ndxs=[0]), network)
    view_b = make_group_view(SimpleNamespace(ball_ndxs=[4]), network)
    group_a = analyze_group_view(view_a)
    group_b = analyze_group_view(view_b)
    identity_a = shared_network_identity(view_a, network)
    identity_b = shared_network_identity(view_b, network)

    assert group_a.system_ids == (0,)
    assert group_b.system_ids == (4,)
    assert group_a.source_network_id == group_b.source_network_id
    assert group_a.solve_universe_id == group_b.solve_universe_id
    assert group_a.selection_id != group_b.selection_id
    assert identity_a.solvent_mode == "MOLECULE_ONLY"
    assert identity_a.solve_universe_id == identity_b.solve_universe_id
    assert identity_a.selection_id != identity_b.selection_id
    assert group_a.volume == pytest.approx(2.0)
    assert group_b.volume == pytest.approx(5.0)
    assert group_a.area == pytest.approx(4.0)
    assert group_b.area == pytest.approx(1.5)
    system_ids = {system_id for _, system_id in context.local_to_system}
    assert system_ids == {0, 2, 4, 6, 8}
    assert not {1, 3, 5, 7}.intersection(system_ids)
    assert group_a.complete and group_b.complete
    assert system_ids_for_group(SimpleNamespace(ball_ndxs=[0]), network) == (0,)
    with pytest.raises(NetworkViewError, match="outside the solved generator universe"):
        make_group_view(SimpleNamespace(ball_ndxs=[1]), network)
    pd.testing.assert_frame_equal(network.balls, balls_before)
    pd.testing.assert_frame_equal(network.surfs, surfaces_before)


def test_molecule_only_views_preserve_metrics_after_archive_round_trip(
    monkeypatch, tmp_path
):
    network = _molecule_only_network()
    network.sys = System(
        file=str(Path("vorpy/data/EDTA.pdb")),
        name="archive-system",
        make_dir=False,
        boundary_config=BoundaryConfig("none"),
    )
    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("archive round-trip rebuilt Network"),
    )

    original = tuple(
        analyze_group_view(
            make_group_view(SimpleNamespace(ball_ndxs=[system_id]), network)
        )
        for system_id in (0, 4)
    )
    archive = tmp_path / "molecule-only.vpy"
    save_network(network, archive)
    restored = load_network(archive)
    restored_results = tuple(
        analyze_group_view(
            make_group_view(SimpleNamespace(ball_ndxs=[system_id]), restored)
        )
        for system_id in (0, 4)
    )

    assert restored.solvent_provenance == network.solvent_provenance
    assert restored.get_solve_context_key() == network.get_solve_context_key()
    assert restored_results[0].solve_universe_id == original[0].solve_universe_id
    assert restored_results[1].solve_universe_id == original[1].solve_universe_id
    assert restored_results[0].selection_id != restored_results[1].selection_id
    for before, after in zip(original, restored_results):
        assert after.system_ids == before.system_ids
        assert after.area == pytest.approx(before.area)
        assert after.volume == pytest.approx(before.volume)


def test_n8_edta_shell_shared_group_results_round_trip(monkeypatch, tmp_path):
    """One F11 EDTA solve feeds two typed group Results without a rebuild."""
    from vorpy.src.results.shared_network import SharedNetworkIdentity

    build_calls = []
    original_build = Network.build

    def counted_build(*args, **kwargs):
        build_calls.append(True)
        return original_build(*args, **kwargs)

    monkeypatch.setattr(Network, "build", counted_build)
    backend = VorPyBackend(
        VorPySolveSettings(
            network_type="aw",
            boundary_mode="shell",
            typed_group_selections={
                "edta-group-0": (0,),
                "edta-group-1": (1,),
            },
            source_network_id="producer:edta-n8-shell",
        )
    )
    result = backend.solve(
        EDTA_PDB,
        progress=lambda *args, **kwargs: None,
        is_cancelled=lambda: False,
    )
    assert len(build_calls) == 1
    network = result.export_group.net
    assert getattr(network, "archive_id", None) is None
    assert network.source_network_id == "producer:edta-n8-shell"
    assert network.solvent_mode == SolventMode.MOLECULE_ONLY.value
    boundary_indices = set(network.boundary_indices)
    physical_indices = [
        index for index in range(len(network.balls)) if index not in boundary_indices
    ]
    assert len(physical_indices) == 32
    assert len(boundary_indices) == 242
    assert all(bool(network.balls.iloc[index]["complete"]) for index in physical_indices)
    assert network.retained_generator_provenance == {
        "policy": "explicit_system_indices",
        "retained_system_indices": (),
        "retained_generator_ids": (),
        "retained_generator_digest": None,
        "boundary_independent": True,
    }

    context = network.get_solve_context_key()
    physical_system_ids = tuple(
        system_id
        for _, system_id in context.local_to_system
        if system_id is not None
    )
    assert physical_system_ids == tuple(range(32))
    assert len(network.generator_membership_evidence) == len(network.balls)
    assert all(
        entry[4] is True
        for entry in network.generator_membership_evidence
        if entry[0] in boundary_indices
    )

    archive = tmp_path / "edta-n8.vpy"
    save_network(network, archive)
    source_balls = network.balls.copy(deep=True)
    source_surfaces = network.surfs.copy(deep=True)
    source_atoms = network.sys.balls.copy(deep=True)

    groups = (
        SimpleNamespace(ball_ndxs=[0], group_id="edta-group-0"),
        SimpleNamespace(ball_ndxs=[1], group_id="edta-group-1"),
    )
    typed_results = result.typed_results
    assert result.typed_results == typed_results
    assert len(typed_results) == 2
    assert result.typed_selection_system_ids == {
        "edta-group-0": (0,), "edta-group-1": (1,)
    }
    assert all(typed.status is Status.CERTIFIED for typed in typed_results)
    for typed in typed_results:
        identity = SharedNetworkIdentity(
            **typed.provenance.details["shared_network_identity"]
        )
        assert not identity.provisional_source_network
        assert typed.identity.source_network_id == "producer:edta-n8-shell"
        assert typed.identity.source_network_id != network.archive_id
    assert typed_results[0].identity.solve_universe_id == typed_results[1].identity.solve_universe_id
    assert typed_results[0].identity.selection_id != typed_results[1].identity.selection_id
    assert tuple(item.identity.group_id for item in typed_results) == (
        "edta-group-0", "edta-group-1"
    )
    assert all(
        typed.provenance.details["source_retained_generator_provenance"]
        == network.retained_generator_provenance
        for typed in typed_results
    )

    system_to_local = {
        system_id: local_id
        for local_id, system_id in context.local_to_system
        if system_id is not None
    }
    analyses = []
    for group in groups:
        selected_local = system_to_local[group.ball_ndxs[0]]
        reference_volume = float(network.balls.iloc[selected_local]["vol"])
        reference_area = 0.0
        reference_boundary_area = 0.0
        reference_physical_area = 0.0
        for _, surface in network.surfs.iterrows():
            generator_indices = tuple(int(value) for value in surface["balls"])
            if selected_local not in generator_indices:
                continue
            generator_system_ids = tuple(
                context.local_to_system[index][1] for index in generator_indices
            )
            competitors = {
                value for value in generator_system_ids
                if value is not None and value != group.ball_ndxs[0]
            }
            if not competitors and None not in generator_system_ids:
                continue
            reference_area += float(surface["sa"])
            if None in generator_system_ids:
                reference_boundary_area += float(surface["sa"])
            else:
                reference_physical_area += float(surface["sa"])

        view = make_group_view(group, network)
        analysis = analyze_group_view(view)
        analyses.append(analysis)
        assert analysis.volume == pytest.approx(reference_volume, abs=1e-12)
        assert analysis.area == pytest.approx(reference_area, abs=1e-12)
        assert analysis.boundary_dependent_area == pytest.approx(
            reference_boundary_area, abs=1e-12
        )
        assert analysis.physical_interface_area == pytest.approx(
            reference_physical_area, abs=1e-12
        )

    expected_metrics = (
        (5.676116886923854, 19.53014188145206, 0.0, 19.53014188145206),
        (
            14.170714484103758,
            46.75203064184723,
            2.202876203619455,
            44.54915443822777,
        ),
    )
    for analysis, expected in zip(analyses, expected_metrics, strict=True):
        assert analysis.volume == pytest.approx(expected[0], abs=1e-12)
        assert analysis.area == pytest.approx(expected[1], abs=1e-12)
        assert analysis.boundary_dependent_area == pytest.approx(
            expected[2], abs=1e-12
        )
        assert analysis.physical_interface_area == pytest.approx(
            expected[3], abs=1e-12
        )

    for typed, analysis in zip(typed_results, analyses, strict=True):
        assert typed.metrics["area"].value == pytest.approx(analysis.area, abs=1e-12)
        assert typed.metrics["volume"].value == pytest.approx(analysis.volume, abs=1e-12)
        assert typed.metrics["area"].units == "\u00c5\u00b2"
        assert typed.metrics["volume"].units == "\u00c5\u00b3"
        assert typed.metrics["area"].scope == "group_view"
        assert typed.provenance.details["area_components"] == {
            "total": analysis.area,
            "boundary_dependent": analysis.boundary_dependent_area,
            "physical_interface": analysis.physical_interface_area,
            "components_proven": True,
            "units": "\u00c5\u00b2",
        }
        assert typed.provenance.details["source_solvent_provenance"] == network.solvent_provenance

    pd.testing.assert_frame_equal(network.balls, source_balls)
    pd.testing.assert_frame_equal(network.surfs, source_surfaces)
    pd.testing.assert_frame_equal(network.sys.balls, source_atoms)

    restored = load_network(archive)
    assert restored.source_network_id == "producer:edta-n8-shell"
    restored_result = AnalysisResult(source=EDTA_PDB, name="EDTA restored")
    restored_typed = attach_group_view_results(
        restored_result,
        restored,
        groups,
        system="EDTA",
        radius_configuration_id="vorpy-backend-radii-v1",
    )
    assert restored.solvent_provenance == network.solvent_provenance
    assert restored.get_solve_context_key() == network.get_solve_context_key()
    for before, after, analysis in zip(typed_results, restored_typed, analyses, strict=True):
        assert after.status is Status.CERTIFIED
        assert after.identity.source_network_id == "producer:edta-n8-shell"
        assert after.identity.source_network_id != restored.archive_id
        assert after.identity.solve_universe_id == before.identity.solve_universe_id
        assert after.identity.selection_id == before.identity.selection_id
        assert after.provenance.details["dependency_digest"] == (
            before.provenance.details["dependency_digest"]
        )
        assert after.metrics["area"].value == pytest.approx(
            analysis.area, abs=1e-12
        )
        assert after.metrics["volume"].value == pytest.approx(
            analysis.volume, abs=1e-12
        )
    assert restored_typed[0].identity.selection_id != restored_typed[1].identity.selection_id
    assert len(build_calls) == 1


def test_n10_retained_mg_backend_views_keep_f12_and_global_ids(monkeypatch, tmp_path):
    """A retained ion is a competitor, not a selected group or a second solve."""
    from vorpy.src.results.shared_network import ResultsDependency, SharedNetworkIdentity
    from vorpy.src.group import shared_view as shared_view_module

    build_calls = []
    original_build = Network.build
    original_attach = shared_view_module.attach_group_view_results
    source_snapshots = []

    def counted_build(*args, **kwargs):
        build_calls.append(True)
        return original_build(*args, **kwargs)

    def checked_attach(result, network, groups, **kwargs):
        before = tuple(
            frame.copy(deep=True)
            for frame in (network.balls, network.verts, network.edges, network.surfs,
                          network.sys.balls)
        )
        typed = original_attach(result, network, groups, **kwargs)
        for current, saved in zip(
            (network.balls, network.verts, network.edges, network.surfs,
             network.sys.balls),
            before,
            strict=True,
        ):
            pd.testing.assert_frame_equal(current, saved)
        source_snapshots.append(True)
        return typed

    monkeypatch.setattr(Network, "build", counted_build)
    monkeypatch.setattr(shared_view_module, "attach_group_view_results", checked_attach)
    result = VorPyBackend(
        VorPySolveSettings(
            network_type="aw",
            boundary_mode="shell",
            retained_system_indices=(32,),
        )
    ).solve(EDTA_PDB, lambda *_args: None, lambda: False)

    network = result.export_group.net
    source_network_id = "producer:edta-retained-mg"
    group_ids = []
    for system_id in range(32):
        group = SimpleNamespace(
            ball_ndxs=(system_id,), group_id=f"edta-atom-{system_id}"
        )
        try:
            make_group_view(group, network, source_network_id=source_network_id)
        except NetworkViewError:
            continue
        group_ids.append(system_id)
        if len(group_ids) == 2:
            break
    assert len(group_ids) == 2, "retained Mg leaves fewer than two complete group views"
    selections = {
        f"edta-atom-{system_id}": (system_id,) for system_id in group_ids
    }
    groups = tuple(
        SimpleNamespace(name=name, group_id=name, ball_ndxs=ids)
        for name, ids in selections.items()
    )
    shared_view_module.attach_group_view_results(
        result,
        network,
        groups,
        system="EDTA",
        radius_configuration_id="vorpy-backend-radii-v1",
        source_network_id=source_network_id,
    )
    result.typed_selection_system_ids = selections
    context = network.get_solve_context_key()
    system_to_local = {
        system_id: local_id
        for local_id, system_id in context.local_to_system
        if system_id is not None
    }
    assert len(build_calls) == 1
    assert getattr(network, "archive_id", None) is None
    assert network.source_network_id == source_network_id
    assert source_snapshots == [True]
    assert len(result.typed_results) == 2
    assert result.typed_selection_system_ids == selections
    assert len(network.sys.balls) == 612  # Original solvent-containing System is retained.
    assert tuple(system_to_local) == tuple(range(33))
    assert network.group == list(range(32))
    assert network.retained_generator_provenance["retained_system_indices"] == (32,)
    assert network.retained_generator_provenance["retained_generator_ids"] == (
        "system:32",
    )
    assert not bool(network.balls.iloc[system_to_local[32]]["complete"])

    for typed, (group_id, system_ids) in zip(
        result.typed_results, selections.items(), strict=True
    ):
        local_id = system_to_local[system_ids[0]]
        assert typed.status is Status.CERTIFIED
        assert typed.identity.source_network_id == source_network_id
        assert typed.identity.source_network_id != getattr(network, "archive_id", None)
        assert not hasattr(typed.identity, "retained_system_indices")
        shared_identity = typed.provenance.details["shared_network_identity"]
        assert shared_identity["retained_system_indices"] == (32,)
        assert shared_identity["retained_generator_ids"] == ("system:32",)
        assert typed.provenance.details["source_retained_generator_provenance"] == (
            network.retained_generator_provenance
        )
        assert typed.identity.group_id == group_id
        assert typed.metrics["volume"].value == pytest.approx(
            float(network.balls.iloc[local_id]["vol"]), abs=1e-12
        )

        reference_area = 0.0
        for _, surface in network.surfs.iterrows():
            local_generators = tuple(int(value) for value in surface["balls"])
            if local_id not in local_generators:
                continue
            system_generators = tuple(
                context.local_to_system[index][1] for index in local_generators
            )
            competitors = {
                value for value in system_generators
                if value is not None and value not in system_ids
            }
            if competitors or None in system_generators:
                reference_area += float(surface["sa"])
        assert typed.metrics["area"].value == pytest.approx(
            reference_area, abs=1e-12
        )

    assert (
        result.typed_results[0].identity.solve_universe_id
        == result.typed_results[1].identity.solve_universe_id
    )
    assert (
        result.typed_results[0].identity.selection_id
        != result.typed_results[1].identity.selection_id
    )

    # The NEXUS cache seam forwards stale lineage to ATLAS without rebuilding.
    stale_dependencies = {
        typed.identity.selection_id: ResultsDependency(
            identity=SharedNetworkIdentity(
                **typed.provenance.details["shared_network_identity"]
            ),
            source_id=(
                f"{typed.identity.source_network_id}:"
                f"{typed.identity.selection_id}"
            ),
            source_revision="older-network-revision",
        )
        for typed in result.typed_results
    }
    cached_result = AnalysisResult(source=EDTA_PDB, name="EDTA stale cache")
    cached = attach_group_view_results(
        cached_result,
        network,
        tuple(
            SimpleNamespace(name=group_id, group_id=group_id, ball_ndxs=ids)
            for group_id, ids in selections.items()
        ),
        system="EDTA",
        radius_configuration_id="vorpy-backend-radii-v1",
        source_network_id=source_network_id,
        cached_dependencies=stale_dependencies,
    )
    assert all(typed.status is Status.STALE for typed in cached)

    archive = save_network(network, tmp_path / "edta-retained-mg.vpy")
    restored = load_network(archive)
    assert restored.source_network_id == source_network_id
    restored_result = AnalysisResult(source=EDTA_PDB, name="EDTA retained reload")
    restored_typed = attach_group_view_results(
        restored_result,
        restored,
        groups,
        system="EDTA",
        radius_configuration_id="vorpy-backend-radii-v1",
    )
    assert restored.retained_generator_provenance == network.retained_generator_provenance
    for before, after in zip(result.typed_results, restored_typed, strict=True):
        assert after.status is Status.CERTIFIED
        assert after.identity.source_network_id == before.identity.source_network_id
        assert after.identity.source_network_id != restored.archive_id
        assert after.identity.solve_universe_id == before.identity.solve_universe_id
        assert after.identity.selection_id == before.identity.selection_id
        assert after.provenance.details["dependency_digest"] == (
            before.provenance.details["dependency_digest"]
        )
        assert after.metrics["area"].value == pytest.approx(
            before.metrics["area"].value, abs=1e-12
        )
        assert after.metrics["volume"].value == pytest.approx(
            before.metrics["volume"].value, abs=1e-12
        )
    assert len(build_calls) == 1


def test_shared_view_analysis_fails_closed_without_complete_view():
    network = _solved_network()
    view = NetworkView.for_group(network, [10])
    view = replace(view, completeness=replace(view.completeness, proven=False))

    with pytest.raises(NetworkViewError, match="completeness"):
        analyze_group_view(view)
    with pytest.raises(NetworkViewError, match="authoritative F10 solvent provenance"):
        shared_network_identity(NetworkView.for_group(network, [10]), network)
    with pytest.raises(NetworkViewError, match="authoritative F10 solvent provenance"):
        attach_group_view_results(
            SimpleNamespace(),
            network,
            [SimpleNamespace(ball_ndxs=[10])],
            system="synthetic",
            radius_configuration_id="synthetic",
        )


def test_shared_view_identity_fails_closed_without_f12_provenance(monkeypatch):
    network = _solved_network()
    network.solvent_mode = SolventMode.MOLECULE_ONLY.value
    network.refresh_solvent_provenance()
    view = NetworkView.for_group(network, [10])
    monkeypatch.setattr(
        Network,
        "retained_generator_provenance",
        property(lambda _network: None),
    )

    with pytest.raises(NetworkViewError, match="F12 retention provenance"):
        shared_network_identity(view, network)


@pytest.mark.skipif(
    not CAMBRIN_ARCHIVE.is_file(),
    reason="local F8 Cambrin archive is not available",
)
def test_f8_cambrin_group_view_reports_competitor_completeness_blocker(monkeypatch):
    """Use the cached molecular solve; do not spend time rebuilding Cambrin."""
    monkeypatch.setattr(
        Network,
        "build",
        lambda *args, **kwargs: pytest.fail("Cambrin regression rebuilt Network"),
    )
    network = load_network(CAMBRIN_ARCHIVE)
    balls_before = network.balls.copy(deep=True)
    surfaces_before = network.surfs.copy(deep=True)

    assert network.settings["net_type"] == "aw"
    assert len(network.group) == 674
    assert network.boundary_indices == ()
    assert network.balls["complete"].value_counts().to_dict() == {
        0: 8516,
        True: 674,
    }

    context = network.get_solve_context_key()
    system_to_local = {
        system_id: local_id
        for local_id, system_id in context.local_to_system
        if system_id is not None
    }
    selected_local = system_to_local[2]
    assert selected_local == 2
    reference_volume = float(network.balls.iloc[selected_local]["vol"])
    surface_ids = [
        int(position)
        for position, surface in network.surfs.iterrows()
        if selected_local in [int(value) for value in surface["balls"]]
    ]
    reference_area = sum(
        float(network.surfs.iloc[position]["sa"])
        for position in surface_ids
    )
    assert reference_volume == pytest.approx(13.572938281053846, abs=1e-12)
    assert reference_area == pytest.approx(38.701579174896224, abs=1e-12)
    assert len(surface_ids) == 21

    with pytest.raises(NetworkViewError, match="completeness") as error:
        make_group_view(SimpleNamespace(ball_ndxs=[2]), network)
    assert "2037" in str(error.value)

    pd.testing.assert_frame_equal(network.balls, balls_before)
    pd.testing.assert_frame_equal(network.surfs, surfaces_before)
