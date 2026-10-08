"""Result-only regressions: saved cache snapshots, no scientific computation."""

import importlib
import json
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.results import (
    HUMAN_SUMMARY_LEVELS,
    METRICS_COLUMNS,
    AtomIdentity,
    Environment,
    GeneratorComplex,
    InteractionLinks,
    MolecularContactSelection,
    MolecularContactSurfaceGeometry,
    MolecularContactSurfaceResult,
    MolecularPatchProvenance,
    MolecularSurfaceConvention,
    GroupResult,
    InterfaceResult,
    NetworkResult,
    Provenance,
    Quantity,
    RepresentationExtension,
    ResultIdentity,
    Status,
    SystemResult,
    Topology,
    adapt_atom_metadata,
    adapt_cached_result,
    adapt_contact_mapping,
    adapt_interface_analysis,
    generator_complex_extension,
    unavailable_interface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def identity(partition="power", **changes):
    data = {
        "result_kind": "interface",
        "system": "2KAI",
        "system_id": "2KAI:H-tier",
        "frame": "static",
        "network_id": "2KAI:matched-probe:" + partition,
        "network_scope": "system",
        "environment": "dry",
        "partition": partition,
        "representation": "voronoi",
        "radius_configuration_id": "2KAI:H-tier:"
        + ("expanded-r+1.4" if partition == "power" else "base-r"),
        "interface_id": "AB_I",
        "group_A": "chains:A+B",
        "group_B": "chain:I",
        "side": "shared",
        "orientation": "group1 -> group2",
        "alpha_method": "power_distance" if partition == "power" else "additive",
        "alpha_value": 0.0 if partition == "power" else 1.4,
        "alpha_units": "A^2" if partition == "power" else "A",
        "probe_radius": 1.4,
        "probe_units": "Å",
    }
    data.update(changes)
    return ResultIdentity(**data)


def snapshot(partition="power"):
    # Test fixtures deserialize a frozen cache, then hand a Python object to the
    # adapter. This is not a live exported-directory import path.
    state = json.loads(
        (FIXTURES / f"2kai_{partition}_cached_analysis.json").read_text()
    )
    return InterfaceGeometryAnalysis(**state)


def sample_quantity(name="area", value=1.0, status=Status.CERTIFIED, units="Å²"):
    return Quantity(name, value, units, status, "cached test measurement", "test_scope")


def test_power_exact_intrinsic_and_completed_values():
    result = adapt_interface_analysis(snapshot(), identity())
    assert isinstance(result, InterfaceResult)
    assert [
        result.topology.quantities[k].value for k in ("vertices", "edges", "faces")
    ] == [605, 955, 350]
    assert result.topology.quantities["euler_characteristic"].value == 0
    assert result.selection["selected_features"].value == 350
    h, k = result.mean_curvature, result.gaussian_curvature
    assert h.integrated_mean_curvature.value == pytest.approx(
        -64.45049624200274, abs=1e-12
    )
    assert h.integrated_mean_curvature.status == Status.CERTIFIED
    assert h.complete_total.value == h.integrated_mean_curvature.value
    assert h.smooth_face_contribution.value == 0
    assert h.edge_contribution.value == h.complete_total.value
    assert h.boundary_contribution.value is None
    assert h.boundary_contribution.status == Status.NOT_SUPPORTED
    assert k.integrated_gaussian_curvature.value == pytest.approx(
        5.188337486733097, abs=1e-12
    )
    assert k.integrated_gaussian_curvature.status == Status.CERTIFIED
    assert k.interior_vertex_defects.value == k.integrated_gaussian_curvature.value
    assert k.boundary_corner_contribution.value == pytest.approx(
        -5.1883374867332295, abs=1e-12
    )
    assert k.gauss_bonnet_completion.value == pytest.approx(0, abs=2e-12)
    assert k.gauss_bonnet_completion.status == Status.CERTIFIED
    assert k.gauss_bonnet_residual.value == pytest.approx(
        -1.3233858453531866e-13, abs=2e-13
    )
    assert k.integrated_gaussian_curvature.value != k.gauss_bonnet_completion.value
    assert result.topology.genus_by_component == (0, 0)


@pytest.mark.parametrize(
    (
        "partition",
        "expected_area",
        "mean_status",
        "smooth_gaussian_status",
        "intrinsic_gaussian_status",
        "completed_gaussian_status",
        "expected_euler",
    ),
    [
        (
            "power",
            898.1185793538328,
            Status.CERTIFIED,
            Status.CERTIFIED,
            Status.CERTIFIED,
            Status.CERTIFIED,
            0,
        ),
        (
            "aw",
            787.0925777403216,
            Status.PARTIAL,
            Status.PARTIAL,
            Status.UNRESOLVED,
            Status.UNRESOLVED,
            -1,
        ),
    ],
)
def test_physical_interface_quantities_keep_units_status_and_provenance(
    partition,
    expected_area,
    mean_status,
    smooth_gaussian_status,
    intrinsic_gaussian_status,
    completed_gaussian_status,
    expected_euler,
):
    provenance = Provenance(f"cached {partition} physical interface")
    source = snapshot(partition)
    result = adapt_interface_analysis(
        source, identity(partition), provenance=provenance
    )

    area = result.metrics["area"]
    mean = result.mean_curvature.integrated_mean_curvature
    smooth_gaussian = result.gaussian_curvature.smooth_face_contribution
    intrinsic_gaussian = result.gaussian_curvature.integrated_gaussian_curvature
    completed_gaussian = result.gaussian_curvature.gauss_bonnet_completion
    euler = result.topology.quantities["euler_characteristic"]

    assert result.identity.partition == partition
    assert area.value == pytest.approx(expected_area)
    assert area.units == "Å²"
    assert area.status is Status.CERTIFIED
    assert mean.units == "Å" and mean.status is mean_status
    assert smooth_gaussian.units == "1"
    assert intrinsic_gaussian.units == "1"
    assert completed_gaussian.units == "1"
    assert smooth_gaussian.status is smooth_gaussian_status
    assert intrinsic_gaussian.status is intrinsic_gaussian_status
    assert completed_gaussian.status is completed_gaussian_status
    assert euler.value == expected_euler and euler.units == "1"
    assert all(
        quantity.provenance is provenance
        for quantity in (
            area,
            mean,
            smooth_gaussian,
            intrinsic_gaussian,
            completed_gaussian,
            euler,
        )
    )
    assert result.metrics["n_contacts"].value is None
    assert result.metrics["n_contacts"].status is Status.NOT_CALCULATED


def test_aw_partial_is_not_full_curvature_or_zero():
    result = adapt_interface_analysis(snapshot("aw"), identity("aw"))
    h, k = result.mean_curvature, result.gaussian_curvature
    assert result.selection["selected_features"].value == 339
    assert h.integrated_mean_curvature.value == pytest.approx(
        -54.150160911300546, abs=1e-10
    )
    assert h.integrated_mean_curvature.status == Status.PARTIAL
    assert h.integrated_mean_curvature.scope == "supported_contributions"
    assert (
        h.complete_total.value is None and h.complete_total.status == Status.UNRESOLVED
    )
    assert (
        h.smooth_face_contribution.support_count == 334
        and h.smooth_face_contribution.total_count == 339
    )
    assert result.metrics["mean_unsupported_faces"].value == 5
    assert result.metrics["mean_unresolved_internal_seams"].value == 14
    assert k.integrated_gaussian_curvature.value is None
    assert k.integrated_gaussian_curvature.status == Status.UNRESOLVED
    assert k.integrated_gaussian_curvature.scope == "full_selected_interface"
    assert k.smooth_face_contribution.scope == "supported_subcomplex"
    assert k.gauss_bonnet_completion.value is None
    assert k.gauss_bonnet_residual.value is None
    assert k.gauss_bonnet_expected.scope == "full_selected_interface"
    subset = result.extensions[0].payload["physical_accounting"][
        "supported_subcomplex_gauss_bonnet"
    ]
    assert subset["residual"] == pytest.approx(-0.0006654952459843599, abs=1e-10)


def test_missing_ratios_counts_and_discrete_aggregate_stay_missing():
    result = adapt_interface_analysis(snapshot(), identity())
    assert result.metrics["mean_curvature_per_area"].value is None
    assert result.metrics["mean_curvature_per_area"].status == Status.NOT_CALCULATED
    assert result.metrics["n_contacts"].value is None
    assert result.gaussian_curvature.intrinsic_discrete_contribution.value is None
    assert result.topology.quantities["orientable"].value is None
    assert result.topology.quantities["nonmanifold_edges"].value is None


def test_model_snapshots_do_not_mutate_cache_and_are_deeply_read_only():
    source = snapshot("aw")
    before = json.dumps(vars(source), sort_keys=True)
    result = adapt_interface_analysis(source, identity("aw"))
    assert json.dumps(vars(source), sort_keys=True) == before
    source.voronoi_side_1["supported_subcomplex_gauss_bonnet"]["residual"] = 123
    assert (
        result.extensions[0].payload["physical_accounting"][
            "supported_subcomplex_gauss_bonnet"
        ]["residual"]
        != 123
    )
    with pytest.raises(TypeError):
        result.metrics["area"] = sample_quantity()
    with pytest.raises(FrozenInstanceError):
        result.status = Status.CERTIFIED
    exported = result.to_dict()
    exported["extensions"][0]["payload"]["physical_accounting"]["area"] = 1
    assert result.metrics["area"].value != 1
    assert json.loads(json.dumps(result.to_dict(), allow_nan=False)) == result.to_dict()


@pytest.mark.parametrize(
    "field,value",
    [
        ("frame", 2),
        ("network_scope", "interface"),
        ("network_id", "another-network"),
        ("alpha_value", 0.5),
        ("alpha_method", "another-convention"),
        ("probe_radius", 1.5),
        ("radius_configuration_id", "another-radii"),
        ("orientation", "group2 -> group1"),
        ("environment", "solvent_competing"),
        ("interface_id", "different-interface"),
        ("system_id", "another-system"),
    ],
)
def test_result_id_covers_scientific_identity(field, value):
    base = identity()
    assert replace(base, **{field: value}).result_id != base.result_id
    assert identity().result_id == base.result_id


def test_identity_and_scope_validation():
    with pytest.raises(ValueError):
        identity(network_scope="full")
    with pytest.raises(ValueError):
        identity(radius_configuration_id="")
    with pytest.raises(ValueError):
        identity(representation="dual_a", side="shared")
    with pytest.raises(ValueError):
        identity(alpha_units=None)
    assert identity(alpha_value=0).result_id == identity(alpha_value=0.0).result_id
    assert identity(alpha_value=-0.0).result_id == identity(alpha_value=0.0).result_id
    with pytest.raises(ValueError):
        identity(alpha_value=True)


def test_environment_and_scope_normalization():
    assert Environment("solvent_competing", True, 100, 5).environment == "solvent_competing"
    assert identity(network_scope="system network").network_scope == "system"
    assert identity(network_scope="interface/dedicated").network_scope == "interface"
    with pytest.raises(ValueError):
        identity(environment="wet")
    with pytest.raises(ValueError):
        identity(network_scope="full")


@pytest.mark.parametrize("status", list(Status))
def test_all_quantity_statuses(status):
    value = 0.0 if status in {Status.CERTIFIED, Status.PARTIAL, Status.STALE} else None
    q = sample_quantity(value=value, status=status)
    assert q.status == status
    if status in {Status.UNRESOLVED, Status.NOT_CALCULATED, Status.NOT_SUPPORTED}:
        with pytest.raises(ValueError):
            sample_quantity(value=0.0, status=status)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"value": float("nan")},
        {"value": float("inf")},
        {"support_count": 3, "total_count": 2},
        {"support_count": -1},
        {"units": ""},
        {"status": "INVENTED"},
        {"value": None, "status": Status.CERTIFIED},
    ],
)
def test_invalid_quantities_rejected(kwargs):
    base = {
        "name": "area",
        "value": 1.0,
        "units": "Å²",
        "status": Status.CERTIFIED,
        "method": "cached",
        "scope": "full",
    }
    base.update(kwargs)
    with pytest.raises(ValueError):
        Quantity(**base)


@pytest.mark.parametrize(
    "kind,cls",
    [("system", SystemResult), ("network", NetworkResult), ("group", GroupResult)],
)
def test_other_result_types_compose_shared_blocks(kind, cls):
    ident = identity(
        result_kind=kind,
        representation="molecular" if kind == "system" else "voronoi",
        group_id="group_A" if kind == "group" else None,
    )
    result = adapt_cached_result(
        ident,
        provenance=Provenance("cached producer"),
        metrics={"area": sample_quantity()},
    )
    assert isinstance(result, cls)
    assert result.to_dict()["identity"]["result_kind"] == kind


@pytest.mark.parametrize("representation,side", [("dual_a", "A"), ("dual_b", "B")])
@pytest.mark.parametrize("partition", ["aw", "power"])
@pytest.mark.parametrize("environment", ["dry", "solvent_competing"])
def test_dual_extension_without_fabricated_geometry(
    representation, side, partition, environment
):
    ident = identity(
        partition,
        representation=representation,
        side=side,
        orientation="producer-defined",
        environment=environment,
    )
    result = unavailable_interface(ident)
    assert result.status == Status.NOT_CALCULATED
    assert result.metrics["area"].value is None
    assert result.topology.quantities["faces"].value is None
    assert (
        unavailable_interface(ident, status=Status.NOT_SUPPORTED).status
        == Status.NOT_SUPPORTED
    )
    with pytest.raises(ValueError):
        adapt_interface_analysis(snapshot(partition), ident)


def test_generator_complex_extension_has_no_surface_requirements():
    def count(dimension, value):
        return Quantity(
            f"simplex_count_dim_{dimension}",
            value,
            "1",
            Status.CERTIFIED,
            "cached generator complex",
            "generator_complex",
        )

    complex_ = GeneratorComplex(
        "AB/contact",
        simplex_counts={0: count(0, 4), 1: count(1, 6), 3: count(3, 1)},
        topology={
            "connected_components": Quantity(
                "connected_components", 1, "1", Status.CERTIFIED,
                "cached generator complex", "generator_complex"
            ),
            "euler_characteristic": Quantity(
                "euler_characteristic", -1, "1", Status.CERTIFIED,
                "cached generator complex", "generator_complex"
            ),
        },
        status=Status.CERTIFIED,
        provenance=Provenance("cached generator complex"),
    )
    extension = generator_complex_extension(complex_)
    ident = identity(
        representation="dual_generator_complex",
        side="AB/contact",
        orientation="producer-defined",
    )
    result = adapt_cached_result(
        ident,
        provenance=Provenance("cached generator complex"),
        status=Status.PARTIAL,
        extensions=(extension,),
    )
    payload = result.to_dict()["extensions"][0]["payload"]
    assert result.metrics == {}
    assert result.topology is None
    assert result.mean_curvature is None
    assert result.gaussian_curvature is None
    assert payload["side"] == "AB/contact"
    assert payload["simplex_counts"]["3"]["value"] == 1
    assert payload["topology"]["euler_characteristic"]["value"] == -1


@pytest.mark.parametrize("side", ["A", "B", "AB/contact"])
def test_generator_complex_sides_are_identity_only(side):
    ident = identity(
        representation="dual_generator_complex",
        side=side,
        orientation="producer-defined",
    )
    assert ident.side == side


def test_dual_separator_is_reserved_without_geometry_assumptions():
    ident = identity(
        representation="dual_separator",
        side="AB/contact",
        orientation="producer-defined",
    )
    assert ident.representation == "dual_separator"


def molecular_identity(side="A", **changes):
    data = {
        "result_kind": "molecular_contact_surface",
        "representation": "molecular_contact_surface",
        "side": side,
        "orientation": "group1 -> group2",
        "partner_id": "chain:I" if side == "A" else "chains:A+B",
        "contact_selection_source": "power_alpha0_interface",
        "contact_selection_id": "2KAI:power:alpha0:AB_I",
        "molecular_surface_convention": "expanded_union_of_balls",
        "definition_version": "mcs-v1",
        "source_network_id": identity().network_id,
        "source_interface_id": "AB_I",
        "interaction_id": "2KAI:AB",
        "orientation_convention": "outward_molecular_normal",
    }
    data.update(changes)
    return identity(**data)


def test_molecular_surface_identity_and_side_reversal():
    a = molecular_identity("A")
    b = molecular_identity("B")
    assert a.molecular_surface_convention == MolecularSurfaceConvention.EXPANDED_UNION_OF_BALLS.value
    assert a.side == "A" and b.side == "B"
    assert a.result_id != b.result_id
    assert molecular_identity("A").result_id == a.result_id
    with pytest.raises(ValueError):
        molecular_identity("AB/contact")


def test_molecular_surface_provenance_and_partial_geometry():
    area = Quantity("area", 12.5, "Ã…Â²", Status.PARTIAL, "cached patch sum", "supported_patches")
    boundary_loops = Quantity(
        "boundary_loops", 0, "1", Status.CERTIFIED, "cached incidence", "selected_surface"
    )
    selection = MolecularContactSelection(
        quantities={
            "selected_ab_contacts": Quantity(
                "selected_ab_contacts", 2, "1", Status.CERTIFIED, "cached selection", "selection"
            ),
            "unresolved_contacts": Quantity(
                "unresolved_contacts", None, "1", Status.UNRESOLVED, "cached selection", "selection"
            ),
        },
        selected_ab_contact_ids=("contact-1", "contact-2"),
        direct_contact_atom_ids=("atom-A", "atom-B"),
        unresolved_contact_ids=("contact-3",),
        status=Status.PARTIAL,
    )
    geometry = MolecularContactSurfaceGeometry(
        quantities={"area": area, "boundary_loops": boundary_loops},
        status=Status.PARTIAL,
    )
    patch = MolecularPatchProvenance(
        "patch-1",
        source_atom_stable_id="atom-A",
        source_sphere_center=(1.0, 2.0, 3.0),
        source_radius=1.7,
        source_radius_units="Ã…",
        partner_contact_atom_ids=("atom-B",),
        parent_selected_contact_ids=("contact-1",),
        clipping_occlusion_atom_ids=("atom-A2",),
        component_id=0,
        boundary_arc_ids=("arc-1",),
        junction_ids=("junction-1",),
        status=Status.PARTIAL,
    )
    result = adapt_cached_result(
        molecular_identity(),
        provenance=Provenance("cached molecular contact surface"),
        status=Status.PARTIAL,
        molecular_selection=selection,
        geometry=geometry,
        patch_provenance=(patch,),
    )
    assert isinstance(result, MolecularContactSurfaceResult)
    assert result.geometry.quantities["area"].value == 12.5
    assert result.geometry.quantities["boundary_loops"].value == 0
    assert result.patch_provenance[0].parent_selected_contact_ids == ("contact-1",)
    assert result.molecular_selection.unresolved_contact_ids == ("contact-3",)


def test_molecular_surface_canonical_area_is_a_result_metric():
    provenance = Provenance("molecular surface area producer")
    result = adapt_cached_result(
        molecular_identity(),
        provenance=provenance,
        status=Status.CERTIFIED,
        metrics={
            "area": Quantity(
                "area", 12.5, "Å²", Status.CERTIFIED,
                "spherical boundary integral", "molecular_contact_surface",
                provenance=provenance,
            )
        },
    )

    assert result.metrics["area"].value == pytest.approx(12.5)
    assert result.metrics["area"].units == "Å²"
    assert result.metrics["area"].provenance is provenance


def test_molecular_contact_counts_are_typed_and_provenance_bearing():
    provenance = Provenance("molecular contact selection")
    selection = MolecularContactSelection(
        quantities={
            "selected_contacts": Quantity(
                "selected_contacts", 3, "1", Status.CERTIFIED,
                "selected-contact cache", "contact_selection",
                provenance=provenance,
            ),
            "unresolved_contacts": Quantity(
                "unresolved_contacts", 1, "1", Status.PARTIAL,
                "unresolved-contact audit", "contact_selection",
                provenance=provenance,
            ),
        },
        selected_ab_contact_ids=("contact-1", "contact-2", "contact-3"),
        unresolved_contact_ids=("contact-4",),
        status=Status.PARTIAL,
    )
    result = adapt_cached_result(
        molecular_identity(),
        provenance=provenance,
        status=Status.PARTIAL,
        molecular_selection=selection,
    )

    selected = result.molecular_selection.quantities["selected_contacts"]
    unresolved = result.molecular_selection.quantities["unresolved_contacts"]
    assert selected.value == 3 and selected.units == "1"
    assert unresolved.value == 1 and unresolved.units == "1"
    assert selected.status is Status.CERTIFIED
    assert unresolved.status is Status.PARTIAL
    assert selected.provenance is provenance
    assert unresolved.provenance is provenance


def test_molecular_surface_unsupported_convention_keeps_values_missing():
    result = adapt_cached_result(
        molecular_identity(molecular_surface_convention="solvent_excluded"),
        provenance=Provenance("unsupported molecular surface"),
        status=Status.NOT_SUPPORTED,
    )
    assert result.status == Status.NOT_SUPPORTED
    assert result.geometry is None
    assert result.mean_curvature is None
    assert result.gaussian_curvature is None


def test_molecular_surface_stale_propagates_without_erasing_values():
    area = Quantity("area", 3.0, "Ã…Â²", Status.CERTIFIED, "cached", "surface")
    result = adapt_cached_result(
        molecular_identity(),
        provenance=Provenance("stale molecular surface"),
        status=Status.STALE,
        geometry=MolecularContactSurfaceGeometry(
            quantities={"area": area}, status=Status.CERTIFIED
        ),
    )
    assert result.geometry.quantities["area"].status == Status.STALE
    assert result.geometry.quantities["area"].value == 3.0


def test_molecular_surface_links_stay_separate_from_other_representations():
    molecular = molecular_identity()
    links = InteractionLinks(
        molecular.interaction_id,
        result_ids={
            "molecular_contact_surface:A": molecular.result_id,
            "physical_interface": "result-v1:physical",
            "molecular_contact_surface:B": "result-v1:molecular-B",
        },
        quantities={
            "area_ratio": Quantity(
                "area_ratio", None, "1", Status.NOT_CALCULATED, "not implemented", "comparison"
            )
        },
    )
    result = adapt_cached_result(
        molecular,
        provenance=Provenance("linked molecular surface"),
        interaction_links=links,
    )
    assert not isinstance(result, InterfaceResult)
    assert not isinstance(result, GeneratorComplex)
    assert result.interaction_links.result_ids["physical_interface"] == "result-v1:physical"
    assert result.interaction_links.quantities["area_ratio"].value is None


def test_genus_requires_justification():
    with pytest.raises(ValueError):
        Topology(genus_by_component=(0,))
    source = snapshot()
    source.voronoi_side_1["manifold_incidence_certified"] = False
    assert (
        adapt_interface_analysis(source, identity()).topology.genus_by_component is None
    )


def atom(atom_name="CA", **changes):
    data = {
        "system_id": "2KAI:H-tier",
        "chain": "A",
        "residue_number": "1",
        "insertion_code": "",
        "residue_name": "GLY",
        "atom_name": atom_name,
    }
    data.update(changes)
    return AtomIdentity(**data)


def test_stable_atom_identity_excludes_generator_and_frame():
    original = atom()
    remapped = replace(original, generator_index=99, network_id="other-network")
    assert original.stable_atom_id == remapped.stable_atom_id
    assert (
        replace(original, insertion_code="A").stable_atom_id != original.stable_atom_id
    )
    assert (
        replace(original, alternate_location="B").stable_atom_id
        != original.stable_atom_id
    )
    assert (
        replace(original, system_id="other-system").stable_atom_id
        != original.stable_atom_id
    )
    assert (
        replace(original, system_atom_id="atom-1").stable_atom_id
        == replace(original, system_atom_id="atom-1", atom_name="CB").stable_atom_id
    )
    with pytest.raises(ValueError):
        AtomIdentity("system", "", "", "", "", "", generator_index=1, network_id="net")


def test_atom_metadata_and_contact_mapping():
    a = adapt_atom_metadata(
        {
            "chain_name": "A",
            "res_seq": 1,
            "res_name": "GLY",
            "name": "CA",
            "stable_id": "A|1|GLY|CA",
            "system_num": 42,
        },
        system_id="2KAI:H-tier",
        generator_index=10,
        network_id=identity().network_id,
    )
    b = atom("CB", chain="I")
    ident = identity()
    record = adapt_contact_mapping(
        {
            "restriction_status": "selected",
            "dual_feature_id": "1:10,13",
            "surface_id": 20,
        },
        ident,
        a,
        b,
    )
    assert record.selection_status == "SELECTED"
    primal_mapping = adapt_contact_mapping(
        {"selection_status": "SELECTED", "dual_edge_id": "1:10,13", "surface_id": 20},
        ident,
        a,
        b,
    )
    assert primal_mapping.dual_feature_id == record.dual_feature_id
    assert record.quantities == {}  # No area or curvature attribution invented.
    assert record.to_dict()["atom_A"]["stable_atom_id"] == a.stable_atom_id
    assert record.to_dict()["identity"]["environment"] == "dry"
    assert (
        replace(record, atom_A=replace(a, generator_index=999)).contact_id
        == record.contact_id
    )
    assert replace(record, physical_feature_id=21).contact_id != record.contact_id
    ident_next = replace(ident, frame=1)
    next_record = replace(record, identity=ident_next, result_id=ident_next.result_id)
    assert next_record.contact_id != record.contact_id
    assert next_record.molecular_pair_id == record.molecular_pair_id
    with pytest.raises(ValueError):
        replace(record, atom_A=replace(a, system_id="different-system"))
    with pytest.raises(ValueError):
        adapt_interface_analysis(snapshot(), ident, contacts=(record, record))


def test_environment_comes_from_competing_generators_not_input_solvent():
    dry = Environment("dry", explicit_solvent_present=True, solvent_generator_count=0)
    assert dry.environment == "dry"
    assert Environment("solvent_competing", True, 100, 5).solvent_interface_competitors == 5
    with pytest.raises(ValueError):
        Environment("dry", True, 100)
    with pytest.raises(ValueError):
        Environment("dry", solvent_interface_competitors=1)
    source = snapshot()
    source.metadata.update(solvent_detected=True, solvent_competing_generator_count=0)
    assert adapt_interface_analysis(source, identity()).environment.environment == "dry"
    with pytest.raises(ValueError):
        adapt_interface_analysis(source, identity(environment="solvent_competing"))


def _explicit_identity_metadata(ident):
    return {
        "system_id": ident.system_id,
        "system": ident.system,
        "network_id": ident.network_id,
        "network_scope": "system network",
        "result_interface_id": ident.interface_id,
        "group_A": ident.group_A,
        "group_B": ident.group_B,
        "frame": ident.frame,
        "radius_configuration_id": ident.radius_configuration_id,
        "alpha_method": ident.alpha_method,
        "alpha_value": ident.alpha_value,
        "alpha_units": ident.alpha_units,
        "probe_radius": ident.probe_radius,
        "probe_units": ident.probe_units,
    }


def test_cached_identity_match_and_old_metadata_gap():
    ident = identity()
    source = snapshot()
    source.metadata.update(_explicit_identity_metadata(ident))
    result = adapt_interface_analysis(source, ident)
    assert result.identity.network_scope == "system"

    # Existing caches omit these identity fields. They remain usable, but those
    # fields are explicitly unavailable to validate.
    assert adapt_interface_analysis(snapshot(), ident).identity == ident


def test_unknown_cached_scope_is_rejected():
    ident = identity()
    source = snapshot()
    source.metadata.update(_explicit_identity_metadata(ident))
    source.metadata["network_scope"] = "whole universe"
    with pytest.raises(ValueError, match="Cached network scope is invalid"):
        adapt_interface_analysis(source, ident)


@pytest.mark.parametrize(
    "field,value",
    [
        ("system_id", "other-system"),
        ("network_id", "other-network"),
        ("network_scope", "interface/dedicated"),
        ("result_interface_id", "other-interface"),
        ("group_A", "other-group"),
        ("frame", 1),
        ("radius_configuration_id", "other-radii"),
        ("alpha_method", "other-alpha"),
        ("alpha_value", 0.5),
        ("alpha_units", "A"),
        ("probe_radius", 2.0),
        ("probe_units", "nm"),
    ],
)
def test_cached_identity_mismatch_is_rejected(field, value):
    ident = identity()
    source = snapshot()
    source.metadata.update(_explicit_identity_metadata(ident))
    source.metadata[field] = value
    with pytest.raises(ValueError, match="Cached identity mismatch"):
        adapt_interface_analysis(source, ident)


def test_equivalent_cached_state_reconstructs_deterministically():
    ident = identity()
    provenance = Provenance("equivalent cached state")
    first = adapt_interface_analysis(snapshot(), ident, provenance=provenance)
    second = adapt_interface_analysis(snapshot(), ident, provenance=provenance)
    assert first.result_id == second.result_id
    assert first.to_dict() == second.to_dict()


def test_stale_invalidates_all_certification_without_erasing_values():
    source = snapshot()
    source.metadata["result_status"] = "STALE"
    result = adapt_interface_analysis(source, identity())
    assert result.status == Status.STALE
    assert result.metrics["area"].status == Status.STALE
    assert result.mean_curvature.complete_total.status == Status.STALE
    assert result.mean_curvature.complete_total.value == pytest.approx(
        -64.45049624200274
    )
    assert result.gaussian_curvature.gauss_bonnet_completion.status == Status.STALE
    assert result.topology.orientation_status == Status.STALE


def test_geometry_profile_marks_skipped_curvature_and_topology_not_calculated():
    source = snapshot()
    source.metadata["curvature_requested"] = False
    source.voronoi_side_1["manifold_incidence_certified"] = False
    source.voronoi_side_2["manifold_incidence_certified"] = False
    for key in (
        "total_integrated_mean_curvature",
        "total_mean_curvature",
        "total_H",
        "partial_integrated_mean_curvature",
        "partial_H",
        "smooth_integrated_mean_curvature",
        "smooth_mean_curvature",
        "smooth_H_partial",
        "smooth_H",
        "internal_seam_integrated_mean_curvature",
        "internal_crease_mean_curvature",
        "edge_H",
        "intrinsic_integrated_gaussian_curvature",
        "total_gaussian_curvature",
        "total_K",
        "smooth_integrated_gaussian_curvature",
        "smooth_gaussian_curvature",
        "smooth_K_partial",
        "intrinsic_discrete_gaussian_curvature",
        "intrinsic_seam_gaussian_curvature",
        "interior_vertex_gaussian_curvature",
        "intrinsic_vertex_gaussian",
        "boundary_geodesic_curvature",
        "boundary_geodesic_gaussian",
        "boundary_corner_turning",
        "boundary_corner_gaussian",
        "gauss_bonnet_lhs",
        "gauss_bonnet_expected",
        "gauss_bonnet_residual",
        "topology_vertices",
        "topology_edges",
        "topology_faces",
        "euler_characteristic",
        "connected_components",
        "components",
        "boundary_loops",
        "boundary_edges",
        "nonmanifold_edges",
        "orientable",
    ):
        source.voronoi_side_1[key] = None
        source.voronoi_side_2[key] = None

    provenance = Provenance("geometry profile cache")
    result = adapt_interface_analysis(source, identity(), provenance=provenance)

    assert result.metrics["area"].value is not None
    assert result.metrics["area"].provenance is provenance
    assert result.mean_curvature.integrated_mean_curvature.status is Status.NOT_CALCULATED
    assert result.mean_curvature.complete_total.status is Status.NOT_CALCULATED
    assert result.mean_curvature.smooth_face_contribution.status is Status.NOT_CALCULATED
    assert result.mean_curvature.edge_contribution.status is Status.NOT_CALCULATED
    assert result.mean_curvature.boundary_contribution.status is Status.NOT_SUPPORTED
    assert result.gaussian_curvature.integrated_gaussian_curvature.status is Status.NOT_CALCULATED
    assert result.gaussian_curvature.smooth_face_contribution.status is Status.NOT_CALCULATED
    assert result.gaussian_curvature.gauss_bonnet_completion.status is Status.NOT_CALCULATED
    assert result.topology.quantities["vertices"].status is Status.NOT_CALCULATED
    assert result.topology.orientation_status is Status.NOT_CALCULATED
    assert result.mean_curvature.complete_total.provenance is provenance
    assert result.gaussian_curvature.gauss_bonnet_completion.provenance is provenance


def test_requested_null_curvature_remains_unresolved():
    source = snapshot()
    source.metadata["curvature_requested"] = True
    for key in (
        "total_integrated_mean_curvature",
        "total_mean_curvature",
        "total_H",
        "partial_integrated_mean_curvature",
        "partial_H",
    ):
        source.voronoi_side_1[key] = None

    result = adapt_interface_analysis(source, identity())

    assert result.mean_curvature.complete_total.value is None
    assert result.mean_curvature.complete_total.status is Status.UNRESOLVED
    assert result.mean_curvature.integrated_mean_curvature.value is None
    assert result.mean_curvature.integrated_mean_curvature.status is Status.UNRESOLVED


def test_extension_and_dict_validation():
    payload = {"mesh_reference": "cache:triangles", "simplex_counts": {"2": 3}}
    ext = RepresentationExtension(
        "vorpy.dual_side", "1", Status.NOT_CALCULATED, payload
    )
    payload["simplex_counts"]["2"] = 9
    assert ext.payload["simplex_counts"]["2"] == 3
    with pytest.raises(ValueError):
        adapt_interface_analysis(snapshot(), identity(), extensions=(ext, ext))
    with pytest.raises(ValueError):
        adapt_cached_result(
            identity(),
            provenance=Provenance("cache"),
            metrics={"wrong": sample_quantity()},
        )
    with pytest.raises(ValueError):
        adapt_cached_result(
            identity(),
            provenance=Provenance("cache"),
            status=Status.NOT_CALCULATED,
            metrics={"area": sample_quantity()},
        )


def test_table_and_human_contract_are_versioned_without_exporter():
    assert len(METRICS_COLUMNS) == len(set(METRICS_COLUMNS))
    for key in (
        "frame",
        "network_scope",
        "mean_curvature_status",
        "mean_curvature_scope",
        "alpha_units",
        "radius_configuration_id",
    ):
        assert key in METRICS_COLUMNS
    assert tuple(HUMAN_SUMMARY_LEVELS) == ("headline", "support", "provenance")


@pytest.mark.parametrize("partition", ["power", "aw"])
def test_exact_dictionary_documentation_examples(partition):
    source = (
        Path(__file__).resolve().parents[3]
        / "docs/development/result_contract_examples"
        / (partition + ".json")
    )
    expected = json.loads(source.read_text(encoding="utf-8"))
    result = adapt_interface_analysis(
        snapshot(partition),
        identity(partition),
        provenance=Provenance("frozen 2KAI InterfaceGeometryAnalysis regression cache"),
    )
    assert result.to_dict() == expected


def test_constructing_and_serializing_results_performs_no_io_or_science(monkeypatch):
    source, ident = snapshot(), identity()

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "Result construction attempted file IO or scientific computation"
        )

    for module, names in (
        ("vorpy.src.network.network", ("Network",)),
        ("vorpy.src.geometry.duals", ("build_dual",)),
        ("vorpy.src.geometry.filtrations.alpha", ("build_alpha_filtration",)),
        (
            "vorpy.src.interface.geometry_analysis",
            ("analyze_interface_geometry", "summarize_physical_interface"),
        ),
        (
            "vorpy.src.analyze.open_interface_curvature",
            (
                "integrated_mean_curvature",
                "polygon_gaussian_terms",
                "gauss_bonnet_accounting",
            ),
        ),
    ):
        mod = importlib.import_module(module)
        for name in names:
            if name == "Network":
                monkeypatch.setattr(mod.Network, "build", forbidden)
            else:
                monkeypatch.setattr(mod, name, forbidden)
    monkeypatch.setattr("builtins.open", forbidden)
    monkeypatch.setattr(Path, "open", forbidden)
    before = json.dumps(vars(source), sort_keys=True)
    result = adapt_interface_analysis(source, ident)
    assert json.dumps(vars(source), sort_keys=True) == before
    assert result.to_dict()["mean_curvature"]["complete_total"][
        "value"
    ] == pytest.approx(-64.45049624200274)
