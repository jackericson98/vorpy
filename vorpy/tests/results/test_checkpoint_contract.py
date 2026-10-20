from dataclasses import replace

import pytest

from vorpy.src.results import (
    CheckpointAtom,
    CheckpointCacheState,
    CheckpointChain,
    CheckpointContractError,
    CheckpointContainer,
    CheckpointDocument,
    CheckpointFeatureIdentity,
    CheckpointGroup,
    CheckpointInterfaceIdentity,
    CheckpointLoadError,
    CheckpointNetworkState,
    CheckpointPayload,
    CheckpointPolicy,
    CheckpointResultState,
    CheckpointResidue,
    CheckpointSaveError,
    CheckpointSource,
    CheckpointUpdateError,
    GeometryLineage,
    GeometryLineageOrigin,
    GroupCExtensionContract,
    Provenance,
    ProvenanceCatalog,
    Status,
    load_checkpoint_container,
    make_checkpoint_dependency,
    prepare_checkpoint_update,
    save_checkpoint_container,
    validate_checkpoint_dependency,
    validate_group_c_extension,
)


def _dependency():
    return make_checkpoint_dependency(
        source_network_id="network:synthetic",
        solve_universe_id="solve:synthetic",
        source_sha256="a" * 64,
        coordinates=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
        radii=(1.5, 1.7),
        solver_settings={"scheme": "pow", "solvent_policy": "none"},
        generator_ids=("g0", "g1"),
        topology={"cells": ("cell:0", "cell:1")},
    )


def _document():
    dependency = _dependency()
    source = CheckpointSource("synthetic.pdb", b"ATOM synthetic\n")
    atoms = (
        CheckpointAtom(
            "atom:0", 0, "CA", "C", "res:A:1", "GLY", 1, "chain:A",
            (0.0, 0.0, 0.0), 1.5, ("group:A",),
        ),
        CheckpointAtom(
            "atom:1", 1, "CA", "C", "res:B:1", "ALA", 1, "chain:A",
            (1.0, 0.0, 0.0), 1.7, ("group:B",),
        ),
    )
    provenance = Provenance(
        "synthetic-checkpoint",
        scientific_state_id="state:partial",
        details={
            "geometry_lineage": GeometryLineage(
                origin=GeometryLineageOrigin.UNRESOLVED,
                reason="synthetic metadata has no geometry solve",
            ).as_dict(),
            "missing_evidence": ("geometry",),
        },
    )
    catalog = ProvenanceCatalog.from_provenances((provenance,))
    network = CheckpointNetworkState(
        "network:synthetic",
        "network:synthetic",
        "system",
        "pow",
        dependency.content_digest,
        ("g0", "g1"),
        ("checkpoint-feature-v1:" + "1" * 64,),
        Status.PARTIAL,
        {"complete": False, "reason": "synthetic partial state"},
        "network:synthetic",
        {"alpha": 1.4, "net_type": "pow"},
        {"mode": "none", "probe_radius": 0.0},
        ("atom:0", "atom:1"),
    )
    return CheckpointDocument(
        system_id="system:synthetic",
        source=source,
        atoms=atoms,
        residues=(
            CheckpointResidue("res:A:1", "chain:A", "GLY", 1, ("atom:0",)),
            CheckpointResidue("res:B:1", "chain:A", "ALA", 1, ("atom:1",)),
        ),
        chains=(CheckpointChain("chain:A", "A", ("res:A:1", "res:B:1"), ("atom:0", "atom:1")),),
        groups=(
            CheckpointGroup("group:A", "A", ("atom:0",), ("g0",), "network:synthetic"),
            CheckpointGroup("group:B", "B", ("atom:1",), ("g1",), "network:synthetic"),
        ),
        interfaces=(CheckpointInterfaceIdentity("system:synthetic", "group:B", "group:A"),),
        dependencies=(dependency,),
        networks=(network,),
        caches=(
            CheckpointCacheState(
                "cache:analysis", "analysis", dependency.content_digest, Status.PARTIAL,
                "payload:analysis", provenance.content_digest,
            ),
        ),
        results=(
            CheckpointResultState(
                "result:interface", "primary", Status.PARTIAL, dependency.content_digest,
                {"metric": "area", "value": 12.5, "units": "A^2"}, provenance.content_digest,
                "network:synthetic",
            ),
        ),
        provenance=catalog,
    )


def _container():
    document = _document()
    dependency = document.dependencies[0]
    network_payload = CheckpointPayload(
        "network-geometry", b"synthetic-network-geometry", dependency.content_digest, Status.PARTIAL
    )
    cache_payload = CheckpointPayload(
        "analysis-cache", b"synthetic-analysis-state", dependency.content_digest, Status.PARTIAL
    )
    networks = (
        replace(document.networks[0], payload_ref=network_payload.payload_id),
        CheckpointNetworkState(
            "network:group", dependency.source_network_id, "group", "pow",
            dependency.content_digest, ("g0",), status=Status.PARTIAL,
            payload_ref=network_payload.payload_id, settings={"net_type": "pow"},
            solvent_policy={"mode": "none"}, competitor_ids=("atom:0",),
        ),
        CheckpointNetworkState(
            "network:interface", dependency.source_network_id, "interface", "pow",
            dependency.content_digest, ("g0", "g1"), status=Status.PARTIAL,
            payload_ref=network_payload.payload_id, settings={"net_type": "pow"},
            solvent_policy={"mode": "none"}, competitor_ids=("atom:0", "atom:1"),
        ),
    )
    results = (
        replace(document.results[0], network_id="network:synthetic"),
        CheckpointResultState(
            "result:group", "group", Status.PARTIAL, dependency.content_digest,
            {"metric": "area", "value": 4.0, "units": "A^2"},
            document.results[0].provenance_ref, "network:group",
        ),
    )
    return CheckpointContainer(
        replace(
            document,
            networks=networks,
            results=results,
            caches=(replace(document.caches[0], payload_ref=cache_payload.payload_id),),
        ),
        (network_payload, cache_payload),
    )


def test_policy_and_source_identity_are_explicit():
    policy = CheckpointPolicy()
    assert policy.filename("1bna") == "1bna_checkpoint.vpy"
    assert policy.path("results", "1bna").as_posix().endswith("results/1bna_checkpoint.vpy")
    source = CheckpointSource("input.pdb", b"ATOM 1\n")
    encoded = source.as_dict()
    assert encoded["source_id"] == source.source_id
    assert CheckpointSource.from_mapping(encoded) == source
    with pytest.raises(CheckpointContractError):
        policy.filename("../unsafe")


def test_checkpoint_round_trip_is_deterministic_and_self_contained():
    document = _document()
    first = document.serialize()
    second = document.serialize()
    restored = CheckpointDocument.deserialize(first)
    assert first == second
    assert restored.as_dict() == document.as_dict()
    assert restored.source.content == b"ATOM synthetic\n"
    assert restored.results[0].typed_result["value"] == 12.5
    assert restored.results[0].status is Status.PARTIAL
    assert restored.networks[0].settings["net_type"] == "pow"
    assert restored.networks[0].solvent_policy["mode"] == "none"
    assert restored.provenance.resolve(restored.results[0].provenance_ref).geometry_lineage.origin is GeometryLineageOrigin.UNRESOLVED


def test_container_bundle_round_trip_preserves_scoped_payloads_and_results():
    container = _container()
    manifest = container.serialize()
    restored = load_checkpoint_container(manifest, container.payload_members())
    assert manifest == restored.serialize()
    assert restored.content_digest == container.content_digest
    assert restored.document.source.content == b"ATOM synthetic\n"
    assert {network.scope for network in restored.document.networks} == {
        "system", "group", "interface"
    }
    assert [result.typed_result["value"] for result in restored.document.results] == [12.5, 4.0]
    assert container.payload_members() == restored.payload_members()
    bundle = save_checkpoint_container(container)
    assert bundle.load().content_digest == container.content_digest


def test_container_load_and_save_failures_are_explicit():
    container = _container()
    with pytest.raises(CheckpointLoadError):
        load_checkpoint_container(container.serialize(), {})
    corrupt_members = container.payload_members()
    first_member = next(iter(corrupt_members))
    corrupt_members[first_member] = b"corrupt"
    with pytest.raises(CheckpointLoadError):
        load_checkpoint_container(container.serialize(), corrupt_members)
    with pytest.raises(CheckpointSaveError):
        save_checkpoint_container(_document())


def test_atomic_update_plan_validates_before_replace(tmp_path):
    current = _container()
    plan = prepare_checkpoint_update(current, current=current)
    assert plan.validate_current(current) is current
    bundle = plan.save_bundle()
    assert plan.validate_written(bundle).content_digest == current.content_digest

    stale = prepare_checkpoint_update(current, expected_current_digest="checkpoint-container-v1:" + "f" * 64)
    with pytest.raises(CheckpointUpdateError):
        stale.validate_current(current)

    destination = tmp_path / "system_checkpoint.vpy"
    previous_bytes = current.serialize().encode("utf-8")
    destination.write_bytes(previous_bytes)
    temporary = destination.with_name(".checkpoint-failed.tmp")
    try:
        temporary.write_bytes(bundle.manifest_bytes)
        plan.validate_written(bundle)
        raise OSError("synthetic write failure before replace")
    except OSError:
        if temporary.exists():
            temporary.unlink()
    assert destination.read_bytes() == previous_bytes


def test_interface_identity_is_shared_by_ab_and_ba():
    ab = CheckpointInterfaceIdentity("system:synthetic", "group:A", "group:B")
    ba = CheckpointInterfaceIdentity(
        "system:synthetic", "group:B", "group:A", "group:B", "group:A"
    )
    assert ab.interface_id == ba.interface_id
    assert ab.reference_orientation == ba.reference_orientation == ("group:A", "group:B")
    assert ab.display_orientation == ("group:A", "group:B")
    assert ba.display_orientation == ("group:B", "group:A")


@pytest.mark.parametrize(
    "field",
    (
        "source_sha256",
        "coordinates_digest",
        "radii_digest",
        "solver_settings_digest",
        "generator_universe_digest",
        "topology_digest",
    ),
)
def test_incompatible_dependency_inputs_invalidate_without_rebuild(field):
    dependency = _dependency()
    changed = replace(dependency, **{field: "b" * 64})
    assert validate_checkpoint_dependency(dependency, changed) is Status.STALE
    assert validate_checkpoint_dependency(dependency, None) is Status.NOT_CALCULATED


def test_group_c_contract_is_identity_only_until_forge_proves_merge():
    document = _document()
    dependency = document.dependencies[0]
    feature = CheckpointFeatureIdentity(
        "surface", dependency.source_network_id, "surface:C:0", ("g2",), dependency.topology_digest
    )
    extension = GroupCExtensionContract(
        document.system_id,
        dependency.source_network_id,
        "group:C",
        ("atom:0",),
        ("g2",),
        (feature,),
        dependency.content_digest,
    )
    assert extension.merge_supported is False
    assert validate_group_c_extension(document, extension) is Status.CERTIFIED
    extended = replace(document, group_extensions=(extension,))
    assert extended.group_extensions == (extension,)
    assert validate_group_c_extension(document, replace(extension, status=Status.CERTIFIED)) is Status.NOT_SUPPORTED


def test_missing_dependency_and_provenance_are_rejected():
    document = _document()
    with pytest.raises(CheckpointContractError):
        CheckpointDocument(**{
            **document.__dict__,
            "networks": (replace(document.networks[0], dependency_digest="checkpoint-dependency-v1:" + "f" * 64),),
        })
    with pytest.raises(CheckpointContractError):
        CheckpointDocument(**{
            **document.__dict__,
            "results": (replace(document.results[0], provenance_ref="provenance-v1:" + "f" * 64),),
        })


def test_invalid_source_and_dependency_digests_fail_closed():
    source = CheckpointSource("input.pdb", b"ATOM 1\n")
    invalid_source = {**source.as_dict(), "sha256": "0" * 64}
    with pytest.raises(CheckpointContractError):
        CheckpointSource.from_mapping(invalid_source)
    with pytest.raises(CheckpointContractError):
        replace(_dependency(), topology_digest="not-a-digest")


def test_prm_certified_typed_result_is_rejected_but_unresolved_is_preserved():
    document = _document()
    dependency = document.dependencies[0]
    prm_network = CheckpointNetworkState(
        "network:prm", dependency.source_network_id, "system", "prm",
        dependency.content_digest, ("g0", "g1"), status=Status.PARTIAL,
        competitor_ids=("atom:0", "atom:1"),
    )
    certified = CheckpointResultState(
        "result:prm", "primary", Status.CERTIFIED, dependency.content_digest,
        {"metric": "area", "value": 1.0}, network_id="network:prm",
    )
    with pytest.raises(CheckpointContractError):
        replace(document, networks=(document.networks[0], prm_network), results=(certified,))
    unresolved = replace(certified, status=Status.UNRESOLVED)
    restored = replace(document, networks=(document.networks[0], prm_network), results=(unresolved,))
    assert restored.results[0].status is Status.UNRESOLVED
