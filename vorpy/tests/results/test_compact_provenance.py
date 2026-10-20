import csv
import io
from types import SimpleNamespace

import pytest

from vorpy.src.results import (
    COMPACT_INTERFACE_LOG_COLUMNS,
    GeometryLineage,
    GeometryLineageOrigin,
    INTERFACE_LOG_QUANTITY_UNITS,
    InterfaceLogRow,
    Provenance,
    ProvenanceCatalog,
    ProvenanceContract,
    ProvenanceDigestError,
    ProvenanceInventory,
    ProvenanceSidecarMissingError,
    Status,
    deserialize_interface_log_provenance,
    serialize_interface_log_provenance,
)


def _lineage():
    return GeometryLineage(
        origin=GeometryLineageOrigin.DIRECT_SOLVE,
        solve_universe_id="solve:1bna",
        geometry_source_id="source:1bna",
        geometry_version="vorpy-aw-v2",
        solver_settings_digest="settings:1bna",
        solver_settings={"scheme": "aw", "alpha": 1.4},
        environment_digest="environment:1bna",
        environment={"solvent_mode": "FULL_SOLVENT"},
        source_network_id="network:1bna",
        source_interface_id="interface:A-B",
        selection_scope={"ordered_groups": ("A", "B")},
        orientation={"convention": "physical_interface_A_to_B"},
        completeness={"completeness_proven": False, "reason": "missing cells"},
    )


def _contract(*, status=Status.PARTIAL, inventory_ref=()):
    return ProvenanceContract(
        source_id="source:1bna",
        source_revision="1bna-db1976",
        source_sha256="a" * 64,
        method="interface-area-v1",
        scheme="aw",
        ordered_selection=("A|15|LYS|NZ", "I|3B|TYR|CA"),
        generator_ids=("atom:A:15:NZ", "atom:I:3B:CA"),
        settings={"alpha": 1.4, "alpha_units": "A", "solvent_mode": "FULL_SOLVENT"},
        dependency_digests=("dependency-v1:source", "dependency-v1:geometry"),
        status=status,
        completeness={"proven": False, "unresolved_simplex_count": 12000},
        geometry_lineage=_lineage(),
        missing_evidence=("simplex:1:1320,4061",),
        unresolved_evidence=("solver cells incomplete",),
        inventory_refs=inventory_ref,
        evidence={
            "unresolved_simplex_ids": tuple(f"1:{index},{index + 1}" for index in range(12000)),
            "source_role": "canonical-1bna-interface",
        },
    )


def _csv_bytes(rows, columns):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)
    return len(stream.getvalue().encode("utf-8"))


def test_synthetic_1bna_metadata_deduplicates_without_changing_metrics():
    inventory = ProvenanceInventory(
        "dual-simplex-incidence",
        {"missing": ["1:1320,4061"], "unresolved": ["2:1,2,3"]},
    )
    contract = _contract(inventory_ref=(inventory.content_digest,))
    small_evidence = {
        **contract.as_dict(),
        "evidence": {
            "unresolved_simplex_ids": tuple(f"1:{index},{index + 1}" for index in range(100)),
            "source_role": "canonical-1bna-interface",
        },
    }
    variants = [
        contract,
        ProvenanceContract.from_mapping(small_evidence),
        ProvenanceContract.from_mapping({**small_evidence, "scheme": "power"}),
    ]
    catalog = ProvenanceCatalog.from_provenances(
        [variants[index % 3] for index in range(106)],
        inventories=(inventory,),
    )
    rows = [
        InterfaceLogRow(
            "A-B",
            "physical_partition_interface",
            "shared",
            "area",
            123.5 + index,
            INTERFACE_LOG_QUANTITY_UNITS["area"],
            Status.CERTIFIED,
            variants[0 if index < 18 else (index % 2) + 1],
        )
        for index in range(106)
    ]
    legacy = _csv_bytes(
        [row.as_csv_row() for row in rows],
        ("interface_id", "representation", "side", "quantity", "value", "units", "status", "provenance"),
    )
    compact = _csv_bytes(
        [row.as_compact_csv_row(catalog) for row in rows],
        COMPACT_INTERFACE_LOG_COLUMNS,
    ) + len(catalog.serialize().encode("utf-8"))

    assert len(catalog.entries) == 3
    assert compact < legacy
    assert compact < legacy // 2
    for legacy_row, compact_row in zip(rows, rows):
        old = legacy_row.as_csv_row()
        new = compact_row.as_compact_csv_row(catalog)
        assert (new["value"], new["units"], new["status"]) == (
            old["value"],
            old["units"],
            old["status"],
        )


def test_catalog_serialization_is_deterministic_and_round_trips_lineage():
    inventory = ProvenanceInventory("simplex", {"ids": ["1:1,2"]})
    contract = _contract(inventory_ref=(inventory.content_digest,))
    catalog = ProvenanceCatalog.from_provenances((contract,), inventories=(inventory,))

    assert catalog.serialize() == catalog.serialize()
    restored = ProvenanceCatalog.deserialize(catalog.serialize())
    resolved = restored.resolve(catalog.reference(contract))

    assert resolved == contract
    assert resolved.geometry_lineage == contract.geometry_lineage
    assert restored.inventories[inventory.content_digest] == inventory


def test_sidecar_missing_and_invalid_digest_fail_closed():
    inventory = ProvenanceInventory("simplex", {"ids": ["1:1,2"]})
    contract = _contract(inventory_ref=(inventory.content_digest,))
    catalog = ProvenanceCatalog.from_provenances((contract,), inventories=(inventory,))
    missing = "provenance-v1:" + "f" * 64

    with pytest.raises(ProvenanceSidecarMissingError):
        catalog.resolve(missing)
    with pytest.raises(ProvenanceDigestError):
        catalog.resolve("bad-digest")
    with pytest.raises(ProvenanceDigestError):
        ProvenanceCatalog(
            entries={contract.content_digest[:-1] + "0": contract},
            inventories={inventory.content_digest: inventory},
        )
    with pytest.raises(ValueError, match="source_sha256"):
        ProvenanceContract(
            source_id="source:1bna",
            method="method",
            scheme="aw",
            source_sha256="not-a-digest",
        )


def test_partial_unresolved_contract_and_legacy_log_remain_readable():
    contract = _contract(status=Status.UNRESOLVED)
    projected = contract.to_provenance()
    assert projected.geometry_lineage == contract.geometry_lineage
    assert projected.details["canonical_provenance"]["missing_evidence"] == (
        "simplex:1:1320,4061",
    )

    legacy = Provenance("legacy", details={"partition": "power"})
    restored = deserialize_interface_log_provenance(serialize_interface_log_provenance(legacy))
    assert isinstance(restored, Provenance)
    assert restored.details["partition"] == "power"
    assert deserialize_interface_log_provenance('{"partition":"power"}') == {
        "partition": "power"
    }
    restored_contract = deserialize_interface_log_provenance(
        serialize_interface_log_provenance(contract)
    )
    assert isinstance(restored_contract, ProvenanceContract)
    assert restored_contract.content_digest == contract.content_digest


def test_export_writes_legacy_and_compact_logs_with_profile(tmp_path, monkeypatch):
    from vorpy.src.output import visualization

    provenance = Provenance(
        "synthetic export",
        details={"unresolved_simplex_ids": tuple(f"1:{i},{i + 1}" for i in range(100))},
    )
    rows = [
        InterfaceLogRow(
            "A-B", "physical_partition_interface", "shared", "area",
            float(index), INTERFACE_LOG_QUANTITY_UNITS["area"], Status.PARTIAL,
            provenance,
        ).as_csv_row()
        for index in range(106)
    ]
    monkeypatch.setattr(
        visualization,
        "_compact_result_rows",
        lambda iface, **kwargs: rows,
    )

    profile = visualization._compact_log_bundle(
        SimpleNamespace(), tmp_path / "logs.csv"
    )

    assert (tmp_path / "logs.csv").is_file()
    assert (tmp_path / "logs.compact.csv").is_file()
    assert (tmp_path / "logs.provenance.json").is_file()
    resolved = visualization.read_compact_log(tmp_path / "logs.compact.csv")
    assert len(resolved) == 106
    assert resolved[0]["provenance"] == provenance
    assert profile["counts"]["unique_provenance"] == 1
    assert profile["counts"]["files_written"] == 3
    assert profile["bytes"]["compact_csv"] < profile["bytes"]["legacy_csv"]
