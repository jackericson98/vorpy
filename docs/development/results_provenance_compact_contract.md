# ATLAS compact Results provenance contract

Status: implementation-ready handoff for VECTOR, 2026-10-19.

## Finding

ATLAS adapters commonly pass the same immutable `Provenance` object to every
`Quantity`. The duplication appears at the canonical-log boundary:
`_compact_result_rows()` places that object on every `InterfaceLogRow`, and the
legacy serializer expands the complete JSON payload independently for every
row. The current path has no digest index or sidecar. Large unresolved simplex
and incidence lists therefore repeat in `logs.csv`; numeric values and geometry
are not the source of this duplication.

## Versioned wire format

The existing `logs.csv` and its eight columns remain unchanged:

```text
interface_id,representation,side,quantity,value,units,status,provenance
```

The new compact format is a separate, explicitly versioned file:

```text
logs.compact.csv
interface_id,representation,side,quantity,value,units,status,provenance_ref
```

Its sibling sidecar is `logs.provenance.json`:

```json
{
  "format": "vorpy.results.provenance",
  "schema_version": "1",
  "entries": {
    "provenance-v1:<sha256>": {
      "record_type": "contract",
      "source_id": "source:1bna",
      "source_revision": "...",
      "source_sha256": "<64 hex characters>",
      "method": "interface-area-v1",
      "scheme": "aw",
      "ordered_selection": ["A", "B"],
      "generator_ids": ["..."],
      "settings": {},
      "dependency_digests": ["..."],
      "status": "PARTIAL",
      "completeness": {},
      "geometry_lineage": {},
      "missing_evidence": ["..."],
      "unresolved_evidence": ["..."],
      "inventory_refs": ["inventory-v1:<sha256>"],
      "evidence": {}
    }
  },
  "inventories": {}
}
```

`provenance_ref` is the complete `provenance-v1:<sha256>` key. The digest is
computed from deterministic UTF-8 JSON with sorted keys, compact separators,
and no non-finite values. An inventory uses the same rule with the
`inventory-v1:` prefix. A catalog stores each unique payload once, so repeated
metric rows carry only the small reference. Inventory payloads are optional,
but a contract that names one must have that inventory in the sidecar.

The explicit contract preserves source ID and SHA-256, method, scheme,
settings, ordered selections, generator IDs, dependency digests, status,
completeness, `GeometryLineage`, missing/unresolved evidence, and arbitrary
producer evidence. Legacy `Provenance` entries are also accepted and stored
verbatim, which keeps old metadata lossless while producers migrate.

The compact typed Results files are additive. Classic scheme geometry logs
remain the authoritative `aw_logs.csv`, `pow_logs.csv`, and `prm_logs.csv`
outputs produced by the existing geometry-log writer; `logs.csv` retains its
legacy eight-column typed Results payload, while `logs.compact.csv` is the
version-2 reference form. Neither typed Results form replaces or reinterprets
the classic physical-network records.

## Reader and failure rules

ATLAS exposes `ProvenanceContract`, `ProvenanceCatalog`,
`ProvenanceReference`, `ProvenanceInventory`, and the compact row helpers in
`vorpy.src.results`. `InterfaceLogRow.as_csv_row()` remains the legacy writer;
`as_compact_csv_row(catalog)` emits the v2 reference column.

Readers must select the format by filename/header. A legacy row is parsed from
its full JSON `provenance` field with `deserialize_interface_log_provenance()`.
A compact row is resolved with `resolve_interface_log_provenance_reference()`.
Missing sidecars or missing referenced entries raise
`ProvenanceSidecarMissingError`; malformed or stale keys/payloads raise
`ProvenanceDigestError`. Readers must not replace either failure with an empty
provenance or a certified status.

Existing `.vpy` archives remain on their current archive schema and are read by
the existing archive loader. This change adds no archive member and does not
change interpretation of any existing column. VECTOR can later emit both
formats during migration; no geometry rebuild or numeric metric conversion is
needed. A compact reader must reject a missing/stale sidecar before exposing
rows as complete or certified; old logs and archives continue through their
existing readers.
