# Shared self-contained checkpoint contract

Status: ATLAS schema/identity groundwork, aligned with the NEXUSA checkpoint
contract approved for implementation review on 2026-10-19. This document does
not implement archive I/O or molecular geometry.

## Public API

`vorpy.src.results.checkpoint` defines the payload semantics:

* `CheckpointPolicy` fixes one output name, `<system>_checkpoint.vpy`, and
  records the atomic-save policy.
* `CheckpointDocument` is a deterministic JSON-safe manifest with format
  `vorpy.system.checkpoint`, schema version `1`.
* `CheckpointSource` embeds the original source bytes, basename, byte count,
  SHA-256, and derived `source-v1:<sha256>` identity.
* `CheckpointAtom`, `CheckpointResidue`, `CheckpointChain`, and
  `CheckpointGroup` preserve coordinates, radii, membership, generator
  identities, and partial status.
* `CheckpointInterfaceIdentity` treats A/B and B/A as one unordered physical
  identity while retaining a deterministic reference orientation and a
  display orientation.
* `CheckpointDependency` hashes source, coordinates, radii, solver settings,
  generator universe, and topology. `validate_checkpoint_dependency` returns
  `CERTIFIED`, `STALE`, `NOT_CALCULATED`, or `UNRESOLVED`; it never rebuilds a
  network.
* `CheckpointNetworkState`, `CheckpointCacheState`, and
  `CheckpointResultState` retain independent network, cache, typed Result,
  dependency, provenance, completeness, and `PARTIAL`/`UNRESOLVED` states.
  Network state records the complete competitor atom universe explicitly;
  certified Results identify their source network. A certified PRM typed
  Result is rejected until its scientific support is approved; unresolved or
  partial PRM records remain preservable.
* `ProvenanceCatalog` is embedded in the manifest, so typed Results retain
  compact references and full evidence even when standalone logs are off.

`CheckpointDocument.serialize()` is canonical, sorted JSON. `deserialize()`
reconstructs and validates every cross-reference; a missing dependency,
provenance record, digest, or source-byte check fails closed with
`CheckpointContractError`.

## VECTOR-facing container API

`CheckpointContainer(document, payloads)` is the one shared container. The
document holds source bytes, identities, metadata, typed Results, and
provenance. `CheckpointPayload` holds opaque lossless bytes for solved network
geometry, analysis, or cache state; its content-addressed ID and SHA-256 are
validated before use.

The public serialization boundary is:

```python
bundle = save_checkpoint_container(container)
bundle.manifest_bytes       # checkpoint/manifest.json
bundle.payload_members      # checkpoint/payloads/<sha256>.bin -> bytes
bundle.container_digest
loaded = load_checkpoint_container(
    bundle.manifest_bytes, bundle.payload_members
)
```

The manifest has format `vorpy.system.checkpoint.container`, container version
`1`, `system_id`, an inner `document` mapping with format
`vorpy.system.checkpoint`, and a `payloads` metadata list. The bundle is
writer-neutral; VECTOR places those members in the ZIP and also writes the
root source copy required by the frozen output contract. No source or payload
bytes are regenerated during load.

For an update, VECTOR calls
`prepare_checkpoint_update(replacement, current=loaded)` (or supplies an
explicit `expected_current_digest`). The plan exposes `validate_current()`,
`save_bundle()`, and `validate_written()`. `CheckpointLoadError` covers bad
JSON, missing members, or digest mismatch; `CheckpointSaveError` covers an
invalid/unserializable bundle; `CheckpointUpdateError` covers stale current
state or failed post-write validation.

## VECTOR handoff

VECTOR owns the archive byte container, encoding/decoding, source-copy layout,
same-directory temporary file, flush/sync, completed-manifest validation, and
atomic replace. The writer should use the policy flags as required behavior:
write the bundle to a same-directory temporary archive, flush and sync it,
reload it through `load_checkpoint_container`, call `validate_written()`, and
only then atomically replace the destination. Any failure before replacement
must delete only the temporary file and leave the prior valid checkpoint in
place. No
`network_archive.py` or output exporter was changed for this milestone.

Existing v3 archives remain readable by the existing loader and are not
automatically rewritten. A v3 archive without embedded source bytes cannot be
silently promoted to a new self-contained checkpoint; it remains a legacy
state and source validation continues to fail closed when its original source
is unavailable.

## FORGE handoff: Group C

`GroupCExtensionContract` records the new group, generator universe, base
source-network identity, dependency digest, and stable feature identities.
`merge_supported` is always `False`; `validate_group_c_extension` only
checks append-only identity compatibility and rejects a claimed certified
merge. FORGE must prove topology/incidence merge semantics before a later
schema revision can enable it.

## Ownership and limitations

This is an additive ATLAS contract. It does not replace classic AW/POW/PRM
geometry logs, mutate a primary network, solve geometry, migrate v3 archives,
or write a checkpoint. VECTOR should integrate the manifest before changing
export timing or archive serialization; FORGE remains responsible for any
future topology merge proof.
