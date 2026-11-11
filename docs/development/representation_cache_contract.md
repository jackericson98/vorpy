# Representation cache contract

This is a persistence extension point, not a geometry implementation.

The current reader rejects unknown object kinds and fields. Therefore v2 is not
an open extension format. New representation-cache fields are reserved by
archive format version 3; the reader still accepts old v1 and v2 archives.
An archive labelled v2 that already contains the reserved representation fields
is rejected rather than silently treated as a v2 document.
Future adapter kinds must be registered before archives containing them are
read.

## Cache layout

An interface may later expose two optional mappings:

```text
representation_caches[key] -> registered cached adapter
representation_links[key]  -> IDs and provenance only
```

The key is a stable representation name, for example:

```text
physical_interface
dual_contact_complex
generator_complex
molecular_contact_surface:A
molecular_contact_surface:B
dual_separator
```

The codec stores registered adapter objects. It does not require every adapter
to contain surfaces, area, or curvature. Link records must use archive IDs, not
Python object IDs.

## Adapter registration

Agent #3's final dataclasses should be added to
`vorpy/src/io/scientific_adapters.py`'s `registry()` with a unique adapter kind.
Nested custom records need their own registry entries. Adapter fields must be
explicit and codec-supported; arbitrary constructors, imports, and pickle are
not allowed. Include an `archive_id` and `archive_provenance` when the record
has independent scientific identity. Network/interface ownership is checked
from those stable provenance IDs.

## Molecular contact surface envelope

Required future metadata:

- definition and schema version
- side (`A` or `B`)
- molecular surface convention
- system and frame identity
- source network and interface IDs
- selected-contact query ID
- contact-selection provenance
- radius model and probe radius
- canonical environment (`dry` or `solvent_competing`)
- source atom stable IDs
- support/status and provenance

Pending geometry fields, exact shape owned by Agent #3:

- spherical patch records
- circular boundary arcs
- junction vertices
- connected-component IDs
- atom/contact provenance on geometric records
- area and topology summaries
- mean/Gaussian curvature summaries
- unresolved/support diagnostics

## Capability and dependencies

Capability reads `representation_links` and cache metadata only. It never
constructs geometry. A link may carry dependency states such as:

```text
atom_centers_radii
stable_atom_ids
selected_ab_contacts
molecular_boundary_convention
spherical_clipping
```

This permits messages such as “selected contacts available” or “clipping
geometry unavailable” while the representation itself is `NOT_CALCULATED`.

## Round trip

After a real molecular contact surface is cached:

```text
save -> fresh process -> load -> inspect/export
```

must use only the cached adapter and links. It must not solve a network, select
contacts, clip spheres, rebuild arrangements, or integrate curvature.

Old archives without these optional mappings load normally. They report the
requested molecular-contact representation as `NOT_CALCULATED`; no cache is
synthesized during load.

## Cache-only visualization

`export_molecular_contact_surface()` and
`export_molecular_contact_bundle()` consume loaded canonical records only:

- canonical: `SphericalPatch`, `CircularArc`, `JunctionVertex`, and available
  `RefinedSeamPiece` records
- derived: rendering triangles, OBJ boundary/junction files, CSV mappings,
  and visualization manifests

The renderer does not reload OBJ/OFF output as scientific state. Missing
representations are reported as `NOT_CALCULATED`; partial topology does not
prevent rendering known patches and arcs. Bundle paths are centralized under
`interfaces/<interaction>/` with independent A/B representation directories
and a `visualization/` manifest directory.

## VECTOR V6 archive identity boundary

The current `.vpy` format is a versioned ZIP archive. Version 3 contains
`metadata.json`, `state.json`, and codec-owned `arrays/` and `tables/` members.
`state.json` is an explicit object table: registered scientific adapter kinds
and their declared fields are encoded with references, while unknown kinds or
fields are rejected. The reader still dispatches archive versions 1, 2, and
3. Version 1/2 projections use their historical field registries; version 3
is the first version that may contain `representation_caches` and
`representation_links`. A legacy archive is never upgraded by inference.

Cached Results are not currently stored as `InterfaceResult`, dual-result, or
alpha-result rows. The archive stores the producer-owned cache objects (for
example `InterfaceGeometryAnalysis`, `DualComplex`, filtration, and registered
representation adapters) plus their explicit metadata and provenance. VECTOR
adapts those objects after load and regenerates canonical logs. Loading and
export must not solve, rebuild geometry, or turn a missing cache into a value.

### Solve-context fields

`SolveContextKey` is not yet an archive field. When FORGE/ATLAS agree on a
first-class context record, these existing fields are safe JSON data after
tuple-to-array conversion:

`schema_version`, `scheme`, `scope_mode`, `network_scope`, `frame_id`,
`partition_scheme`, `boundary_mode`, all `*_digest` fields,
`physical_generator_ids`, `local_to_system`, `boundary_generator_ids`,
`membership_digest`, `resolution`, `solver_version`, and `digest`.

`solve_universe_id`, `context_id`, and `competitor_digest` are derived from
the persisted fields and must be recomputed, not stored as independent
authority. Runtime `geometry_cache_key` values contain Python table/object
identities and must never be archived. Paths, timestamps, metric values, and
row positions are not identity evidence. A positional generator fallback is
local evidence only and cannot certify a shared-network membership.

### Shared-network round-trip contract

- `Network.archive_id` identifies the source-network namespace. It is
  distinct from `solve_universe_id`, which identifies solve inputs.
- A future `selection_id` must identify the exact ordered group/interface view,
  membership evidence, orientation, and contact-selection definition. The
  current archive preserves existing `interface_id`, group IDs, and explicit
  local membership fields, but has no first-class NetworkView selection record
  and must not infer one.
- `local_to_system` and stable generator IDs preserve the local-to-system
  namespace. Boundary generators retain an explicit null system mapping.
- Cached dependencies must retain source network/interface IDs, source
  revision or digest, solve-context identity, selection identity, producer
  schema, and completeness/status evidence. Until a typed manifest exists,
  `Provenance.details["cache_dependencies"]` is the reserved location.

On restore, a missing required identity is unavailable or unresolved; a known
changed source geometry, selection, revision, or dependency digest is
`STALE`; unsupported representations remain `NOT_SUPPORTED`. A matching but
incomplete producer snapshot remains `PARTIAL`. No status is promoted by a
parent network, matching values, or the ability to rebuild geometry.

### Migration rule

This milestone makes no production schema change. Existing v1/v2 readers and
the v3 writer remain unchanged. A future context/dependency manifest should
follow an expand/verify/contract migration: add an explicitly versioned,
optional record or introduce a new archive version only when the fields are
agreed; read old archives unchanged; validate new records fail-closed; and
retain a writer path that emits the pre-manifest form for rollback. No old
`network_id` or provenance field may be reinterpreted as a solve key or view
selection ID.

## VECTOR V7 view-aware boundary

The current FORGE view API is a read-only snapshot over an already-solved
source network. It exposes three distinct identities that must remain distinct
after a future archive integration:

- `source_network_id` is the stable namespace owning local topology rows.
- `solve_universe_id` is the source solve-context identity.
- `selection_id` is the deterministic digest of the view kind, ordered A/B
  system-generator selection, source network, and source solve context.

`NetworkView` itself is deliberately not an archive object in V7. After a
round trip, VECTOR should reconstruct a view from the restored source network
and the canonical selection record, then verify all three identities. It must
not persist row snapshots, Python object IDs, or a view under a new unversioned
field. The restored view must preserve `local_to_system`, stable generator IDs,
boundary mappings, explicit feature incidence, and the distinction between
`A -> B` and `B -> A`; reversing the groups creates a different selection ID.

Cached dual memberships remain producer evidence, not view inference. A/B
membership IDs must be checked against the restored source-network namespace,
solve context, ordered selection identity, producer revision, and completeness
before cached dual Results are adapted. A mismatch in any known dependency is
`STALE`; missing identity is unavailable or unresolved; unsupported A/B output
remains `NOT_SUPPORTED`. VECTOR must never rebuild geometry or infer a missing
membership from coordinates or row order.

The NetworkView/archive selection schema remains unchanged until ATLAS
finalizes canonical identity fields and FORGE validates view semantics. The
safe interim contract is
test-only: prove that a restored source network can recreate the same view
identity and cached memberships without a `network_view` archive record. A
future manifest may add explicit source/view/dependency records through the
V6 expand/verify/contract migration path, with v1/v2/v3 readers preserved.

## VECTOR V9 boundary metadata

v3 now has one optional network field, `boundary_generators`, using the
existing registered `boundary_generator` adapter. This preserves a
subordinate network's explicit boundary generator positions/radii when they
differ from the parent System. The field was added after the v1/v2 field
snapshots, so v1/v2 compatibility is unchanged and older v3 archives remain
readable with the field absent.

An archive without this field does not synthesize network-owned boundary
generators from topology rows or coordinates. `boundary_indices` can still
identify a boundary cell (for example, `boundary:4`), but the missing
generator metadata produces a different boundary-generator digest and is not
equivalent to the certified round-trip identity. New v3 archives must retain
the explicit adapter records; loaders never rebuild geometry to recover them.

Boundary configuration and generator adapter records are object-table
references, not scientific identity records. Archive preparation must not add
runtime `archive_id` values to these frozen values because that would alter
the boundary and environment digests merely by saving an otherwise unchanged
network. Older archives that already contain such IDs remain readable.

## VECTOR integration checkpoint

F10 producer metadata is now carried by optional v3 fields. The System stores
the explicit solvent mode; the Network stores the mode, local-to-System
generator indices and stable IDs, participating solvent IDs/digest, and
solvent-cell completeness. These fields are added after the v1/v2 field
snapshots, so existing archive versions and older v3 archives remain
readable. Loading a legacy archive without them sets the optional provenance
to explicitly unresolved (`None`/empty membership); it does not infer a mode
or membership from atom labels, coordinates, row counts, or settings.

The archived producer evidence is sufficient for a no-rebuild round trip, but
it does not create a NetworkView or certify a new analysis. View selection and
cached Results still require their own source/solve/selection identity and
completeness checks.

## VECTOR V11 archive/export integration

Archives persist the solved Network and cached producer objects, not runtime
typed Results rows. Canonical export adapts those caches after load into the
same deterministic Results identities and log rows, including separate dual
contact, A, and B representations. A matching cache may be exported without
geometry rebuild; a missing, stale, incomplete, or unsupported dependency
remains unavailable or non-certified according to the Results contract.

Two selections may therefore share one restored Network and solve identity
while retaining distinct selection identities. Export tests must compare the
source namespace, solve context, selection, partition, provenance, and row
bytes rather than relying on Python object identity.

### V11 integration handoff

The authoritative optional v3 archive fields for F10 are:

- `system.solvent_mode` — producer-normalized `MOLECULE_ONLY`,
  `SOLVENT_CONTEXT`, or `FULL_SOLVENT`.
- `network.generator_system_indices` and `network.generator_stable_ids` —
  local-generator to System/global membership, including explicit boundary
  `None` mappings where present.
- `network.solvent_local_indices` and `network.solvent_generator_ids` —
  explicit participating solvent membership only.
- `network.solvent_mode`, `network.solvent_cells_complete`, and
  `network.solvent_generator_digest` — producer evidence and completeness;
  these are never inferred during load.
- `network.retained_system_indices`, `network.retained_generator_ids`, and
  `network.retained_generator_digest`, plus `system.retained_system_indices` —
  F12 explicit retained molecular competitors. Retention is distinct from
  solvent membership and must remain in solve identity and provenance.
- `network.boundary_generators` — the V9 explicit boundary-generator adapter
  records. Boundary adapter object IDs are not scientific identity.

NEXUS should reuse `network.archive_id` as the source namespace, the loaded
`get_solve_context_key()` for solve identity and mapping, and `make_view()` for
selection identity only when completeness is proven. ATLAS should consume the
cached dual's explicit `metadata["explicit_group_membership"]` through
`adapt_dual_contact_complex()` and `adapt_dual_generator_complex()`; A/B
membership is never reconstructed from simplex rows. PRISM may display only
the typed Result status/provenance and the archived solvent evidence.

The canonical export entry points are `_cached_dual_result_bundle()` for
`dual_contact`, `dual_a`, and `dual_b` typed Results,
`_cached_alpha_selection_result()` for cached alpha state, and `_compact_log()`
for deterministic eight-column interface rows plus serialized provenance.
Missing or unsupported cache dependencies remain `NOT_CALCULATED`,
`NOT_SUPPORTED`, `PARTIAL`, or `STALE` as supplied by the Results adapters.

Reusable V11 archive assertions are: pre-save versus post-load solve-key
equality; boundary/environment digest equality; source archive ID equality;
local-to-System and stable-generator mapping equality; solvent provenance and
completeness equality; explicit cache dependency/provenance equality; absence
of serialized `NetworkView` or typed Result objects; deterministic repeated
`logs.csv` bytes; and `Network.build()`/analysis rebuild guards. The current
real EDTA view regression intentionally expects `NetworkViewError` when the
selected closure contains an incomplete generator. V8 synthetic archive-view
coverage remains the positive reconstruction case.

## VECTOR V13 group-layer archive contract

New v3 archives stamp `metadata.producer_contracts["group.layer_atoms"]` as
`system_global_nonnegative_v1`. Writers and readers then require every layer
atom to be a valid nonnegative System row ID. Boundary generators remain in
Network topology mappings such as `boundary_indices`, `generator_system_indices`,
and `layer_net_atoms`.

Marker-free v1/v2/v3 archives retain the historical reader rule that permits
the legacy `-1` layer sentinel. The format version remains unchanged; absence
of the producer contract is the compatibility signal. Re-saving a legacy
archive writes the current contract and therefore rejects any still-present
sentinel instead of rewriting it.
