# A9/A11: Canonical shared-network, view, and solvent identity

Status: Results contract finalized against the current FORGE F6/F7 APIs, with
A11 solvent-mode provenance added as an additive Results boundary.
These milestones add Results-owned identity/dependency checks and focused tests;
it does not edit FORGE-owned solve/view code, geometry builders, archives,
visualization, or exporters.

## Goal

Represent several group and interface analyses derived from one solved global
network without making them look like independent solves. The contract must
also prevent an isolated-group solve from being confused with a view over a
shared network, or a cache from being accepted after its solve context
changes.

The contract is an identity and evidence boundary. It does not build a network,
select geometry, infer memberships from coordinates, or promote a cached value
to `CERTIFIED`.

## A9 canonical decisions

The following names are canonical at the Results boundary:

| Concept | Canonical name and value |
| --- | --- |
| Source network namespace | `source_network_id`: the explicit producer network namespace, preserved through archive reload. It is opaque and distinct from solve identity. |
| Source identity authority | `source_identity_status`: `AUTHORITATIVE`, `PROVISIONAL`, or `UNRESOLVED`; namespace presence alone is never certification evidence. |
| Solve universe | `solve_universe_id`: the exact F6 `SolveContextKey.solve_universe_id`, currently `solve:<context-digest>`. F6's raw `context_id` is not a Results alias. |
| Selection/view | `selection_id`: the exact F7 `NetworkView.selection_id`, currently `view:<selection-digest>`. There is no second `view_id`. |
| Frame | `frame_id`: the producer's canonical non-empty string, including the F6-normalized numeric frame. FORGE exposes `STATIC_FRAME_ID = "static"` for explicit static-structure attestation; `default` remains the legacy unknown/provisional sentinel. |
| Results analysis scope | `analysis_context_kind`: `GLOBAL_NETWORK` or `ISOLATED_GROUP`. |
| F6 solve scope | `scope_mode`: F6's `full`, `group`, or `interface`; it is not a view selection and not the Results `network_scope`. |
| Solvent scientific context | `solvent_mode`: `MOLECULE_ONLY`, `SOLVENT_CONTEXT`, or `FULL_SOLVENT`, with explicit solvent-generator and cell-completeness evidence where applicable. |
| A/B order | `group_A_selection_id`, `group_B_selection_id`, and `orientation`, in that order. Reversal is a different selection. |

`source_network_id` must be supplied explicitly for new certified Results. F7's
fallback `network:<solve-context-digest>` is a provisional namespace, not a
canonical source-network ID. Results may retain that value as legacy evidence,
but it is `UNRESOLVED` for identity-required validation and can never support
`CERTIFIED` or `EXACT_REUSE`. It is not silently rewritten. A future migration
must carry an explicit alias-to-canonical mapping, the exact
`solve_universe_id`, producer revision, and an attestation that both names
refer to the same source network.

`SharedNetworkIdentity.source_identity_status` is the explicit authority
boundary. It defaults to `UNRESOLVED` for compatibility with legacy callers,
but `validate_results_dependency()` and solve-context validation require
`AUTHORITATIVE`. `PROVISIONAL` is retained for non-certifiable internal
namespaces; `UNRESOLVED` covers missing or contradictory producer evidence.
An `AUTHORITATIVE` status paired with the solve-derived `network:<digest>`
namespace is rejected as contradictory. `archive_id`, `source_revision`, and
`source_digest` identify cache artifacts or dependencies only; none can supply
source-network authority.

## Existing Results to reuse

| Existing type | A9 canonical reuse |
| --- | --- |
| `ResultIdentity` | Per-result system, frame, network, current canonical partition scheme, group/interface, side, and representation identity. |
| `Provenance` | Producer/source revision, immutable context details, diagnostics, and (until a typed field exists) serialized cache dependency records. |
| `GeometryLineage` | Direct solve versus materialized view versus exact reuse, source network/view, selection scope, and completeness evidence. |
| `Status` | Independent context, result, quantity, and dependency status. |
| `Quantity` and existing result subclasses | Measurements remain producer-owned and status-bearing. |
| `adapt_cached_result()` / `adapt_interface_analysis()` | Continue adapting already-produced blocks; adapters validate identity and copy cached values only. |
| `adapt_atom_metadata()` | Continue resolving stable atom identity; it must not use a generator index as an atom identity. |

No new CSV columns are needed. VECTOR continues to own serialization of the
existing eight-column interface log. The fields below are a proposed typed
contract, not a request to change that export.

The concrete reuse mapping is:

- `ResultIdentity.system_id`, `system`, `frame`, `network_id`,
  `network_scope`, `partition`, `radius_configuration_id`, `group_id`,
  `interface_id`, `group_A`, `group_B`, `side`, and `orientation` remain the
  per-result projection. `source_network_id`, `source_interface_id`, and
  `contact_selection_source`/`contact_selection_id` carry parent and contact
  selection references where already supported.
- `GeometryLineage.solve_universe_id`, `geometry_source_id`,
  `geometry_version`, `solver_settings_digest`, `environment_digest`,
  `source_network_id`, `source_view_id`, `selection_scope`, `compatibility`,
  and `completeness` are the existing compact geometry dependency boundary.
  A9 requires them to agree with the typed context; it does not replace them.
- `Provenance.source`, `source_revision`, `scientific_state_id`, `details`,
  and `diagnostics` remain the provenance envelope. A typed
  `Provenance.dependencies` field is a later additive option; until then,
  only `details["cache_dependencies"]` is the reserved serialized location.
- Existing `Status`, `Quantity`, `NetworkResult`, `GroupResult`, and
  `InterfaceResult` continue to carry quantity-level evidence. Context status
  must not overwrite a quantity status.

## 1. Two analysis contexts

`network_scope` remains the existing Results ownership vocabulary (`system` or
`interface`). It is not a solve-context kind. The current FORGE
`SolveUniverseKey.scope_mode` (`full`, `interface`, or `group`) is also not a
semantic substitute for analysis context: it describes the solve's competitor
scope, while the field below describes how the Result relates to that solve.

```python
AnalysisContextKind = Literal[
    "GLOBAL_NETWORK",       # one solve contains the shared universe
    "ISOLATED_GROUP",       # the group is the solve universe
]
```

### `GLOBAL_NETWORK`

- `source_network_id` identifies the stable network namespace and is distinct from a
  Python object ID, archive path, or the solve-key digest.
- All group and interface records are selections/views over that network.
- A group or interface analysis must retain the same solve-key identity,
  network namespace, frame, geometry scheme, partition scheme, resolution,
  environment, solver settings, and competitor universe as the global solve.
- A view may be `MATERIALIZED_VIEW` or `EXACT_REUSE` in geometry lineage, but
  that lineage never raises its scientific status.

### `ISOLATED_GROUP`

- The selected group is the solve universe, not merely a view. It has its own
  solve-key identity, even if its values happen to match a global view.
- The result must use the isolated network identity and an explicit selection
  policy describing excluded competitors/boundary generators.
- If it originated from a global network, `source_network_id` may point to the
  global network, but the result must not claim `EXACT_REUSE` without the
  existing positive compatibility and completeness evidence.
- It must not share a `result_id` with a global-network view merely because
  group names and measurements match.

Proposed typed context envelope:

```python
@dataclass(frozen=True)
class AnalysisContext:
    analysis_context_id: str          # typed envelope digest, not F6 context_id
    kind: AnalysisContextKind
    system_id: str
    frame_id: str                    # canonical non-empty FORGE frame key
    source_network_id: str            # explicit source-network namespace
    network_scope: Literal["system", "interface"]
    geometry_scheme: str              # SolveUniverseKey.scheme
    partition_scheme: str             # SolveUniverseKey.partition_scheme
    scope_mode: str                   # SolveUniverseKey.scope_mode
    radius_configuration_id: str
    solve_universe_id: str             # SolveContextKey.solve_universe_id
    solve_key_schema_version: int
    solver_version: str
    resolution: object
    solver_settings_digest: str
    environment_digest: str
    generator_digest: str
    geometry_digest: str
    competitor_digest: str
    membership_digest: str
    solvent_mode: Literal["MOLECULE_ONLY", "SOLVENT_CONTEXT", "FULL_SOLVENT"]
    solvent_generator_digest: str | None
    solvent_cells_complete: bool | None
    status: Status
    completeness: "CompletenessEvidence"
```

`solve_universe_id` is the exact F6 `SolveContextKey.solve_universe_id`
(`solve:<digest>`), not a second independently generated ID. `network_id` is
the legacy per-result field; new shared-network provenance uses the explicit
`source_network_id` namespace above. The Results contract never derives a
source network namespace from a solve digest.

`geometry_scheme` and `partition_scheme` are separate. The former is the
geometric producer scheme (`SolveUniverseKey.scheme`); the latter is the
partition vocabulary currently projected by Results as `ResultIdentity.partition`.
`ResultIdentity.partition` must equal `partition_scheme`, not
`geometry_scheme`. A future PRM or non-AW producer requires an explicit
vocabulary/schema migration; it must not be silently aliased into `power`.

`frame_id` is the canonical, non-empty string that participates in the solve
key. An integer legacy `ResultIdentity.frame` may be projected into it only by
an explicit adapter. FORGE's `STATIC_FRAME_ID` (`"static"`) is the explicit
sentinel for a producer-attested static structure. The legacy `default`
sentinel means unknown/unresolved unless the producer separately attests its
static meaning, and cannot support `CERTIFIED` shared references by itself.

The solve identity digest is over canonical, sorted, versioned solve data only:
no metric values, paths, timestamps, row positions, or Python object IDs. The
solve-key component digests (`generator_digest`, `geometry_digest`,
`environment_digest`, `solver_settings_digest`, `competitor_digest`),
resolution, solver version, frame, schemes, scope mode, and schema version are
all part of validation. A change in any one creates a different solve
identity.

Any future Results `analysis_context_id` is a typed envelope ID and must not
alias F6's raw `context_id`. It hashes identity fields only: no metric values,
paths, timestamps, or Python object IDs. Results stores the canonical
`solve_universe_id` for solve identity.

### Additive Results identity fields

The existing fields remain the source of truth for legacy results. New A9
results may add these optional fields to `ResultIdentity` and require them when a
shared/isolated context is known:

```python
analysis_context_id: str | None
analysis_context_kind: AnalysisContextKind | None
selection_id: str | None
solvent_mode: Literal["MOLECULE_ONLY", "SOLVENT_CONTEXT", "FULL_SOLVENT"] | None
```

`selection_id` is the stable digest/ID of the exact group or interface
selection, including membership evidence and ordered A/B identity. It is not a
directory name or a Python object ID. `source_network_id` continues to identify
the parent global network when an isolated result was derived from one. These
fields are additive and required for new shared-network results; legacy
results remain valid under their existing identity rules.

### F6/F7 reconciliation and migration boundary

The current FORGE implementation provides `SolveUniverseKey` with
`scheme`, `scope_mode`, `frame_id`, `partition_scheme`, environment/generator/
geometry/solver digests, resolution, solver version, and a final `digest`.
`SolveContextKey` remains the compatibility API name for that type. Its
`solve_universe_id` property is the canonical `solve:<digest>` Results value.
F7 `NetworkView` carries the explicit source namespace, the same solve ID, and
the canonical `view:<digest>` selection ID.
F7's `NetworkView.network_id` property is a compatibility alias for
`source_network_id`; Results does not use it as a second identity. Likewise,
F6's raw `context_id` is not persisted as a competing solve ID.

The A6 draft's ambiguous `scheme`, nullable frame, and solve/network/view
conflation are resolved here. Results uses `geometry_scheme` for F6 `scheme`,
`partition_scheme` for F6 `partition_scheme`, `scope_mode` for F6 solve scope,
`source_network_id` for the explicit source namespace, and `selection_id` for
F7 view identity. No FORGE field is renamed in this milestone. FORGE still
needs to replace its fallback `network:<solve-context-digest>` with an explicit
source namespace; until then Results treats that fallback as provisional.

## A11 solvent-mode provenance

The default scientific context is `MOLECULE_ONLY`. Solvent atoms may remain in
the loaded molecular structure, but their presence alone does not make them
part of the generator universe and must not change a Results identity. The
producer must explicitly attest the participating generator evidence:

| Mode | Required evidence | Certification meaning |
| --- | --- | --- |
| `MOLECULE_ONLY` | No participating solvent-generator digest; solvent cells are not complete. | Solvent is excluded from the solve universe. |
| `SOLVENT_CONTEXT` | A stable digest of the participating solvent generator IDs; `solvent_cells_complete` is false or unknown. | Solvent generators influence molecular geometry; complete solvent cells are not required. |
| `FULL_SOLVENT` | A stable solvent-generator digest and `solvent_cells_complete=True`. | Molecular and solvent cells were solved; other completeness gates still apply. |

This is provenance, not a request to change network construction. Results never
classifies solvent from coordinates, atom names, loaded-structure presence, or
the legacy `Environment` values (`dry` and `solvent_competing`). That existing
environment vocabulary and archive identity remain unchanged.

F6 already hashes the actual generator records through `generator_digest`,
`geometry_digest`, and `membership_digest`. Therefore a molecule-only solve and
a solvent-participating solve are distinguishable when the producer's solved
network actually excludes or includes different generators. F6 currently does
not classify which generators are solvent and does not encode solvent-cell
completeness; its opaque `environment_digest` is insufficient to distinguish
`SOLVENT_CONTEXT` from `FULL_SOLVENT`. The three fields above are therefore
required in the Results identity/provenance envelope. They change the
`ResultsDependency` digest even when the supplied F6 solve ID happens to remain
the same, so cached geometry-dependent Results become `STALE` on a mode,
solvent-generator, or cell-completeness change.

If a producer cannot attest the participating solvent IDs or cell completeness,
the adapter must remain `UNRESOLVED` (or preserve `STALE` for an invalidated
cached result); it must not infer a mode or promote a result to `CERTIFIED`.
No geometry is rebuilt during this validation. F10 now provides the explicit
classification and evidence through `Network.solvent_provenance`; ATLAS
consumes that mapping without recomputing membership or completeness.

The Results adapter accepts an optional authoritative `source_solve_context`.
When supplied, `validate_solve_context_identity()` compares the exact F6
`solve_universe_id`, `frame_id`, `scheme`, `partition_scheme`, and `scope_mode`
against the typed Results identity. Missing or contradictory fields are
`UNRESOLVED`. This validates solve lineage; it does not derive solvent mode
from F6's current opaque environment or generator fields. If the supplied
producer context exposes the canonical `solvent_mode`,
`solvent_generator_digest`, and `solvent_cells_complete` fields, all three are
also compared; a partial or contradictory set is `UNRESOLVED`. F10 provides
those three fields explicitly; current callers must pass its mapping as
`source_solvent_provenance` before solvent-aware certification.
The A10 group adapter also accepts F10's `source_solvent_provenance` mapping
(the existing `Network.solvent_provenance` value) and compares the same three
canonical fields. Additional membership evidence is retained as provenance but
is not re-derived by ATLAS.

## 2. Explicit membership evidence

Atom membership and generator membership are different namespaces and must be
retained separately. Network-local generator indices are valid only together
with the `source_network_id` that owns them.

```python
@dataclass(frozen=True)
class MembershipEvidence:
    membership_id: str
    atom_ids: tuple[str, ...]                    # stable system atom IDs
    generator_indices: tuple[int, ...]           # local to source_network_id; supplemental
    stable_generator_ids: tuple[str, ...]        # producer-stable IDs
    atom_to_generator: tuple[tuple[str, int], ...] = ()
    generator_namespace: str                    # must equal context.source_network_id
    source: str                                  # explicit producer evidence
    digest: str
    status: Status
    completeness: "CompletenessEvidence"
```

Rules:

- FORGE supplies the membership sets and mapping; ATLAS never derives them
  from geometry, row order, distances, or group names.
- `stable_generator_ids` are required for a shared reference that is to be
  `CERTIFIED`. A positional fallback such as `ball:{row}` is local evidence
  only and forces `PARTIAL` or `UNRESOLVED`; it cannot certify cross-network
  membership.
- `atom_to_generator` is a declared mapping, not a promise that every atom has
  exactly one generator. Any multiplicity or omitted mapping is explicit.
- `generator_indices` are meaningful only with `generator_namespace`; they do
  not replace stable generator IDs.
- A group can have a valid atom selection while generator membership is
  `NOT_CALCULATED`; the resulting group/interface result is never certified
  and is partial or unavailable as appropriate.
- `membership_digest` is the canonical digest of the ordered stable atom IDs,
  stable generator IDs, namespace, explicit mapping, producer schema, and
  completeness evidence. Changing one selected atom/generator creates a
  different identity.

## 3. Group and interface selections

The shared context contains immutable selection records. These are references
and evidence, not geometry results.

```python
@dataclass(frozen=True)
class GroupSelection:
    selection_id: str
    group_id: str
    label: str
    analysis_context_id: str
    membership: MembershipEvidence
    status: Status


@dataclass(frozen=True)
class InterfaceSelection:
    selection_id: str
    interface_id: str
    analysis_context_id: str
    group_A_selection_id: str
    group_B_selection_id: str
    orientation: str                 # ordered A -> B convention
    contact_selection_source: str | None = None
    contact_selection_id: str | None = None
    status: Status = Status.NOT_CALCULATED
```

`group_A_selection_id` and `group_B_selection_id` are ordered references to
the exact group selections, not merely group labels or group definitions.
Reversing them creates a distinct interface selection and must not be
normalized away. The interface selection digest includes both ordered group
selection IDs, orientation, interface identity, and any contact-selection
definition/source revision. The existing
`ResultIdentity.group_id`, `group_A`, `group_B`, `interface_id`, `side`, and
`orientation` fields remain the per-result projection of these records.

### Selection does not change the solve universe

F7 views are created from an already-solved source network. Therefore changing
the selected group, changing the ordered A/B pair, or changing a contact
selection changes `selection_id` and dependent Results digests only; it does
not change `solve_universe_id`, `frame_id`, or the solve-key component set.

This rule applies to view selections, not to the solve itself. Changing the
F6 competitor universe or solve membership (`scope_mode`, source generator
membership, or boundary membership) is a solve-context change and must produce
a different `solve_universe_id`. Likewise, generator coordinates/radii or
stable IDs, frame, geometry scheme, partition scheme, environment/boundary
configuration, resolution, build options, solver settings, or solver version
are geometry/solve identity inputs and must change the solve ID.

Global and isolated analyses are therefore distinct even when their selected
atoms or measurements match: a global view reuses the global solve ID, while
an isolated group owns a separate solve ID and `ISOLATED_GROUP` context.

Proposed context envelope addition:

```python
@dataclass(frozen=True)
class SharedNetworkAnalysis:
    context: AnalysisContext
    groups: tuple[GroupSelection, ...]
    interfaces: tuple[InterfaceSelection, ...]
    dependencies: tuple["CacheDependency", ...]
    provenance: Provenance
    status: Status
```

This is a linkage envelope, not a replacement for `NetworkResult`,
`GroupResult`, or `InterfaceResult`. Each scientific result remains an
independent typed snapshot with its own status and quantities.

## 4. Cached geometry dependencies

Shared analysis must be able to say which cached objects it used without
rebuilding them during adaptation.

```python
@dataclass(frozen=True)
class CacheDependency:
    dependency_id: str
    kind: Literal[
        "network", "dual", "filtration", "physical_interface",
        "molecular_contact", "topology",
    ]
    source_id: str
    source_revision: str | None
    source_digest: str | None
    source_network_id: str
    solve_universe_id: str
    frame_id: str
    analysis_context_kind: AnalysisContextKind
    geometry_scheme: str
    partition_scheme: str
    scope_mode: str
    selection_id: str | None
    required: bool
    dependency_digest: str
    status: Status
    completeness: "CompletenessEvidence"
```

`dependency_digest` is the Results-owned `dependency-v1:<sha256>` digest over
the canonical identity payload: source network ID, solve-universe ID, frame,
analysis context kind, geometry/partition schemes, solve scope, selection/view
ID, source ID, source revision or content digest, required flag, and sorted
upstream dependency digests. Status and measured values are not digest inputs.
The implementation is `ResultsDependency.dependency_digest`.

Dependency validation is an identity check, not a cache-miss recovery path.
For every required dependency, the adapter compares the expected context
identity with the cached entry: source network namespace, solve-key digest, frame,
geometry scheme, partition scheme, scope mode, selection ID, source ID and
revision/digest, producer schema, and completeness evidence. A missing
required field is not a match. A mismatch against a previously valid snapshot
is `STALE`; contradictory fields within one record are `UNRESOLVED`. No
network solve or geometry rebuild is allowed during adaptation.

Dependency IDs and the manifest digest use canonical, sorted, versioned data.
The manifest includes the context digest, every dependency identity, required
flags, upstream dependency digests, and producer revisions. It excludes paths,
timestamps, object IDs, and measured values. Reordering independent records
must not change the manifest digest; changing any identity-bearing field must.
An archive restore is valid only when the manifest digest and all required
entries validate recursively. Optional dependency failure makes only the
dependent representation unavailable/partial; it never disappears silently.

`GeometryLineage` remains the compact per-result lineage record. Its
`solve_universe_id`, `source_network_id`, `source_view_id`, `selection_scope`,
and `completeness` fields should reference the same context and dependency
evidence. `EXACT_REUSE` keeps its existing strict compatibility requirements.

If dependency records need first-class access rather than only serialized
details, the additive proposal is `Provenance.dependencies:
tuple[CacheDependency, ...] = ()`. Until that field exists, the canonical
serialized location is `Provenance.details["cache_dependencies"]`; adapters
must validate it as typed data and must not treat arbitrary detail keys as
dependencies.

## 5. Completeness and certification

```python
@dataclass(frozen=True)
class CompletenessEvidence:
    scope: str
    completeness_proven: bool
    expected_count: int | None = None
    observed_count: int | None = None
    missing_ids: tuple[str, ...] = ()
    unresolved_ids: tuple[str, ...] = ()
    method: str = ""
    reason: str | None = None
```

Status rules:

| Status | Meaning in the shared-network contract |
| --- | --- |
| `CERTIFIED` | The producer explicitly certified the requested result and all required context, membership, selection, and geometry dependencies are complete. |
| `PARTIAL` | Some requested subset or dependency is valid and measured, but completeness is not proven. Values retain their narrower scope. |
| `UNRESOLVED` | Evidence is contradictory or a required mapping/identity cannot be resolved. |
| `NOT_CALCULATED` | The requested analysis/dependency was not run or is absent from the cache. |
| `NOT_SUPPORTED` | The representation or context mode is not supported by the producer. |
| `STALE` | A cached snapshot exists, but its context/dependency identity no longer matches. |

Fail-closed precedence for a requested dependency is:

1. contradictory identity or membership evidence -> `UNRESOLVED`;
2. unsupported scheme/context/representation -> `NOT_SUPPORTED`;
3. no prior cache/result -> `NOT_CALCULATED`;
4. known prior-valid cache with changed identity/revision -> `STALE`;
5. matching identity with incomplete producer evidence -> `PARTIAL`;
6. matching identity, stable memberships, complete required dependencies, and
   explicit producer certification -> eligible for `CERTIFIED`.

Missing identity evidence is not a mismatch: it is `UNRESOLVED`. A status is
never promoted because a parent context is certified, because values happen
to agree, or because geometry could be rebuilt. There is no zero-fill fallback.

Context status never promotes a result. A `CERTIFIED` global network does not
certify every group/interface selection, and `EXACT_REUSE` does not promote a
quantity. Result and quantity status must be preserved independently.

## 6. Proposed Results API surface

The minimal API addition is a context validator/factory layer around existing
adapters:

```python
def make_result_identity(
    context: AnalysisContext,
    *,
    result_kind: str,
    representation: str,
    group: GroupSelection | None = None,
    interface: InterfaceSelection | None = None,
    side: str | None = None,
    **representation_fields,
) -> ResultIdentity:
    """Project one validated selection into the existing identity model."""


def make_context_provenance(
    context: AnalysisContext,
    *,
    selection: GroupSelection | InterfaceSelection | None = None,
    dependencies: tuple[CacheDependency, ...] = (),
    source: str,
    diagnostics: tuple[str, ...] = (),
) -> Provenance:
    """Create immutable provenance; never compute or refresh geometry."""


def validate_cache_dependencies(
    context: AnalysisContext,
    dependencies: tuple[CacheDependency, ...],
) -> Status:
    """Return dependency validity without rebuilding missing geometry."""


def adapt_shared_network_results(
    analysis: SharedNetworkAnalysis,
    *,
    cached_results: tuple[ResultModel, ...],
) -> tuple[ResultModel, ...]:
    """Validate linkage and return existing typed result snapshots."""
```

`adapt_shared_network_results()` should delegate actual blocks to existing
`adapt_cached_result()`, `adapt_interface_analysis()`, and dual/molecular
adapters. It must reject a result whose identity does not match the context;
it must not aggregate group/interface measurements into a global value.

The API should also accept an explicitly unavailable request and return the
existing `unavailable_interface()`/`unavailable_dual_*()` forms with
`NOT_CALCULATED` or `NOT_SUPPORTED`, rather than zero-filled placeholders.

## 7. FORGE F5 conceptual handoff

FORGE owns the producer-side evidence and supplies:

1. One stable solve-key identity, network namespace, canonical frame identity,
   geometry scheme, partition scheme, scope mode, resolution, solver version,
   solver-settings digest, environment digest, generator/geometry digests, and
   competitor digest. The canonical mapping is F6
   `SolveContextKey.solve_universe_id` -> Results `solve_universe_id`.
   `source_network_id` must be explicit and must not be the solve-derived
   fallback namespace.
2. Explicit atom and network-generator membership records with namespace and
   membership digests. Stable generator IDs are required for certification;
   row-position fallbacks are not.
3. Ordered group/interface selections and selection digests. Interfaces must
   reference the exact ordered group selections.
4. Read-only cached network/dual/filtration/interface objects with source IDs,
   revisions, and completeness/status evidence.
5. A dependency manifest that can be checked after archive restore without
   invoking geometry.

ATLAS owns the typed identity projection, provenance validation, status rules,
and adapter composition. GAUSS-owned measurements remain producer-supplied
typed `Quantity`/curvature/topology blocks. VECTOR owns presentation and
serialization only; this design does not alter its exports.

### Review basis and current limitations

The current F6/F7/F10 implementation is the available implementation evidence
in this checkout. F6 now exposes the compatibility name `SolveUniverseKey` and
the canonical `SolveContextKey` alias, with `solve:<digest>` IDs. F7 exposes
read-only `NetworkView` snapshots with deterministic `view:<digest>` selection
IDs, ordered A/B groups, fail-closed source-network identity validation, and F10
solvent provenance. Any provisional source namespace and positional generator
fallback remain insufficient for certified shared Results.

The Results group adapter validates the supplied solve-key fields and F10
solvent provenance, but broader solve-key component, membership, and dependency
manifest validation remains an implementation gap. That is not permission to
rebuild geometry during adaptation.

### Phase 2 migration requirements

- FORGE must publish an explicit source-network namespace and stop emitting
  `network:<solve-context-digest>` as the only source ID. No automatic aliasing
  is permitted; use an attested migration mapping if old caches must survive.
- Keep the public F6/F7 names stable: `solve_universe_id`, `source_network_id`,
  `selection_id`, and `frame_id`; retain `SolveUniverseKey` only as the
  compatibility type alias until callers migrate to `SolveContextKey`.
- FORGE must document the static-frame sentinel and distinguish it from an
  unknown/default frame. Unknown frames cannot be certified.
- Stable generator IDs remain required for certification; positional fallbacks
  remain local evidence only.
- F10 callers must pass `Network.solvent_provenance` through the Results adapter;
  a raw F6 environment value or loaded molecular structure is not sufficient to
  derive the A11 mode.
- FORGE and ATLAS must version the dependency manifest and agree on source
  revision/content-digest guarantees and explicit global/isolated evidence.

## 8. Acceptance gates before implementation

- A synthetic one-network/two-group case produces distinct group/interface
  result IDs while sharing one `source_network_id` and `solve_universe_id`.
- The same group analyzed in global and isolated contexts produces distinct
  identities and explicit lineage, even when values match.
- Atom membership, generator membership, and A/B ordering round-trip without
  coordinate-based inference.
- Missing or mismatched dual/filtration/physical caches become
  `NOT_CALCULATED`, `NOT_SUPPORTED`, `UNRESOLVED`, or `STALE` as appropriate;
  no adapter rebuilds geometry.
- A partial group or interface selection remains partial even when the global
  network is certified.
- Existing `ResultIdentity.result_id`, `GeometryLineage`, eight-column log
  schema, and current adapters remain backward compatible for legacy caches.

Focused validation remains limited to existing, implemented contracts: the
current solve-key tests cover deterministic identity, frame/scope changes,
and fail-closed unsupported/incomplete sources; the geometry-lineage tests
cover strict exact-reuse compatibility and completeness. A9 adds focused
Results identity/dependency tests without depending on unfinished FORGE
callers.

Out of scope: changes to existing scientific result row schemas, archive
registry changes, visualization/export changes, new log columns, geometry
algorithms, numerical tolerance changes, and broad benchmark validation. The
small Results-owned identity/dependency contract and its focused tests are the
only production-facing additions in A9.

## A14 retained-generator provenance and multi-result handoff

F12's authoritative retention mapping is consumed exactly as emitted by
`Network.retained_generator_provenance`:

```python
{
    "policy": "explicit_system_indices",
    "retained_system_indices": tuple[int, ...],
    "retained_generator_ids": tuple[str, ...],
    "retained_generator_digest": str | None,
    "boundary_independent": True,
}
```

The two ordered membership tuples must have equal length and the digest must
match F12's canonical SHA-256 payload. A default molecule-only solve has empty
tuples and a null digest; that explicit empty mapping is valid. Missing,
partial, reordered, duplicate, contradictory, or digest-mismatched evidence is
`UNRESOLVED`. Results never reconstruct it from coordinates, atom labels, or
generator positions.

`SharedNetworkIdentity` carries the three retention values when explicit
retention exists. They are included in `context_payload()` and therefore in
`ResultsDependency.dependency_digest`; changing retained Mg membership cannot
reuse a prior cached result. When no authoritative retention evidence exists,
the fields are omitted from the payload so legacy dependency IDs remain stable.
The group adapter accepts `source_retained_generator_provenance`, preserves the
raw mapping in `Provenance.details`, and requires that mapping for an explicitly
retained identity. A legacy identity with no retention evidence remains
compatible, but it does not claim F12 retention certification.

`SharedResultsCollection` and `collect_shared_results()` are the stable
Results-owned handoff for NEXUS/PRISM. The collection is immutable, requires
one explicit source network and solve-universe ID, and rejects entries with
different lineage or duplicate `selection_id` values. It contains typed
`GroupResult`/other `ResultModel` snapshots only; it performs no geometry or
aggregation.

VECTOR coordination: no archive schema or exporter changes are made here.
v3 archive persistence must continue to round-trip F12's four retention fields
and preserve their absence in legacy archives. A reloaded network should pass
its retained-generator mapping to the Results adapter before certification.

## A15 production certification requirements

NEXUS is the producer-side handoff owner. For every source network it must
forward the exact `Network.retained_generator_provenance` mapping and project
the same values into `SharedNetworkIdentity`:

```python
retention = source_network.retained_generator_provenance
identity = SharedNetworkIdentity(
    ...,
    retained_system_indices=retention["retained_system_indices"],
    retained_generator_ids=retention["retained_generator_ids"],
    retained_generator_digest=retention["retained_generator_digest"],
)
adapt_group_view_analysis(
    ..., source_retained_generator_provenance=retention
)
```

A retained-generator Result may be `CERTIFIED` only when all of the following
are true: source network ID is explicit and non-provisional; the F6 solve
context, frame, schemes, and scope match; solvent provenance is complete; the
three retention identity fields match the producer mapping; the retention
policy is `explicit_system_indices`; `boundary_independent` is true; the
ordered index/ID tuples are equal-length, unique, and digest-consistent; the
view selection and dependency revision match; and the producer reports the
requested group analysis complete. Missing or contradictory evidence is
`UNRESOLVED`; a dependency mismatch is `STALE`. No membership is inferred.

After VECTOR archive reload, certification must repeat these checks against the
reloaded network's retained mapping and solve context. Archive persistence is
outside ATLAS; the corrected boundary-marker validation permits legacy
synthetic boundary rows with `system_num=None` while still rejecting missing
physical generator IDs.

## Phase III Round 3 source-identity rule

`source_network_id` is authoritative only when it is explicitly supplied by the
producer and survives archive reload as that field. `archive_id` identifies a
serialized artifact/revision, not the source-network namespace, and must never
be copied into `source_network_id`. A missing or blank producer ID is therefore
an unresolved source identity: NetworkView may retain a provisional bookkeeping
namespace, but a Results adapter must not certify or reuse the resulting
analysis. Solve-context digests and provisional `network:<digest>` names
likewise remain dependency evidence only; they are not source-network
identities.

## A16 authoritative source identity

The canonical Results API is:

```python
SharedNetworkIdentity(
    source_network_id="producer:<opaque-id>",
    source_identity_status="AUTHORITATIVE",
    ...,
)
```

`source_identity_status` is carried into `ResultIdentity` and the serialized
`Provenance.details["shared_network_identity"]` payload. New producer handoffs
must set `AUTHORITATIVE` only from an explicit FORGE/NEXUS source-network ID;
the existing F7 `network:<solve-digest>` fallback is `PROVISIONAL`. Omitting
the field is backward-compatible construction but is `UNRESOLVED` at every
Results certification boundary.

`ResultsDependency.dependency_digest` includes the authority status, so a
selection can change its dependency while retaining solve identity, and an
authority or source-identity change cannot reuse a certified cache. A cached
dependency with contradictory authority status is `UNRESOLVED`; ordinary
revision, frame, geometry, solvent, or retained-generator mismatches remain
`STALE` where the existing dependency rules require it. `archive_id`,
`source_revision`, and `source_digest` remain artifact/dependency evidence and
cannot substitute for producer authority.

FORGE F15 must publish the explicit source-network ID and its authority
attestation; NEXUS N14 must pass that attestation into
`SharedNetworkIdentity.source_identity_status`. Legacy archive reloads may
retain a nonempty internal namespace for inspection, but without the
attestation their Results remain `UNRESOLVED`.
