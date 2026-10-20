# Shared-network group analysis prototype

Status: bounded, opt-in prototype. `Group.build()` does not use this path.

`vorpy.src.group.shared_view.analyze_group_view()` consumes an existing,
complete `NetworkView` and returns only group-owned volume plus exposed surface
area. Two views can therefore read one solved full-system network while
retaining the source network namespace, solve-universe identity, global system
generator IDs, competitor IDs, and completeness evidence. The function reads
immutable view snapshots; it does not call `Network.build()` or mutate the
source network.

## Legacy adapter boundary

`Group.get_info()` is not currently a safe direct consumer of `NetworkView`. It
expects a mutable `Network` shape with pandas `balls`, `surfs`, `edges`, and
`verts` tables, local topology indices, and cached layer/curvature state. It
also writes `Group.layer_*` and curvature result fields. Remapping immutable
view rows into a synthetic `Network` would risk losing full-system competitor
incidence, generator ownership, boundary rows, and edge/vertex curvature
dependencies.

The minimal future adapter contract is therefore an explicit FORGE-owned
view/materialization protocol that provides, without rebuilding geometry:

1. source network and solve-context identity;
2. global-to-local generator ownership, including boundary generators;
3. complete cell/surface/edge/vertex incidence and source IDs; and
4. the producer-owned curvature/edge caches required by the legacy pipeline.

NEXUS should not implement that materialization locally. The next integration
step is to coordinate this contract with FORGE, then add an explicit opt-in
adapter for the legacy pipeline after its numerical outputs are compared with
an isolated `Group.get_info()` run.

## N2 Cambrin checkpoint

The available F8 artifact is `output/cambrin_33/cambrin.vpy`. It loads without
calling `Network.build()` and contains AW geometry with no boundary generators,
674 complete protein cells, 7,674 vertices, 14,503 edges, and 7,495 surfaces.
For global generator [2], an independent raw-table reference is:

- volume: `13.572938281053846`
- exposed area: `38.701579174896224` across 21 incident surfaces

The current `NetworkView` fails closed before the shared-view adapter can run:
the source has 8,516 water competitor rows with `complete=0`, and the selected
view reaches competitor generator `2037`. The cached solve therefore proves
geometry and source immutability, but it does not satisfy the current view
completeness contract for competitor-only rows. The focused regression records
this exact blocker rather than marking water cells complete or copying the
network. FORGE must first agree how competitor-only generators are represented
in view completeness; then the same test can compare the adapter to the two
`[2]` references at `1e-12`.

## N3 molecule-only integration audit

The minimal integration path is already present at the view boundary:

1. FORGE constructs one solved `Network` whose `balls` and incidence tables
   contain only the participating molecular generators. Input solvent may stay
   in `System.balls`, but its global IDs must not be inserted into the solve
   tables or remapped over. Each retained row keeps its input-global
   `system_num` and stable generator ID.
2. NEXUS obtains each group selection as those explicit global IDs and calls
   `Network.make_view(group_ids=...)` once per group. The solve context maps
   local topology indices back to the retained global IDs; absent solvent IDs
   remain absent. `Group.ball_ndxs` is a System-level selection, not a compact
   solve-row position, so this conversion must use the producer's retained
   `system_num`/stable-ID map rather than `.iloc` positions.
3. NEXUS passes each read-only view to `analyze_group_view`. No `Group.build()`
   or second `Network.build()` is required, and the source tables stay shared.

The opt-in `system_ids_for_group()` / `make_group_view()` helpers now enforce
the System-ID-to-solved-universe check at this seam and reject an omitted
solvent ID instead of filtering or remapping it implicitly.
`shared_network_identity()` additionally requires F10 provenance before
projecting a view into the Results identity contract; legacy sources remain
unresolved rather than receiving an inferred solvent mode.

F10 now makes `Group.make_net()` producer-compatible: it uses
`prepare_network_geometry()` and attaches `Network.solvent_provenance`, compact
local-to-System membership, stable generator IDs, and cell-completeness
evidence. The shared-view prototype remains separate from `Group.build()`
analysis routing. `NetworkView` should continue to require complete source rows
and must not infer or filter water competitors.

The current checkout exposes the authoritative F10 fields on `Network`. The v3
archive field registry does not list separate solvent fields, but a bounded
round trip restores the provenance through the serialized settings and ball
membership columns. This remains provisional until VECTOR validates the same
behavior on a real F10 archive and confirms that legacy archives remain
explicitly unresolved rather than inferred.

### Bounded real-molecule checkpoint

Using the existing `vorpy/data/EDTA.pdb` with F10's default producer mode
completed one small solve without a second build: 32 molecular generators,
`solvent_mode=MOLECULE_ONLY`, no participating solvent IDs, and 8/32 complete
cells. `make_group_view()` mapped System ID `[0]` correctly, then
`NetworkView` rejected the view at incomplete competitor `5`. This is real
producer/mapping evidence, but not a certified group area/volume result.

### N8 bounded EDTA shell acceptance

One AW solve with the explicit F11 shell configuration produced 274 generators:
32 molecular generators and 242 virtual boundary generators. All 32 molecular
cells were complete. Two immutable views from that one source produced ATLAS
`GroupResult` values through `AnalysisResult.typed_results`, using an explicit
producer source ID without creating `archive_id`; the source balls, surfaces,
and System atoms were unchanged and no second `Network.build()` was called.
The separate Workbench regression omits that source ID and verifies the
provisional/`UNRESOLVED` path.

Independent sums over the source cell/surface tables matched the view results
to `1e-12`:

| selection | volume (Å³) | total exposed area (Å²) | boundary-dependent area (Å²) | physical interface area (Å²) |
| --- | ---: | ---: | ---: | ---: |
| System `[0]` | 5.676116886923854 | 19.53014188145206 | 0.0 | 19.53014188145206 |
| System `[1]` | 14.170714484103758 | 46.75203064184723 | 2.202876203619455 | 44.54915443822777 |

The typed area metric is explicitly the exposed group-view area; the split
above prevents shell-clipped area from being reported as a solvent-free
physical molecular interface. The opt-in `attach_group_view_results()` seam
also preserves one shared solve identity and distinct selection identities.

The shared checkout's archive repair now permits reload of this real EDTA
network. The focused regression passes the same explicit producer source ID
to both analyses and compares solve identity, selection identity, provenance,
area, and volume; the build counter remains exactly one across solve, views,
save, and reload. The archive's separate ID is assigned only during save.

### N9–N11 production and identity contract

`VorPySolveSettings.typed_group_selections` is the smallest producer-side
handoff. It accepts a named mapping of System-global ID tuples, or the
backward-compatible tuple form such as `((0,), (1,))`. The default is empty,
so legacy backend execution does not construct typed Results. When configured,
the backend calls `attach_group_view_results()` once after the existing solve
and stores its tuple of ATLAS Results in `AnalysisResult.typed_results`.
`AnalysisResult.typed_selection_system_ids` preserves the label-to-System-ID
mapping. Group analysis creates immutable views and never calls
`Network.build()`.

Identity has three distinct roles:

- `source_network_id` is a producer namespace. An explicit
  `Network.source_network_id` or view argument is authoritative and may be
  used on a fresh Network without assigning `archive_id`.
- `solve_universe_id` is the canonical F6 solve-context digest. It is available
  on a fresh Network and is shared by all selections of that solve.
- `archive_id` is the persisted object key assigned by the archive layer. It
  is not a substitute for manufacturing source identity during a live solve.

If no producer source ID or existing archive ID exists, `NetworkView` uses
`network:<solve-context-digest>` as an explicitly provisional source namespace.
Group metrics are still calculated, but typed Results remain `UNRESOLVED`; the
producer does not invent an archive UUID or certify incomplete identity
evidence. A fresh Network with an explicit producer source ID can certify when
F10/F12 and view completeness validate. If a separate source ID was not
archived, a reloaded Network may resolve its source namespace from its
persisted `archive_id`; this changes source identity and must not reuse cache
lineage from the provisional fresh solve. Cache dependencies record the
solve-context digest as `source_digest`; `archive_id` is not used as the
default source revision.

Each typed result carries the same source/solve identity, its own
`NetworkView.selection_id`, F10 solvent provenance, and F12 retained-System
indices, stable generator IDs, and digest. Missing F12 provenance fails
closed. Retained-generator fields belong to `SharedNetworkIdentity`, recorded
in `typed_result.provenance.details["shared_network_identity"]`; they are not
fields of `ResultIdentity`. The full F12 producer mapping is retained under
`source_retained_generator_provenance`. Cached dependencies passed to
`attach_group_view_results()` are keyed by `selection_id`; ATLAS validates
their lineage without rebuilding geometry, though this NEXUS adapter still
recalculates the selected metrics.

PRISM contract: pass the named Workbench selection mapping as
`VorPySolveSettings(typed_group_selections={label: (system_global_id, ...)})`.
When an authoritative producer namespace is available, also pass
`source_network_id="<producer-id>"`; it is optional and defaults to `None`.
NEXUS never generates a source or archive ID, and the existing
`Network.source_network_id` is honored when the setting is omitted.
The backend returns the single authoritative payload at
`AnalysisResult.typed_results: tuple[object, ...] | None`: `None` means no
typed producer handoff, otherwise it contains one ATLAS `GroupResult` per
selection, in mapping order. The keys become `group_id`; values remain
System-global IDs and are also recorded in
`AnalysisResult.typed_selection_system_ids`. Consume the tuple read-only and
preserve each result's status and provenance. Do not rebuild Results in
presentation code. An absent payload preserves legacy display behavior. If
neither the setting nor the fresh Network supplies a producer source ID, its
provisional status is expected; PRISM must not promote it.

VECTOR follow-up: v3 archive fields currently persist `archive_id` but not a
separate `source_network_id`. To preserve an explicit producer namespace across
save/reload, persist that optional field independently. When absent, preserve
current fallback behavior: use an existing persisted `archive_id`, otherwise
remain provisional. NEXUS does not change archive ownership or infer the
missing producer field.

Cache note: the ATLAS adapter validates cached dependency identity and status;
the NEXUS helper still computes group metrics when called. Reusing cached
metric values without analysis needs an ATLAS-owned cached-Result contract and
is not implied by a matching dependency digest.
