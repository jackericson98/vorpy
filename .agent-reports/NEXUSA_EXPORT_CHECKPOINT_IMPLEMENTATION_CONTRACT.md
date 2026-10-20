# NEXUSA — VorPy Export and Checkpoint Implementation Contract

**Status:** FROZEN FOR IMPLEMENTATION REVIEW  
**Owner:** NEXUSA  
**Date:** 2026-10-19  
**Scope:** export layout, checkpoint persistence, CLI export planning, and
cross-agent ownership. No solver or scientific-formula change is authorized by
this document.

This is the authoritative follow-up to
`.agent-reports/NEXUSA_EXPORT_ARCHITECTURE_AUDIT.md`. Where an older report
describes a different default, this contract wins. In particular:

* every normal tier, including `tiny`, updates one shared checkpoint;
* default interface output does not write standalone Results/provenance CSVs
  or classic interface scheme logs;
* `all` does not silently enable alpha, dual, detailed-water, experimental, or
  obsolete exports;
* export presets select files only. They never request geometry or analysis.

No implementation, molecular solve, output migration, Git operation, or
existing-output rewrite was performed to publish this contract.

## 1. Non-negotiable invariants

1. Calculation scope and export scope are separate.
2. A requested interface never implies a system or standalone group solve.
3. A requested export never implies a calculation, optional analysis, or
   completeness promotion.
4. There is one canonical checkpoint per molecular-system output root, never a
   checkpoint per group or interface.
5. Physical A/B interface identity is unordered. A stable reference orientation
   is stored once; display orientation is derived and does not duplicate
   geometry.
6. Only direct physical A/B facets contribute to primary interface area.
   Internal and supporting surfaces remain separately labelled.
7. No-contact is a valid empty interface, not a geometry error. Unavailable
   curvature, topology, water, or alpha values remain unavailable rather than
   becoming zero.
8. Solvent membership and solvent-cell completeness are different facts.
   `SOLVENT_CONTEXT` may support locally certified direct contact geometry.
   `FULL_SOLVENT` remains fail-closed for complete solvent-cell-dependent
   quantities.
9. Existing source identity, generator identity, competitors, facet/edge/
   vertex incidence, contact identity, and completeness evidence are never
   discarded to make an export smaller.
10. The deferred requirement for reconstructed dual-contact surfaces—complete
    DNA-side and ligand-side triangulated patches with natural outer
    boundaries—is not implemented or assigned to a default tier in this cycle.

## 2. Authoritative output layout

For a source with system name `<system>` and original basename `<input-name>`:

```text
<system>/
  <input-name>                         # byte-for-byte source copy
  <system>_checkpoint.vpy              # one shared checkpoint
  system/
  groups/
    <stable-group-id>/
  interfaces/
    <stable-interface-id>/
```

`<input-name>` retains its actual extension and bytes. The implementation must
record its size and SHA-256 in checkpoint provenance and must not rename a CIF,
MOL2, GRO, or other supported source to `.pdb`. If a run is loaded solely from
an existing checkpoint and the original source file is unavailable, the
self-contained checkpoint remains sufficient for reload; no fabricated source
copy is created.

The root `system/` directory is created for system metadata. `groups/` is
created only when a group or interface request needs group representations.
`interfaces/` is created only when an interface is requested. Empty top-level
scope directories are not scientific outputs and need not be created.

The checkpoint is root-owned. Scope directories contain derived, scope-owned
files only. A scope directory must not receive a copy of a full-system network
merely because it was selected as a view.

### Scope-to-calculation rules

| Request | Geometry requested | Scope directories and ownership |
|---|---|---|
| no `-g`, no `-i` | complete system network | `system/` gets system network outputs |
| `-g A` | standalone A network | `groups/A/` gets a group network |
| `-i -g A -g B` | direct A/B interface only | `system/` metadata, minimal `groups/A/` and `groups/B/`, and one interface folder; no standalone A/B networks |
| `-i -g A -g B -c groups` | direct interface plus standalone A and B networks | interface folder plus calculated `groups/A/` and `groups/B/` |
| `-i -g A -g B -c system` | one full system network plus direct interface view | `system/` network, selected group outputs only when valid source-derived group data is explicitly required, and one interface folder |
| `-c system` | complete system network | `system/` network only unless other scopes are explicitly requested |
| `-c groups` | requested standalone group networks | only requested `groups/<id>/` folders |

For interface-only runs, `groups/A/` and `groups/B/` are placeholders, not
network results. They contain `A.pdb` and `B.pdb` with original atom
identities. The interface folder contains `interface_atoms.pdb`, the union of
atoms participating in direct A/B contact. A no-contact interface writes a
valid empty `interface_atoms.pdb` and explicit empty status.

`-i -c system` may reuse valid complete system-derived group views. It must not
rebuild equivalent geometry or claim a group result from incomplete source
records. `-i -c groups` is the explicit request for independent group solves.

## 3. Scope ownership and folder conditions

### System

`system/` owns system-level settings, source provenance, network metadata,
system-network meshes and system-network classic scheme logs. It does not own
interface contact tables, group surroundings, or derived dual/alpha products.

The original source copy is root-owned, not duplicated into `system/` by
default. A system network is written only when the calculation plan requests
one or an already valid cached system network is explicitly exported.

### Groups

`groups/<stable-group-id>/` owns a requested standalone group network and its
group-specific representations. A group folder created only to support an
interface is limited to the source-identity PDB required by the interface
contract. It must not contain `surfs.off`, atom-cell meshes, group logs, or
standalone analysis unless a group calculation was requested.

Group network logs are retained for solved group networks at every tier that
exports a solved network. This is separate from interface default behavior.

### Interfaces

`interfaces/<stable-interface-id>/` owns direct A/B physical geometry, contact
atom representations, interface information, and explicitly requested
interface analysis. The identity is derived from the system identity and the
unordered pair of stable group identities. A/B and B/A therefore resolve to
the same folder and geometry. The stored reference orientation determines the
sign of mean curvature; reversing display perspective changes the displayed
normal and mean-curvature sign only. Gaussian curvature is invariant.

Supporting geometry, alpha/dual representations, shell products, and detailed
water products are separate components. None is folded into direct interface
area.

## 4. File contract by tier

The tables use these marks:

* **R** — required when the scope exists and the source/result is available;
* **C** — required only when that scope has a calculated network;
* **O** — explicit optional component, never implied by the preset;
* **L** — legacy compatibility output, opt-in;
* **—** — not written by the preset.

An incomplete or unsupported result still gets its scope info/status record.
Geometry files contain only records whose own completeness/evidence gates pass.
The checkpoint retains records that are absent from a separate export.

### 4.1 System files

| Tier | Required system files |
|---|---|
| `tiny` | `system_info.txt` (R). If a system network is calculated: `network_info.txt`, `edges.off`, and the method log `<net_type>_logs.csv` (C). |
| `small` | `tiny` plus `surfs.off` for a calculated system network (C). System classic log remains present for a calculated system network. |
| `medium` | `small` plus `verts.off` and the existing source/PyMOL representations (`set_atoms.pml` and any selected source representation) when requested by the system writer (C). |
| `large` | `medium` plus supported shell meshes (`shell_surfs.off`, `shell_edges.off`, `shell_verts.off`) and complete-cell atom meshes under `atoms/` when explicitly available (C/O by feature). |
| `all` | `large` plus supported separate primary meshes, external/surrounding representations, and diagnostics that are explicitly part of the primary system writer. It does not add `logs.csv`, alpha, dual, detailed water, or experimental products automatically. |

The root source copy and shared checkpoint apply to every row and are listed
separately in Section 5. If no system network was requested, no system mesh or
system classic log is created merely because an interface is being exported.

### 4.2 Calculated group files

| Tier | Required group files |
|---|---|
| `tiny` | `group_info.txt` (R/C; current `network_info.txt` remains a compatibility alias), `edges.off`, and `<net_type>_logs.csv` for a calculated group network (C). |
| `small` | `tiny` plus `surfs.off` and `group_atoms.pdb` for a calculated group network (C). |
| `medium` | `small` plus `verts.off`, `surr_atoms.pdb`, and `surr_resids.pdb` when those validated selections exist (C). |
| `large` | `medium` plus supported shell meshes and complete-cell atom-cell meshes: `atoms/ball_<atom>_<net_type>.off`, with their supported `_edges.off` and `_verts.off` companions (C/O by completeness). |
| `all` | `large` plus supported separate meshes (`sep_surfs`, `sep_edges`, `sep_verts`), `ext_atoms.pdb`, and supported layer products. Legacy `logs.csv`, old vertex text files, obsolete `_net.csv`, and experimental/nonpolar products remain explicit opt-ins. |

The interface-only placeholders are an exception to the calculated-group
table: they contain only `A.pdb` or `B.pdb` and a provenance/status record if
the implementation needs one. They never inherit group tier files.

### 4.3 Interface files

| Tier | Required interface files |
|---|---|
| `tiny` | `interface_info.txt`, `interface_edges.off`, and `interface_atoms.pdb`. The info file records canonical identity, reference orientation, contact count, direct-area status, solvent policy, and completeness evidence. |
| `small` | `tiny` plus `interface_surfs.off`. Direct physical facets only; supporting surfaces are not silently merged. |
| `medium` | `small` plus `interface_verts.off`, validated direct/support contact tables, and `surr_atoms.pdb`/`surr_residues.pdb` when explicitly requested and available. These are physical/evidence tables, not standalone Results/provenance logs. |
| `large` | `medium` plus supported curvature/topology tables and explicitly requested shell/supporting geometry. Each metric carries its own status and dependency evidence. |
| `all` | `large` plus optional legacy interface scheme logs and typed Results/provenance CSV bundles only when the `logs`/`results` component is explicitly selected. Alpha, dual, reconstructed contact surfaces, detailed water, and experimental products remain explicit add-ons even with `all`. |

The default interface tiers do **not** write `logs.csv`,
`logs.compact.csv`, `logs.provenance.json`, or `<net_type>_logs.csv`. The
checkpoint still stores the Results, provenance catalog, dependencies, and
analysis state needed to produce them later. This prevents a large repeated
provenance CSV from becoming an interface default while preserving scientific
state.

`interface_atoms.pdb` is not a copy of all A/B atoms: it is the source-identity
union of atoms participating in accepted direct A/B contact. The selected
`A.pdb`/`B.pdb` files in the minimal group folders preserve the full selected
group representations and are distinct from the contact union.

## 5. Root source and checkpoint policy

### Required behavior

* Every normal export tier (`tiny` through `all`) creates or atomically updates
  `<system>_checkpoint.vpy`, unless checkpoint writing is explicitly disabled.
* The checkpoint is self-contained for reload without the source input. It
  includes original generator/atom/residue/chain identities, geometry settings,
  source provenance, solver/method settings, network and scope identity,
  partial geometry, Results, dependency provenance, analysis caches, and
  completeness/status evidence where available.
* Checkpoint update uses a same-directory temporary file, flush/sync, validation
  of the completed archive, and atomic replace. A failed update leaves the
  previous valid checkpoint intact.
* A checkpoint save never creates a second archive per group or interface.
* Normal automatic saves do not emit a size warning. `all` may report the
  projected number of files and bytes for derived output. In noninteractive
  mode the warning is informational and the requested export continues; it
  never silently changes the preset or calculation.
* Older archive versions remain loadable. No format migration or rewrite is
  performed automatically. An explicit, separately approved migration may
  produce a new checkpoint after validation.
* Incrementally solving a previously unsolved Group C is a later feature. The
  checkpoint must preserve partial records now, but this contract does not
  authorize an incremental-merge solver.

### Checkpoint contents that are not separate files

The checkpoint is the source of truth for typed Results, full provenance,
cache/dependency identities, partial records, and analysis state. Omitting a
standalone log or table from a tier must not remove it from the checkpoint. A
missing dependency, stale digest, unresolved geometry, or incomplete solvent
cell is represented explicitly and causes dependent reload/export requests to
fail closed or remain unavailable.

## 6. Canonical component vocabulary and aliases

The planner operates on canonical component names. Aliases are accepted only
at the CLI boundary and normalize before precedence is applied.

| Canonical component | Existing/accepted aliases | Scope and meaning |
|---|---|---|
| `info` | `group_info`, `interface_info`, `system_info` when scope-qualified | concise scope metadata/status |
| `edges` | `edge`, `interface_edges` | primary edge visualization/records |
| `surfaces` | `surfs`, `surface`, `interface_surfs` | primary surface visualization/records |
| `vertices` | `verts`, `vertex`, `interface_verts` | primary vertex visualization/records |
| `shell` | `shell_surfs`, `shell_surfaces` | shell surface representation |
| `shell_edges` | `shell_edge` | shell edge representation |
| `shell_vertices` | `shell_verts`, `shell_vertex` | shell vertex representation |
| `group_atoms` | `atoms` in group context | selected group atom PDB/records |
| `interface_atoms` | `contact_atoms` | direct A/B contact atom union |
| `surrounding_atoms` | `surr_atoms` | validated surrounding atom representation |
| `surrounding_residues` | `surr_resids`, `surr_residues`, `surrounding_resids` | validated surrounding residue representation |
| `external_atoms` | `ext_atoms` | explicit external atom representation |
| `separate_surfaces` | `sep_surfs` | separate primary surface files |
| `separate_edges` | `sep_edges` | separate primary edge files |
| `separate_vertices` | `sep_verts` | separate primary vertex files |
| `atom_surfaces` | `atom_surfs` | atom-cell surfaces; complete-cell gated |
| `atom_edges` | `atom_edges` | atom-cell edges; complete-cell gated |
| `atom_vertices` | `atom_verts` | atom-cell vertices; complete-cell gated |
| `logs` | `scheme_logs`, `classic_logs`, `lgs` | classic AW/POW/PRM log; never default for interface |
| `results` | `canonical_results`, `typed_results` | typed Results CSV export; explicit |
| `provenance` | `provenance_sidecar` | compact Results provenance sidecar; explicit with `results` |
| `legacy_results` | `legacy_logs`, `logs_csv` | full legacy canonical `logs.csv`; explicit only |
| `contact_tables` | `contacts`, `contact_table` | direct contact/evidence tables |
| `curvature` | `curv` | curvature tables/results; analysis must already exist |
| `topology` | `topo` | topology/incidence tables; status-gated |
| `water` | `solvent`, `buried_water` | solvent/buried-water products; strict dependency gates |
| `dual` | `apollonius`, `alpha_dual` | derived dual-complex products; explicit |
| `alpha` | `alpha_complex` | derived alpha products; explicit |
| `visualization` | `viz`, `pymol` | supported visualization helpers |
| `source_structure` | `pdb` in system context | selected molecular/source representation; the unchanged root source copy remains canonical |
| `pymol_atoms` | `set_atoms` | PyMOL atom-selection helper |
| `mol` | `mol2` | legacy/API molecular writer, explicit only where a working writer exists |
| `cif` | `mmcif` | legacy/API structure writer, explicit only where a working writer exists |
| `xyz` | — | legacy/API coordinate writer, explicit only where a working writer exists |
| `txt` | `text` | legacy/API text writer, explicit only where a working writer exists |
| `checkpoint` | `archive`, `network_archive`, `save_network` | one shared root checkpoint |
| `source` | `input`, `pdb` in legacy contexts | root byte-preserved source copy; not a geometry request |
| `timing` | `timing_report`, `profile` | stage timing metadata; no calculation request |

Scope-qualified names are preferred where `info`, `atoms`, `edges`, or
`surfaces` could otherwise be ambiguous. The old `network_info.txt` filename
remains a compatibility alias during migration; it is not a second canonical
scientific artifact.

The following are not canonical default components: `nonpolar`, `mol`, `cif`,
`xyz`, `txt`, old `aw_verts.txt`/`pow_verts.txt` text dumps, `_net.csv`,
experimental contact OBJ products, and reconstructed dual-contact surfaces.
They may remain explicit legacy/API paths only where an existing writer is
known to support them.

## 7. CLI planning grammar and precedence

The calculation planner and export planner are separate state machines.

### Calculation grammar

The existing selection syntax is preserved:

```text
-g <group-selector> [-g <group-selector> ...]
-i
-c system|groups|interface
```

The ordered `-g` selections remain ordered for display/reference orientation.
The canonical interface identity is normalized as an unordered pair after
selection. `-i` without `-c system` is interface-only. `-c system` is the
explicit full-system opt-in. `-c groups` is the explicit standalone-group
opt-in. No export token can alter this calculation plan.

### Export grammar

Repeated `-e` blocks are accepted. Each block is a left-to-right sequence of
operations separated by the existing `and` convention where that convention
is already accepted:

```text
-e <preset-or-component> [add|no|remove|only <component> ...]
-e <same grammar>                         # repeatable
-e dir <path>
-e format <format>
-e checkpoint|no_checkpoint
```

Canonical presets are `tiny`, `small`, `medium`, `large`, and `all`.
`micro` is an internal/legacy alias for `tiny`. `no_archive` and
`no-network` remain accepted compatibility aliases for `no_checkpoint` only
where the current parser already accepts them; they must not mean “do not
solve the network.”

Precedence is deterministic:

1. Start with the default `small` selection for the requested scope.
2. Expand each preset at the position where it occurs.
3. Apply `only` by clearing the current export set for that scope, then add its
   named components.
4. Apply `add`, `remove`, and `no` in token order.
5. Process repeated `-e` blocks in command-line order.
6. For each component, the last explicit operation wins. A later preset is an
   explicit reset/expansion and a later `no`/`remove` overrides it.

Examples:

```text
-e small add surfs
-e large no shell
-e tiny add edges
```

The final state contains the `tiny` baseline and `edges`, with the preceding
explicit operations resolved in order; it never triggers a new solve. A
component unavailable from the cached result is reported as unavailable, not
fabricated or silently replaced with a larger calculation.

`dir`, `format`, `checkpoint`, `no_checkpoint`, and timing controls are
planning controls, not scope selectors. A single `dir` must resolve to the one
output root; it must not create scope-local archives.

## 8. Legacy migration and compatibility

Migration is reader-first and additive.

| Existing artifact/behavior | Contract treatment |
|---|---|
| Existing `.vpy` archive versions | Load read-only through the current reader; no automatic rewrite. Explicit migration is a separate approved operation. |
| Existing `.vorpy/checkpoints/*.geometry.vpy` | Keep readable; do not create another such file for new output. Explicitly validate and re-save into the shared root checkpoint when migration is authorized. |
| Existing flat `network_info.txt` | Continue reading and optionally emit as a compatibility alias; canonical new names are scope-specific info files. |
| Existing flat group/interface directories | Do not delete, move, or overwrite. New runs use the root/scope hierarchy. A manifest/migration report may map old paths to new scope IDs. |
| Existing `logs.csv` | Readable and writable only through explicit `legacy_results`/compatibility selection. It is never a default interface file. |
| Existing `<net_type>_logs.csv` | Retained for solved system and group networks; interface use is explicit. |
| Existing compact Results/provenance pair | Read/write only when `results` and `provenance` are selected; digest/reference failures remain fail-closed. |
| Existing alpha/dual/water/nonpolar outputs | Readable as derived products where supported; never automatically regenerated by a preset. |

No migration may change geometry, atom identities, source competitor policy,
interface identity, completeness flags, or scientific formulas. A future
manifest may record old-to-new paths, but it is not permission to copy every
legacy duplicate into every scope.

## 9. ATLAS/VECTOR checkpoint boundary

The shared boundary is a typed checkpoint document and a versioned export plan.
Neither agent may change the other’s contract silently.

### ATLAS owns

* canonical object identities for system, generator, group, interface, and
  Results records;
* typed Results fields, status vocabulary, dependencies, provenance catalogs,
  source/solve/frame/selection identity, and cache validity;
* archive payload schema, format versions, compatibility rules, and explicit
  migration validation;
* preservation of partial records and fail-closed dependency checks;
* the scientific rule that omitted files do not imply omitted Results.

ATLAS does not own CLI grammar, output directory selection, mesh writers,
classic log layout, or primary geometry construction.

### VECTOR owns

* scope-directory creation and filename mapping from an accepted export plan;
* OFF/PDB/table/classic-log/compact-log serialization;
* source-copy emission, compatibility aliases, and output file accounting;
* archive byte encoding/decoding and atomic temporary-file replacement, using
  the ATLAS schema and validators;
* export-stage and archive-save timing labels, without double-counting mesh or
  archive work.

VECTOR must not invent Results fields, mark geometry complete, drop competitors,
or decide whether a metric is scientifically certified. `network_archive.py`
is a shared seam: ATLAS approves schema/validation changes, and VECTOR owns
the I/O implementation only after that approval.

PRISM requests the plan and reports progress. FORGE supplies geometry and
evidence. GAUSS defines metric certification. NEXUSA owns the cross-agent
acceptance gate.

## 10. Status, completeness, and no-contact semantics

The export contract preserves both geometry evidence and analysis status.

### Geometry evidence

`COMPLETE` means the requested local construction evidence passes for the
stated scope. `INCOMPLETE_CONSTRUCTION` means required records or incidence
are missing. `UNBOUNDED` is a proven geometric condition, not a construction
failure. `UNSUPPORTED_ANALYTIC_GEOMETRY` means the analytic branch cannot
certify the requested feature. A producer flag alone cannot promote a failed
independent audit.

### Analysis status

`CERTIFIED` is supported, finite, and complete for the named metric.
`PARTIAL` contains a named supported subset and coverage/unresolved evidence.
`UNRESOLVED` means the request was attempted but required evidence failed.
`NOT_CALCULATED` means it was not requested or not attempted.
`EMPTY_NO_CONTACT` is a certified empty direct A/B interface when the contact
test found no accepted physical facets.

For no contact:

* `interface_info.txt` records the stable identity, selected groups, zero
  direct-facet/contact counts, and `EMPTY_NO_CONTACT` status;
* `interface_atoms.pdb` and primary interface OFF files are valid empty
  artifacts;
* direct area is a certified finite zero only when the no-contact decision is
  itself certified;
* curvature, topology, signed normal, and solvent-cell metrics are
  `NOT_CALCULATED` or `UNRESOLVED` as appropriate, never numeric zero;
* no supporting or invented contact surface is emitted.

For wet runs, direct contact may use `SOLVENT_CONTEXT` evidence while the
producer retains the complete `FULL_SOLVENT` competitor universe. Complete
solvent-cell metrics, buried-water products, closed-cell volume/topology, and
AW results requiring complete participating cells retain strict gates.

## 11. Acceptance and test matrix

This contract is accepted only after bounded tests cover the following. No
real EDTA, 1BNA, p53tet, or large molecular solve is part of this gate.

| Area | Required bounded acceptance |
|---|---|
| Routing | `-i` produces only direct interface geometry; `-i -c groups` adds standalone group solves; `-i -c system` builds/reuses one system network and the interface; no implicit system/group solve. |
| Ordering | A/B and B/A share identity and geometry; selected order remains the display/reference orientation. |
| Folder scope | Interface-only creates minimal `groups/A`, `groups/B`, `interface_atoms.pdb`, and no unrelated group/system meshes/logs. Explicit system/groups create only requested scope outputs. |
| Tiers | Each tier has its exact scope inventory; every calculated scope has at least one visual representation; default interface excludes standalone Results/provenance/classic logs. |
| CLI modifiers | Repeated `-e`, `add`, `remove`/`no`, `only`, `-e small add surfs`, and `-e large no shell` resolve with last-explicit-operation precedence. |
| Calculation isolation | Preset changes do not call geometry, optional analysis, or full-system planning. Cached export does not call `Network.build()` or analysis. |
| Source copy | Synthetic PDB/CIF/MOL2/GRO-like inputs retain exact basename, extension, bytes, and digest; source files are never overwritten. |
| Checkpoint | All five tiers create/update one shared checkpoint; reload works without the original source; atomic failure preserves the prior valid checkpoint; old archive fixtures remain readable. |
| Identity/provenance | Typed Results, dependencies, partial records, source identity, and cache digests survive checkpoint round-trip even when standalone logs are omitted. |
| Solvent | Explicit water competitors remain in `SOLVENT_CONTEXT` and `FULL_SOLVENT`; direct contact geometry is invariant to scope; solvent-cell-dependent metrics fail closed when cells are incomplete. |
| Methods | Bounded AW/POW/PRM direct geometry, area, contacts, topology/evidence and scope equivalence; AW complete-cell gate remains distinct; PRM typed Results limitation remains a tracked unresolved status. |
| Completeness | Complete selected facets export; incomplete unrelated solvent cells do not invalidate independent direct planar evidence; incomplete required local features remain partial/unresolved. |
| No contact | Empty interface artifacts and statuses are explicit and contain no invented facets. |
| Logs | Solved system/group networks retain classic logs; default interface does not; explicit `logs`, `results`, `provenance`, and `legacy_results` produce only the selected forms. |
| Progress | PRISM reports the actual requested calculation and separates load, solve, analysis, meshing, export, and archive-save timing without double counting. |
| Preservation | Existing output directories and input files remain unchanged; no destructive Git operation is used. |

The current bounded scientific evidence supports the status rules above. The
PRM exact geometry producer is supported, but the typed PRM Results handoff is
still a tracked blocker and must not be hidden by an export preset.

## 12. Estimated output counts and size expectations

Counts are planning ranges, not scientific requirements; per-cell atom meshes,
diagnostics, and optional tables are variable.

| Scope/tier | Typical fixed files | Variable additions | Expected size class |
|---|---:|---|---|
| system/tiny | 2–4 including info and network edge/log when solved | source/checkpoint | KB to MB plus checkpoint |
| system/small | 3–6 | surface mesh and network log | MB to tens of MB |
| system/medium | 5–9 | source/PyMOL representations | MB to tens of MB |
| system/large | 8–15 | shell and complete-cell atom meshes | tens to hundreds of MB |
| system/all | 10+ | supported separate/legacy primary products | unbounded by preset |
| group/tiny | 2–4 | network log when solved | KB to MB |
| group/small | 4–6 | atom PDB and surface mesh | MB-scale |
| group/medium | 8–12 | surrounding representations | MB to tens of MB |
| group/large | 12+ | approximately three meshes per complete atom cell | high |
| group/all | 16+ | layers and separate primary products | very high |
| interface/tiny | 3 | empty/direct edges and contact atom PDB | KB to low MB |
| interface/small | 4 | direct surface mesh | MB-scale |
| interface/medium | 8–12 | contact tables, vertices, surrounding files | MB to tens of MB |
| interface/large | 12+ | certified curvature/topology/support products | high |
| interface/all | 14+ | explicit log/results/legacy products | potentially dominated by CSVs |
| shared root | 1 checkpoint plus source copy | optional timing/manifest metadata | KB to hundreds of MB |

The old bounded evidence explains the range: a preserved interface directory
had 14 files and about 51 MB, with one `logs.csv` accounting for about 47 MB;
the physical/contact subset without that CSV was about 3.85 MB. The contract
therefore makes the large legacy canonical CSV opt-in.

## 13. Specialist file ownership and change gates

No specialist may edit another specialist’s implementation area without an
explicit seam review. The following ownership is the implementation map for
this contract:

| Owner | Files/areas | Allowed responsibility |
|---|---|---|
| **PRISM** | `vorpy/src/command/vpy_cmnd.py`, `vorpy/src/command/export.py`, calculation-plan/CLI tests | calculation target planning, ordered group selection, `-i` routing, export grammar/precedence, progress, timing aggregation, all-preset warning |
| **VECTOR** | `vorpy/src/output/**`, output layout and visualization writers, archive byte I/O in `vorpy/src/io/network_archive.py` after ATLAS schema approval | scope folders, source copy, OFF/PDB/table/log serialization, compatibility aliases, export/archive timing, atomic checkpoint write plumbing |
| **ATLAS** | `vorpy/src/results/**`, scientific adapter/identity and archive-field contracts, schema review for `vorpy/src/io/network_archive.py` | typed Results, provenance, dependency/cache validity, checkpoint payload/version/migration rules, partial/status preservation |
| **FORGE** | network construction/topology and interface extraction modules | primary geometry, source/context competitor universe, direct interface views, facets/edges/vertices, completeness and solvent evidence; no export filenames |
| **GAUSS** | scientific analysis/status modules and acceptance fixtures | area/contact/curvature/topology definitions, method-specific gates, certified/partial/unresolved semantics; no CLI or writer behavior |
| **NEXUSA** | `.agent-reports/NEXUSA_*`, root `results.txt`, integration acceptance tests | contract, cross-agent coordination, bounded end-to-end tests, migration acceptance, consolidated reports |

The shared `network_archive.py` seam requires both ATLAS and VECTOR review:
ATLAS approves payload semantics and validation; VECTOR implements encoding,
atomic replacement, and path plumbing. PRISM calls the API but does not embed
archive schema. FORGE and GAUSS do not edit export serialization.

## 14. GO/NO-GO handoff

**GO for implementation planning:** PRISM may implement the frozen calculation/
export plan grammar and progress contract; VECTOR may implement scope-aware
serialization and the five-tier file plan after the ATLAS checkpoint seam is
reviewed; ATLAS may define the checkpoint schema/version and validation; GAUSS
may review each status/completeness acceptance row.

**NO-GO for molecular validation or broad migration:** do not launch 1BNA,
EDTA, p53tet, or another full solve; do not rewrite existing output trees; do
not auto-migrate old archives; do not implement reconstructed dual-contact
surfaces; do not mark PRM typed Results as resolved.

The next integration gate is a bounded synthetic export/checkpoint matrix,
followed by user approval of one isolated molecular pilot command. This
contract itself authorizes no pilot.

