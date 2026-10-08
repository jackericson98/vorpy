# Scientific result contract, version 1.0

The `vorpy.src.results` package is an isolated read-only adapter layer. It owns
no scientific calculations, live file imports, export paths, session persistence,
GUI behavior, presets, or dual-side construction. Importing the package imports
only the standard library. Scientific producers remain authoritative.

## Composition and consumer API

`SystemResult`, `NetworkResult`, `GroupResult`, and `InterfaceResult` are thin,
frozen result types over the same `ResultModel` fields. Their only specialization
is validation of identity kind. They compose `ResultIdentity`, `Provenance`,
`Environment`, `Quantity`, `Topology`, `MeanCurvature`, `GaussianCurvature`,
`ContactRecord`, selection quantities, and versioned `RepresentationExtension`
records. Collections are immutable snapshots, including nested extension maps.

`adapt_interface_analysis(analysis, identity)` consumes an existing Python
`InterfaceGeometryAnalysis` cache. It reads its dictionaries, not Network
properties. It does not call `analyze_geometry`, `summarize_physical_interface`,
filtration queries, selection helpers, mesh/dual builders, or curvature kernels.
It never reverses normals; it reads the explicitly requested cached side.
When cached identity metadata is present, the adapter validates it against the
supplied identity, including normalized network scope, partition, alpha, probe,
and solvent context. Older caches may omit identity fields; omitted fields are
reported as unavailable to validate, while known mismatches still fail.
`adapt_cached_result` accepts already produced typed blocks for the other result
types and future dual-side interfaces. `unavailable_interface` explicitly models
NOT_CALCULATED or NOT_SUPPORTED representations without fabricated geometry.

`to_dict()` returns detached JSON-safe Python state with IDs and schema version.
No production JSON/CSV/TXT writer or dictionary deserializer is included.
Snapshots are not updated when the source cache changes: reconstruct the adapter
or invalidate it as STALE. Default result-level PARTIAL means the complete
requested descriptor set is not explicitly certified by its producer. Consumers
must inspect each quantity's status and scope; a result-level status does not
override a missing descriptor or certify every metric.

## Identity and provenance

`ResultIdentity` requires kind, system, frame (possibly null for unknown/static),
network ID, network scope, environment, canonical partition, representation, and
radius-configuration ID. Canonical environments are `dry` and
`solvent_competing`; legacy boolean solvent-context metadata is normalized at the
adapter boundary. Interface identity additionally requires interface ID,
ordered group identities, side, and orientation. Alpha values require an explicit
method/convention and units. Probe values require units. System namespace can be
specified separately from the human-facing system name.

Network scope is **system** or **interface**, regardless of directory names or
whether a complex is called full. Partition values are **aw**, **power**, or
**power**; the physical adapter normalizes existing `pow` metadata for validation.
Persisted scope aliases `system network` and `interface/dedicated` normalize to
the canonical scope values. Interface representations are **voronoi**, **dual_a**, **dual_b**,
with sides **shared**, **A**, **B** respectively. Future dual orientation strings
are producer-defined. Physical cache orientations must be `group1 -> group2`
or `group2 -> group1`.

`result-v1:<sha256>` hashes canonical JSON of every identity field. It excludes
metric values, statuses, paths, and timestamps. Consequently, scope, frame,
partition, representation, normal direction, alpha convention/value/units, probe,
radii, and network changes distinguish results. A radius-configuration ID must
refer to exact scientific configuration (base radii, expansion/application policy,
and any overrides), preferably a persisted fingerprint supplied by the producer.
An unknown radius configuration must have an explicit stable unknown-state ID,
not an arbitrary default that merges distinct analyses.

`Provenance` records source, optional source revision/scientific-state ID,
JSON-safe details, and diagnostics references. These identify the snapshot and
numerical implementation independently of its scientific analysis identity.

### Geometry lineage

Geometry reuse is provenance, not a scientific `Status`. When present,
`Provenance.details["geometry_lineage"]` is a frozen `GeometryLineage` record
with `schema_version: 1` and one of these origins:

| Origin | Meaning |
| --- | --- |
| `DIRECT_SOLVE` | Geometry independently constructed for the reported solve universe. |
| `EXACT_REUSE` | Geometry obtained from a compatible certified solve universe and source view. |
| `MATERIALIZED_VIEW` | A representation assembled from authoritative stored primitives. |
| `UNRESOLVED` | Geometry origin or reuse compatibility cannot be certified. |

Resolved origins require deterministic `solve_universe_id`,
`geometry_source_id`, `geometry_version`, solver-settings and environment
digests, the corresponding canonical metadata or immutable resolvable
references, source network identity, explicit selection/orientation scope,
and completeness evidence. `EXACT_REUSE` additionally requires a source view
identity, positive generator/competitor/solve-context compatibility evidence,
and `completeness_proven: true`. A materialized view may retain
`completeness_proven: false`; it is not thereby certified.

`GeometryLineage` never changes a quantity or result status. `CERTIFIED`,
`PARTIAL`, `UNRESOLVED`, `NOT_CALCULATED`, `NOT_SUPPORTED`, and `STALE`
remain independent. Legacy provenance with no `geometry_lineage` key remains
absent; loaders and adapters must not infer `DIRECT_SOLVE` from an archive,
directory, or available geometry. If a caller explicitly normalizes legacy
data for this contract, it may report `UNRESOLVED` with null unknown IDs and
an explicit reason.

## Quantities and unavailable state

Every `Quantity` has `name`, `value`, `units`, `status`, `method`, and `scope`.
Optional fields are `support_count`, `total_count`, and `provenance`. Counts must
be nonnegative integers and support cannot exceed total. Values must be finite.
Units are explicit; dimensionless quantities use `1`.

Statuses are CERTIFIED, PARTIAL, UNRESOLVED, NOT_CALCULATED, NOT_SUPPORTED, STALE.
CERTIFIED/PARTIAL require a value. UNRESOLVED/NOT_CALCULATED/NOT_SUPPORTED require
null. STALE can retain a previous value for inspection. Marking a result STALE
propagates STALE to its typed quantities and contact selection, invalidating
certification without erasing cached measurements. Raw extension caches retain
historical producer flags under an extension whose status is STALE.

The adapter copies only cached values. It does not sum terms, divide by area,
count atoms/contacts/residues, infer rejection counts, derive Euler characteristic,
or synthesize planar zeros. Cached planar zeros remain valid values. Missing
per-area metrics remain NOT_CALCULATED. Producers can supply precomputed
`cached_metrics` Quantity records, including independently certified counts.

## Mean curvature

`MeanCurvature.integrated_mean_curvature` is the headline available measurement.
For a certified complete physical interface it is the full total. For AW with
missing faces/seams it is the cached supported sum, PARTIAL with scope
`supported_contributions`. `complete_total` remains a separate UNRESOLVED/null
quantity. This makes the partial number useful without suggesting a full total.

Decomposition fields are `smooth_face_contribution`, `edge_contribution`, and
`boundary_contribution`. The physical convention is:

    C_H = sum(face smooth H) + 1/2 sum(internal seam beta ds)

A free boundary has no automatic scalar mean-curvature term. Although the legacy
physical cache encodes that omission as boundary_H=0, the contract exposes it as
NOT_SUPPORTED/null under the explicit omission method, not a measured zero.
The original value remains in the versioned cache extension for provenance.
Dual producers supply their own validated discrete-edge convention and boundary
policy; this adapter makes no assumption about their construction.

## Gaussian curvature

`integrated_gaussian_curvature` is the intrinsic surface integral, **not** the
Gauss-Bonnet completion. Decomposition fields retain smooth faces, an optional
already aggregated intrinsic/discrete contribution, intrinsic seams, interior
vertex defects, free-boundary geodesic and corner terms. Completion, expected
value, and residual remain separate named quantities. Extrinsic fold beta is
never used as Gaussian curvature.

The cached Power example has intrinsic integral +5.188337486733097 and corner
turning -5.1883374867332295; completion and residual are approximately zero.
For AW, full intrinsic curvature/completion/residual are unresolved. Available
smooth/seam/vertex/boundary terms carry `supported_subcomplex` scope. The raw
supported-subcomplex Gauss-Bonnet audit remains in the extension. Its terms and
residual cannot replace the full-interface quantities. The expected full GB value
retains full-interface scope independently of the subcomplex terms.

### Authoritative curvature field map

The typed curvature blocks are the handoff contract to geometry producers:

| Results field | Meaning | Units |
| --- | --- | --- |
| `MeanCurvature.smooth_face_contribution` | smooth integrated H over curved faces | Å |
| `MeanCurvature.edge_contribution` | signed physical-seam contribution, `1/2 sum(beta ds)` | Å |
| `MeanCurvature.integrated_mean_curvature` | total integrated H, smooth plus physical seams when complete | Å |
| `MeanCurvature.complete_total` | separately certified complete total; otherwise null/status-bearing | Å |
| `MeanCurvature.boundary_contribution` | measured scalar boundary-H term only | Å |
| `GaussianCurvature.smooth_face_contribution` | smooth integrated K | 1 |
| `GaussianCurvature.intrinsic_seam_contribution` | intrinsic seam term | 1 |
| `GaussianCurvature.interior_vertex_defects` | summed interior junction/vertex defects | 1 |
| `GaussianCurvature.boundary_geodesic_contribution` | exterior boundary geodesic curvature | 1 |
| `GaussianCurvature.boundary_corner_contribution` | boundary corner turning | 1 |
| `GaussianCurvature.integrated_gaussian_curvature` | intrinsic total before exterior boundary completion; explicit alias `intrinsic_total_gaussian_curvature` | 1 |
| `GaussianCurvature.gauss_bonnet_completion` | completed Gauss--Bonnet total | 1 |
| `GaussianCurvature.gauss_bonnet_expected` | topology-derived expected `2 pi chi` value | 1 |
| `GaussianCurvature.gauss_bonnet_residual` | completed total minus expected total | 1 |

`intrinsic_discrete_contribution` is an optional aggregate audit field for the
interior-junction terms. If both it and `interior_vertex_defects` are present,
they must describe the same contribution and must not be added twice. The
completed total is the intrinsic total plus exterior boundary geodesic and
corner terms. Therefore smooth K and completed Gauss--Bonnet are intentionally
different quantities; a producer must never copy one into the other.

For free boundaries, an omitted scalar mean-curvature term is represented by
`boundary_contribution = null` with `NOT_SUPPORTED` (or another explicit
unavailable status). A numerical zero is valid only when the producer has
actually measured or proved that term to be zero.

GAUSS conversion requirements for molecular contact surfaces are therefore:

1. Preserve H as `smooth_face_contribution + edge_contribution`, with the
   physical seam convention and orientation stated in provenance.
2. Put smooth K only in `smooth_face_contribution` (also exposed as the
   `smooth_integrated_gaussian_curvature` semantic alias).
3. Put the intrinsic pre-boundary total in
   `integrated_gaussian_curvature`; never use the completed GB value for this
   field.
4. Keep seam, interior-junction, exterior-boundary, corner, completion,
   expected, and residual terms independently status-bearing.
5. Exclude topology-only cuts from every physical curvature contribution.

## Consolidated interface log contract

VECTOR owns writing `interfaces/logs.csv`; ATLAS owns the row vocabulary and
validation. The stable eight-column header is:

    interface_id,representation,side,quantity,value,units,status,provenance

The first seven columns are the scientific row payload. `provenance` is always
present as the eighth, nullable JSON column, so the header does not change
between runs. `INTERFACE_LOG_BASE_COLUMNS` is an internal seven-column prefix;
it is not a valid `logs.csv` header and must not be used to discard provenance.
Empty value and provenance fields are serialized as empty CSV fields, never as
zero or the string `None`.

Canonical `representation` values are:

    physical_partition_interface
    molecular_contact_surface
    dual_contact_complex
    dual_generator_complex
    dual_separator
    alpha_selection
    water_analysis
    timing

The Results identity mapping is `voronoi -> physical_partition_interface`,
`dual -> dual_contact_complex`, `dual_a/dual_b -> dual_generator_complex`, and
the corresponding same-name mappings for molecular surfaces and reserved dual
representations. The AW/Power partition is not silently discarded: when it is
not already part of the interface identity, it must be included in the
provenance JSON as `partition: "aw"` or `partition: "power"`.

Canonical log quantity names and units are exported as
`INTERFACE_LOG_QUANTITY_UNITS`; the Results curvature-field mapping is exported
as `RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG`; metric aliases are exported as
`RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG`; molecular-selection aliases are
exported as `RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG`. They include
`area` (`Å²`);
contact counts;
geometric counts; topology counts and flags (`1`); explicit curvature names
such as `smooth_integrated_mean_curvature` (`Å`),
`physical_seam_integrated_mean_curvature` (`Å`),
`total_integrated_mean_curvature` (`Å`),
`complete_integrated_mean_curvature` (`Å`),
`intrinsic_integrated_gaussian_curvature` (`1`),
`smooth_integrated_gaussian_curvature` (`1`),
`intrinsic_seam_gaussian_curvature` (`1`),
`exterior_boundary_geodesic_curvature` (`1`),
`interior_junction_gaussian_defect` (`1`), `boundary_corner_turning` (`1`),
`completed_gauss_bonnet` (`1`), `gauss_bonnet_expected` (`1`), and
`gauss_bonnet_residual` (`1`), `mean_curvature_per_area` (`Å^-1`), and
`gaussian_curvature_per_area` (`Å^-2`). Alpha counts and selected-pair structures use
`alpha_*` and dimensionless units; `alpha_value` retains its explicit native
units. Water quantities use `water_*` and dimensionless units. Timing rows use
`timing_<event>` and seconds (`s`).

Existing physical-interface metric aliases are `metrics.n_contacts ->
contact_count` and `metrics.n_residue_contacts -> residue_contact_count`, both
dimensionless. Other legacy `n_*` metric names remain future-table fields until
a producer supplies an explicit long-log mapping; they are not renamed by
guesswork.

Structured values are compact deterministic JSON in the `value` field:
UTF-8, sorted object keys, no insignificant whitespace, and JSON `true`/`false`
for booleans. Scalar strings remain plain strings; numeric values remain JSON
numbers. `Status` values are exactly `CERTIFIED`, `PARTIAL`, `UNRESOLVED`,
`NOT_CALCULATED`, `NOT_SUPPORTED`, or `STALE`. `AVAILABLE` is not a scientific
status.

`CERTIFIED` and `PARTIAL` rows require values. `UNRESOLVED`, `NOT_CALCULATED`,
and `NOT_SUPPORTED` rows require an empty value. `STALE` may retain a previous
value for inspection, or be empty if no value was retained. A measured zero is
serialized as `0` and remains distinct from every unavailable status.

When provenance is needed, it is deterministic JSON containing the source
snapshot identity and, for a `Quantity` row, its `method`, `scope`, support and
total counts, source revision/scientific-state ID, details, and diagnostics.
Rows are ordered deterministically by `(interface_id, representation, side,
quantity)` with stable producer order only as a final tie-breaker. Only
calculated quantities are emitted as values; unavailable quantities are emitted
as status-bearing empty rows so missing does not become zero.

The typed curvature aliases are semantic accessors only. They are not dataclass
fields and must not produce extra rows. Every actual mean/Gaussian quantity
field has exactly one canonical log quantity name; the aggregate
`intrinsic_discrete_contribution` remains distinct from the primary
`interior_vertex_defects` field and must never be added twice.

### Current Results coverage audit

The current typed coverage is intentionally uneven:

| Representation | Results coverage | Current gap or handoff |
| --- | --- | --- |
| `physical_partition_interface` | `InterfaceResult`, AW/Power identity, metrics, topology, H/K, selection | The current writer emits all eight columns, preserves provenance, and sorts by `(interface_id, representation, side, quantity)`. It must preserve `partition: "aw"` or `"power"` whenever identity is not available to the row consumer. |
| `molecular_contact_surface` A/B | `MolecularContactSurfaceResult`, geometric counts, topology, H/K, patch provenance | The GAUSS conversion currently does not place surface `area` in `metrics`; the current writer has the canonical molecular-selection mapping, so it will log those aliases once the typed surface result is populated. |
| `dual_contact_complex` | Canonical representation name and extension channel | No generic typed dual-result row adapter exists; no measurements may be copied from the primal interface. |
| `dual_generator_complex` | `GeneratorComplex` and `generator_complex_extension()` with typed simplex/topology quantities | The extension is not traversed by the current compact log path; VECTOR/GAUSS must provide explicit side and quantity rows if dual counts are to be logged. |
| `dual_separator` | Reserved canonical representation only | Mathematics and measurements are not implemented; emit only explicit unavailable rows when requested. |
| `alpha_selection` | Alpha identity fields and physical-interface selection quantities | There is no standalone typed `AlphaSelectionResult`; selection rows must use the canonical `alpha_*` names and retain native units for `alpha_value`. |
| `water_analysis` | Output metadata only | No typed Results water block currently exists; do not promote raw water metadata to certified science without a producer contract. |
| `timing` | Output metadata only | No typed Results timing block currently exists; timings remain diagnostic rows with seconds and provenance. |

Legacy curvature is represented honestly: an absent source key becomes
`NOT_CALCULATED`, while a present source key with a null value becomes
`UNRESOLVED`; neither becomes a numerical zero. A retained cached value whose
scientific state is invalidated is `STALE` and may keep its value.

### Execution-profile coverage

The command profiles control production stages; they do not create alternate
Result schemas. The Results adapter consumes the resulting
`InterfaceGeometryAnalysis` cache and preserves the distinction between a
stage that was not requested and a requested stage that failed to resolve.
Selecting a profile does not itself select interface mode; physical-interface
rows below apply when the interface workflow is also requested.

| Profile | Producer stages | Results state when adapted from the cache |
| --- | --- | --- |
| `geometry` | Network geometry only; no curvature analysis. In ordinary group mode no physical-interface cache is produced. | Geometry/topology/curvature are absent unless another producer supplies them. If a no-curvature interface cache is adapted, selected area and alpha selection are retained; skipped topology/H/K quantities are `NOT_CALCULATED`. |
| `interface` | Physical interface and alpha selection; curvature/topology and water analysis are skipped. | Area and alpha-selection quantities are retained with producer certification; skipped topology/H/K quantities are `NOT_CALCULATED`. |
| `analysis` | Interface/network curvature and topology analysis; water analysis is skipped. | Cached area, topology, H, K, and alpha-selection quantities retain their producer statuses and provenance. Requested-but-null terms are `UNRESOLVED`; unsupported boundary scalar H remains `NOT_SUPPORTED`. |
| `full` | Analysis stages plus interface-water analysis. | Same typed physical Results coverage as `analysis`; water and execution timing remain metadata/extension concerns because no typed Results blocks exist for them. They must not be promoted to certified quantities by the adapter. |

For every profile, `adapt_interface_analysis` is cache-only: it does not
recalculate geometry or analysis. `curvature_requested: false` in the cache is
the authoritative signal for deliberate omission. A missing source key is
`NOT_CALCULATED`; a present-null value from a requested stage is `UNRESOLVED`;
unsupported science is `NOT_SUPPORTED`; retained invalidated cache values are
`STALE`. Cached values preserve the adapter provenance on each Quantity.
Profile names are orchestration metadata and are not representation names:
the physical result remains `physical_partition_interface`, with the AW/Power
partition carried by identity/provenance.

## Topology and chemistry

`Topology.quantities` standardizes vertices, edges, faces, Euler characteristic,
connected components, boundary edges/loops, nonmanifold edges, and orientability.
Each has its own support status. An absent orientability flag is not invented
from a broader cached manifold certificate. `orientation_status` describes the
cached orientation certificate, independently of whether a boolean orientable
quantity was exposed by the producer. Existing certified component genus is
retained with a justification; no genus is calculated or inferred.

`AtomIdentity` contains a stable system namespace, chain, residue number,
insertion code, residue name, atom name, and alternate location. A stable parent
system-atom ID takes precedence when available; otherwise the molecular tuple
defines `atom-v1:<sha256>`. The existing stable-ID string is retained as an alias.
Generator index requires a network ID and never defines the stable atom key.
The metadata adapter accepts `system_num` as the established parent-system ID,
never a generic network `num` field. Producers must preserve parent-system ID
semantics and distinguish alternate locations/duplicate molecular identities.
Stable atom identities deliberately exclude frame, radii, and representation.

`ContactRecord` stores ordered atom_A/atom_B chemistry, owning result identity,
dual/physical feature IDs, selection status, local Quantity records, provenance.
Local area/curvature attribution must come from its producer; shared seam values
are not apportioned automatically. Contact ID hashes result ID, ordered stable
atom IDs, and feature IDs, permitting multiple mapped features for one atom
pair. `molecular_pair_id` hashes only the ordered stable atom pair, supporting
cross-frame/representation comparisons. Selection statuses distinguish SELECTED,
REJECTED, UNRESOLVED, NOT_CALCULATED, NOT_SUPPORTED, STALE. Unknown selection
strings fail validation. Chemistry side order is supplied explicitly, never
inferred from sorted generator indices. Protein-water records are not converted
into A/B contacts by the adapter.

## Environment

`Environment` records dry/solvent-competing, explicit solvent present, solvent generator count,
interface competitor count when available, and provenance. It validates the
canonical `dry`/`solvent_competing` values against a supplied participating-
generator count; legacy `wet` is rejected. The physical adapter copies
the cached solvent-context and competing-generator fields and rejects identity
disagreement. Solvent present in input alone does not imply solvent competition. A missing
interface-specific competitor count remains null; total solvent input atom count
is not relabeled as generator or interface-competitor count.

## Representation extensions and future tables

Extensions have namespace, schema version, status, immutable JSON-safe payload,
and provenance. `vorpy.interface_geometry_analysis` retains original physical
accounting, selection, coverage, unresolved reasons, and solvent comparison,
including information with no universal mapping. The adapter retains
`excluded_pairs` separately: existing exclusions can overlap unresolved pairs
and therefore cannot be renamed certified rejections.

Agent #3 can supply a `vorpy.dual_side` extension with side identity, cached mesh
reference or arrays converted to immutable JSON state, exact selected simplex
IDs, incidence mapping, support/selection policy, and AW/regular provenance.
Typed topology and discrete curvature blocks enter through `adapt_cached_result`.
No existing primal measurements are copied to a requested dual representation.

For non-surface duals, `GeneratorComplex` and
`generator_complex_extension()` provide a dimension-counted simplicial/generator
complex with side `A`, `B`, or `AB/contact`, plus arbitrary typed topology
quantities. Area, faces, boundary loops, and curvature are optional result
blocks and are not fabricated. `dual_separator` is reserved as a separate
representation; its mathematics is not implemented here.

`METRICS_COLUMNS` defines the ordered future tidy-table contract. One row is one
analysis representation/side/frame. It includes identity and units/probe/radii,
chemistry counts, area, H/K and ratios, status/scope, complete H, separate GB
completion/status, topology, selection, and solvent fields. All fields are
scalar or null. Future table units: lengths Å, area Å², H Å, H/area Å^-1,
K dimensionless, K/area Å^-2; alpha retains its explicit native units. Unknown
values must become empty CSV fields, not zero. No table writer is implemented.

`HUMAN_SUMMARY_LEVELS` defines shared GUI/text hierarchy: headline identity and
area/H/K/chi/contacts; support topology/decomposition/selection/status and scope;
then network/frame/alpha/radii/probe/numerical provenance and diagnostics. Human
renderers must display status/scope beside numbers, especially PARTIAL and STALE.

## Agent #2: persistence requirements (prefer scientific-state persistence)

Choose **B: persist underlying scientific state and cheaply rebuild adapters**.
Result objects are snapshots, not the scientific authority. Do not persist
duplicated derived Result dictionaries as authoritative session objects.

For deterministic reconstruction, retain:

- Stable system namespace, display name, frame and group/interface identities;
  ordered side membership, parent-system atom IDs and complete molecular metadata.
- Stable network ID and explicit system/interface scope; generator-to-parent atom
  mappings and boundary-generator identity distinct from physical atoms.
- Exact base/effective radii, radius configuration fingerprint, probe application
  policy/value/units; scheme, alpha method/convention/value/native units.
- Relevant solve's competing-solvent context/count, input-solvent presence, and
  interface competitors if actually measured; never infer wet from input alone.
- InterfaceGeometryAnalysis metadata, alpha_selection, both orientation caches,
  coverage, unresolved reasons, dry_solvent_comparison. Retain nulls, statuses,
  complete versus partial H, all Gaussian terms, full versus subcomplex scopes,
  topology/manifold/genus certificates, and boundary policy.
- Cached physical representation records and existing mapping/contact-selection
  records, stable atom references, local contributions and attribution status.
- Already built dual/filtration scientific state if available, including feature
  IDs, support/completeness, unresolved records, method/options and incidence;
  reconstruction of a Result must not build either scientific object.
- Producer numerical methods/options/revisions, scientific-state revision or
  fingerprint and invalidation state. STALE must survive reload. Preserve
  absent-versus-present-null fields, since these distinguish NOT_CALCULATED from
  UNRESOLVED in the adapter.
- Precomputed per-area metrics and chemistry/count descriptors if the session
  expects them immediately; this layer will not calculate missing aggregates.

If older state lacks mandatory identity/provenance, report that absence or use an
explicit unknown configuration identity. Do not infer scope from a directory or
fabricate a certificate. No archive/session files were changed for this layer.

## Molecular contact surface extension

The interaction triplet remains three separate scientific objects: the A
and B `MolecularContactSurfaceResult` objects and the physical partition
`InterfaceResult`. `InteractionLinks` associates their result IDs with one
interaction ID and the selected-contact source and identity.

The contract reserves these interaction-level quantities without calculating
them: `area_A_contact`, `area_partition_interface`, `area_B_contact`,
`area_ratio_A_to_partition`, `area_ratio_B_to_partition`,
`area_ratio_A_to_B`, and eventual integrated mean- and Gaussian-curvature
quantities for A, the partition, and B. They are optional typed quantities;
missing values remain `NOT_CALCULATED`, `NOT_SUPPORTED`, or another explicit
status. No curvature complementarity subtraction is defined.

Orientation is explicit. Molecular contact surfaces use outward molecular
normals (`outward_molecular_normal`). A physical partition uses its own
A-to-B orientation (`physical_interface_A_to_B`). Signed curvature
comparisons must not be made until both conventions and any required
orientation transform are explicit.

Use the canonical terms `probe-expanded union-of-balls boundary`, `selected
partner-restricted union-of-balls patch`, `physical Power interface`,
`physical AW interface`, `dual contact complex`, and `A/B generator complex`.
The molecular boundary convention and contact-selection provenance are
separate fields; the curved molecular result is not a partition or dual
representation.

### Geometry records versus topology cells

Topology-only cuts are bookkeeping edges introduced solely to regularize the
cell attachment model. They are separate from geometric edges and must not
contribute to length, area, mean curvature, Gaussian curvature, or physical
boundary length. They are never passed to geometric-measure or curvature
integration as physical seams.

The certified dry 2KAI Power topology regression is:

    A: geometric edges 784, topological edges 786, cuts 2,
       chi -5, beta (1, 6, 0), genus 1
    B: geometric edges 754, topological edges 757, cuts 3,
       chi -3, beta (1, 4, 0), genus 0

Patch attachment-record counts remain separate from both spherical-patch
records and topological face counts, even where their numerical values match.
`topological_edge_count` is the explicit topology-layer name; any legacy
`edge_count` field is not automatically substituted for geometric edges.

The molecular result keeps geometric records separate from any topological
cell complex:

    geometric representation
        carrier atoms, spherical patches, circular arcs, junction records,
        refined seams, and geometric boundary arcs
            ↓
    canonical cell decomposition
        vertices, edges, faces, connected components, and boundary components
            ↓
    topology invariants
        Euler characteristic, orientability, genus, Betti numbers, and
        manifold/singular diagnostics

`SphericalPatch` is the canonical curved geometry and provenance record. It
is not automatically a topological 2-cell. A patch may later be shown to
contain multiple retained regions, boundary cycles, or holes.

The contract uses `boundary_arc_count` for geometric boundary-arc records
and `boundary_component_count` for topological connected boundary
components. They are independent quantities.

Euler characteristic is meaningful only when a valid cell decomposition is
available. Genus additionally requires certified manifold topology,
component and boundary-component accounting, and an orientability result.
Raw `V-E+F` is not by itself a genus calculation. Unknown orientability is
represented by a null/status-bearing quantity, never by `False`; a future
nonorientable result has a separate `nonorientable_genus` quantity.

Cell-complex quantities such as `chain_complex_valid`,
`boundary_of_boundary_zero`, `beta0`, `beta1`, and `beta2` remain optional
and may be `NOT_CALCULATED` or `PARTIAL`. Every quantity carries its own
status, so certified area and local geometry can coexist with unresolved
topology. Legacy `patch_count` values are retained as source geometry
metadata and are never promoted to `face_count` without explicit source
evidence.

The same rule applies across `InterfaceResult`, generator complexes, and dual
contact complexes: geometric records, a canonical cell decomposition, and
topology invariants are separate layers. Matching numbers are not evidence
of equivalence.

`MolecularContactSurfaceResult` is separate from physical partition
`InterfaceResult` and `GeneratorComplex`. Its identity records side A/B,
surface convention, definition version, contact-selection source and ID,
partner, source network/interface IDs, and shared interaction ID. Geometry and
curvature blocks are optional; missing values remain status-bearing nulls.
`MolecularPatchProvenance` records source atoms, spheres, contacts, clipping
atoms, patches, components, arcs, and junctions. `InteractionLinks` associates
the A surface, physical interface, and B surface while leaving comparison
quantities such as area ratio uncalculated.
Reserved conventions are `expanded_union_of_balls`, `solvent_accessible`, and
`solvent_excluded`; the contract does not claim implementation of any of them.

## Examples and tests

Exact `to_dict()` examples are provided in `result_contract_examples/power.json`
and `result_contract_examples/aw.json`. They use self-contained copies of the
already saved InterfaceGeometryAnalysis regression caches, materialized as Python
records before invoking the adapter. Runtime adapters do not read those files.
The numeric source snapshots come from the existing open-interface-integration
2KAI cache; these examples do not run `fixture_analysis` or curvature kernels.

Focused tests cover precise Power/AW science, statuses/nulls/scopes, immutability,
stable IDs, environment, unsupported dual combinations, all result types,
extension validation, and blocked file/scientific entry points. Production
exporters, persistence, Workbench, presets, and output directories are unchanged.
