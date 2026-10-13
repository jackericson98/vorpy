# 2KAI input-scope acceptance criteria

This is an input-audit and test-design record only. It does not authorize a
PRM, POW, or AW molecular solve, and it does not change frozen values or
tolerances.

## Requested comparison

The intended biological partition is the deposited 2KAI assignment:

- group A: kallikrein chains A+B;
- group B: inhibitor chain I;
- orientation: A+B -> I.

The canonical comparison inputs are the new chain-restored solvated
`vorpy/data/2KAI_MD_centered_chainrestored_solvated.pdb` and the derived dry
`vorpy/data/2KAI_MD_centered_no_sol_corrected.pdb`. The original
`vorpy/data/2KAI_MD_centered.pdb` remains the immutable source, and the supplied
`vorpy/data/2KAI_MD_centered_no_sol.pdb` remains a rejected incomplete fixture.
The raw MD chain fields are blank in the supplied source, so the strings `A`,
`B`, and `I` must never be inferred from that column. A comparison is eligible
only after an independent identity/sequence audit maps every retained protein
atom to the deposited reference `vorpy/data/2KAI.pdb`, with no unresolved
duplicates, alternate locations, or missing atoms in A, B, or I.

## Current input audit (no geometry computed)

| input | atom records | protein atoms | water atoms | Na/Cl atoms | raw chains | scope finding |
|---|---:|---:|---:|---:|---|---|
| `2KAI_MD_centered.pdb` | 44,520 | 4,340 | 40,086 | 94 (Na 53, Cl 41) | blank only | solvent-present; full A+B+I protein ordering is present, but chain reconstruction is not independently resolved by the current blank-chain preflight |
| `2KAI_MD_centered_no_sol.pdb` | 3,101 | 3,101 | 0 | 0 | blank only | solvent-absent, but not full A+B+I: it is exactly the solvated protein records with serials 1240–4340 |

The solvated file has 13,362 contiguous water residue occurrences under the
comparative preflight definition; the older key-based audit reports 10,000
because wrapped residue numbers collide. The acceptance audit must use an
explicit residue-occurrence definition, not silently mix these counts.

The no-solvent file starts at solvated serial 1240 and matches the solvated
records through serial 4340, including coordinates. The omitted serials 1–1239
are the leading protein segment; in the deposited/reference order, the
solvated prefix begins with chain A and serial 1240 begins with chain B. Thus
the current no-solvent file is not an A+B versus I input. It cannot be used for
a scientifically equivalent two-file metric comparison. This is a blocker,
not evidence for changing a frozen result.

## Corrected derived protein-only input

The supplied files remain untouched. The new derived input is
`vorpy/data/2KAI_MD_centered_no_sol_corrected.pdb`. It contains all 4,340
protein atoms from the solvated source, in original serial/order, with atom
names, residue names/numbers, elements, and coordinates preserved. The only
PDB identity fields changed are the blank chain and insertion-code columns;
both are restored from the authoritative deposited `2KAI.pdb` ATOM topology.

The mapping is an exact ordered match of 288 N-terminated residue blocks by
`(residue_number, residue_name)`, not an inference from blank chain IDs. The
resulting atom membership is A=1,239, B=2,243, I=858, so A+B=3,482 and I=858.
Restoring insertion codes removes the duplicated stable identities around
deposited insertion residues (zero duplicates remain). The per-atom evidence
is in `vorpy/data/2KAI_MD_centered_no_sol_corrected_atom_audit.csv`, with the
summary and hashes in
`vorpy/data/2KAI_MD_centered_no_sol_corrected_audit.json`.

The corrected no-solvent view retains no nonprotein species. The source ion
inventory is Na=53 and Cl=41, with no Mg; 40,086 SOL atoms are also omitted.
Those species are candidates for an explicitly labeled `SOLVENT_CONTEXT` or
`FULL_SOLVENT` run, not for this no-solvent input. No molecular solve was run.

Under the same `MOLECULE_ONLY` policy, `prepare_network_geometry` produces
exactly matching protein generators for the solvated source and corrected
input: system indices, stable IDs, locations, radii, and masses all agree, and
both retain zero participating solvent generators. This is the canonical
protein-only control; it does not make the solvent-context runs scope
equivalent to the dry runs.

## Chain-restored solvated input

The new `vorpy/data/2KAI_MD_centered_chainrestored_solvated.pdb` retains all
44,520 source records in their original order: 4,340 protein atoms, 40,086
SOL water atoms, and 94 ions (53 NA and 41 CL). Its SHA-256 is
`2e869bee8330cffd4ddffe62f1e6063726fd1f8ebee59f5834bf04411d275e6a`.
Coordinates, serials, residue names/numbers, atom names, elements, occupancy/
B-factor fields, and every nonprotein line are preserved byte-for-byte at the
record level. Only protein chain and insertion-code provenance fields differ
from the blank-chain source.

Labels were established from the already validated corrected protein mapping:
an exact ordered match of 288 N-terminated protein residue blocks by
`(residue_number, residue_name)` to authoritative `vorpy/data/2KAI.pdb`
ATOM/DBREF/COMPND topology. This yields A=1,239, B=2,243, and I=858. The
solvent and ion records were not assigned inferred chains or new topology;
they are retained physical context only. The new file is scientifically
acceptable for selecting protein groups A+B versus I, but it is not evidence
of solvent-cell completeness or solvent/ion participation.

## Excluded trajectory frame

`vorpy/data/2KAI_Frames/frame_0100ps.pdb` is not a canonical comparison source.
Its SHA-256 is
`d15ce739150bc47d08d444211b8546af2f9da956840aec9c69d0f933f06185cd`. The
`TRJCONV` header identifies `t=100.00000 step=50000`; it has the same 44,520
record categories and protein atom order as the solvated source when chain
labels are ignored, but uses raw protein chains A/B/C rather than blank/A/B/I
provenance. Its box is 85.754 Å instead of 85.632 Å, and the best-fit rigid
alignment leaves 1.514 Å heavy-protein RMSD over 2,238 atoms. It is therefore
a distinct MD conformation, not a rigidly transformed copy. No PRM, POW, or AW
comparison may substitute this frame without a new scientific equivalence
audit.

## Solvent presence versus participation

These are separate required fields in every future result row:

| field | meaning | current pair |
|---|---|---|
| `solvent_present` | recognized water/ion records exist in the input PDB | `true` / `false` |
| `solvent_mode` | generator policy used to build the network | must be recorded explicitly |
| `solvent_participating` | solvent generators were retained in the solved generator universe | cannot be inferred from `solvent_present` |
| `solvent_cells_complete` | complete solvent-cell evidence exists | required for `FULL_SOLVENT`; not certified by `SOLVENT_CONTEXT` |

The accepted meanings are:

- `MOLECULE_ONLY`: solvent may be present in the source file, but it is not a
  participating generator and must not affect geometry.
- `SOLVENT_CONTEXT`: solvent generators participate as context; this does not
  certify complete solvent cells.
- `FULL_SOLVENT`: solvent participates and complete solvent-cell evidence is
  required before claiming full-cell completeness.

For this pair, the no-solvent file has no solvent to participate. The solvated
file must be run in a separately named mode; a solvent-present file run as
`MOLECULE_ONLY` is not a solvent-participation comparison.

## Six-run interface acceptance matrix

The planned comparison is three methods by two explicit competitor scopes:

| run family | canonical input | solvent mode | physical generators | retained competitors |
|---|---|---|---:|---|
| PRM dry | `2KAI_MD_centered_no_sol_corrected.pdb` | `MOLECULE_ONLY` | 4,340 protein | none |
| PRM with competitors | `2KAI_MD_centered_chainrestored_solvated.pdb` | `SOLVENT_CONTEXT` | 44,520 total | 40,086 SOL + 53 NA + 41 CL |
| POW dry | `2KAI_MD_centered_no_sol_corrected.pdb` | `MOLECULE_ONLY` | 4,340 protein | none |
| POW with competitors | `2KAI_MD_centered_chainrestored_solvated.pdb` | `SOLVENT_CONTEXT` | 44,520 total | 40,086 SOL + 53 NA + 41 CL |
| AW dry | `2KAI_MD_centered_no_sol_corrected.pdb` | `MOLECULE_ONLY` | 4,340 protein | none |
| AW with competitors | `2KAI_MD_centered_chainrestored_solvated.pdb` | `SOLVENT_CONTEXT` | 44,520 total | 40,086 SOL + 53 NA + 41 CL |

Every row must record the canonical input fingerprint, A+B/I membership,
method, solvent presence, solvent mode, participating-generator count, and
solvent-cell status. The metric payload must separately report interface
surface area (Å²), contact count and definitions, edge count, vertex count,
curvature terms with their method-native units/sign convention, completeness
status, and numerical uncertainty status. Counts must be integers; geometric
values must be finite when a result is claimed, and missing values remain
missing rather than zero.

The dry and solvent-context rows are intentionally not an equivalent
two-input metric comparison: the latter has 40,180 additional participating
solvent/ion generators. A solvated `MOLECULE_ONLY` control is the exact
protein-generator equivalence check against the corrected dry input. The
`FULL_SOLVENT` mode is a separate future configuration with the same explicit
retention as `SOLVENT_CONTEXT`, but it is not accepted until complete
solvent-cell evidence is independently supplied. PRM/POW rows may be
`COMPLETE` only when their selected interface features and topology evidence
are accounted for. AW rows may reach at most
`LOCAL_COMPLETE_GLOBAL_UNRESOLVED` for this study; no AW row may emit
`GLOBAL_COMPLETE`.

## Metric-equivalence gates

Metrics may be compared across the two inputs only when all of the following
match: mapped A+B/I atom identities, coordinates for the compared atoms,
orientation, atom radii/weights, boundary model, probe and alpha/selection
settings, surface resolution, maximum-vertex policy, and the explicit solvent
mode/generator universe. Otherwise report the metric as an individual
diagnostic with `scope_equivalence=NON_EQUIVALENT`, never as a solvent effect.

The method-specific acceptance scope is:

- **PRM:** compare ordinary Voronoi interface geometry only within equivalent
  molecule and participation scopes. Area, component/topology counts, and any
  polygonal curvature terms must retain their units and selection status.
- **POW:** apply the same gate, with identical power weights/radii. A POW
  value is not interchangeable with a PRM value; a matched-radius POW control
  is a method contrast, not an input equivalence proof.
- **AW:** apply the same gate and retain the analytic per-surface/per-edge
  provenance. AW curvature and Gaussian terms may be compared only as local or
  partial diagnostics when their support is complete for the reported subset.
  No result from these files may claim global AW curvature certification.

Reversing A+B and I is an orientation check: signed mean-curvature quantities
must change sign, while unsigned quantities, area, and interface membership
must remain invariant. This check does not repair missing input scope.

## Required statuses

The comparison harness should keep the following statuses distinct:

| status | acceptance meaning |
|---|---|
| `INPUT_AUDIT_BLOCKED` | raw identity, chain mapping, or required group atoms are unresolved |
| `SCOPE_NON_EQUIVALENT` | files differ in molecule or participating-generator scope |
| `COMPLETE` | all required features/cells for the stated PRM/POW view are accounted for |
| `INCOMPLETE` | required features/cells or topology evidence are missing |
| `LOCAL_COMPLETE_GLOBAL_UNRESOLVED` | local AW incidence/contributions are complete for a stated subset, but global AW closure is not proven |
| `LOCAL_UNRESOLVED` | local AW feature or numerical evidence is missing |
| `GLOBAL_COMPLETE` | reserved for an independently proven global certificate; do not emit for the current AW molecular comparison |

Numerical uncertainty is separate from completeness:

- `CONVERGED`: analytic/quadrature or other numerical convergence evidence is
  present at the predeclared tolerance;
- `MAXIMUM_ORDER_REACHED`: a finite estimate exists but convergence was not
  demonstrated;
- `UNRESOLVED`: no scientifically defensible estimate is available.

Missing values remain missing; they are not replaced with zero. A small global
residual cannot promote local AW evidence to `GLOBAL_COMPLETE`.

## Focused test design

Before any molecular solve, focused tests should assert:

1. the exact record/category counts above and blank raw chain fields;
2. exact no-solvent suffix identity against the solvated protein records;
3. failure-closed A+B/I eligibility when the chain map is blank or chain A is
   absent;
4. solvent presence and solvent participation are reported as different facts;
5. PRM/POW/AW comparison rows carry method, scope, solvent mode,
   completeness, uncertainty, and units;
6. AW rows never promote local evidence to global certification.

The first two checks are implemented by the focused no-solve input-audit test.
The remaining method-output checks belong beside the shared Results contract
and require coordination with ATLAS before any schema change.
