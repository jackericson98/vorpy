# F5: polymer solvation and Voronoi scheme comparison

This pipeline analyzes existing AW, Power, and Primitive logs. It never runs or
changes the VorPy solver, edits logs, or writes into the organized source dataset.
All plots use NumPy/pandas/Matplotlib; no new dependencies are required.

## Run from the repository root

All F5 analysis scripts also support PyCharm's **Run** button and direct file
execution. The scripts locate the repository from their own file paths, so the
PyCharm working directory can be anywhere. Your settings at the top of each
script still apply. For example:

```powershell
.venv/Scripts/python.exe vorpy/src/analyze/Part_II_Molecular_Analysis/F5/F5G_exposed_vs_buried.py
```

PowerShell examples (use your Python environment if `.venv` differs):

```powershell
# All 12 models, expected 11 frames, and all six figure modules:
.venv/Scripts/python.exe -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all

# Representative real-data validation and all figures:
.venv/Scripts/python.exe -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all --models PEG3000 dp24 --frames 1 --output-root vorpy/src/analyze/Part_II_Molecular_Analysis/F5/representative_output

# A family, selected lengths and frames, and one figure:
.venv/Scripts/python.exe -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all --polymer PEG --models PEG1000 PEG3000 PEG6000 --frames 1-11 --figures E

# A figure module can also be run directly with the same arguments:
.venv/Scripts/python.exe -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5A_terminal_effects --models PEG3000 --frames 1 2 3

# Portable input/output overrides:
.venv/Scripts/python.exe -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all --dataset-root 'D:/OrganizedModels' --output-root 'D:/F5_results' --dpi 300 --width 13 --height 9

# Tests, including real PEG3000/dp24 checks when the external dataset is present:
.venv/Scripts/python.exe -m pytest vorpy/src/analyze/Part_II_Molecular_Analysis/F5/tests -q
```

Other options are `--schemes aw pow prm`, `--include-incomplete`,
`--exclude-hydrogens`, `--no-cache`, `--show`, `--area-rtol`, and `--area-atol`.
Frame selectors accept `1 2 3`, `1-11`, `1,3,5`, and `frame_0001`.
Inspect `--help` for the complete CLI. Other settings, identifiers, thresholds,
colors, PNG/SVG export and font sizes live in `F5_config.py`. Programmatic runs
accept `Config` and `F5_run_all.run(config, ['A', ...])`.

AW is blue, Power red, Primitive green. Chemical classes have separate colors
in F5F, with scheme-colored panel titles. Every generated plot has PNG and SVG
versions by default. Display is off unless requested.

## Reader and source conventions confirmed with real data

The existing reader is `vorpy/src/analyze/tools/compare/read_logs2.py`.
`read_logs2(path, all_=False, balls=True, surfs=True)` returns a dictionary with
`data`, `group data`, `atoms`, `surfs`, `edges`, and `verts`. The last four entries
are DataFrames; group/build entries are dictionaries. Only atoms and surfaces
are requested, avoiding repeated parsing of geometry not used here.

- Atomic `Index` is zero-based PDB **atom-record order**, not the one-based PDB
  serial. Surface `Balls` contains those same indices. Row order is never used
  to associate atom and surface records. PDB serials and original coordinates
  are retained. Log/PDB names, residue IDs, chains and coordinates are checked
  for every logged atom; a mismatch rejects that log instead of guessing a shift.
- These are polymer-scoped logs. Their atom sections contain polymer atoms;
  their surfaces also reference water and ions absent from the atom sections.
  The full frame PDB supplies endpoint identities. Complete polymer cells have
  all their incident areas represented, as checked against their recorded SA.
  Water cells and water–water surfaces generally are **not** fully logged, so
  global interface totals describe the logged scope, not the full solvent box.
- The tested PDBs use `SOL` with `OW`, `HW1`, and `HW2`. Hydrogen-water surfaces
  are included. A known water residue **and** a recognized site name are required;
  polymer oxygens are not classified as water. Water residue occurrence IDs,
  including residue-number wraps and TER segments, distinguish molecules.
  Water compositions require exactly one recognized oxygen and either zero
  (an oxygen-only representation) or two recognized hydrogens. Unsupported
  sites or compositions receive identity-validity flags.
- Heparin's four-character PDB names `AGLC` and `AIDO` appear as `AGL` and `AID`
  in the three-column VorPy log field. Both names are preserved and that precise
  truncation is validated. Sodium `SOD` is classified as an ion.
- Corresponding atoms are independently checked against the same frame PDB
  in every scheme, then matched by polymer/model/frame/index. Paired differences
  additionally check names, residues, chains and coordinates. Repeat-level
  matching checks its full set of atom identities and exact member coordinates;
  this is more conservative than the configurable atom-coordinate tolerance.
- Surface records in the representative real logs have unique unordered atom
  pairs. Repeated identical pair records are counted once even if reversed.
  Conflicting areas for one pair, or a surface ID identifying different pairs,
  reject the log. Separate disconnected faces cannot be reconstructed from such
  ambiguous records; the pipeline does not invent their interpretation.
- Surface `Surface Area` and `Contact Area` differ numerically. Only **Surface
  Area** enters hydration fractions. Contact area is retained for audit and
  interface totals; it is never substituted for missing geometric area.
- The reader omits missing optional atom numerics, but historically converts a
  missing surface area to zero and malformed numerics to defaults. F5 adds a
  narrow raw integrity pass: malformed/nonfinite numeric strings and missing
  endpoint identities reject a log; blank surface areas are restored to NaN.
  This pass is an audit, not a replacement general-purpose VorPy CSV parser.

## Repeating units and chemical classes

`residues.csv` retains original PDB residues. `units.csv` holds the repeat units
used for positional analysis and the exposure reference. Both identifiers are
present in `atoms.csv`.

PEG is one PDB residue per molecule. Its backbone is accepted only if the
observed heavy-atom names and elements form `O1–C2–C3–O4–...`, i.e. `O(CCO)n`,
and every successive bond distance lies in 0.9–1.85 Å. Each hydrogen is assigned
to its closest polymer heavy atom only if its distance is 0.6–1.35 Å. The unit
assignment is deterministic: `C,C,O` is one unit, and the first terminal oxygen
belongs to the first unit. Hydrogens follow their assigned heavy atoms. Internal
oxygens are `PEG ether oxygen`; end oxygens are `PEG terminal oxygen`. These
unit boundaries are a structural bookkeeping convention, not separate PDB residues.

Observed frame-1 repeat counts (also validated separately for every loaded frame):

| Model | EG | PEG300 | PEG1000 | PEG2000 | PEG3000 | PEG4000 | PEG6000 |
|---|---:|---:|---:|---:|---:|---:|---:|
| O(CCO)n units | 1 | 6 | 22 | 45 | 68 | 90 | 136 |

Heparin units are the monosaccharide PDB residues. Their ordering is discovered
from inter-residue C1–O4 distances in 0.9–1.9 Å, requiring one connected,
unbranched path with two ends. The dp labels are not used as ordering evidence.
No periodic unwrapping is guessed. Discontinuous, cyclic or ambiguous paths
retain solvation data but have unavailable end-distance/length measurements.
No convergence model is fitted.

## Measurements, validity and statistical units

For each polymer atom, water/polymer/other interface areas sum all unique
incident surface areas. A polymer–polymer face contributes to both atoms'
atomic properties and once to global logged-interface totals. `other` includes
ions and other non-water/non-polymer endpoints.

`water_fraction = water_area / recorded_atomic_surface_area` (and similarly
for polymer/other). Validity requires a positive finite denominator, a complete
interface sum, a match to recorded SA within `atol + rtol * SA`, valid endpoint
identities, and a complete cell by default. Invalid fractions are **retained,
unclipped**, and excluded from primary summaries. `fraction_valid` always
requires completeness; `primary_valid` can relax just completeness with
`--include-incomplete`, never a failed area partition. Missing numerics remain
NaN. Explicit zero-water area is valid when the total partition is valid.

`--exclude-hydrogens` filters target **polymer** atoms. All incident solvent
surfaces, including water-hydrogen surfaces, remain in a heavy atom's area
partition. Excluding water hydrogens from that partition would change its
denominator and is not supported.

Residue/repeat fractions are ratios of summed water areas to summed atomic SAs,
not means of atomic fractions. These are sums of atomic surfaces, including
internal polymer interfaces, **not** the exposed surface of a union of atoms.
All members must be primary-valid for a residue/repeat to enter primary summaries.
Residue volumes and mean atomic volumes propagate missing members as NaN.

Water coordination counts unique molecule-occurrence IDs, so OW/HW1/HW2 from
one molecule contribute one water. Atomic polymer coordination counts unique
polymer atoms. Residue/repeat polymer coordination counts unique polymer atoms
outside that residue/repeat; polymer-residue coordination similarly excludes
its own residue. Neighbor-set Jaccard uses all neighboring atom IDs and defines
empty versus empty as 1. Zero-count paired differences are valid. Area per water
neighbor is undefined (NaN) for zero water coordination, never divided by zero.

Atomic same-residue/repeat and adjacent-repeat descriptors count incident
polymer-atom contacts. Their residue/repeat sums remain **atom-contact
incidences** (an internal contact contributes twice), while nonlocal contacts
at residue/repeat level count distinct nonlocal repeat IDs. Global physical
interface counts always come from the unique interface table. Adjacent versus
nonlocal uses validated repeat positions and chain identity; unresolved ordering
is unavailable, not a fabricated zero. No contact is interpreted as a hydrogen bond.

Exposure bins are provisional: Buried `[0, 0.10]`, Intermediate `(0.10, 0.40]`,
Exposed `(0.40, 1]`. F5B and F5E environments are assigned using **AW repeat-unit**
fractions and carried unchanged to the other schemes. Zero-water units are
reported separately. F5C's buried fraction intentionally uses each scheme's
own provisional bin to show scheme dependence; its caption makes this explicit.

- F5A averages repeat measurements at each end distance within each frame,
  then averages frames. Bands are frame SD. Curves permit plateau exploration;
  they do not establish convergence or fit a plateau model.
- F5B boxes contain model/frame means, with optional chemical-class subsets.
  They do not treat individual atoms as simulation replicates.
- F5C uses actual repeat counts. Per-frame water fraction is area weighted;
  mean atomic volume is total valid volume divided by valid atom count;
  coordination is the mean per valid unit. Error bars are frame SD and are
  unavailable for a single frame. PEG and Heparin are plotted separately.
- F5E reports signed deviation, mean absolute deviation and median absolute
  deviation **from AW** for matched atom, repeat, and original PDB-residue
  observations. Original PDB residues receive their own AW reference bins;
  these differ from repeat bins for PEG, which occupies one PDB residue. The latter statistic is
  median `abs(S-AW)`, not median deviation about the sample median. Fractions
  use percentage points, other metrics native units. Binned plots average
  within model/frame, then over frames, then give each model equal weight.
  Bands describe between-model SD, not an uncertainty interval for one polymer.
- F5F uses chemical-class bins to avoid overplotting. Its curves average
  model/frame means with equal model weighting. Bin centers are plotted;
  accompanying tables give counts. The secondary figure uses total water area.

The 11 frames may be correlated. Means, SDs and distributions are descriptive;
there are no independent-molecular-replicate claims, confidence intervals,
significance tests, hydrogen-bond estimates, or frame-persistence analysis.

## Topology and internal water (F5D)

Logged Euler characteristic, integrated Gaussian curvature, Gauss–Bonnet
expected/error fields, boundary completeness/closed/manifold/orientable flags,
and component counts are preserved. Blank diagnostic values remain unavailable.
Several representative long-polymer boundaries are logged as non-orientable;
their Gauss–Bonnet expected/error values are absent. The script does not repair,
reinterpret or replace these geometry results. A separate `conditional_2pi_chi`
is calculated only when all required boundary flags are explicitly true.

There is no validated molecular envelope, probe-based accessibility computation,
cavity segmentation or solvent-connectivity graph in these inputs. Internal
water occupancy is therefore **unavailable**, with a status row for each
model/frame/scheme. No probe radius or sensitivity result is invented. Proximity
to a center of mass is not used. Enclosed cavities, open channels, interior
versus exterior water, penetration and diffusion cannot reliably be inferred
here. Euler characteristic is never converted to a cavity or water count.

## Outputs and validation

Default output root: `F5/output/`.

```text
output/
  PEG/F5A/ ... F5F/          PNG, SVG, and relevant CSV companions
  Heparin/F5A/ ... F5F/
  summaries/
    atoms.csv               Original atom identities, metrics, flags, neighbors
    residues.csv            Original PDB residue aggregation
    units.csv               Validated repeat aggregation and AW environments
    atom_pairs.csv          Matched/unmatched atom rows and paired differences
    unit_pairs.csv          Matched/unmatched repeat rows and differences
    residue_pairs.csv       Matched original PDB residues and differences
    interfaces/<family>/<model>/<frame>/<scheme>.csv
    interface_shards.csv    Index of reusable interface-level CSVs
    interface_totals.csv    Unique logged interfaces, by interface class
    F5A_*.csv ... F5F_*.csv  Frame counts, counts of atoms/units, figure summaries
  validation/
    run_config.json
    file_inventory.csv      Every selected/expected PDB-log combination
    failures.csv            Missing/failed scheme/frame records and reasons
    validation.csv          Per-log counts and exclusions
    pdb_residue_inventory.csv
    scheme_matching.csv     Matched, unmatched, identity-mismatched, valid counts
    independent_area_checks.csv
    figure_failures.csv
    run_report.json
    validation_report.md
    scientific_findings.md   Descriptive solvation/terminal/scheme results
  cache/                    Per-log normalized tables, not scientific source data
```

Set-valued neighbor columns are JSON arrays in CSVs. Surface rows retain both
endpoints' indices, identities and coordinates. Interfaces are sharded by
frame/scheme to avoid retaining all surface tables in RAM. Each log is parsed
once per uncached run, and all figure modules reuse the normalized tables.
Failed logs do not prevent other frames/models from running. Figures with no
valid observations are explicitly marked rather than filled with zeros.

Cache keys include source path, size, nanosecond mtime, full configuration, and
SHA-256 hashes of the F5 loader/config and existing reader code. Changing source
metadata or analysis settings invalidates the cache. Cache files are generated
locally; `--no-cache` bypasses them. The original source files are never overwritten.
Generated outputs/caches are ignored by F5's local `.gitignore` and can be regenerated.
Use a different output root for a different scientific cohort. Unrequested
figure files from earlier runs are preserved; the run report and findings only
describe the figures selected for the current run.

The tests cover synthetic classification, shuffled/non-contiguous indices,
water molecule counting, area partitions and residue weighting, zero water,
incomplete/missing faces, duplicate records, AW reference bins, paired identity
and coordinate checks, percentage-point differences, empty Jaccard sets, end
distance and backbone validation, missing topology, cache invalidation, and
continuation after a failed scheme. Real PEG3000/dp24 tests use the external
dataset when available and otherwise explicitly skip.

For frame 1 of PEG3000 and dp24, five selected atoms per scheme receive an
independent CSV-level water-area sum from the raw surface records: 30 checks
across the six representative logs. This validates the aggregation independently
of the main reconstruction path. It does **not** scientifically revalidate
the VorPy solver, dynamics, or topology construction. Full-run counts,
exclusions and limits are recorded in the generated validation report.
