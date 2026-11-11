# Solvent-layer investigation — 2026-11-11

Historical findings below describe the state before implementation. The normal
CLI now supports `python vorpy EDTA -c l 4 -e large`; see
[current CLI documentation](../cli/overview.md) and [layer exports](../cli/exports.md#solvent-layers).

The existing layer traversal is close to the requested molecular definition, but
there is no working end-to-end CLI command for four water-only layers in this
checkout. A fresh EDTA solve and a read-only audit of an existing solvated EDTA
frame produced inspectable PDBs and tables. Production solver/CLI code was not
changed. The investigation exporter is `scripts/investigate_solvent_layers.py`.

## Definition and existing code

The intended layer is the shortest path through water molecules joined by
positive-area Voronoi faces: the selected solute is distance zero, waters sharing
a face with it are layer 1, and each later layer contains previously unassigned
waters touching the previous layer. All atoms of a water participate, not only
its oxygen. Ions are neither layer members nor bridges between layers.

`vorpy/src/group/sort.py` computes group statistics and requests the first shell.
The traversal itself is `vorpy/src/group/layers.py:get_layers` (line 131).
It follows surface adjacency and, by default, expands each neighbor to its whole
residue. It handles system/topology index mapping explicitly. Those parts can
be reused.

On both real EDTA tests, legacy traversal assigned all observed waters to the
same layers as an independent water-only breadth-first traversal. It also
included Mg in layer 1; the larger frame included Na in layers 3 and 4.
**Neither real case showed repeated atoms.** The initial suspicion about
backtracking was not reproduced; absence of a global visited set alone is not
evidence that these test results are wrong. Water-only traversal should still
express the shortest-path/visited rule directly and prevent ion-mediated paths.

## CLI and export gaps

1. `-c layers 4` is advertised by the old help, but the current calculation
   planner rejects it: `Unknown calculation target: layers`.
2. `-e layers` is also rejected: `Unknown export type: 'layers'`. Export parsing
   passes only the export command name onward, with no layer-count handling.
3. `-c s` currently produces an empty `[['f']]` group on EDTA. The full-group
   shortcut in `command/group.py` compares a nested command list with strings.
   The successful test used an explicit atom range instead. `-g full` reaches
   the same shortcut problem. This needs a small group-command fix and regression
   tests, not a new solver.
4. The Python API `group.exports(layers=4)` exists, but its depth check can reuse
   a previously cached one-layer result instead of computing four. Its legacy
   directory is `layers`, filenames start at zero, and it writes PDBs and meshes
   without per-layer native logs. Its loop includes trailing empty layer slots,
   and the information export reads atom counts at a different offset from the
   volumes. It does not implement the requested `solvent/layer_01/...` bundle.
5. A group solve ordinarily solves the selected cells, retaining other atoms as
   geometric context. Context atoms appearing in a network or PDB do **not** mean
   their water cells and water-water faces have been solved. Calling
   `get_layers(max_layers=4)` cannot create those missing faces.

## EDTA tests and outputs

Root: `output/solvent_layers_edta_investigation/`.

### Fresh small model

`vorpy/data/EDTA.pdb` contains 612 atoms: 32 EDTA atoms, one Mg, and 193 three-site
waters. A new AW calculation explicitly targeted all 612 atoms, using maximum
vertex radius 5 Å and the existing default surface resolution/box multiplier.
The network build took about 102 seconds and marked 455 cells complete.

The fresh archive is `full_aw.vpy`; native full-network logs and the complete
input PDB are in `full_aw/EDTA/`. Audited results are in `fresh_EDTA_checked/`.

| Layer | Water molecules | Water atoms | Complete water molecules |
|---|---:|---:|---:|
| 1 | 49 | 147 | 49 |
| 2 | 113 | 339 | 46 |
| 3 | 31 | 93 | 0 |
| 4 | 0 | 0 | 0 |

The small box has no fourth layer in its recorded face graph. Empty layers have
CSV headers and an empty directory cohort, but no PDB because the existing PDB
writer skips empty selections.

### Existing larger EDTA frame

There are ten saved AW frame archives under `output/EDTA_frames/frames/`.
Frame 1 has 1,511 atoms: 32 solute atoms named UNK, one Mg, two Na, and 492 waters.
Its existing data are:

- `output/EDTA_frames/frames/frame_0001/EDTA_frames.vpy`
- `output/EDTA_frames/frames/frame_0001/EDTA_frames.pdb`
- `output/EDTA_frames/frames/frame_0001/aw_logs.csv`

Normal archive loading currently raises `ArchiveError: Unknown fields for
network`. The saved network includes fields absent from the current schema,
including `solvent_mode`, `solvent_cells_complete`, retained-generator metadata,
and provenance. This is an archive compatibility issue, not proof of file
corruption. The root CSV has 1,223 atom rows, rather than all 1,511, so it is not
by itself a complete replacement for the archived atom table.

For this investigation, the existing archive codec read only the four numeric
geometry tables. Every generator mapped once to the saved PDB and its coordinates
matched within 0.002 Å. No archived object graph was migrated or saved back.
This table-only inspection is not a repaired `.vpy` reload for the normal CLI.

Results are in `cached_frame_0001_checked/`:

| Layer | Water molecules | Water atoms | Complete water molecules |
|---|---:|---:|---:|
| 1 | 49 | 147 | 49 |
| 2 | 121 | 363 | 120 |
| 3 | 220 | 660 | 138 |
| 4 | 101 | 303 | 9 |

One water is unassigned within four steps. One nonpositive/nonfinite-area face
was excluded. These are **observed-graph cohorts**: incomplete cells in inner
layers mean some missing contacts could change later assignments. Both test
reports explicitly use `status: observed_graph_only`. The outer-layer volume
and area statistics must exclude incomplete molecules; a PDB containing all
three atoms does not establish complete Voronoi geometry.

Each audited result contains:

```text
seed.pdb
view_layers.pml                 PyMOL loader with four shell colors
report.json
legacy_comparison.csv
solvent/
  summary.csv
  layer_01/                    also layer_02, layer_03, layer_04
    water.pdb
    atoms.csv                  recorded cell volume, area, completeness, identity
    incident_surfaces.csv       all recorded positive-area faces touching this layer
```

The layer CSVs are flat audit tables, **not** the multi-section native VorPy log
format. Incident faces can appear in both neighboring layer exports, and their
endpoints retain original system indices. Use the source PDB/archive for context;
re-solving an isolated layer PDB would change its geometry.

## Commands that work now

From the repository root, this is the tested full-solvent small-EDTA solve. Choose
new output names if preserving another run:

```powershell
.venv/Scripts/python.exe -m vorpy vorpy/data/EDTA.pdb --boundary none -g a 0-611 -s nt aw and mv 5 -e dir output/edta_full_water_new and only logs pdb --save-network output/edta_full_water_new.vpy
```

The explicit atom range includes solvent. `mv 5` is the setting used for this
bounded test, not a universal guarantee of geometric completeness. Inspect
complete-cell flags and retain the full solvent box as geometric context.

To audit the newly saved geometry into a new output folder:

```powershell
.venv/Scripts/python.exe scripts/investigate_solvent_layers.py --network output/edta_full_water_new.vpy --seed-residue MOL --layers 4 --output output/edta_layers_new
```

To inspect the existing larger EDTA frame without solving it again:

```powershell
.venv/Scripts/python.exe scripts/investigate_solvent_layers.py --archive-tables output/EDTA_frames/frames/frame_0001/EDTA_frames.vpy --structure output/EDTA_frames/frames/frame_0001/EDTA_frames.pdb --seed-residue UNK --layers 4 --output output/edta_saved_frame_layers_new
```

Open the generated `view_layers.pml` with PyMOL, or read `solvent/summary.csv` and
the layer atom/face tables directly. The script refuses to reuse an existing
output directory.

## PEG inventory and rebuild scope

`OrganizedModels/PEG` has seven models (EG, PEG300, PEG1000, PEG2000, PEG3000,
PEG4000, PEG6000), with eleven frame directories each. The first-frame AW log
for every model contains only polymer atom rows. Water cell geometry needed for
layers 2–4 cannot be recovered from those polymer tables alone.

The checked `output/PEG3000/PEG3000.vpy` and `output/PEG3000_1/PEG3000.vpy` retain
6,071 atom records but target only 479 polymer cells. The latter marks exactly
those 479 cells complete. These are not solved full-water networks.

For the current local 6,071-atom `vorpy/data/PEG3000.pdb`, the full-system baseline
would use `-g a 0-6070` and the model's validated solver settings, followed by an
audit with `--seed-residue PEG`. That full PEG solve was not run. Other models
require their own atom range; do not reuse this bound blindly. `--all-frames`
and `--parallel-frames N` address trajectories, not a directory of unrelated
PDB files. The 77 already split model/frame inputs would need a batch loop.

## Smallest useful production improvements

1. Fix full-system group selection and restore compatible loading of the saved
   solvent metadata, with focused regressions. Preserve existing archives.
2. Add one explicit CLI layer-count option, propagate it to frame workers, and
   document water-only, molecule-level shortest-path membership. A spelling such
   as `--solvent-layers 4` would be new functionality; it does not exist today.
3. Initially use a correctly solved full-water network for classification. For
   larger PEG systems, progressively solve each discovered water frontier while
   retaining all original generators. Solve the requested outer layer too if its
   volume/area is being exported. Do not substitute a distance cutoff or solve
   each layer in isolation. A larger solvent box may be needed for complete
   fourth-layer cells; changing export settings cannot repair missing geometry.
4. Export stable `solvent/layer_01` through `layer_04` selections with source-index
   provenance, completeness counts, PDBs, and correctly scoped native logs.
   Add meshes/cells according to the existing export preset. Fix cached-depth
   handling and trailing-empty-layer/off-by-one export behavior at the same time.
5. Gate a full PEG rerun on a complete-layer EDTA acceptance case and one PEG
   pilot. Reuse intact full-solvent archives after schema compatibility is fixed;
   rebuild missing water cells and boundary-truncated geometry where necessary.

Validation performed: synthetic cycle/depth/ion-bridge self-checks, independent
queue and set-frontier traversals on both real graphs, legacy comparisons,
disjoint water cohorts, and exported PDB atom counts. All passed. Power/Primitive
schemes, the other nine saved EDTA frames, all 77 PEG logs, and a full PEG rebuild
were not run. Incomplete outer geometry remains the main scientific limitation.
