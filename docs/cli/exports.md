# CLI Exports

VorPy supports **presets** and **individual export components**.

Standard presets include:
```text
tiny
small
medium
large
all
```

Examples:
```bash
python vorpy example.pdb -e small
python vorpy example.pdb -e shell
python vorpy example.pdb -e small and shell
```

Presets are convenience bundles rather than mutually exclusive modes. Individual outputs can be added to a preset.

The standard system summary is written as `system_info.txt`; each solved group
or interface network writes `network_info.txt`. Older `info.txt` exports remain
readable by the Workbench. The normal `large`/`all` presets retain physical AW
surfaces, edges, vertices, contacts, and curvature outputs but do not create
the experimental `nonpolar_interface_geometry/` tables. Request those tables
explicitly with `-e nonpolar`.

Every preset also saves `<system-name>.vpy` in the output directory. This is a
solved-network ZIP archive containing JSON metadata/state and NumPy arrays; it
can be loaded without rebuilding vertices, edges, or surface triangulations.
Multiple presets in one command save the archive only once.

When an interface export includes `logs`, VorPy writes both the legacy
sectioned interface CSV and `interface_<A>_<B>/logs.csv`. The latter is the
canonical eight-column Results log (`interface_id`, representation, side,
quantity, value, units, status, provenance); it serializes cached state only
and does not rebuild geometry. The `small` preset also writes this canonical
log for its compact interface bundle.

Use `-e none` to disable optional files; a completed geometry solve still
writes its mandatory canonical checkpoint. The compatibility modifiers
`-e no_archive` and `-e no checkpoint` do not suppress that durability
checkpoint. An explicit `--save-network path.vpy` still saves the requested
archive. `-e only logs` remains an explicit analysis export and therefore
requires an analysis-capable execution profile.

## Execution-profile dependencies

`--profile geometry` is deliberately geometry-only: it does not populate cell
volume, curvature, or the full analysis dictionary. Geometry-only components
such as edges, surfaces, and vertices can be requested explicitly with
`-e only ...` or `-e add ...`. The bundled `small`, `medium`, `large`, and
`all` presets include analysis-dependent summaries/logs for system and group
scopes, so VorPy rejects those combinations with an actionable message after
the mandatory checkpoint has been verified. Use `--profile full` for those
bundles, or use `-e none` for checkpoint-only geometry. Interface scope keeps
its compact physical-interface summary available under the geometry profile;
unsupported curvature/topology fields remain explicitly unavailable.

For a cache-only archive run, add `--profile-cache-export`. The opt-in
`cache_export_profile.json` beside the normal export directory records load,
export, and explicit-save wall/CPU time, peak RSS when supported, archive
member counts, cached interface state, and existing export phases. It accepts
`.vpy` input only and does not rebuild geometry.

Archive timing remains available in regular logs and archive provenance. The
standalone `<system>.vpy.timing.json` sidecar is emitted only for explicit
verbose, diagnostic, or profiling runs.

Atom cell meshes combine their already solved surface points and triangles,
offsetting triangle indices for each surface. Exporting a cell does not run
triangulation again. OFF retains geometry and colors; VTP also retains scientific
face fields.

AW interface exports support `-e dual` (alias `-e apollonius`), also included in
large/all presets. See [Apollonius interface visualization](apollonius_interface.md)
for the PyMOL launcher, selected/rejected/unresolved layers, and contact mapping.

**TODO:** Add the complete list of export components and the exact contents of every preset from the export implementation.

## Solvent layers

```bash
python vorpy EDTA -c l 4 -e large
```

The run folder retains the usual system PDB, information, and `EDTA_checkpoint.vpy`
checkpoint. The ordinary `groups/EDTA/` folder keeps its group exports. The
sibling `solvent/` applies the same preset to the selected water layers:

```text
<run>/
  EDTA_checkpoint.vpy
  EDTA.pdb
  system/                     system information
  groups/EDTA/                ordinary group exports
  solvent/
    aw_logs.csv               one native log for all requested water layers
    aw_verts.txt              combined solvent vertex checkpoint
    group_info.txt            combined solvent group information
    layers.csv                atom/water/completeness and boundary counts by layer
    README.txt
    layer_01.pdb              also layer_02.pdb through layer_04.pdb
    layer_01_surfs.off         also layer_02_surfs.off through layer_04_surfs.off
    layer_01_edges.off         also layer_02_edges.off through layer_04_edges.off
    layer_01_verts.off         also layer_02_verts.off through layer_04_verts.off
    atoms/                    per-atom geometry from the large preset
```

Surface set 1 is the solute–layer-1 interface. Set 2 is the layer-1–layer-2
interface, and so on. Edges/vertices belong to the corresponding interface.
Curvature colors use the inner side of each interface. The selected geometry
format (`off`, `ply`, or `vtp`) and edge/vertex sizing settings apply normally.

Each water appears in exactly one layer PDB. An empty layer still has an explicit
empty PDB and empty geometry files when requested by the preset. There are no
duplicated native logs per layer. The combined atom log adds `Solvent Layer`,
which the normal log reader preserves as an integer. Inspect `layers.csv`
before using volumes, areas, or deeper memberships: missing faces can affect
later layers. A larger solvent box or different solver limits may be needed.

Presets retain their normal intensity: small writes information, layer surfaces,
and the combined log; medium/large add the layer PDBs, edges, and vertices; large
also adds atom geometry. Surrounding-atom PDBs are omitted
from `solvent/` because it already contains the explicitly selected cohorts.

Only consecutive complete water layers, starting at layer 1, enter the native
solvent log. A layer must be nonempty and every water atom must pass the existing
complete-cell check. At the first incomplete or empty layer, that layer and all
later layers are excluded from log atom rows, outer-only geometry rows, and group
totals. Shared geometry incident to retained cells remains as necessary context.
VorPy prints a warning giving the requested count, complete count, and excluded
range. `layers.csv` records `complete` and `included_in_logs` for each layer.
Parallel-frame runs also relay these warnings from each frame's `run.log` to
the main terminal.
PDBs, meshes, and checkpoint membership remain available for inspection.

The root `.vpy` stores both networks and layer membership/boundary indices.
Re-export it with `python vorpy <checkpoint.vpy> -e large`; this does not solve
geometry again. To calculate layers from a checkpoint, use:

```powershell
python vorpy "path/to/dp12_checkpoint.vpy" -c l 4 -e small -e dir "path/to/rerun"
```

The `.vpy` must be the actual input argument. A comment containing its path does
not make a command that opens a `.pdb` reuse saved geometry; PDB input starts a
fresh solute calculation.

Checkpoint layer calculations preserve the saved solute network and input
checkpoint. Matching saved layers and already computed solvent cells are reused.
Missing water frontiers and the combined solvent export network are calculated
as needed. If the saved network excluded solvent generators, its dry-cell faces
cannot identify water contacts: VorPy explicitly reports a separate wet-context
seed calculation, using the saved atom coordinates, radii, and solver settings.
This does not overwrite or rebuild the original dry network. Explicit settings
and interface recovery remain separate checkpoint operations.

For many independent checkpoints, use a CSV manifest with `checkpoint,output`
columns. Paths are relative to the manifest's folder, so a Dropbox manifest can
move between computers without retaining drive letters:

```powershell
python scripts/run_checkpoint_layers.py "path/to/vorpy_layer4_jobs.csv" --workers 20 --dry-run
python scripts/run_checkpoint_layers.py "path/to/vorpy_layer4_jobs.csv" --workers 20
```

The runner validates every checkpoint and output destination before launching,
sorts by decreasing checkpoint bytes as an estimate of runtime, and immediately
fills free worker slots. Each child uses one numerical-library thread. This is
a concurrency limit, not CPU affinity; memory use grows with active jobs. Reduce
`--workers` if RAM is limiting. Console logs and a flushed completion summary
are written under `vorpy_layer4_runlogs` beside the manifest. Nonzero exits are
recorded and the remaining jobs continue; the runner exits nonzero if any fail.
Runs do not automatically skip existing outputs or retry failed jobs.
