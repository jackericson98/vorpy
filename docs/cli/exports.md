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

Every preset also saves `<system-name>.vpy` in the output directory. This is a
solved-network ZIP archive containing JSON metadata/state and NumPy arrays; it
can be loaded without rebuilding vertices, edges, or surface triangulations.
Multiple presets in one command save the archive only once.

Use `-e large -e no_archive` to retain the preset outputs without the automatic
archive. `-e none` disables automatic exports, and `-e none -e logs` or
`-e only logs` exports only the requested component without an archive.
An explicit `--save-network path.vpy` still saves the requested archive.

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

The run folder retains the usual system PDB, information, and `EDTA.vpy`
checkpoint. The original `EDTA/` group folder keeps its ordinary exports. Its
sibling `solvent/` applies the same preset to the selected water layers:

```text
<run>/
  EDTA.vpy
  EDTA.pdb
  EDTA/                       ordinary group exports
  solvent/
    aw_logs.csv               one native log for all requested water layers
    aw_verts.txt              combined solvent vertex checkpoint
    info.txt                  combined solvent group information
    layers.csv                atom/water/completeness and boundary counts by layer
    README.txt
    layer_01.pdb              also layer_02.pdb through layer_04.pdb
    layer_01_surfs.off         also layer_02_surfs.off through layer_04_surfs.off
    layer_01_edges.off         also layer_02_edges.off through layer_04_edges.off
    layer_01_verts.off         also layer_02_verts.off through layer_04_verts.off
    atoms/                    per-atom geometry from the large preset
    nonpolar_interface_geometry/
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
also adds atom geometry and analysis tables. Surrounding-atom PDBs are omitted
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
geometry again. New layer calculations (`-c l N`) require the structure input.
