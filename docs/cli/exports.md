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

When an interface export includes `logs`, VorPy writes both the legacy
sectioned interface CSV and `interface_<A>_<B>/logs.csv`. The latter is the
canonical eight-column Results log (`interface_id`, representation, side,
quantity, value, units, status, provenance); it serializes cached state only
and does not rebuild geometry. The `small` preset also writes this canonical
log for its compact interface bundle.

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
