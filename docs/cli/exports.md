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

Use `-e large -e no_archive` to retain the preset outputs without the automatic
archive. `-e none` disables automatic exports, and `-e none -e logs` or
`-e only logs` exports only the requested component without an archive.
An explicit `--save-network path.vpy` still saves the requested archive.

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
