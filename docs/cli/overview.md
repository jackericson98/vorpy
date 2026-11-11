# Command-Line Interface

The CLI is intended for reproducible calculations, batch workflows, large systems, and automated analysis.

```bash
python vorpy <input_file> [options]
```

Examples:
```bash
python vorpy example.pdb
python vorpy example.pdb -s nt pow
python vorpy example.pdb -g a 0-100
python vorpy example.pdb -e small and shell
python vorpy EDTA -c l 4 -e large
```

`-c l N`, `-c layer N`, and `-c layers N` calculate N water-molecule layers
around the selected solute (the normal non-solvent group by default). Layer 1
shares a positive-area Voronoi face with the solute. Each later layer contains
previously unassigned waters sharing a face with the previous layer. All atoms
of each water are retained; ions do not participate in the water-layer graph.
The count must be a positive integer. Multiple selected groups act as one seed.

VorPy solves each newly reached water frontier against the full input geometry,
then builds one combined network for the requested solvent cohort. The ordinary
solute group remains separate. Choose a single network scheme; combining this
calculation with interface/comparison mode is not supported.

Layer calculations work with `--all-frames` and `--parallel-frames N`. Every
frame uses the same export preset and contains its own checkpoint, ordinary
group folder, and sibling `solvent/` folder. An explicit `--save-network` path is
propagated to every worker/frame, with `_frame_NNNN` suffixes.

See [layer exports](exports.md#solvent-layers) for the files and completeness flags.

See [Boundaries](boundaries.md), [Apollonius analysis](apollonius.md), [Power interface analysis](power_interface.md), [Settings](settings.md), [Selections](selections.md), [Exports](exports.md), and [Examples](examples.md).

**TODO:** Generate a complete authoritative option table from the current CLI parser.

## Analytic AW edge diagnostics

Use the argumentless `--diagnose-edges` flag to audit constructed AW edges
after network construction and before export:

```powershell
python vorpy edta -s mv 5 --diagnose-edges
```

The terminal report includes analytic curved, turning-conic, and straight-edge
matches, coverage, maximum and RMS sample reconstruction residuals, maximum
absolute length difference, status counts, and a bounded list of problem
edges.

This option is read-only. It does not replace sampled edge points or change
analysis fields, logs, surfaces, or exported geometry. POW and PRM networks are
reported as unsupported and are skipped.
