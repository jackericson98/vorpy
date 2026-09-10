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
```

See [Settings](settings.md), [Selections](selections.md), [Exports](exports.md), and [Examples](examples.md).

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
