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
