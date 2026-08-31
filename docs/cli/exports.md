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

**TODO:** Add the complete list of export components and the exact contents of every preset from the export implementation.
