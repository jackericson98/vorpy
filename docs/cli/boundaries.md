# Molecular boundaries

The default is **No boundary** (`none`): use only the loaded generators,
without adding a virtual shell. Loaded solvent remains part of the geometry.
For an unsolvated structure in an interactive terminal, VorPy prompts for
`0` (none, the default when pressing Enter) or `1` (virtual solvent shell).
Non-interactive runs use no boundary without prompting. An explicit
`--boundary none` or `--boundary shell` also skips the prompt.
Exterior cells can remain unbounded when no boundary is added.

```bash
python vorpy molecule.pdb
python vorpy molecule.pdb --boundary none
python vorpy molecule.pdb --boundary shell
python vorpy solvated.pdb --boundary explicit
```

`explicit` requires recognized solvent or ion residues. `shell` may be used
with either solvated or unsolvated structures. The shell is a deterministic,
box-like arrangement of geometry-only generators; it is not an SES, SAS, or
solvent model. Shell generators remain separate from molecular atom and
residue tables.

The defaults are a 1.4 Å probe radius, 2.5 Å shell spacing, 2.8 Å shell offset,
and 1.5 Å generator radius. The spacing is slightly below the 3.0 Å diameter
of an explicit-water-sized generator to avoid exactly tangent neighbors. The
offset is twice the probe radius from the molecular union after accounting for
atomic radii.

```bash
python vorpy molecule.pdb --boundary shell \
    --probe-radius 1.4 --shell-spacing 2.0 --shell-offset 2.8
```

The shell parameters are included in boundary metadata. `ses`, `sas`, and
`probe` are reserved for future implementations and currently report that
they are planned but unavailable.
