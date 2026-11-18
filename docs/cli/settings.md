# CLI Settings

Settings are modified with `-s`.

```bash
python vorpy example.pdb -s nt pow and sr 0.05 and mv 80
```

| Setting | Meaning | Example |
|---|---|---|
| `nt` | Network type | `aw`, `pow`, `prm` |
| `sr` | Surface resolution | `0.05` |
| `mv` | Vertex search/build limit | `5`, `10`, `40` |
| `bm` | Bounding-box multiplier | `1.25` |
| `ss` | Surface coloring scheme | scheme-dependent |
| `sc` | Colormap | Matplotlib colormap |

For PDB input with a valid `CRYST1` record, `bm` scales the lattice's centered
Wigner-Seitz cell and candidate Voronoi vertices outside that convex cell are
discarded. Lengths are interpreted in Å and angles in degrees. This supports
orthogonal and triclinic cells, including rhombic-dodecahedron representations.
Without `CRYST1`, VorPy retains the atom-coordinate axis-aligned bounding box.
The filter limits vertex acceptance; it does not add periodic images or create
new boundary faces. Multipliers below `1.0` contract the cell about its center;
for example, `-s bm 0.95` moves every cell plane inward by 5% and retains
approximately 85.7% of the original cell volume.

**TODO:** Verify every setting, alias, default, unit, and accepted range against the implementation.
