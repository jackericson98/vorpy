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

**TODO:** Verify every setting, alias, default, unit, and accepted range against the implementation.
