# Power/Laguerre interface analysis

VorPy's Power interface mode is an optional analysis of a weighted regular
complex. It uses Cazals/Chothia atomic radii, expands each radius by the probe,
and assigns power weight \(w_i=(r_i+probe)^2\). It is separate from the
ordinary VorPy additively weighted (AW) solve and does not change AW geometry.

The `cazals` preset is explicit: it selects the Cazals/Chothia radius
assignment and defaults to probe 1.40 Å, alpha 0, and condition-b threshold
M=5. Water and other HETATM records are excluded from this protein AB model.
Radius fallbacks are reported in the summary.

Use the existing VorPy group syntax. Multiple selectors joined by `and` form
one group:

```powershell
python vorpy vorpy/data/2KAI.pdb --power-interface cazals `
  -g c A and c B -g c I
```

The output defaults to `output/2KAI`. The normal files are
`power_interface_summary.txt`, `power_interface_atoms.csv`,
`power_interface_facets.csv`, `power_interface_edges.csv`, and
`power_interface_components.csv`. An export root can be selected with the
existing `-e dir <directory>` syntax. The optional settings are
`--power-probe-radius`, `--power-alpha`, and `--power-M`.

```powershell
python vorpy vorpy/data/2KAI.pdb --power-interface cazals `
  --power-probe-radius 1.4 --power-alpha 0 --power-M 5 `
  -g c A and c B -g c I -e dir output/cazals
```

The analyzer reports Cazals' raw turning-angle sum
\(C_{raw}=\sum l\beta\) separately from the conventional
\(H=(k_1+k_2)/2\) edge measure \(C_H=C_{raw}/2\). It does not substitute
the latter into the published Cazals statistic.

## Cazals validation

The focused validation suite is separate from routine solves. It writes
inventory, curvature-distribution, condition-b, facet, and forensic outputs:

```powershell
python -m vorpy.src.analyze.cazals_validation --system 2KAI
python -m vorpy.src.analyze.cazals_validation --system 1UDI --output-dir output/cazals_validation
python -m vorpy.src.analyze.cazals_validation --system all --output-dir output/cazals_validation
```

The frozen observations are recorded in
[the validation baseline](../cazals_power_validation_baseline.md). Synthetic
geometry is exact; local singleton-angle checks agree for 945/945 1UDI edges
and 864/864 2KAI edges. The 2KAI interface atom ratio and significant
component count reproduce the reported observations, but VorPy has **not**
reproduced the published whole-interface \(+17^\circ\) value: the current
magnitude is \(7.634^\circ\). The remaining source/methodology difference is
unresolved, so this is not described as an exact reproduction of all Cazals
results.

## Power versus AW

To run the existing AW interface analysis and the independent Power analysis
side by side, use the existing two-group `-i` interface workflow plus
`--power-vs-aw`:

```powershell
python vorpy vorpy/data/2KAI.pdb -g c A and c B -g c I -i `
  -s mv 5 --power-vs-aw -e dir output/2KAI-comparison
```

The command writes `power_vs_aw_summary.txt` and
`power_vs_aw_metrics.csv`, along with the model-specific Power and AW exports.
It reports areas, components, interface localization, and curvature terms
under separate model labels. Power has planar facets and discrete edge
turning. AW has curved surfaces, edge turning, and smooth surface curvature.
No combined Power/AW total is emitted. The AW edge orientation and surface
normal sign conventions have not been reconciled for that sum; Power's
normalized Cazals statistic is also not interchangeable with AW's smooth
integrated mean curvature. Lys15-associated features are reported as a
spatial localization diagnostic, not as a quantitative literature match.
The first 2KAI comparison currently reports zero supported AW surfaces: all
480 candidate pair surfaces are excluded because the constructed AW interface
network has incomplete generator cells. The report records this instead of
treating unsupported geometry as a valid comparison. A bounded, validated AW
network is required before the shared interface-area and curvature comparison
is scientifically interpretable.
