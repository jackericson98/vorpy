# Figure 5

The active panels are:

- **F5A_exposed_vs_buried.py** (formerly F5G): PEG atomic geometry differences
  between Power and AW, grouped by atom type and AW solvent exposure, plus
  volume differences versus exposure and water-contact count.
- **F5B_water_burial.py** (formerly F5H): water-molecule geometry differences
  between Power and AW, grouped by AW polymer-contact fraction.

Both support direct-file execution, including PyCharm's Run button:

```powershell
.venv/Scripts/python.exe vorpy/src/analyze/Part_II_Molecular_Analysis/F5/F5A_exposed_vs_buried.py
.venv/Scripts/python.exe vorpy/src/analyze/Part_II_Molecular_Analysis/F5/F5B_water_burial.py
```

A reads `output/summaries/atoms.csv`. B reads the dataset configured by
`DATASET_ROOT` in its script and requires water-inclusive logs with complete
water cells and incident faces. Polymer-only logs produce coverage and summary
CSVs with unavailable water geometry rather than fabricated measurements.

To prepare the shared tables and run both panels:

```powershell
.venv/Scripts/python.exe -m vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all --models PEG3000 dp24 --frames 1
```

The runner retains dataset/model/frame selection, caching, validation, and
`--output-root`. Use `--figures A` or `--figures B` to select one panel, or
`--figures` alone to prepare shared tables without plotting. Run `--help` for
all loader options. Standalone scripts use their own plotting settings.

PNG/SVG plots and CSV companions are written to `output/F5A/` and `output/F5B/`.
Existing G/H products have been renamed to these prefixes. Shared tables and
caches remain in `output/summaries/` and `output/cache/`.

A uses an AW atomic water-area fraction strictly greater than 0.01 for exposed
atoms. Absolute differences are calculated per atom before averaging; signed
plots retain the sign. Error bars describe sample SD of pooled atom/frame
observations, so larger models contribute more weight. Missing values and zero
AW denominators are excluded. Terminal O-H hydrogens remain in the atom-type
audit and are excluded from plotted types.

B sums complete water cells, removes intramolecular faces from the external
area, and counts distinct neighboring molecules. AW polymer-contact fractions
at or below 0.10 are exterior-like; fractions at or above 0.90 are
polymer-surrounded. Polymer-surrounded waters are internal-hydration candidates;
these measurements do not establish enclosure or solvent accessibility.

`F5_common.py` and `F5_config.py` retain shared input checks, identity matching,
area validation, and configuration. Frames are correlated observations, not
independent molecular replicates.

Retired A-F scripts, their plots and summaries, old run reports, the empty
`SolventExposedDiff.py` placeholder, and the previous README are preserved in
`archive/`. Archived scripts are historical copies, not active entry points.
The original scientific dataset is unchanged.

Run the focused checks from the repository root:

```powershell
.venv/Scripts/python.exe -m pytest vorpy/src/analyze/Part_II_Molecular_Analysis/F5/tests --override-ini addopts= -q
```
