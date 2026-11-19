# Figure 5

## Screened PEG reruns (2026-11-12)

`F5A_PEG_chemistry_vs_hydration.py` and `F5C_water_contacts.py` use
`output/peg_one_layer_investigation_20261112/one_layer_run_status.csv` at the
repository root. Only the 69 model/frame sets with usable AW/Power/Primitive
rows are selected. Each row's `analysis_directory` routes verified re-exports,
including recovered AW PEG3000 frame 5. Both polymer group-name layouts are
accepted. Missing or invalid screened inputs stop generation instead of silently
reducing the cohort. Original logs and solver geometry are read only.

Run the scripts directly to regenerate `output/PEG/F5A` and `output/PEG/F5C`.
F5A retains its Power-versus-AW atomic calculations and plot styling. F5C now
includes all seven PEG models and only identity-checked whole layer-1 waters;
native logs pooling two layers are filtered using `solvent/layer_01.pdb`.
All three complete atomic cells, incident-area/neighbor agreement and external
SA = sum(atomic SA) - 2 sum(intramolecular SA) are required at existing tolerances.
Counts are atom/water **frame observations**, not independent replicates.

The PEG branch of `F5B_water_burial.py` uses the same screened layer membership
and intersects validated full PDB water identities across all three schemes
before its existing AW/Power comparison. Unmatched waters are excluded. Its
non-PEG loader is unchanged. F5C remains an AW-only descriptive plot, so it
uses all valid AW layer-1 waters within the screened frames.

Pre-regeneration outputs are preserved in
`output/PEG/before_screened_20261112/`. F5C had no saved numerical summary;
its previous statistics are reconstructed from the current unchanged inputs
and former configuration, not claimed as an exact historical plot snapshot.

## Earlier shared-table panels

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
