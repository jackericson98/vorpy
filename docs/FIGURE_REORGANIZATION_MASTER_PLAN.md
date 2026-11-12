# VorPy Figure Reorganization Master Plan

Checkpoint date: 2026-11-12 (America/Chicago)
Owner: CODY — source-code inventory and future migration
Artwork owner: FIGY — Dropbox artwork coordination
Repository branch at checkpoint: `master`
Parent commit before this documentation checkpoint: `a7011bd65971001d1cbfd7b732bcbc79a4d54459`

## Checkpoint status and authorization

**Inventory COMPLETE. CODY migration NOT STARTED. This is a planning checkpoint, not an implementation or scientific-results milestone.**

The prior read-only audit is available in this session. It covered all 321 Python files under `vorpy/src/analyze/`, including ignored Python source, and relevant external generators, exporters, tests, and documentation. All 321 analysis files parsed successfully. A second content-hash scan found no changes to those files during the audit. There are 25 package initializer files among the 321.

No generators, molecular analyses, renders, or pytest tests were executed during that audit. No scientific input, numerical tolerance, reference value, output, code path, or dependency was changed. Top-level optional-module discovery used the verified repository `.venv/Scripts/python.exe` (Python 3.10.11); this was not a scientific runtime validation.

The working tree already contained extensive unrelated changes and apparent partial relocations before CODY's checkpoint. Thus "migration not started" means **CODY has not begun this planned migration**; it does not claim that the working tree matches the old committed layout. Those pre-existing changes must be preserved and coordinated with their owner.

The user authorized creating, staging, and committing ONLY this document, with commit message `docs: checkpoint figure reorganization plan F1-F7`. No push, analysis-code rename/move, dependency change, cleanup, or figure regeneration is authorized by this checkpoint. The Git index was empty before document creation.

## 1. Approved figure registry

This registry supersedes the provisional F1–F5 classifications from the earlier audit. Existing filenames and historical figure/panel prefixes are provenance, not authoritative assignments.

| Collection | Approved scope | Status |
|---|---|---|
| F1 | Voronoi methodology and geometry | Approved figure scope; panels unresolved |
| F2 | Molecular dynamics frame persistence | Approved figure scope; contact-persistence implementation incomplete |
| F3 | Pure-water GROMACS baseline and PEG hydration | Approved scope; pure-water baseline PLANNED, NOT COMPLETED |
| F4 | Heparin | Approved figure scope; panels unresolved |
| F5 | DNA and Hairpin | Approved figure scope; panels unresolved |
| F6 | Complex proteins, ligands, interfaces, and all multiscale analysis | Approved scope; panels unresolved |
| F7 | Coarse-graining | TENTATIVE collection; boundary with F6 unresolved |
| Curvature | Separate research collection | Separate from main figures; cross-references allowed |
| SI | To be determined | No final membership or panel assignments approved |

Classification rules for future review:

- PEG chemistry/hydration belongs to F3; Heparin chemistry/hydration belongs to F4.
- Dedicated DNA/Hairpin analysis belongs to F5. Mixed protein/nucleic multiscale code belongs to F6 until a deliberate split is approved.
- Complex proteins, ligand studies, molecular interfaces and atom-to-residue-to-molecule analyses belong to F6.
- Actual comparisons of coarse-grained representations are F7 candidates. Aggregating atomistic measurements by residue does not itself establish a coarse-grained model.
- Curvature-focused code belongs to the separate collection, even when it also concerns frames, proteins or coarse-graining; shared builders retain stable APIs.
- Statistical-ensemble methodology may support F1; supplemental statistical analyses remain SI candidates, not approved SI content.
- Shared code receives ownership/use labels, not an automatic move into a figure directory.
- No final panel numbers are assigned by this plan.

## 2. Existing directory structure

Paths are relative to `vorpy/src/analyze/`. This is the complete source-directory structure from the audit. The appendix inventories every Python file. Repeated generated model/frame branches are expressed as patterns rather than listing data files.

```text
analyze/
├── [root analysis modules]
├── nonpolar_interface/
├── tools/
│   ├── plot_templates/
│   ├── compare/
│   ├── CleanData/
│   │   └── neighbor_validation/                 [reports]
│   └── batch/
├── Part_I_Statistical_Ensembles/
│   ├── Fig1/
│   │   ├── AWVd_vs_Power/
│   │   ├── AWVd_vs_Power1/
│   │   ├── Cells/
│   │   ├── Edges/
│   │   ├── Surfaces/
│   │   ├── Vertices/
│   │   └── Building_Blocks/{Edges,Surfaces,Vertices}/
│   ├── Fig3/
│   ├── Fig4/
│   ├── Fig5/
│   ├── Fig6/
│   ├── Fig7/Other/
│   ├── Other/
│   ├── SupplementaryMaterial/
│   ├── VerifyVorPy/
│   └── old/
├── Part_II_Molecular_Analysis/
│   ├── F1_schemes_curvature_models/
│   │   └── 1C_curvature_representations/
│   ├── F2_frame_comparison/
│   ├── F3_SimpleModels/
│   │   ├── tests/
│   │   ├── output/
│   │   └── archive/{output,representative_output}/
│   ├── F4_atomic_deviations/
│   │   ├── Extra/
│   │   ├── Other/
│   │   ├── F2C_All_Data_Sets_Clustered_full/
│   │   └── F3_old/
│   ├── F5_Coarse_Grained_Models_and_Totals/
│   │   └── F7_Coarse_grained_residues/
│   │       └── F8_CG_Curvature/
│   │           ├── Old/
│   │           └── figure8_outputs/plots/
│   ├── Curvature/
│   │   ├── F4_atomic_curvature/
│   │   │   ├── 4A_curvature_assignments/
│   │   │   └── Extra/
│   │   ├── 5c2A_alanine_dna_facing/
│   │   ├── alanine_curvature_analysis/
│   │   ├── alanine_curvature_by_chain/
│   │   ├── alanine_dna_facing/
│   │   ├── curvature_panel_c2/
│   │   └── curvature_panel_violin/
│   ├── Network_Dynamics/
│   ├── SuplimentalInformation/
│   │   ├── comparative_figures/
│   │   └── Curvature/{alpha_figures,criteria_figures}/
│   └── extra/
└── Part_III_ML_Processes/
```

Current simple-model products include:

- `F3_SimpleModels/output/F5A/` and `output/F5B/`.
- `output/PEG/F5A/`, `output/PEG/F5C/`, and `output/Heparin/F5B/`.
- `output/PEG/before_screened_20261112/{F5A,F5C}/`.
- `output/summaries/interfaces/{PEG,Heparin}/MODEL/frame_N/` and `output/cache/{PEG,Heparin}/MODEL/frame_N/`.
- Archived A–F products under `archive/output/{PEG,Heparin}/F5[A-F]/`, with summaries and validation reports.
- `archive/representative_output/` and its `single_scheme_smoke/` subtree, with family/panel products, caches, interface shards, summaries and validation.

Generated CSVs, PNGs, SVGs, PDFs, reports and caches are mixed with source. Their existence does not prove that current imports can reproduce them. No artifact content was rewritten.

Relevant external source areas: `scripts/`, `src/geometry/visualization/`, `src/geometry/{validation_shapes,topology_diagnostic}.py`, `src/interface/{selected_export,geometry_analysis}.py`, `src/output/`, and `src/visualize/`. Relevant tests are under `vorpy/tests/analyze/`, embedded simple-model tests, and the Figure3B recovery test. Documentation includes `docs/analysis/`, simple-model READMEs and `SCRIPT_MAP.txt`.

No notebooks, R, MATLAB or Gnuplot scripts were found in the searched source inventory. This is not a search of all external Dropbox material or generated Python scene files outside source.

## 3. Principal generators and new classifications

In this section `P/` means `vorpy/src/analyze/Part_II_Molecular_Analysis/`. Output labels below are existing filenames, NOT newly assigned panels.

| Current generator/family | Inputs, purpose, entry | Existing products | New collection |
|---|---|---|---|
| `P/F1_schemes_curvature_models/1A_sphere_schemes.py` | Empty sphere-scheme slot | None | F1 draft |
| `P/F1_schemes_curvature_models/1C_curvature_representations/*.py` | Synthetic two-ball ratio/distance and curvature illustrations; scripts/helpers | Display; `2_1_ratio_*_dist.{png,svg}` cartoons | Curvature; F1 cross-reference |
| `P/F2_frame_comparison/FrameAnalysis.py` | Selected frames; AW/Power/Primitive volume comparisons; main | Display; no automatic figure filename | F2; scientific issue below |
| `P/F2_frame_comparison/CurvatureComparison{,_updated}.py` | AW group metrics versus frame; main | `frame_variation_analysis/frame_variation_{data,summary}.csv`, `frame_variation_<metric>.{png,svg}` | Curvature; F2 cross-reference |
| `P/F2_frame_comparison/AtomicCorrelation{,_updated,_multi_plots}.py` | Selected atomic-log descriptors; main | `atomic_correlation_analysis/`: combined-data/matrix CSVs, matrix and pair PNG/SVG, pair statistics CSV | Needs Review; curvature-heavy descriptor analysis is not persistence |
| `P/Network_Dynamics/residue_flexibilty_using_frame_dynamics.py` | Intended contact persistence, volatility, contact-area variability | None; docstring only | F2 draft |
| `P/F3_SimpleModels/F5A_PEG_chemistry_vs_hydration.py` | EG and PEG300/1000/2000/3000/4000/6000; paired Power/AW atoms in screened cohort; main | `output/PEG/F5A/F5A_PEG_{volume,surface_area}_vs_water_contacts.{png,svg}`, chemistry/coverage CSV; legacy alias | F3 |
| `P/F3_SimpleModels/F5C_water_contacts.py` | Same PEG models; complete AW water/frame observations; main | `output/PEG/F5C/F5C_water_contacts.{png,svg,csv}`, bins and coverage CSV | F3 |
| `P/F3_SimpleModels/F5A_exposed_vs_buried.py` | Cached PEG atomic table; AW/Power; main | `output/F5A/F5A_<metric>[_signed].{png,svg}`, exposure/contact plots and CSVs | F3 |
| `P/F3_SimpleModels/F5B_heparin_chemistry.py` | dp2/12/18/24/30; paired AW/Power chemistry, coordination and chain position; main | `output/Heparin/F5B/F5B_<metric>_<comparison>.{png,svg}`, chemistry summary CSV | F4 |
| `P/F3_SimpleModels/F5B_water_burial.py` | PEG/Heparin water geometry; AW/Power, validated PEG three-scheme membership; main | `output/F5B/`: metric PNG/SVG, contact histogram/counts, coverage/summary CSVs | Shared F3/F4; split decision pending |
| `P/F3_SimpleModels/F5_run_all.py` | PDB/scheme logs, cache, validation; main | Summary/interface CSVs, PKL caches, validation CSV/JSON/MD | Shared F3/F4; dispatches only older A/B panels |
| `P/F3_SimpleModels/{F5_config,F5_common,F5_peg_dataset}.py` | Identity, completeness, area, matching and screened layer membership helpers | Tables/configuration; screened status CSV consumed | Shared F3/F4; PEG-specific helper F3 |
| `P/F3_SimpleModels/F5C_rerun_missing.py` | Preview or explicitly run missing AW checkpoint exports | `output/PEG/F5C/reruns/` | F3 preparation; do not run during migration |
| `P/F3_SimpleModels/archive/F5[A-F]*.py` | Retired terminal/environment/length/topology/scheme/coordination analyses | Historical PNG/SVG and summary/validation tables | F3/F4 historical; topology SI candidate |
| `P/F4_atomic_deviations/F3A_multiscale_summary.py` | Atom/residue/molecule AW/Power summary; main | `F3A_multiscale_summary.{png,svg,pdf,csv}` | F6 |
| Same folder: `F3B_atomic_groups.py` | Protein/nucleic atomic grouping; main | `F3B_atomic_groups.{png,svg}`, caches, summaries/outliers/candidates CSV | F6 |
| Same folder: `F3C_residue_groups.py` | Protein/nucleic residue distributions; main | `F3C_*` plots/cache/summary/candidate CSV | F6 |
| Same folder: `F3E_solvent_exposure.py` | Atomic/residue deviations versus AW solvent-facing SA; main | `F3E_exposure_multiscale_ranges.{png,svg}`, slopes/cutoffs/cache/summary CSV | F6 |
| Same folder: `F3_common.py` | Shared matching, metric definitions, aggregation, discovery | In-memory data/tables | F6 shared helper; retain consistent definitions |
| Same folder: `F3D_geometry_examples.py`, `F3_build_figure.py` | Intended selections/compositor | None; docstring/comment-only | F6 drafts |
| Same folder: `F3_atomic_umap.py` | Standardized AW/Power descriptors, UMAP/clustering; main | Assignments/standardized/cluster-summary CSV, clusters/metrics PNG, report JSON | F6 exploratory; possible SI placement unresolved |
| Same folder: `F3_curvature_pca.py` | AW/Power atomic descriptors, correlation/PCA; main | `F3_atomic_curvature_pca/`: descriptor/correlation/PCA CSV, five PNG, report JSON | Curvature |
| Same folder: `F2C_All_Data_Sets_Clustered_full/{config,pipeline,run_single,run_batch}.py` | Split legacy grouping pipeline; single/batch main | PNG/SVG and statistics CSV | F6; canonical alternative unresolved |
| Same folder: `Extra/`, `Other/`, `F3_old/` | Atomic, residue, molecule, exposure, ligand, NCP and p53 analyses; see appendix | Many display-only; diagnostic/training/summary CSVs, selected PNG/SVG | F6 mixed/multiscale; frame utilities F2; curvature separate |
| `P/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme{,2,3}.py` | Compare atomistic/CG volume/SA; AW/Power; main | Display; variant 2 writes `residue_data.csv` | F7 tentative |
| Same folder: `color_residues.py`, `diff_by_res_vor_pow.py` | Residue coloring/instance comparison; main | Colored PDB; three `*_read_logs2.csv` diagnostics/summary | Needs Review: F6 residue aggregation versus F7 representation comparison |
| Nested `F7_Coarse_grained_residues/coarse_grained_residue_volume_difference.py` | Protein/nucleic CG residue grid; AW/Power/Primitive; main | Display by default; optional specified PNG | F7 tentative |
| Nested `F8_CG_Curvature/f8_build_surfaces.py` | Inter-residue/interface surface table; main/helper | `surf_df.csv` | Shared F6/F7 builder |
| Nested `F8_CG_Curvature/f8b_raincloud_interfaces.py` | AW-focused interface curvature raincloud; main | `figure8_outputs/plots/F8B_raincloud_interfaces*.png` | Curvature; F7 cross-reference |
| Nested F8 `f8a_plot_visualize_example.py`, `f8d_interface_heatmaps.py` | Empty slots | None | Drafts; assignment unresolved |
| Nested F8 `Old/*.py` | CG curvature exploration and 1BNA coloring | Display, colored PDB | Curvature; coloring F7 tentative |
| `P/Curvature/*.py` | Sequence localization, concentration, violin, alanine/area/chain/DNA-facing analyses | PNG/PDF/SVG, residue tables and summaries in existing generated subdirectories | Curvature |
| `P/Curvature/F4_atomic_curvature/` | Bin/peak assignment and interaction-energy analyses | `*_curvature_bins.csv`, `*_curvature_peaks.csv`, `combined_peaks.csv`, energy/curvature PNG/SVG/CSV | Curvature |
| `P/SuplimentalInformation/comparative_figures/F5_alpha_selection.py` | Persisted static comparison CSV; main | `F5a_alpha_selection_absolute`, `F5b_alpha_selection_retention`, PNG/PDF/SVG | F1 |
| `P/SuplimentalInformation/Curvature/{DulaunayDual,visualize_interface_criteria}.py` | Synthetic weighted Delaunay/Power/alpha and condition cartoons; main | `alpha_figures/01_*…05_*.png`, `criteria_figures/01_*…04_*.png` | F1 |
| Same folder: `CubeVertexGaussianCurvature.py` | Cube/flat-vertex Gaussian illustration; top-level | `{cube,cube_vertex,flat_vertex}.{png,svg}` | Curvature |
| `P/SuplimentalInformation/cazals*.py` | PDB/radii adapters, Power-alpha implementation, diagnostic and independent 2KAI audits | NPZ, filtration/facet/edge CSV, stage audit JSON/TXT; details below | Shared core adapters stay stable; diagnostics Curvature/SI candidate |
| `P/extra/*.py` | Mixed-system element comparisons, feature-table preparation and ML experiments | Display; `atomic_comparisons*.csv`, `pow_model_training.csv`, clustering CSV | F6 exploratory / SI candidate; mixed simple cohorts need explicit selection |
| `Part_III_ML_Processes/Surfaces_Training.py` | Incomplete synthetic-ensemble AW/Power descriptor exporter | Caller-selected CSV | SI candidate; experimental |
| `Part_I_Statistical_Ensembles/Fig1/` | Synthetic geometry/building-block illustrations | Mostly display; proximity GIF/optional PNG | F1 |
| Statistical `Fig3/`, `Fig5/`, `SupplementaryMaterial/`, `old/`, `Other/` | Radius/coordinate/overlap/density distributions and experiments | Mostly display; coordinate/overlap/density CSVs | SI candidates; method-comparison subset F1 |
| Statistical `Fig4/`, `Fig6/`, `Fig7/` | Scheme differences, neighborhood/sphericity/radial comparisons | Mostly display; colored PDB/TXT, selected PNG/SVG | F1; supplemental alternatives SI candidates |
| Statistical `VerifyVorPy/` | Independent volume/vertex comparisons | PNG/SVG and `all_info.csv` | SI validation candidate |

Historical 1BNA/DNA and Hairpin content exists in mixed-system scripts, but the simple-model folder has no equivalently organized dedicated DNA/Hairpin pipeline. F5 entry points must be selected deliberately, not inferred from an old "F5" prefix.

### Root-level analysis and publication outputs

- `comparative_geometry_study.py`: audit/compute/plots-only entry point for 2KAI, 1AK4 and 1BRC, comparing AW, Cazals alpha-0 Power, full Power and radius-matched Power. Writes static/trajectory/solvent comparison CSVs, provenance and case geometry.
- `comparative_study_figures.py`: called by that runner; reads persisted CSVs and writes `figures/figure_A_static_interface_area`, `figure_B_conventional_curvature`, `figure_C_signed_unsigned_curvature`, `figure_D_aw_decomposition`, `figure_E_power_aw_parity`, `figure_F_2kai_trajectory`, `figure_G_finite_solvent_effect`, `figure_H_method_effect_decomposition`, each PNG/PDF/SVG. Classification spans F1/F2/F6/Curvature. These letters remain historical output labels.
- `curvature_showcase_2kai.py`: saved static/solvated 2KAI AW/Power evidence; `static_interface_projection`, `alpha_dual_contacts_2d`, `finite_solvent_context_2d` PNG/SVG, curvature/GB CSV, PML/Python scenes and reports.
- `matched_probe_power_aw.py`, `open_interface_curvature.py`: matched geometry/topology/curvature CSVs, identity audits, cache/report/PML outputs; Curvature collection, shared API constraints apply.
- `dual_restriction_audit.py`, `cazals_validation.py`, `cazals_edge_forensics.py`, `aw_pipeline_audit.py`, `finite_solvent_boundary.py`: diagnostics and validation; SI candidates, not final SI assignments.
- `apollonius.py`, `apollonius_export.py`, `condition_b_geometry.py`, `power_interface_selection.py`: F1 mathematical/selection support.
- `aw_interface_curvature.py`, `power_interface_curvature.py`, `power_interface_cli.py`, `power_vs_aw.py`: shared Curvature APIs/exporters; not disposable figure scripts.
- `nonpolar_interface/{classification,table,explore,export}.py`: F6 interface table/export support used by production output code.

### Relevant external scripts and helpers

| Current source | Role / products | Proposed classification |
|---|---|---|
| `scripts/{plan_figure3b_recovery,validate_figure3b_recovery,smoke_figure3b_pymol}.py` | Recovery plan/validation tables and PyMOL PNG smoke render | F6 preparation |
| `scripts/render_pocket_preview.py` | Saved diagnostic mesh → `pocket_intact.png`, `pocket_cutaway.png`; optional VTK | F1/SI candidate |
| `scripts/{export_cached_2kai_interface_dual,export_cached_interface,export_interface_geometry_example}.py` | Saved network/cache → dual/interface geometry, PDB/CSV/JSON/PML | F1/F6 |
| `scripts/export_open_interface_geometry_examples.py` | Frozen matched-probe fixtures → geometry/curvature tables | Curvature |
| `scripts/{investigate_solvent_layers,investigate_peg_first_layer,audit_peg_boundary_patches,probe_peg_missing_water}.py` | Layer/geometry diagnostics, reports and PML; some modes build diagnostic geometry | F3 support / SI candidate |
| `scripts/{reexport_peg_layer1,validate_peg_layer1_bundle,summarize_peg_layer1_investigation}.py` | Recover/validate/index PEG layer inputs; status CSV | F3 preparation |
| `scripts/{audit_cached_interface_contacts,verify_2kai_pow_dry_pilot}.py` | Existing-export verification CSV/JSON/TXT | SI candidate |
| `scripts/{run_analysis_frames,run_checkpoint_layers}.py` | Launch jobs and exports | F2/F3 preparation; not read-only figure generators |
| `src/geometry/visualization/demo_2kai.py` | Loads AW and builds matched Power geometry; writes visualization checkpoint | F1/F6; not a cache-only renderer |
| `src/geometry/visualization/{export,interface,molecular_contact}.py` | Shared dual/interface/contact PML/Python/mesh/CSV/JSON exporters | Stable shared support |
| `src/geometry/{validation_shapes,topology_diagnostic}.py` | Synthetic shape/topology examples and presentation meshes | F1/SI support |
| `src/interface/{geometry_analysis,selected_export}.py` | Runtime interface analysis and scene export | Stable F6/Curvature support |
| `src/visualize/mpl_visualize.py` | Shared Matplotlib 3D primitives | Stable shared support |
| `src/visualize/VorpyIcon.py` | Lattice/icon experiment | Needs Review; branding rather than manuscript evidence |
| `src/output/{visualization,curvature_colors,color_tris,set_pymol_atoms}.py` | Runtime scene/color support | Stable shared support |

Benchmarks, performance instrumentation and NEXUS queue runners remain operational tooling outside the publication migration.

## 4. Dependencies and migration-sensitive references

### Data/import relationships

```text
PDB + AW/Power/Primitive logs
  → tools.compare.read_logs2
    → F3_common → multiscale/atomic/residue/exposure/UMAP/PCA
      → current F3B/F3E also import legacy F3_old clustered implementation
        → F2C_Helper_Auto_Cluster
    → F5_config + F5_common
      → F5_run_all → old exposed/buried A + water-burial B
      → PEG chemistry / Heparin chemistry
      → F5_peg_dataset ↔ lazy water-contact/water-burial validation imports
    → F8 surface builder → surf_df.csv → raincloud generator
    → legacy statistical/molecular tools

PDB/radii/group definitions
  → SI cazals_pdb + cazals_power_interface
    → weighted Power/alpha geometry and curvature
      → diagnostics, matched-probe study, comparative study
        → saved comparison CSVs → publication plots
```

`F3_common.discover_systems` expects immediate `<single-letter>_<molecule>` children with `aw/aw_logs.csv` and `pow/pow_logs.csv`; Primitive is optional. Several scripts exclude A/B/C by default. A renamed root or changed folder cohort can silently alter the scientific sample unless validated.

The screened PEG helper consumes `output/peg_one_layer_investigation_20261112/one_layer_run_status.csv`, per-row analysis directories and identity-checked `solvent/layer_01.pdb`. The prior README describes 69 usable three-scheme model/frame sets. That is documented cohort provenance, not a newly computed result in this checkpoint.

### Confirmed broken and fragile imports

| Reference/problem | Affected source |
|---|---|
| Removed `F3_Multiscale_Scheme_Differences.F3_common` | Current F3B, F3C, F3E, atomic UMAP and curvature PCA |
| Removed `F3_Multiscale_Scheme_Differences.F3_old.F2C_All_Data_Sets_Clustered` | Current F3B and F3E |
| Removed historical UMAP import | F3B `--umap` route |
| `src.analyze....F3_Multiscale_Scheme_Differences...F2C_Helper_Auto_Cluster` | Legacy monolithic clustered generator; wrong historical package and alternate package root |
| Removed top-level `Part_II_Molecular_Analysis.F8_CG_Curvature.f8_build_surfaces` | Current nested F8 raincloud |
| Removed `Part_II_Molecular_Analysis.F5.*` | Four embedded simple-model test modules and their internal imports |
| Explicit `__package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'` | Seven current simple-model modules; breaks direct-file bootstrap after folder move |
| Relative `.F5_config/.F5_common/.F5_run_all` inside `archive/` | Archived A–F generators resolve to absent archive siblings under normal package imports |
| Removed `Part_I_Statistical_Ensembles.Review1.Other.NeighborsAverageRadii` | `Fig7/Other/LocalDensityPolydispersity.py` |
| Absent historical `Data.Analyze.*`, `System.sys_funcs.*` | Statistical SI, volume verification and legacy density tools |
| Absent `build_edge_old` symbol | Building-block `Edges.py` and `Surfaces.py` |
| Absent `vorpy.src.analyze.read_logs2` | Fallback in both frame-curvature scripts; primary modern import exists |
| Bare `from F3_common import ...` | F3A; direct-file/package invocation differ |

### Shared code that must not move casually

| Shared analysis API | Production consumer |
|---|---|
| `tools/compare/read_logs2.py` parsing functions | `src/inputs/logs.py` |
| `nonpolar_interface` exports | `src/output/output.py` |
| `apollonius.py` | Dual builders, AW alpha code, CLI and scientific adapters |
| AW/Power interface curvature modules | Geometry interface builders, CLI and result adapters |
| `open_interface_curvature.py` | `src/interface/geometry_analysis.py` |
| `SuplimentalInformation/cazals_pdb.py` | Geometry interface diagnostics |
| `SuplimentalInformation/cazals_power_interface.py` | Geometry alpha-filtration implementation |
| `cazals_validation.power_facet_polygon_area` | AW/Power comparison implementation |

Coordinate shared-interface changes with FORGE, GAUSS, ATLAS and VECTOR owners before incompatible changes. Moving serializable public classes or altering module-name strings can also affect saved-result loading; archives were not opened merely to test that risk.

### Paths, execution and outputs

- Fixed-depth `sys.path` adjustments are already wrong in many moved/nested folders, especially `F3_old`, `Extra`, atomic curvature and F7/F8.
- Multiple statistical/batch tools use `os.chdir('../../../..')`; some modules launch Tkinter or plotting at import time.
- Hardcoded inputs include `E:\Molecular`, E-drive manuscript/OneDrive roots, current-user Dropbox paths, `C:/Users/jacke/Documents/data1/`, and `/home/jack/PycharmProjects/foam_gen/...`.
- Legacy batch ledgers use `Data/Analyze/tools/batch/loggy_woggys*.txt`; one radius script uses `../../tools/batch/loggy_woggys_old.txt`.
- Outputs may be source-relative, selected by dialog, caller-supplied, or current-directory relative. `CubeVertexGaussianCurvature.py` writes in CWD; `DulaunayDual.py` uses relative `alpha_figures/`; criteria cartoons anchor beside their source.
- Moving simple-model code changes its default output/cache roots. Cache signatures include source/input metadata and hashes; do not assume cache reuse across relocation.
- Many old plotting scripts only display figures. `plot_templates/scatter.py` accepts `save` but does not save; line/histogram/table/inline helpers likewise do not establish automatic exported filenames.
- `CleanData` and `batch/wipe_user_data.py` include move/delete/rename operations. They are not safe discovery or import-smoke targets.
- The current simple-model direct-execution test globs every root Python file. Only two generators receive non-generating self-check execution; newer generators can run normally. Do not run that test unchanged during a no-regeneration migration.

### Dependencies and available inputs

The prior environment check found scikit-learn, seaborn, GUDHI, Pillow and Tkinter. It did not find `umap`, `hdbscan`, `kmodes`, Plotly, matplotlib-venn or PyMOL in the project interpreter. Some are optional branches; absence does not mean all related scripts fail. Optional VTK rendering was not availability-tested. Module discoverability does not prove ABI/version compatibility.

Several analysis dependencies are absent from core requirements; GUDHI is an optional project extra. Installing/changing dependencies requires separate approval.

The PEG status CSV and representative 2KAI, 1AK4, 1BRC, 1UDI and solvated-2KAI PDBs exist. The sampled default `output/comparative_study/cases/2KAI/H/full_aw.vpy` and matched-probe Power cache path were absent in this checkout. External Dropbox datasets were not fully validated.

### Tests and documentation with old paths

- `vorpy/tests/analyze/test_f3_atomic_umap.py`: removed multiscale imports at lines 5 and 90.
- `docs/analysis/figure3_atomic_umap.md:14`: removed generator path.
- Embedded simple-model classification/loading/scheme/area tests: removed F5 package.
- Active and archived simple-model READMEs: old F5 commands and old figure identity.
- `F4_atomic_deviations/SCRIPT_MAP.txt`: old intended `F3_multiscale_differences/` location and `Figure_3` asset names.
- Root pytest configuration targets `tests`; relevant tests are under `vorpy/tests` and embedded analysis folders. Use explicit focused paths after approval.

## 5. Duplicate, obsolete and scientific-review candidates

Seven statistical Fig7 pairs are identical after text decoding/newline normalization:

1. `3_Radius_UnsharedNeighbors.py` / `B_Radius_UnsharedNeighbors.py`
2. `7_RadialDistance_VolDiff.py` / `F_RadialDistance_VolDiff.py`
3. `ExclusiveNeighbors_SphericityDifference.py` / `D_ExclusiveNeighbors_SphericityDifference.py`
4. `6_Overlap_Volume.py` / `H_Overlap_Volume.py`
5. `5_pizza_pie.py` / `G_pizza_pie.py`
6. `2_NeighborDistanceSphericity.py` / `Other/NeighborDistanceSphericity.py`
7. `1_SphericityDifference_NerighborRadCVDiff.py` / `Other/SphericityDifference_NerighborRadCVDiff.py`

Fig1 `Surfaces/Fill_Surface_points.py` and its `Building_Blocks` copy are identical; the two vertex `doublets.py` copies have identical ASTs. Other building-block copies are overlapping alternatives, not all identical.

Compare before selecting a canonical implementation: three atomic-correlation scripts, two frame-curvature scripts, three CG-total scripts, monolithic/split/new atomic-grouping code, residue-deviation generations, alanine variants, and `CombineLogs{,2}.py`. Duplication is not deletion authorization.

Retired A–F polymer generators are explicitly documented as retired. Empty/docstring-only slots include sphere schemes, F8 example/heatmap, archived `SolventExposedDiff.py`, F3 geometry/compositor, frame-network dynamics and `PurgeLogs.py`. `CreateFinalFoamData.py` has a `pass` function; foam verification is a scaffold. `Surfaces_Training.py` uses its atom writer after file closure and has only an empty surface header. These are unfinished, not completed publication analyses.

Scientific issues to resolve separately from mechanical moves:

- `FrameAnalysis.py:149–150` uses current AW volume in Power/Primitive frame-deviation numerators while dividing by each scheme's baseline. Review before declaring a canonical F2 implementation.
- Its difference sign differs from the newer `100 * (Power - AW) / AW` definition.
- Preserve units: atomic volume Å³, surface area Å², integrated mean-curvature quantities Å, and raw angle-weighted edge quantities Å·rad as applicable. Do not combine them as one "curvature" measure.
- Preserve signed versus absolute deviations, missing/zero-denominator handling, identity matching, completeness and area-partition checks.
- Legacy `vor/del` versus `aw/prm` labels need source verification, not blind text replacement.
- Atomic/residue sums, CG cells and external molecular union surfaces are distinct measurements.
- Pooled atom/water/frame observations are correlated descriptive observations, not independent replicates or persistence evidence.
- Polymer-surrounded water is not proof of enclosure/accessibility; archived internal-water occupancy is explicitly unavailable.
- No scientific reference values or tolerances are changed by this plan.

Naming issues: `SuplimentalInformation`, `DulaunayDual.py`, `residue_flexibilty...`, `52ca_...` versus `5c2a_...`, `Nerighbor`, numeric-leading import names, inconsistent `Extra/extra/Other/Old/old/archive`, and old figure prefixes inside newly numbered folders. Correcting spelling can still break imports, commands and saved provenance.

## 6. Planned FIGY coordination — Dropbox artwork

No FIGY communication or Dropbox edit has been performed or implied by this checkpoint.

CODY owns code provenance, generator classification, dependency repair and future migration validation. FIGY owns locating, reviewing and organizing Dropbox artwork/source design files. Both use the approved figure registry above; neither assigns F3–F7 panel letters independently.

Before any move, prepare a jointly reviewed mapping with:

- Existing code path and current artwork path/name.
- Existing historical figure/panel label, approved destination figure, and panel status (`UNASSIGNED` until agreed).
- Molecular system, scheme/radii/probe definitions, frame/cohort, input evidence and units.
- Generator commit, data/output hashes when available, rendered export filename/format and editable artwork source.
- Whether an artwork file is a direct generator product, edited composite, or historical/unverified variant.
- Owner, dependencies, preservation copy/location and approval status.

FIGY should inventory editable originals, composites, linked/raster assets, fonts and export versions without overwriting originals. CODY should supply code/output provenance and identify stale aliases. Resolve artwork/code mismatches before filenames change. A future approved artwork migration needs verified preservation copies/version history and checked linked assets; no Dropbox cleanup is authorized now.

## 7. Pure-water GROMACS baseline — PLANNED, NOT COMPLETED

F3 reserves scope for a pure-water baseline paired with PEG hydration. The previous audit did not establish a completed pure-water GROMACS pipeline, validated baseline dataset, plots, or a canonical generator. Do not label this baseline as reproduced or complete. The audit did not comprehensively search external simulation stores.

Before implementation, agree on:

1. Water model, force field, simulation provenance, ensemble, temperature/pressure, box size and periodic-boundary treatment.
2. Equilibration/production windows, selected frames, sampling interval and correlation/replicate interpretation.
3. Coordinate units and conversions; atom/water identity and wrapping/centering policy.
4. AW/Power/Primitive scope, radii, hydrogen representation, periodic or finite-box geometry, completeness and boundary exclusions.
5. Measurements needed to support PEG comparisons, with identical definitions and denominators where comparisons are claimed.
6. Appropriate validation: identity/coverage, nonfinite/missing values, geometric area accounting, periodic-volume accounting if mathematically applicable, existing tolerance policy and descriptive statistical checks.
7. Approved compute budget, data location, responsible owner and reproducible command record.

A baseline pilot and validation are future authorized work. No simulation, molecular solve, dataset download or figure regeneration is part of this checkpoint.

## 8. Future organization and migration strategy

Proposed destination names below are PLANNING ONLY. F7 and SI remain tentative. The existing shared APIs listed above stay in place initially.

```text
vorpy/src/analyze/publication/
├── F1_methodology_geometry/
├── F2_frame_persistence/
├── F3_water_baseline_peg/
│   ├── water_baseline/          [planned work; no completed results implied]
│   └── peg/
├── F4_heparin/
├── F5_dna_hairpin/
├── F6_complex_multiscale/
├── F7_coarse_graining/          [tentative]
├── curvature/
├── SI/                         [membership TBD]
└── archive/                    [preserved history; no automatic deletion]
```

### Migration sequence and gates

1. **Resume safely:** read this document, check current branch/index/status, identify owners of pre-existing relocations, and compare current files with this audit snapshot. Do not rebuild another agent's unfinished work.
2. **Approve mapping:** settle canonical generators, F3–F7 boundaries and artwork mapping. No final panel letters yet.
3. **Repair current breakage first:** agree a bounded milestone for old imports/direct-file bootstraps/tests/docs at their present locations; preserve numerical behavior and datasets.
4. **Protect shared APIs:** explicitly list public imports, serialized-result expectations, CLI and GUI consumers. Coordinate incompatible FORGE/GAUSS/ATLAS/VECTOR changes.
5. **Record a baseline:** capture approved source/data/artifact hashes, invocation, environment versions, cohort membership and output inventory. Use existing outputs; do not regenerate merely to prepare migration.
6. **Move one family per approved milestone:** update all callers, test paths, docs, direct-file/package execution and explicit output configuration together. Keep historical output filenames initially unless separately approved.
7. **Validate:** static AST/import-reference checks first, then carefully scoped tests using the repository interpreter. Do not import top-level mutating/dialog/generating modules indiscriminately.
8. **Scientific comparison:** if a bounded rerun is authorized, compare identical inputs/settings/cohorts before/after; check identity, completeness, scheme/status semantics, units, tables and topology. Image appearance alone is not validation.
9. **FIGY review:** verify rendered-to-artwork provenance and linked assets without changing scientific values.
10. **Checkpoint/review:** one focused commit only with approval; record paths changed, checks, limitations and rollback instructions. Do not begin the next family without authorization.
11. **Output renaming and retirement later:** separate approved milestones for output-root changes, final panel naming and obsolete-code disposition.

Use the required `.venv/bin/python -m pytest -q <explicit focused paths>` environment convention where available. On this Windows checkout the verified interpreter exists at `.venv/Scripts/python.exe`; any future Windows invocation should explicitly use that verified executable under the repository's environment rules. No bare pytest/Conda assumptions and no dependency installation without approval.

Validation must preserve scientific checks rather than weakening tolerances. Safe tests include bounded synthetic geometry and existing focused molecular regression cases; full p53tet, expensive benchmarks and large exports require explicit approval. The direct-execution test must be reviewed/restricted before running because it can regenerate newer panels.

### Preservation and rollback

- Preserve all tracked/untracked work, original logs, PDBs, checkpoints, caches, raw exports, publication images, editable artwork and reports.
- Never use `git reset --hard`, `git clean -fd`, forced push, broad staging or recursive deletion to implement or undo this migration.
- Before a future approved move, agree a recoverable source/artifact snapshot with checksums. Creating backups outside this repository requires approval.
- Keep original output directories and scientific artifacts unchanged during the initial code migration; choose a separate approved output directory for validation reruns.
- Record an explicit old→new path manifest and every caller/doc/test change. Avoid moving public serializable classes unless a compatibility plan is validated.
- Roll back a committed migration via a reviewed inverse change or approved revert commit, preserving concurrent edits. Git reverts and filesystem restores require the applicable authorization; do not execute them automatically.
- Roll back artwork through FIGY's approved preserved versions/history, not by overwriting unrelated Dropbox changes.
- Confirm restored hashes and focused validations; retain failed-migration evidence for diagnosis.
- This checkpoint commit contains only this plan and cannot restore pre-existing uncommitted relocations or scientific data. Do not treat it as a repository/data backup.

## 9. Unresolved decisions and next action

- F3–F7 panel assignments, counts, layouts and canonical artwork are UNASSIGNED.
- Pure-water baseline design, input availability, ownership, GROMACS provenance and validation are unresolved; baseline is planned.
- Choose dedicated F5 DNA/Hairpin entry points from mixed historical analyses without silently changing cohorts.
- Decide F6 residue/multiscale versus F7 true-CG boundaries and whether F7 remains a main figure.
- Decide how mixed runners/panels are cross-referenced or split; no automatic duplication of scientific calculations.
- Select authoritative atomic-correlation, frame-curvature, atomic grouping, CG-total and alanine variants.
- Resolve the frame-deviation formula/sign issue through scientific review, separately from relocation.
- Decide historical output-name compatibility, cache treatment and source-independent output roots.
- Determine SI membership and the relationship of statistical-ensemble figures to this manuscript.
- Determine whether retired scripts must remain runnable or are historical source only.
- Coordinate the already-present working-tree moves with their owner before CODY touches them.

**Exact next step:** after this checkpoint is safely committed, obtain authorization for one bounded import/path reconciliation milestone and jointly approve the figure-to-generator/artwork mapping with FIGY. Do not start migration, baseline implementation, panel numbering or artifact regeneration solely because this plan exists.

## 10. Validation and limitations of this checkpoint

This is a documentation-only change. Validate document presence/required registry and status statements, stage only this path, inspect the staged name list, commit with the requested message, and verify commit hash and file membership. Record branch and unrelated status in the user-facing completion report. Do not push.

Prior audit validation: 321 successful AST parses, read-only module/path/reference checks, optional-module discoverability, duplicate comparisons and unchanged analysis-source hashes. No pytest, generators, geometry recomputation, figure rendering, exhaustive external-data validation or Dropbox inspection. Missing external context and dynamic imports remain limitations; static evidence is not scientific reproducibility.

## Appendix A. Exhaustive analysis Python-source inventory

The following table is derived from the prior audit's retained static metadata, not from generator execution. All paths are relative to `vorpy/src/analyze/`. "Source summary" uses the module docstring or declared function/class names when no docstring exists; it is evidence of implemented/intended scope, not a claim that the code runs. "File literals" includes input and output references and is deliberately not mislabeled as a list of generated results. The principal inventory above identifies known outputs. "Main" means an explicit `__main__` guard; "module" can mean helper or top-level script and is not safe-to-import certification. Collection labels are proposals under the approved registry; SI remains a candidate collection.


| Current relative path | Entry | Source summary / declared functions | System tokens observed | Scheme tokens observed | File literals (inputs and/or outputs) | Proposed collection |
|---|---|---|---|---|---|---|
| `__init__.py` | module | Package initializer / exports | — | — | — | Shared/support; review |
| `apollonius_export.py` | module | CSV and text exports for Apollonius and AW alpha complexes. | — | aw | _simplices.csv; _summary.txt; apollonius_edges.csv; apollonius_summary.txt; apollonius_tetrahedra.csv; apollonius_triangles.csv; apollonius_vertices.csv | F1 shared (retain) |
| `apollonius.py` | module | Combinatorial straight-line Apollonius dual of a solved VorPy network. | — | aw, voronoi | — | F1 shared (retain) |
| `aw_interface_curvature.py` | module | Bicolor interface curvature on existing additively weighted networks. | — | aw | aw_interface_edges.csv; aw_interface_summary.txt; aw_interface_surfaces.csv | Curvature (shared APIs retained) |
| `aw_pipeline_audit.py` | module | Forensic comparison of full and interface-filtered AW networks. | — | aw | — | SI candidate / shared |
| `cazals_edge_forensics.py` | module | Forensic, non-production audits of Cazals power-interface edge sets. | 2kai, protein | power | 2KAI_13_alpha_triangle_discrepant_edges.csv; 2KAI_condition_b_audit.csv; 2KAI_curvature_forensic_report.txt; 2KAI_edge_set_comparison.csv; 2KAI_incidence_audit.csv; 2KAI_interface_edge_forensics.csv; 2KAI_pattern_curvature.csv; 2KAI_sign_convention_audit.txt | SI candidate / shared |
| `cazals_validation.py` | Main | Compact validation suite for the weighted Power/Cazals interface analysis. | 2kai, protein | aw, power | .pdb; _beta_distribution.csv; _beta_distribution.png; _beta_histogram.csv; input_inventory.json; lys15_curvature_edges.csv; lys15_interface_facets.csv; power_interface_components.csv; power_interface_edge_incidence.csv; power_interface_facets.csv; validation_metrics.json; validation_summary.txt | SI candidate / shared |
| `comparative_geometry_study.py` | Main | Reproducible comparative-study preflight and atom-model construction. | 2kai, protein | aw, power | *.pdb; */H/metrics.csv; .json; .pdb; H.pdb; H/metrics.csv; H/selection_membership.csv; HS_context.pdb; alpha0_vs_full_power.csv; atom_identity_matching.csv; atom_matching_details.json; atom_model_manifest.csv; comparative_study_summary.txt; complete.json; excluded_geometry.csv; excluded_geometry_audit.csv; facet_selection.csv; frame_0100ps.pdb; full_aw.vpy; full_aw_cache.json; group_definitions.json; input_audit.csv; lys15_localized_region.json; measured_edges.csv; method_definition_audit.json; method_effect_decomposition.csv; metrics.csv; paired_interfaces.pml; phase_status.json; power_interface_*.csv; power_vs_aw_effects.csv; protein_heavy.pdb; provenance.json; radii_audit.csv; selection_membership.csv; solvent_boundary_audit.csv; solvent_excluded_geometry.csv; solvent_paired_case_metrics.csv; solvent_paired_effects.csv; static_method_comparison.csv; trajectory_metrics.csv; trajectory_summary.csv; validation.csv; viz_radii_audit.csv | Needs Review: F1/F2/F6/Curvature |
| `comparative_study_figures.py` | module | Publication plots generated only from persisted comparative-study CSVs. | 2kai | aw, power | .pdf; .png; .svg; figure_*.png; solvent_paired_effects.csv; static_method_comparison.csv; trajectory_metrics.csv; trajectory_summary.csv | Needs Review: F1/F2/F6/Curvature |
| `condition_b_geometry.py` | module | Independent mathematical diagnostics for Cazals condition beta. | — | power | — | F1 shared (retain) |
| `curvature_showcase_2kai.py` | Main | Generate an evidence-led 2KAI interface curvature/visualization checkpoint. | 2kai, protein | aw, power, voronoi | alpha_bicolor_contact_dual_edges.csv; alpha_contact_patch_status.txt; alpha_dual_contacts_2d.png; alpha_dual_contacts_2d.svg; aw_edge_curvature_components.csv; aw_surface_curvature_components.csv; bicolor_pair_comparison.csv; boundary_curvature_components.csv; camera_metadata.json; comparison_summary.json; curvature_comparison.csv; curvature_comparison_report.txt; finite_solvent_context.pml; finite_solvent_context_2d.png; finite_solvent_context_2d.svg; gauss_bonnet_audit.csv; gaussian_components_cache.pkl; input_identity_audit.csv; local_example.json; local_example.pml; output/comparative_study/cases/2KAI/H/aw/aw_interface_edges.csv; output/comparative_study/cases/2KAI/H/aw/aw_interface_surfaces.csv; output/comparative_study/cases/2KAI/H/full_aw.vpy; power_edge_curvature_components.csv; power_facet_area_cache.json; power_geometry_cache.pkl; power_surface_curvature_components.csv; shared_pair_surfaces.csv; showcase_layers.pml; showcase_metadata.json; solvated_atom_audit.csv; solvated_input_audit.json; static_aw_alpha.pml; static_interface_projection.png; static_interface_projection.svg; static_overlay.pml; static_power_alpha.pml; total_interface_metrics.csv; vertex_gaussian_components.csv; vorpy/data/2KAI.pdb; vorpy/data/2KAI_MD_centered.pdb | Curvature (shared APIs retained) |
| `dual_restriction_audit.py` | Main | Read-only audit of molecular filtering stages on dual connections. | 2kai | aw, power | 2KAI_aw_full_connections.csv; 2KAI_aw_restriction.pml; 2KAI_aw_restriction_metadata.json; 2KAI_dual_restriction_stages.pml; 2KAI_dual_restriction_summary.json; 2KAI_dual_restriction_summary.txt; 2KAI_longest_connections_by_stage.csv; 2KAI_longest_full_regular_connections.csv; 2KAI_power_connections.csv; 2KAI_power_distance_summary.csv; 2KAI_probe_radius_control.csv; 2KAI_probe_radius_control_summary.json; radii_audit.csv; visualization_metadata.json; vorpy/data/2KAI.pdb | SI candidate / shared |
| `finite_solvent_boundary.py` | module | Finite explicit-solvent hull-distance audit for sampled interface geometry. | — | aw | — | SI candidate / shared |
| `matched_probe_power_aw.py` | Main | Controlled 2KAI matched-probe Power/AW restriction audit. | 2kai, protein | aw, power | aw_pair_birth_checkpoint.csv; bicolor_pair_comparison.csv; comparison_report.txt; comparison_summary.json; connection_distances.csv; disagreement_pairs.csv; full_dual_pairs.csv; input_identity_audit.csv; matched_probe_layers.pml; output/comparative_study/cases/2KAI/H/full_aw.vpy; output/comparative_study/cases/2KAI/H/radii_audit.csv; power_facet_area_cache.json; power_geometry_cache.pkl; restricted_pairs.csv; shared_pair_curvature.csv; shared_pair_surfaces.csv; total_interface_metrics.csv; visualization_metadata.json; vorpy/data/2KAI.pdb | Curvature (shared APIs retained) |
| `nonpolar_interface/__init__.py` | module | Package initializer / exports | — | — | — | F6 shared (retain) |
| `nonpolar_interface/classification.py` | module | Explicit, versioned chemical and interface labels for analysis tables. | — | — | — | F6 shared (retain) |
| `nonpolar_interface/explore.py` | module | Descriptive, system-aware summaries for the first geometry comparison. | — | — | — | F6 shared (retain) |
| `nonpolar_interface/export.py` | module | Stable file export for nonpolar/interface geometry analysis tables. | — | — | .csv; _metadata.json | F6 shared (retain) |
| `nonpolar_interface/table.py` | module | Build local geometry analysis tables without mutating VorPy networks. | — | — | — | F6 shared (retain) |
| `open_interface_curvature.py` | Main | Topology and curvature accounting for open, face-based interfaces. | 2kai | aw, power | aw_curvature_components.csv; aw_edge_curvature_components.csv; aw_supported_face_GB_audit.csv; aw_surface_curvature_components.csv; aw_topology_audit.csv; bicolor_pair_comparison.csv; input_identity_audit.csv; open_interface_curvature_report.txt; orientation_reversal_audit.csv; output/comparative_study/cases/2KAI/H/full_aw.vpy; power_curvature_components.csv; power_edge_curvature_components.csv; power_facet_area_cache.json; power_geometry_cache.pkl; power_independent_crease_audit.csv; power_topology_audit.csv; synthetic_validation.csv; total_interface_metrics.csv | Curvature (shared APIs retained) |
| `Part_I_Statistical_Ensembles/Fig1/AWVd_vs_Power/Surfaces.py` | module | Top-level script; see principal family inventory | — | pow | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/AWVd_vs_Power1/Surfaces.py` | module | Top-level script; see principal family inventory | — | pow | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/Edges.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/ellipse_edge.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/Inscribed_Circle.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/loop_edge.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Encapsulation.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Fill_Surface_points.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Plane_Between_Centers.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Surfaces.py` | module | Top-level script; see principal family inventory | — | aw, pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Vertices/doublets.py` | module | Top-level script; see principal family inventory | — | aw | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Vertices/TangentialCircle.py` | module | Declared: plot_circles | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Vertices/Vertices.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Cells/cells.py` | module | Top-level script; see principal family inventory | — | aw | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Edges/Edges.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Edges/ellipse_edge.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Edges/Inscribed_Circle.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Edges/loop_edge.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/proximity.py` | module | Declared: plot_proximity | — | aw, pow, power, prm | .png; proximity.gif | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/SortBalls.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Encapsulation.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Fill_Surface_points.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Plane_Between_Centers.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Surfaces.py` | module | Top-level script; see principal family inventory | — | aw, pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Vertices/doublets.py` | module | Top-level script; see principal family inventory | — | aw | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig1/Vertices/Vertices.py` | module | Top-level script; see principal family inventory | — | — | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/DistributionOfXYZ.py` | module | Top-level script; see principal family inventory | — | — | location_data.csv | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/histofradiiCompare.py` | Main | Declared: calculate_cdf, inverse_transform_sampling, get_hist_of_radii | — | aw, pow | ../../tools/batch/loggy_woggys_old.txt | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/HistofRadiiRaw.py` | module | Declared: calculate_cdf, inverse_transform_sampling, gamma | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/NumOlaps_csv.py` | Main | Declared: parse_data, plot_overlapping_line_plots | — | power | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/NumOverlaps.py` | Main | Declared: num_olaps, plot_olap_percent | — | aw, pow | pbc_olap1.csv | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/OneOlapPlot.py` | Main | Declared: format_ticks, distribution_of_overlaps | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/OneXYZPlot.py` | module | Declared: convert_point | — | — | location_data1.0.csv | SI candidate |
| `Part_I_Statistical_Ensembles/Fig3/ranging_pdfs.py` | module | Declared: lognormal, gamma, weibull, physical_DeVries, physical_Ranadive_Lemlich, physical_GalOr_Hoelscher, plot_function1, plot_function | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig4/color_by_diff.py` | module | Top-level script; see principal family inventory | — | aw, pow, power | _set_diff.txt; _vor_diff_colored.pdb | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig4/NetNeighbors.py` | Main | Declared: get_net_neighbors, plot_data | — | aw, pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig4/overlaid_line_plots.py` | module | Top-level script; see principal family inventory | — | pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig4/plot_distributions_complete_cells.py` | module | Declared: lognormal, gamma, weibull | — | aw, pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig4/RadiusVsPercentDifference.py` | module | Top-level script; see principal family inventory | — | aw, pow, power | _vor_diff_colored.pdb | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig5/overlaid_line_plots_physical.py` | module | Top-level script; see principal family inventory | — | pow, power | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig5/physical_pdfs.py` | module | Declared: normalize_pdf, p1_unnormalized, p2_unnormalized, p3_unnormalized, calculate_cdf, inverse_transform_sampling | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig5/plot_hist_of_real_dists.py` | module | Declared: p1, p2, p3, gamma, calculate_cdf, inverse_transform_sampling | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig5/plot_physical_mpl.py` | module | Top-level script; see principal family inventory | — | pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/Fig6/DivergenceVsSphericity.py` | Main | Declared: select_folder, get_files, read_logs_data, parse_pdb_file, plot_spcy_results, plot_radii_vs_sphericity, display_metrics, process_all_subfolders … | — | aw, pow, power | .csv; .pdb | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig6/overlaid_line_plots_extra_data.py` | module | Top-level script; see principal family inventory | — | aw, pow | *.png; *.svg; .png | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/1_SphericityDifference_NerighborRadCVDiff.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/2_NeighborDistanceSphericity.py` | Main | This file plots the sphericity difference as a function of std neighbor distance and neighbor polydispersity | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/3_Radius_UnsharedNeighbors.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/4_MissingNeighborDistance.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/5_pizza_pie.py` | Main | This script calculates the volume difference between the additively weighted and power surfaces as a function of radial distance | — | aw, pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/6_Overlap_Volume.py` | Main | Calculates the volume of overlap between the power and the additively weighted surfaces | — | aw, pow, power | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/7_RadialDistance_VolDiff.py` | Main | This script calculates the volume difference between the additively weighted and power surfaces as a function of radial distance | — | aw, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/B_Radius_UnsharedNeighbors.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/C_MissingNeighborDistance.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/D_ExclusiveNeighbors_SphericityDifference.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/ExclusiveNeighbors_SphericityDifference.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/F_RadialDistance_VolDiff.py` | Main | This script calculates the volume difference between the additively weighted and power surfaces as a function of radial distance | — | aw, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/G_pizza_pie.py` | Main | This script calculates the volume difference between the additively weighted and power surfaces as a function of radial distance | — | aw, pow, power | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/H_Overlap_Volume.py` | Main | Calculates the volume of overlap between the power and the additively weighted surfaces | — | aw, pow, power | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/Other/BallRadius_AbsoluteVolumeDifference.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/Other/LocalDensityPolydispersity.py` | Main | Declared: sort_balls, calculate_local_density, calculate_local_polydispersity, calculate_std_neighbor_distance | — | aw | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/Other/NeighborDistanceSphericity.py` | Main | This file plots the sphericity difference as a function of std neighbor distance and neighbor polydispersity | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/Other/NeighborsAverageRadii.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | — | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/Other/Radius_NeighborCount.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/Fig7/Other/SphericityDifference_NerighborRadCVDiff.py` | Main | Declared: select_file, get_data, plot_data | — | aw, pow | /aw/aw_logs.csv; /balls.pdb; /pow/pow_logs.csv | F1 / SI candidate |
| `Part_I_Statistical_Ensembles/old/Check_Max_Vertex_Parameters.py` | module | Top-level script; see principal family inventory | — | aw, pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/old/check_num_foams.py` | module | Top-level script; see principal family inventory | — | — | /home/jack/PycharmProjects/foam_gen/Data/user_data | SI candidate |
| `Part_I_Statistical_Ensembles/old/density_deviation.py` | module | Declared: sqrt_model | — | — | /density_adjustments.txt | SI candidate |
| `Part_I_Statistical_Ensembles/old/Open_density_factor.py` | module | Declared: get_box_vols | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/old/plot_physical_plotly.py` | module | Top-level script; see principal family inventory | — | pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/old/plot_surf_graphy_plotly.py` | module | Top-level script; see principal family inventory | — | aw, pow, power | — | SI candidate |
| `Part_I_Statistical_Ensembles/old/StandardDeviationFoamGens.py` | module | Top-level script; see principal family inventory | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Other/check_deviation.py` | module | Declared: calculate_standard_error, calculate_standard_deviation, coefficient_of_variation | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/Other/DensityVoxelization.py` | Main | Declared: estimate_bounding_box_overlap, find_sub_box_density | — | aw, pow | densy_wensities.csv | SI candidate |
| `Part_I_Statistical_Ensembles/Other/Overlaid_Line_Plots_Density_Voxelization.py` | module | Top-level script; see principal family inventory | — | power | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/MissingBallsHeatMap.py` | Main | Declared: num_missing_balls, plot_heatmap | — | aw, pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/MonteCarloPlacement.py` | module | Declared: is_point_inside_circles, can_place_circle | — | — | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/NeighborsMatrix.py` | Main | Declared: get_neighbors_data, plot_neighbor_data | — | aw, pow, power | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/OverlapsMatrix.py` | Main | Declared: num_olaps_folder, parse_data, distribution_of_overlaps | — | aw, pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/RadiusDistributionMatrix.py` | Main | Declared: gamma, get_syses1, get_syses, distribution_of_overlaps | — | aw, pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/SurfaceOverlapsMatrix.py` | Main | Declared: get_overlap_data, distribution_of_overlaps | — | aw, pow | overlaps_1.csv | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/XYZMatrix.py` | Main | Declared: get_syses, distribution_of_overlaps | — | aw, pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/XYZMatrixFake.py` | Main | Declared: get_syses, distribution_of_overlaps | — | aw, pow | — | SI candidate |
| `Part_I_Statistical_Ensembles/VerifyVorPy/VerifyVerts.py` | Main | Declared: calc_dist, get_information, compare_vertices, plot_comparison, run_it_all | — | aw | .pdb; .png; .svg; all_info.csv; aw_verts.txt; vertices.txt | SI candidate |
| `Part_I_Statistical_Ensembles/VerifyVorPy/VerifyVerts1.py` | Main | Declared: calc_dist, get_information, compare_vertices, plot_comparison, run_it_all | — | aw | .pdb; .png; .svg; V_Vertices.txt; all_info.csv; aw_verts.txt; vertices.txt | SI candidate |
| `Part_I_Statistical_Ensembles/VerifyVorPy/VerifyVolumes.py` | Main | Declared: get_information, compare_vertices, plot_comparison, run_it_all | — | aw | .pdb; .png; .svg; V_Vertices.txt; aw_verts.txt; vertices.txt | SI candidate |
| `Part_II_Molecular_Analysis/__init__.py` | module | Empty file | — | — | — | Shared/support; review |
| `Part_II_Molecular_Analysis/Curvature/52ca_alanine_curvature_vs_area.py` | Main | Declared: Config, choose_logs_file, build_residue_instances, plot_residue_scatter, main | — | — | *.csv; alanine_curvature_vs_area.csv; alanine_curvature_vs_area.pdf; alanine_curvature_vs_area.png | Curvature |
| `Part_II_Molecular_Analysis/Curvature/52ca_alanine_dna_facing_violin.py` | Main | Declared: Config, choose_logs_file, sanitize_filename, parse_neighbors_cell, build_residue_instances_from_atoms, label_residue_dna_facing, plot_violin_two_groups, main | dna, rna | — | *.csv; _dna_facing_table.csv; _dna_facing_violin.pdf; _dna_facing_violin.png | Curvature |
| `Part_II_Molecular_Analysis/Curvature/5b_color_pdbs_by_curvature.py` | Main | Declared: get_atom_dict, color_pdb_by_curvature, write_pymol_script | — | — | _colored_by_curvature.pdb; _set_colors.pml | Curvature |
| `Part_II_Molecular_Analysis/Curvature/5c1_curvature_localization_along_sequence.py` | Main | Declared: CurvConfig, residue_to_one_letter, choose_logs_file, moving_average, percentile_clip, build_residue_df_from_atoms, aggregate_across_chains, expand_to_full_residue_range … | protein | — | *.csv; _panel_c1.pdf; _panel_c1.png; _panel_c1.svg; _panel_c1_residue_table.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/5c2_curvature_concentration.py` | Main | Declared: residue_to_one_letter, CurvConfig, choose_logs_file, sanitize_filename, _rankdata, spearman_rho, build_residue_df_from_atoms, aggregate_across_chains … | 5c2a, dna, protein, rna | — | *.csv; _panel_c2A_area_coupling.csv; _panel_c2A_area_coupling.pdf; _panel_c2A_area_coupling.png; _panel_c2A_area_coupling.svg; _panel_c2_lorenz.pdf; _panel_c2_lorenz.png | Curvature |
| `Part_II_Molecular_Analysis/Curvature/5c2a_alanine_curvature_vs_area_by_chain.py` | Main | Declared: Config, choose_logs_file, sanitize_filename, build_residue_instances, per_chain_summary, plot_residue_scatter_by_chain, main | — | — | *.csv; _chain_summary.csv; _curv_vs_area_by_chain.pdf; _curv_vs_area_by_chain.png; _points.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/5c2a_alanine_dna_facing.py` | Main | Declared: choose_logs_file, build_atom_index_maps, compute_residue_curvature, label_dna_facing, plot_violin, main | dna, rna | — | *.csv; alanine_dna_facing_table.csv; alanine_dna_facing_violin.pdf; alanine_dna_facing_violin.png | Curvature |
| `Part_II_Molecular_Analysis/Curvature/5c3_violin.py` | Main | Declared: residue_to_token, token_class, print_token_outliers_iqr, residue_to_one_letter, ViolinConfig, choose_logs_file, sanitize_filename, build_residue_instances_from_atoms … | dna, protein, rna | — | *.csv; .pdf; .png; .svg; _residue_instances.csv; outliers_token_A.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/__init__.py` | module | Empty file | — | — | — | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/__init__.py` | module | Empty file | — | — | — | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/curvature_assignments.py` | module | Top-level script; see principal family inventory | — | — | — | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/curvature_name_bins.py` | Main | Top-level script; see principal family inventory | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow | — | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/CurvatureByInterface.py` | Main | Declared: classify_surf | 1bna, cambrin, edta, hairpin, ligand, ncp, p53tet, protein, streptavidin | — | _vor_logs.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/CurvatureDistributions.py` | Main | Top-level script; see principal family inventory | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow | — | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_Peaks.py` | Main | Declared: _parse_ball_indices, _load_pdb_index_map, main | bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | aw | .pdb; _curvature_bins.csv; aw_logs.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4B_Peak_Assignment.py` | Main | Declared: _peak_signature, _cosine_similarity, merge_consecutive_similar_peaks, merge_consecutive_almost_identical_peaks, compute_non_overlapping_peak_bounds, find_peak_bounds_by_inflection, find_peaks, write_peaks_txt … | bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | — | *.csv; .csv; .txt; Select *_curvature_bins.csv; _curvature_bins.csv; _curvature_peaks.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4C_Combine_Curvature_Peak_Files.py` | Main | combine_curvature_peaks_folder.py | hairpin | — | Output CSV path. Default is <selected_root>/combined_peaks.csv; _curvature_bins.csv; combined_peaks.csv | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/Extra/surface_energy.py` | module | Declared: surface_tension, molar_free_surface_energy, stokes_theorem, basic_surface_energy | — | — | — | Curvature |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/SurfaceEnergyCurvature.py` | Main | Declared: parse_balls, normalize_pair, infer_element_from_atom_name, get_element_column, assign_energy_rank, add_surface_pair_columns, add_atom_pair_types, classify_pair_energy … | protein | aw | Choose model folder containing aw/aw_logs.csv; all_surfaces_curvature_by_energy_class.png; all_surfaces_curvature_by_pairtype.png; all_surfaces_curvature_energy_classes.csv; aw_logs.csv; bonded_surfaces_curvature_by_energy_class.png; bonded_surfaces_curvature_energy_classes.csv; curvature_by_energy_class_summary.csv; energy_rank_curvature_summary.csv; energy_rank_vs_curvature.png; literature_vs_calculated_surface_energy.png; literature_vs_calculated_surface_energy.svg; literature_vs_calculated_surface_energy_data.csv; nonbonded_surfaces_curvature_by_energy_class.png; nonbonded_surfaces_curvature_energy_classes.csv; pairtype_curvature_bar.png; pairtype_curvature_bar.svg; pairtype_curvature_summary.csv | Curvature |
| `Part_II_Molecular_Analysis/extra/__init__.py` | module | Empty file | — | — | — | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/AtomAssignmentClustering.py` | module | Top-level script; see principal family inventory | — | aw, pow | — | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/AtomAssignmentFileMaker.py` | Main | Declared: clustering_measures, estimate_volume, average_pairwise_distance, get_group_info, get_atoms_info, get_rad_vals, write_files, bool_assign … | bsa, protein | aw, pow, power, vdw | /atomic_comparisons.csv; /atomic_comparisons_simple.csv; /pow_model_training.csv | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/AvgAtomDiff.py` | Main | Declared: atoms_per_diff | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow, power, primitive | — | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/ElementBinningDiffs.py` | Main | Declared: atoms_per_diff | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow, power, primitive | — | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/ElementBinningTots.py` | Main | Declared: atom_tots | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow, power, primitive | — | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/KMeansClustering.py` | module | Top-level script; see principal family inventory | — | pow, power | /atomic_comparisons_with_clusters.csv; atomic_comparisons_with_clusters.csv | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/extra/TrainAWVolModel.py` | module | Declared: mean_absolute_percentage_error, load_data_and_train, predict_from_loaded_data, print_results, plot_results | — | aw, pow, vdw | *.csv | F6 exploratory / SI candidate |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/__init__.py` | module | Empty file | — | — | — | F1 draft/support |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1A_sphere_schemes.py` | module | Empty file | — | — | — | F1 draft/support |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/__init__.py` | module | Empty file | — | — | — | Curvature (F1 link) |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/curvature_range_cartoons.py` | module | Declared: midpoint_between_surfaces_point, get_vor_points, plot_voronoi_partition | — | voronoi | .png; .svg | Curvature (F1 link) |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/curvature_range_plot.py` | Main | Declared: curvature_distance_ratio_sweep, plot_curvature_distance_ratio_lines | — | — | — | Curvature (F1 link) |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/Curvature_Shapes_mpl.py` | Main | Declared: offset_point, normalize_vectors, set_axes_equal, style_axis, create_figure, offset_curve, plot_double_sided_surface, plot_principal_line … | — | — | — | Curvature (F1 link) |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/Show_Ratio_Dist.py` | Main | Declared: plot_circle, plot_circle_rads, plot_line_between, main | — | — | — | Curvature (F1 link) |
| `Part_II_Molecular_Analysis/F2_frame_comparison/__init__.py` | module | Empty file | — | — | — | Needs Review: descriptors |
| `Part_II_Molecular_Analysis/F2_frame_comparison/AtomicCorrelation_multi_plots.py` | Main | VorPy atom-level correlation analysis. | — | — | *.csv; .png; .svg; atomic_correlation_combined_data.csv; atomic_correlation_matrix.csv; atomic_correlation_matrix.png; atomic_correlation_matrix.svg; atomic_selected_pair_statistics.csv; aw_logs.csv | Needs Review: descriptors |
| `Part_II_Molecular_Analysis/F2_frame_comparison/AtomicCorrelation_updated.py` | Main | VorPy atom-level correlation analysis. | — | — | *.csv; .png; .svg; atomic_correlation_combined_data.csv; atomic_correlation_matrix.csv; atomic_correlation_matrix.png; atomic_correlation_matrix.svg; atomic_selected_pair_statistics.csv; aw_logs.csv | Needs Review: descriptors |
| `Part_II_Molecular_Analysis/F2_frame_comparison/AtomicCorrelation.py` | Main | VorPy atom-level correlation analysis. | — | — | *.csv; .png; .svg; atomic_correlation_combined_data.csv; atomic_correlation_matrix.csv; atomic_correlation_matrix.png; atomic_correlation_matrix.svg; atomic_selected_pair_statistics.csv | Needs Review: descriptors |
| `Part_II_Molecular_Analysis/F2_frame_comparison/CurvatureComparison_updated.py` | Main | VorPy Frame Variation Analysis ============================== | — | aw | f#/aw/aw/aw_logs.csv; .png; .svg; aw_logs.csv; frame_variation_data.csv; frame_variation_summary.csv | Curvature (F2 link) |
| `Part_II_Molecular_Analysis/F2_frame_comparison/CurvatureComparison.py` | Main | VorPy Frame Variation Analysis ============================== | — | aw | f#/aw/aw/aw_logs.csv; .png; .svg; aw_logs.csv; frame_variation_data.csv; frame_variation_summary.csv | Curvature (F2 link) |
| `Part_II_Molecular_Analysis/F2_frame_comparison/FrameAnalysis.py` | Main | Declared: compare_by_atom, get_frames_by_group, plot_frame_diffs, plot_frame_deviations | — | aw, pow, prm | aw/aw_logs.csv; pow/pow_logs.csv; prm/prm_logs.csv | F2 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5A_terminal_effects.py` | Main | Terminal profiles of validated repeat units (no convergence-model fitting). | — | — | summaries/F5A_frame_profiles.csv; summaries/F5A_terminal_interior_frame_contrasts.csv; summaries/F5A_terminal_profiles.csv | F3/F4 historical; review |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5B_solvation_environments.py` | Main | AW repeat exposure defines provisional environments for every scheme. | peg | aw | summaries/F5B_chemical_environment_frame_means.csv; summaries/F5B_environment_frame_means.csv; summaries/F5B_zero_water_units.csv | F3/F4 historical; review |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5C_polymer_length.py` | Main | Length dependence; error bars are temporal SD, not independent-replicate CIs. | — | — | summaries/F5C_length_frame_means.csv; summaries/F5C_length_summary.csv | F3 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5D_internal_water.py` | Main | Available boundary diagnostics only; internal occupancy is explicitly unavailable. | — | — | F5D/boundary_flags.csv; summaries/F5D_internal_occupancy_unavailable.csv; summaries/F5D_topology_diagnostics.csv | F3/F4 historical; review |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5E_scheme_differences.py` | Main | Paired scheme deviations against AW, in native units or percentage points. | — | aw, pow, prm | summaries/F5E_environment_deviations.csv; summaries/F5E_exposure_binned_deviations.csv; summaries/F5E_model_deviation_summary.csv | F3/F4 historical; review |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5F_water_coordination.py` | Main | Binned hydration geometry by chemical class; no hydrogen-bond inference. | peg | — | F5F/water_area_coordination_frame_bins.csv; summaries/F5F_coordination_frame_bins.csv | F3/F4 historical; review |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/SolventExposedDiff.py` | module | Empty file | — | — | — | F3/F4 historical; review |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_common.py` | module | F5 scientific tables. Source logs and solver state are never modified. | heparin, peg | aw, pow, prm | — | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_config.py` | module | Configuration shared by the F5 loader, runner, and figures. | dna, heparin, peg | aw, pow, power, primitive, prm, voronoi | C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\OrganizedModels | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_peg_dataset.py` | module | Screened PEG paths and identity-checked first-layer membership; no solver calls. | peg, peg3000 | aw, pow, prm | .pdb; _logs.csv; layer_01.pdb; output/peg_one_layer_investigation_20261112/one_layer_run_status.csv | F3 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_run_all.py` | Main | Run F5 analysis directly or as a module (see README). | dp12, dp18, dp2, dp24, dp30, heparin, peg, peg1000, peg2000, peg300, peg3000, peg4000, peg6000 | aw, pow, prm | .csv; .pkl; _logs.csv; frame.pdb; summaries/interface_shards.csv; summaries/interface_totals.csv; validation/failures.csv; validation/figure_failures.csv; validation/file_inventory.csv; validation/independent_area_checks.csv; validation/pdb_residue_inventory.csv; validation/run_config.json; validation/run_report.json; validation/scheme_matching.csv | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5A_exposed_vs_buried.py` | Main | Compare matched PEG Power and AW cells using the cached F5 atom table. | peg | aw, pow, power, voronoi | F5A_atom_types.csv; F5A_summary.csv; F5A_summary_signed.csv; F5A_volume_vs_exposure_atoms.csv; F5A_volume_vs_exposure_atoms_signed.csv; F5A_volume_vs_water_contacts_atoms_signed.csv; output/summaries/atoms.csv; summaries/atoms.csv | F3 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5A_PEG_chemistry_vs_hydration.py` | Main | Plot matched PEG atom volume and surface-area changes against AW contacts. | dna, peg, peg1000, peg2000, peg300, peg3000, peg4000, peg6000 | aw, pow, power, voronoi | C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\P2\NewModels\PEG\Reruns\layers_4; F5A_PEG_chemistry_vs_hydration.csv; F5A_coverage.csv | F3 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5B_heparin_chemistry.py` | Main | Figure 5B: Heparin chemistry, AW water coordination, and chain position. | dna, dp12, dp18, dp2, dp24, dp30, heparin | aw, pow, power, voronoi | C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\OrganizedModels\Heparin; F5B_heparin_chemistry_summary.csv; _logs.csv; frame.pdb | F4 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5B_water_burial.py` | Main | Compare AW and Power water-molecule geometry from water-inclusive F5 logs. | heparin, peg | aw, pow, power | F5B_coverage.csv; F5B_summary.csv; F5B_surrounded_counts.csv | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5C_rerun_missing.py` | Main | Preview six missing/empty AW solvent exports; pass --run to rerun sequentially. | dna, peg, peg2000, peg300, peg4000 | aw, voronoi | C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\P2\NewModels\PEG\Checkpoints\aw; _checkpoint.vpy | F3 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5C_water_contacts.py` | Main | Figure 5C from AW solvent logs; one observation per complete water/frame. | dna, peg, peg1000, peg2000, peg300, peg3000, peg4000, peg6000 | aw, voronoi | C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack\Ericson Voronoi DNA\P2\NewModels\PEG\Reruns\layers_4\aw; F5C_coverage.csv; F5C_water_contacts.csv; F5C_water_contacts_bins.csv | F3 |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/__init__.py` | module | Empty file | — | — | — | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_classification.py` | module | Declared: test_provisional_exposure_boundaries_and_invalid_values, test_environment_is_aw_repeat_reference_for_other_schemes | — | pow | — | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_direct_execution.py` | module | PyCharm-style file execution must not depend on the working directory. | — | — | — | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_loading.py` | module | Declared: case, test_classification_requires_residue_and_water_name, test_missing_topology_remains_unavailable, test_end_distances, test_peg_repeat_order_requires_named_connected_backbone, test_reader_preserves_missing_surface_area_and_rejects_bad_numeric, tes | dp24, heparin, peg, peg3000 | aw, pow, prm | F5A/F5A_summary.csv; F5A/F5A_volume.png; F5B/F5B_coverage.csv; F5B/F5B_summary.csv; _logs.csv; aw_logs.csv; frame.pdb; log.csv; pow_logs.csv; validation/pdb_residue_inventory.csv | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_scheme_comparison.py` | module | Declared: paired, test_paired_differences_percentage_points_and_zero_counts, test_scheme_identity_mismatch_cannot_generate_differences, test_unmatched_atoms_reported_and_empty_neighbor_sets_have_jaccard_one, test_repeat_matches_verify_member_coordinates, test_ | — | pow | — | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_surface_areas.py` | module | Declared: test_index_matching_area_partition_and_unique_water_count, test_residue_fraction_uses_area_weights_and_deduplicates_molecules, test_zero_water_contacts_are_valid, test_duplicate_records_count_once_but_conflicts_fail, test_invalid_cells_retained_and_f | — | — | — | F3/F4 shared/support |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/__init__.py` | module | Empty file | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/__init__.py` | module | Empty file | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/AllinaRow.py` | Main | Top-level script; see principal family inventory | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/compare_frames.py` | module | Top-level script; see principal family inventory | — | — | — | F2 experimental |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/Distribution_of_element_radii.py` | module | Top-level script; see principal family inventory | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/Distribution_Types.py` | module | Top-level script; see principal family inventory | 1bna, cambrin, hairpin, ncp, p53tet, protein | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/frame_analysis.py` | module | Top-level script; see principal family inventory | cambrin | aw, pow, power, primitive, prm | .pdb; aw.csv; pow.csv; prm.csv | F2 experimental |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/num_atoms_log.py` | Main | Top-level script; see principal family inventory | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/Plot_by_system.py` | Main | Declared: get_vals, plot_distribution, plot_vals | — | pow | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/ResidueBinningDiffs.py` | Main | Declared: residue_per_diff | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow, power, primitive | /res_data/{}_del_res.csv; /res_data/{}_pow_res.csv; /res_data/{}_vor_res.csv; _{}_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/ResidueBinningTots.py` | Main | Declared: residue_per_diff | 1bna, bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow, power, primitive | /res_data/{}_del_res.csv; /res_data/{}_pow_res.csv; /res_data/{}_vor_res.csv; _{}_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/vis_table_all_models.py` | Main | Declared: get_system_data | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/__init__.py` | module | Empty file | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/config.py` | module | Top-level script; see principal family inventory | dna, edta, protein, rna | aw, pow, voronoi | E:\Molecular; E:\OneDrive - Georgia State University\GSU NSC\Manuscripts\Ericson Voronoi DNA\P2\fig2_atomic_level_scheme_deviations\2C_Full_Plots_Clustered | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/pipeline.py` | module | Declared: run_pca_analysis, summarize_sol_facing, resolve_group_folders, canonicalize_atom_name, infer_molecule_class, _read_scheme_logs, select_folders_multi, collect_atom_volume_points … | dna, protein, rna | aw, pow, prm | aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/run_batch.py` | Main | Declared: run_batch_cluster_plots | dna, protein, rna | — | _failed_runs.txt | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/run_single.py` | Main | Top-level script; see principal family inventory | dna | aw, pow | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_atomic_umap.py` | Main | Joint AW/Power atom-descriptor UMAP and exploratory HDBSCAN clustering. | — | aw, power | F3_atomic_umap_assignments.csv; F3_atomic_umap_cluster_summary.csv; F3_atomic_umap_clusters.png; F3_atomic_umap_metrics.png; F3_atomic_umap_report.json; F3_atomic_umap_standardized.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_build_figure.py` | module | Assemble final Figure 3 from already-rendered panel assets. | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_common.py` | module | Shared data utilities for Figure 3. | — | aw, pow, power, primitive, prm, voronoi | aw\aw_logs.csv; aw_logs.csv; pow\pow_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_curvature_pca.py` | Main | Figure 3 atomic geometry/curvature correlations and PCA. | — | aw, pow, power, vdw | .xlsx; analysis_report.json; atomic_descriptors_merged.csv; atomic_descriptors_standardized.csv; correlation_pearson.csv; correlation_pearson.png; correlation_spearman.csv; correlation_spearman.png; pca_explained_variance.csv; pca_loadings.csv; pca_loadings.png; pca_scores.csv; pca_scores_pc1_pc2.png; pca_scree.png | Curvature |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/__init__.py` | module | Empty file | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3A_atom_diffs.py` | Main | Declared: _stats, get_atom_sa_neighbors, get_atom_data, plot_data | — | aw, pow, prm | aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3A_residue_diffs.py` | Main | Declared: get_res_sa, get_res_data, plot_data | — | aw, pow, power, primitive, prm | aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3B_molecule_totals.py` | Main | Declared: plot_data | — | aw, pow, prm | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3B_residue_diffs_updated.py` | Main | Declared: _res_key, _stats, get_res_sa_neighbors, get_res_data, plot_data | — | aw, pow, prm | aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3C_molecule_diffs_updated.py` | Main | Declared: _res_key, get_sol_neighbor_count, _pct, get_molecule_data, plot_data | — | aw, pow, prm | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3C_residue_sphericity_vs_deviation.py` | Main | Declared: get_res_sa, build_residue_sphericity_vs_deviation_points, plot_residue_sphericity_vs_deviation | — | aw, pow, prm | aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3C_Residue_type_diffs.py` | Main | Declared: normalize_residue_label, pick_folder, safe_sem, get_matching_atom_by_index, check_index_uniqueness, build_atom_lookup, get_matching_atom, get_residue_key … | — | aw, pow, power, primitive, prm | aw_logs.csv; figure_3c_residue_type_deviation_distribution.png; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3D_molecule_sphericity_vs_deviation.py` | module | Declared: plot_data | — | aw, pow, prm | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2_AnalysisTable.py` | Main | Declared: infer_element, normalize_na_residue, is_nucleic_acid_res, classify_atom, infer_element_from_atom_name, neighbor_count, load_atoms_from_logs, summarize_by_atom_type … | dna, protein, rna | aw, power, vdw | *.csv; /atom_type_summary.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2A_AtomicVolumeDifferences.py` | Main | Declared: add_outside_element_list, apply_2c_axis_style, plot_vols | dna, protein, rna | aw, pow, prm, voronoi | E:/Molecular; E:/OneDrive - Georgia State University/GSU NSC/Manuscripts/Ericson Voronoi DNA/P2/fig2_atomic_level_scheme_deviations/2A_Full_Plots; _element_volume_plot.png; _element_volume_plot.svg; aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2Ac_atomic_deviations.py` | Main | Declared: MoleculeInputs, _safe_float_series, _extract_element, _key_cols, _normalize_atoms_df, _merge_pow_aw, _abs_percent_diff, _choose_curvature_col … | p53tet | aw, pow | _exemplar_candidates.csv; _pow_vs_aw_atoms_merged.csv; _summary.json; all_summaries.json; hammerhead_aw_logs.csv; hammerhead_pow_logs.csv; p53tet_aw_logs.csv; p53tet_pow_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_All_Data_Sets_Clustered.py` | Main | Declared: sort_plot_groups, get_chem_class, is_bonded, build_bonds, classify_carbon, classify_nitrogen, classify_oxygen, classify_hydrogen … | dna, edta, protein, rna | aw, pow, prm, voronoi | .csv; .pdb; .png; .svg; E:/Molecular; E:/OneDrive - Georgia State University/GSU NSC/Manuscripts/Ericson Voronoi DNA/P2/fig2_atomic_level_scheme_deviations/2C_Full_Plots_Clustered; E:\Molecular; E:\OneDrive - Georgia State University\GSU NSC\Manuscripts\Ericson Voronoi DNA\P2\fig2_atomic_level_scheme_deviations\2C_Full_Plots_Clustered; _failed_runs.txt; _group_key.csv; _plot_points.csv; _stats.csv; aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper_Auto_Cluster.py` | module | Declared: ensure_required_columns, make_sol_facing_binary, build_feature_matrix, build_kprototypes_matrix, run_kmeans, run_agglomerative, run_dbscan, run_hdbscan … | — | aw, pow | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper_Protein_Assignment.py` | Main | Declared: pick_folders, get_output_directory, safe_system_name, load_logs_from_folder, parse_group_string, build_atom_group_mapping, standardize_columns, ensure_required_columns … | protein | aw, pow, prm | .png; aw_logs.csv; pow_logs.csv; prm_logs.csv; protein_atom_grouping_progress.csv; protein_atom_grouping_progress.json | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper.py` | Main | Declared: choose_folder, find_first_existing, require_path, load_patch_atoms_csv, parse_pdb_atoms, normalize_atom_df, _norm_str, _coord_key … | — | aw, pow | .pdb; aw_logs.csv; input.pdb; model.pdb; patch_atom_aw_neighbor_pairs.csv; patch_atom_aw_pow_pair_compare.csv; patch_atom_feature_correlations.csv; patch_atom_outlier_summary_by_atom_name.csv; patch_atom_outlier_summary_by_neighbor_type.csv; patch_atom_pair_feature_correlations.csv; patch_atom_pow_neighbor_pairs.csv; patch_atom_volume_neighbor_summary.csv; patch_atoms.csv; pdb_mapping_diagnostics.csv; pdb_mapping_summary.csv; pow_logs.csv; structure.pdb | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper2.py` | Main | Declared: choose_csv, classify_neighbor_type, shannon_entropy, safe_std, safe_mean, safe_sum, unique_count, build_patch_level_features … | — | aw, pow | *.csv; Select patch_atom_aw_pow_pair_compare.csv; patch_level_feature_correlations.csv; patch_level_feature_heatmap.png; patch_level_feature_scatter_grid.png; patch_level_feature_summary.txt; patch_level_features.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Ligand_Interface_Residues.py` | Main | Declared: get_ligand_residue_interfaces, print_interfaces, plot_interfaces, main | ligand | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Mean_Vol_SA.py` | Main | Declared: get_res_sa, get_res_vol, main | — | — | *.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/NCP_Histone_DNA_Interfaces.py` | Main | Declared: get_dna_contact_profile, print_profile, plot_profile, main | dna, protein | voronoi | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/NCP_Stacked_DNA_Contact_Profile.py` | Main | Declared: get_histone_contact_profile, plot_stacked_profile, print_histone_totals, main | dna, ncp | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/NCP_System_Analysis.py` | Main | Declared: get_group_stats, get_surface_areas, main | dna, protein | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/p53tet_Chain_Equivalency_deviation.py` | Main | Declared: get_chain_data, print_data, plot_data, main | p53tet | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Res_Vol_SA.py` | Main | Declared: get_res_vol, get_res_sa, plot_residue_data, main | — | aw, pow, power, primitive, prm | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Tester_for_power_residues.py` | Main | Declared: pick_folder, normalize_residue_label, get_first_present, normalize_chain_value, normalize_atom_name, normalize_index, normalize_resseq, get_atom_identity … | — | aw, pow, prm | _pow_atom_comparisons.csv; _pow_residue_summary.csv; _pow_system_summary.csv; _prm_atom_comparisons.csv; _prm_residue_summary.csv; _prm_system_summary.csv; _unmatched_aw_to_pow.csv; _unmatched_aw_to_prm.csv; aw_logs.csv; pow_logs.csv; prm_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3A_multiscale_summary.py` | Main | Figure 3A — multiscale AW vs Power deviation summary. | — | aw, pow, power | F3A_multiscale_summary.csv; F3A_multiscale_summary.pdf; F3A_multiscale_summary.png; F3A_multiscale_summary.svg | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3B_atomic_groups.py` | Main | Figure 3B — Atomic AW vs Power comparisons. | dna, protein, rna | aw, pow, power | F3B_atomic_group_outliers.csv; F3B_atomic_groups.png; F3B_atomic_groups.svg; F3B_atomic_instance_cache.csv; F3B_atomic_instance_outliers.csv; F3B_atomic_top_candidates.csv; F3B_nucleic_group_summary.csv; F3B_protein_group_summary.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3C_residue_groups.py` | Main | Figure 3C — Residue-type AW vs Power deviations. | dna, protein, rna | aw, power, voronoi | F3C_residue_instance_cache.csv; F3C_residue_instance_outliers.csv; F3C_residue_top_candidates.csv; F3C_residue_type_outliers.csv; F3C_residue_type_summary.csv; F3C_residue_type_violins_legend_counts.png; F3C_residue_type_violins_legend_counts.svg | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3D_geometry_examples.py` | module | Figure 3D — representative atomic and residue geometry examples. | — | aw, power | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3E_solvent_exposure.py` | Main | Figure 3E — Direct solvent-exposure dependence of AW/Power deviations. | dna, protein, rna | aw, power | .pdb; F3E_atom_exposure_cache.csv; F3E_atom_type_controlled_slopes.csv; F3E_exposure_group_cutoffs.csv; F3E_exposure_multiscale_ranges.png; F3E_exposure_multiscale_ranges.svg; F3E_exposure_multiscale_summary.csv; F3E_residue_exposure_cache.csv; F3E_residue_type_controlled_slopes.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/__init__.py` | module | Empty file | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/AvgSurfDiff.py` | Main | Declared: per_diffs_surfs, surfs_per_diff | bsa, cambrin, edta, hairpin, ncp, p53tet, streptavidin | aw, pow, power, primitive, prm, voronoi | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/InsideOutside.py` | Main | Declared: inside, inside_out | cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow, power, primitive | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/OtherSADiffs.py` | Main | Declared: classify_surf | 1bna, cambrin, edta, hairpin, ligand, ncp, p53tet, protein, streptavidin | power, primitive | _del_logs.csv; _pow_logs.csv; _vor_logs.csv | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/SolventAccessible.py` | Main | Declared: draw_circle, get_sol_acc_points | — | — | — | F6 (legacy/support as applicable) |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/color_residues.py` | Main | Declared: color_pdb_by_res | — | — | _colorized.pdb; _residue_data.csv | Needs Review: F6/F7 |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/diff_by_res_vor_pow.py` | Main | Declared: add_vorpy_root_to_path, import_read_logs2, normalize_col_name, find_column, get_schema_columns, get_atoms_df_from_logs, clean_value, make_residue_key … | — | aw, pow, voronoi | _logs.csv; residue_instance_percent_difference_diagnostic_read_logs2.csv; residue_instance_raw_volume_diagnostic_read_logs2.csv; residue_percent_difference_summary_read_logs2.csv | Needs Review: F6/F7 |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/coarse_grained_residue_volume_difference.py` | Main | Declared: normalize_residue_label, compute_aa_rank_groups, choose_root_folder, iter_logs_files, build_residue_volume_df, _quantiles, plot_figure_7a_grid, plot_figure_7a … | rna | aw, pow, prm | _logs.csv | F7 tentative |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/__init__.py` | module | Empty file | — | — | — | F7 tentative |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/f8_build_surfaces.py` | Main | Declared: classify_residue, iter_log_files, build_atom_tables, _empty_surface_df, _classify_ball, _canonicalize_surface_pair, _pair_class, build_surface_table … | dna, ncp, protein | aw, pow, prm | _logs.csv; surf_df.csv | F6/F7 shared |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/f8a_plot_visualize_example.py` | module | Empty file | — | — | — | Needs Review: draft |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/f8b_raincloud_interfaces.py` | Main | Declared: add_regime_shading, compute_dual_modes, plot_dual_mode_markers, adjust_color, _trim_to_y_range, _square_grid_shape, _ensure_surf_df, _cg_order_key … | ncp | aw | .png; F8B_raincloud_interfaces_AW.png; surf_df.csv | Curvature |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/f8d_interface_heatmaps.py` | module | Empty file | — | — | — | Needs Review: draft |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/__init__.py` | module | Empty file | — | — | — | Curvature |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/color_code_residues.py` | Main | Declared: color_pdb_by_res | 1bna | — | 1BNA.pdb; 1BNA_coarse_ad.pdb; 1BNA_coarse_ncap.pdb; 1BNA_coarse_scbb_ad.pdb; 1BNA_coarse_scbb_ncap.pdb; 1BNA_residue_data.csv; C:/Users/jacke/Documents/data1/; _colorized.pdb | F7 tentative |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/curvature_assignments.py` | Main | Top-level script; see principal family inventory | — | — | — | Curvature |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/curvature_by_radius.py` | Main | Top-level script; see principal family inventory | — | aw | .csv; .pdb; _a.pdb; a.pdb | Curvature |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/curvature_distribution.py` | Main | Top-level script; see principal family inventory | cambrin, edta, hairpin, ncp, p53tet, streptavidin | pow | — | Curvature |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/max_avg_curv.py` | Main | Top-level script; see principal family inventory | — | — | .csv; .pdb; _a.pdb; a.pdb | Curvature |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme.py` | Main | Top-level script; see principal family inventory | — | aw, pow, power | .csv; .pdb; _a.pdb; a.pdb | F7 tentative |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme2.py` | Main | Declared: percent_diff | — | power | .csv; .pdb; /residue_data.csv; _a.pdb; a.pdb | F7 tentative |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme3.py` | Main | Declared: percent_diff | — | aw, pow, power | — | F7 tentative |
| `Part_II_Molecular_Analysis/Network_Dynamics/__init__.py` | module | Empty file | — | — | — | F2 draft/support |
| `Part_II_Molecular_Analysis/Network_Dynamics/residue_flexibilty_using_frame_dynamics.py` | module | We are going to input frame logs files for a molecule and study the dynamics of the chemical relevance of the network graph through comparisons of 1. Contact persistence, Contact area variability, Residue network volatility, and Weighted strength fluctuations | — | — | — | F2 draft/support |
| `Part_II_Molecular_Analysis/SuplimentalInformation/__init__.py` | module | Package initializer / exports | — | — | — | Curvature / SI candidate |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_2kai_curvature_forensics.py` | Main | Read-only forensic exports for the frozen 2KAI Cazals Power interface. | 2kai, protein | power | 2KAI_CURVATURE_FORENSIC_REPORT.txt; 2kai_M5_facet_forensics.csv; 2kai_M5_summary.txt; 2kai_curvature_by_pattern.csv; 2kai_curvature_by_singleton_atom.csv; 2kai_curvature_by_singleton_chain.csv; 2kai_curvature_by_singleton_residue.csv; 2kai_curvature_edges.csv; 2kai_lys15_radial_curvature.csv; 2kai_sign_counterfactuals.csv; vorpy/data/2KAI.pdb | Curvature / SI candidate |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_2kai_independent_audit.py` | Main | Independent stage-by-stage 2KAI Cazals reproduction audit. | 2kai, protein | power | 2KAI.pdb; audit_report.txt; audit_result.json; reproduction_metadata.json; stage_01_atoms_and_radii.csv; stage_01_census.csv; stage_02_gudhi_regular_and_alpha_simplices.csv; stage_03_independent_power_vertices.csv; stage_04_bicolor_m5_facets.csv; stage_04_facet_orthogonal_ball_candidates.csv; stage_05_independent_interface_edges.csv; stage_07_final_sH.json; stage_counts_and_comparisons.csv | Curvature / SI candidate |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_alpha0_audit.py` | Main | cazals_alpha0_audit.py | 2kai, protein | power, voronoi | .npz; .pdb; .pdb or .npz; Input must be .pdb or .npz; bicolor_facets.csv; interface_edge_candidates.csv; simplex_filtration.csv | Curvature / SI candidate |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_pdb.py` | module | cazals_pdb.py | 2kai, protein | power | — | Curvature shared (retain) |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_power_interface_diagnostic.py` | Main | Diagnostic comparison of unrestricted and Cazals-style 2KAI interface curvature. | 2kai | power, voronoi | .npz; .pdb; .pdb or .npz; Input must be .pdb or .npz | Curvature / SI candidate |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_power_interface.py` | Main | cazals_power_interface.py | 2kai | power | .npz; .pdb; Input must be .pdb or .npz | Curvature shared (retain) |
| `Part_II_Molecular_Analysis/SuplimentalInformation/comparative_figures/__init__.py` | module | Empty file | — | — | — | F1 support |
| `Part_II_Molecular_Analysis/SuplimentalInformation/comparative_figures/F5_alpha_selection.py` | Main | Declared: _extract_progression, _plot_absolute, _plot_retention, _print_summary, main | — | power | — | F1 support |
| `Part_II_Molecular_Analysis/SuplimentalInformation/comparative_figures/plot_common.py` | module | Declared: select_csv, select_output_directory, load_static_comparison, static_protein_only, ordered_systems, require_methods, method_label, configure_matplotlib … | 2kai, protein | aw, power | *.csv; .pdf; .png; .svg; Select static_method_comparison.csv | F1 support |
| `Part_II_Molecular_Analysis/SuplimentalInformation/Curvature/CubeVertexGaussianCurvature.py` | module | CubeVertexGaussianCurvature.py | — | — | .png; .svg | Curvature |
| `Part_II_Molecular_Analysis/SuplimentalInformation/Curvature/DulaunayDual.py` | Main | Two-group weighted Delaunay, power diagram, and alpha selection. | protein | power | 01_circles.png; 02_weighted_delaunay.png; 03_power_diagram.png; 04_alpha_selected_power.png; 05_cazals_selected_power.png | F1 |
| `Part_II_Molecular_Analysis/SuplimentalInformation/Curvature/visualize_interface_criteria.py` | Main | Four simple, geometrically consistent 2D diagrams of interface criteria. | protein | power, voronoi | 01_condition_a_accepted.png; 02_condition_a_rejected.png; 03_condition_b_accepted.png; 04_condition_b_rejected.png | F1 |
| `Part_III_ML_Processes/Surfaces_Training.py` | module | Declared: gather_file_names, get_atoms_info, get_surfs_info | — | aw, pow | aw_logs.csv; balls.txt; pow_logs.csv | SI candidate |
| `power_interface_cli.py` | module | User-facing Cazals-compatible Power/Laguerre interface runner. | protein | aw, laguerre, power | .pdb; power_interface_atoms.csv; power_interface_components.csv; power_interface_edges.csv; power_interface_facets.csv; power_interface_summary.txt | Shared/support; review |
| `power_interface_curvature.py` | module | Cazals-style bicolor power-interface curvature on a weighted regular complex. | — | aw, power | power_interface_edges.csv; power_interface_summary.txt; power_interface_surfaces.csv; power_simplex_filtration.csv | Curvature (shared APIs retained) |
| `power_interface_selection.py` | module | Selection views of solved Power geometry; no curvature/geometry formulas. | — | power | — | F1 shared (retain) |
| `power_vs_aw.py` | module | Transparent side-by-side reporting for Power and AW interface results. | — | aw, power | power_vs_aw_metrics.csv; power_vs_aw_summary.txt | Curvature (shared APIs retained) |
| `tools/__init__.py` | module | Package initializer / exports | — | — | — | Shared / SI candidate |
| `tools/batch/__init__.py` | module | Package initializer / exports | — | — | — | Shared / SI candidate |
| `tools/batch/compile_logs.py` | Main | Declared: get_logs_and_pdbs, make_new_logs | — | aw, pow, power | /Data/Analyze/tools/batch/loggy_woggys.txt; logs_pdbs.txt; new_log_info.csv | Shared / SI candidate |
| `tools/batch/compile_new_logs.py` | Main | Declared: get_logs_and_pdbs, make_new_logs | — | aw, pow | /Data/Analyze/tools/batch/loggy_woggys1.txt; loggy_woggys1.txt; logs_pdbs.txt; new_log_info.csv | Shared / SI candidate |
| `tools/batch/generate_batch.py` | module | Top-level script; see principal family inventory | cambrin | pow | — | Shared / SI candidate |
| `tools/batch/get_files.py` | Main | Declared: get_files, get_all_files | — | aw, pow, prm | .csv; .pdb; _logs.csv | Shared / SI candidate |
| `tools/batch/InDepthLogAnalysis.py` | Main | Declared: get_logs_and_pdbs, make_new_logs | — | aw, pow | /Data/Analyze/tools/batch/loggy_woggys1.txt; loggy_woggys1.txt; logs_pdbs.txt; new_log_info.csv | Shared / SI candidate |
| `tools/batch/IsoMetricQuotient.py` | Main | Declared: get_logs_and_pdbs, make_new_logs | — | aw, pow | /Data/Analyze/tools/batch/loggy_woggys1.txt; Isoq.csv; loggy_woggys1.txt; logs_pdbs.txt | Shared / SI candidate |
| `tools/batch/make_batch_fixed_missing_atoms.py` | Main | Declared: print_run_statements | — | aw, pow, prm | .pdb | Shared / SI candidate |
| `tools/batch/make_foam_runs.py` | module | Top-level script; see principal family inventory | — | aw | .csv; /aw/aw_logs.csv; /verts.txt; balls.pdb | Shared / SI candidate |
| `tools/batch/wipe_user_data.py` | Main | Top-level script; see principal family inventory | — | aw, pow | info.txt | Shared / SI candidate |
| `tools/calcs.py` | module | Declared: bin_data, sort_lists | — | — | — | Shared / SI candidate |
| `tools/CleanData/__init__.py` | module | Package initializer / exports | — | — | — | Shared / SI candidate |
| `tools/CleanData/AdjustAfterRuns.py` | Main | Declared: move_contents, integrate_folders, adjust_after_runs | — | aw, pow | /aw/aw_logs.csv; /aw/aw_verts.txt; /aw_verts.txt; /balls_Network_aw/balls_Network_aw_logs.csv; /balls_Network_pow/balls_Network_pow_logs.csv; /pow/pow_logs.csv; /pow/pow_verts.txt; /pow_verts.txt | Shared / SI candidate |
| `tools/CleanData/check_mol_data.py` | Main | Declared: _extract_system_name, _find_pdb_file, _load_complete_indices, _missing_in_order, check_mol_data | — | aw, pow, prm | .pdb; aw_logs.csv; pow_logs.csv; prm_logs.csv | Shared / SI candidate |
| `tools/CleanData/check_neighbors.py` | Main | Declared: Config, choose_logs_file, sanitize_filename, build_surface_neighbor_map, build_neighbors_map_from_atoms, infer_neighbors_column, compare_and_print_neighbors, main | — | — | *.csv; _neighbors_vs_surfaces.csv | Shared / SI candidate |
| `tools/CleanData/CheckFoamFolders.py` | Main | Declared: completeness_check | — | aw, pow | /aw/aw_logs.csv; /aw/aw_logs.pdb; /pow/pow_logs.csv; /pow/pow_logs.pdb | Shared / SI candidate |
| `tools/CleanData/CleanFoamData.py` | Main | Declared: clean_foam_data | — | — | /all_foam_data.csv | Shared / SI candidate |
| `tools/CleanData/combine_surfs.py` | Main | Declared: read_off, combine_off_files, main | — | — | — | Shared / SI candidate |
| `tools/CleanData/CombineLogs.py` | Main | Declared: parse_string_lists, combine_build_information, combine_group_information, combine_atoms_lines, combine_surface_lines, combine_edges_lines, combine_vertex_lines, warn_on_missing_atoms … | — | aw, vdw | /Total_logs.csv | Shared / SI candidate |
| `tools/CleanData/CombineLogs2.py` | Main | Declared: read_logs_sections, write_logs_sections, parse_string_lists, combine_build_information, combine_group_information, combine_atoms_lines, combine_surface_lines, combine_edges_lines … | — | vdw | /Total_logs.csv | Shared / SI candidate |
| `tools/CleanData/CreateFinalFoamData.py` | module | Declared: create_final_foam_data_file | — | — | — | Shared / SI candidate |
| `tools/CleanData/CreateVertsFromLogs.py` | Main | Declared: write_log_to_vert, create_verts_from_logs | — | aw, pow | /aw/aw_logs.csv; /aw/aw_verts.txt; /aw_logs.csv; /pow/pow_logs.csv; /pow/pow_verts.txt; /pow_logs.csv | Shared / SI candidate |
| `tools/CleanData/ExtractLogs.py` | Main | Declared: parse_ndxs_spec, _safe_int, _safe_float, _filter_atoms, _filter_surfaces, _filter_edges, _filter_vertices, _recompute_group_information … | dna, protein | aw, pow, prm, vdw | E:/Molecular/K_NCP/DNA/frames/f; E:/Molecular/K_NCP/Protein/frames/f; Extracted_logs.csv; H:/Data/Molecular/K_NCP/frames/f; _logs.csv | Shared / SI candidate |
| `tools/CleanData/FinalCleaning.py` | Main | Declared: get_txt_from_pdb, clean_folder | — | aw, pow | /aw/aw_logs.csv; /aw/aw_verts.txt; /set_atoms.pml; /set_balls.pml; aw_logs.csv; aw_verts.txt; balls.pdb; balls.txt; pow_logs.csv; pow_verts.txt; set_balls.pml | Shared / SI candidate |
| `tools/CleanData/Fix_Logs_SA.py` | Main | Declared: get_correct_SA, replace_csv_cell_by_index, correct_logs_SA, correct_multiple_logs_SA | — | aw, pow, prm | .pdb; aw/aw_logs.csv; aw/aw_logs_new.csv; pow/pow_logs.csv; pow/pow_logs_new.csv; prm/prm_logs.csv; prm/prm_logs_new.csv | Shared / SI candidate |
| `tools/CleanData/GetBasicFiles.py` | module | Top-level script; see principal family inventory | — | aw, pow | — | Shared / SI candidate |
| `tools/CleanData/MakeFoamRuns.py` | Main | Declared: make_foam_runs | — | aw | .csv; .pdb; /verts.txt | Shared / SI candidate |
| `tools/CleanData/PurgeLogs.py` | module | Purges a logs file of a | — | — | — | Shared / SI candidate |
| `tools/CleanData/recolor_surfs.py` | Main | Declared: recolor_off, main | — | — | — | Shared / SI candidate |
| `tools/CleanData/RenameFiles.py` | Main | Declared: move_contents, rename_pdb, get_nets_in_order | — | aw, pow | /aw/aw_logs.csv; /aw/aw_verts.txt; /aw_verts.txt; /pow/pow_logs.csv; /pow/pow_verts.txt; /pow_verts.txt; aw_logs.csv; balls.pdb; group_balls.pdb; info.txt; pow_logs.csv; set_balls.pml; surrounding_balls.pdb | Shared / SI candidate |
| `tools/CleanData/RenameFolders.py` | Main | Declared: get_the_best_folders, move_contents, copy_contents, rename_folders | — | aw, pow | /file_move_ledger.txt | Shared / SI candidate |
| `tools/CleanData/verify_frames.py` | Main | Declared: verify_frames | — | aw, pow, prm | /aw/aw_logs.csv; /pow/pow_logs.csv; /prm/prm_logs.csv | Shared / SI candidate |
| `tools/CleanData/VerifyFiles.py` | Main | Declared: verify_pdb, verify_logs | — | — | — | Shared / SI candidate |
| `tools/CleanData/VerifyFoamData.py` | module | Top-level script; see principal family inventory | — | — | — | Shared / SI candidate |
| `tools/compare/__init__.py` | module | Package initializer / exports | — | — | — | Shared / SI candidate |
| `tools/compare/color_by_residue.py` | module | Declared: color_pdb_by_res | — | — | _colorized.pdb | Shared / SI candidate |
| `tools/compare/compare_files.py` | module | Declared: compare_files, compare_files2 | — | — | — | Shared / SI candidate |
| `tools/compare/compare_verts.py` | Main | Declared: check_verts | — | — | C:/Users/jacke/PycharmProjects/vorpy/Data/user_data/cambrin_10/sys/cambrin_logs.csv; C:/Users/jacke/PycharmProjects/vorpy/Data/user_data/cambrin_11/sys/cambrin_logs.csv | Shared / SI candidate |
| `tools/compare/get_res_data.py` | module | Declared: residue_data | — | — | — | Shared / SI candidate |
| `tools/compare/in_out.py` | module | Declared: in_out_data | — | — | — | Shared / SI candidate |
| `tools/compare/pdb_names.py` | module | Top-level script; see principal family inventory | — | — | — | Shared / SI candidate |
| `tools/compare/read_logs.py` | module | Declared: read_atom, read_surf, read_edge, read_vert, _read_legacy_logs, read_logs | — | — | — | Shared / SI candidate |
| `tools/compare/read_logs2.py` | module | Backward-compatible VorPy log reader with modern header-driven parsing. | — | vdw | — | Shared / SI candidate |
| `tools/compare/sol_facing_percentage.py` | module | Declared: _find_first_matching_column, _coerce_index_value, _normalize_residue_name, _parse_neighbor_list, _build_solvent_lookup, _resolve_surface_columns, _build_surface_index, _get_neighbor_map_from_atom_df … | — | — | — | Shared / SI candidate |
| `tools/compare/sort_verts.py` | module | Top-level script; see principal family inventory | — | — | C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/Outputs/EDTA_Mg/EDTA_Mg_verts.txt; C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/vdos results/EDTA_Mg_generic_atoms/EDTA_Mg_generic.txt | Shared / SI candidate |
| `tools/plot_templates/__init__.py` | module | Package initializer / exports | — | — | — | Shared / SI candidate |
| `tools/plot_templates/bar.py` | module | Declared: bar | — | — | — | Shared / SI candidate |
| `tools/plot_templates/blank_axes.py` | Main | Declared: create_graph | — | — | — | Shared / SI candidate |
| `tools/plot_templates/blank_legend.py` | Main | Declared: legend_plot | — | — | — | Shared / SI candidate |
| `tools/plot_templates/box_whisker.py` | module | Declared: box_whisker | — | — | — | Shared / SI candidate |
| `tools/plot_templates/histogram.py` | module | Declared: histogram | — | — | — | Shared / SI candidate |
| `tools/plot_templates/inline_plot.py` | module | Declared: inline_plot | — | — | — | Shared / SI candidate |
| `tools/plot_templates/line.py` | module | Declared: line_plot | — | — | — | Shared / SI candidate |
| `tools/plot_templates/scatter.py` | module | Declared: scatter | — | — | — | Shared / SI candidate |
| `tools/plot_templates/starting_point.py` | Main | Top-level script; see principal family inventory | — | — | — | Shared / SI candidate |
| `tools/plot_templates/table.py` | module | Declared: color_cells, table | cambrin, edta, hairpin, ncp, p53tet, streptavidin | — | — | Shared / SI candidate |

System/scheme columns are lexical observations, not cohort assertions: words can occur in examples, docstrings or directory names. Detailed cohort claims are confined to the principal inventory. There are 321 entries, including 25 package initializers.

## Appendix B. Per-script import dependencies from the audit

This lists static import statements relevant to each analysis source file. Relative imports are written with their original dot level. Nested/optional/fallback imports are included without asserting that the branch executes. Dynamic imports and generated code require separate review. The broken-reference table above distinguishes confirmed missing targets from valid conditional imports.

| Current relative path | Internal, relative, or sibling imports | Other imported module roots |
|---|---|---|
| `__init__.py` | vorpy.src.analyze.apollonius [AlphaComplex, ApolloniusComplex, ApolloniusSimplex, PrimalFeatureRef, UnsupportedFeature]; vorpy.src.analyze.aw_interface_curvature [AWInterfaceCurvatureResult, InterfaceEdgeRecord, InterfaceSurfaceRecord, analyze_aw_interface_curvature] | — |
| `apollonius_export.py` | — | csv, json, pathlib |
| `apollonius.py` | vorpy.src.calculations.edge_geometry [AdditivelyWeightedTrisectorBranch, AdditivelyWeightedTrisectorConic, AffineEdgeGeometry, LineEdgeGeometry]; vorpy.src.network.edge_geometry_diagnostics [resolve_aw_network_edge] | dataclasses, itertools, math, numpy |
| `aw_interface_curvature.py` | vorpy.src.analyze.apollonius [ApolloniusComplex]; vorpy.src.calculations.aw_interface_orientation [aw_interface_edge_orientation]; vorpy.src.calculations.edge_geometry [AdditivelyWeightedTrisectorBranch, match_aw_trisector_conic_to_samples]; vorpy.src.network.edge_geometry_diagnostics [resolve_aw_network_edge]; vorpy.src.calculations.surface_mean_curvature [calculate_aw_network_surface_mean_curvatures] | collections, csv, dataclasses, math, numpy, pathlib |
| `aw_pipeline_audit.py` | vorpy.src.network.edge_geometry_diagnostics [resolve_aw_network_edge] | __future__, csv, math, numpy, pandas, pathlib |
| `cazals_edge_forensics.py` | — | __future__, collections, csv, itertools, math, numpy, pathlib |
| `cazals_validation.py` | vorpy.src.analyze.power_interface_curvature [analyze_power_interface_curvature, analyze_power_regular_complex, cazals_signed_beta]; vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb [prepare_cazals_pdb]; vorpy.src.analyze.cazals_edge_forensics [run_2kai_forensics] | __future__, argparse, collections, csv, hashlib, itertools, json, math, matplotlib, numpy, pathlib |
| `comparative_geometry_study.py` | vorpy.src.analyze.cazals_validation [interface_components, interface_atom_ids]; vorpy.src.analyze.power_vs_aw [_aw_surface_components]; vorpy.src.system [System]; vorpy.src.group [Group]; vorpy.src.io [load_network, save_network]; vorpy.src.analyze.power_interface_curvature [analyze_power_regular_complex]; vorpy.src.analyze.power_interface_selection [power_selection_views]; vorpy.src.analyze.aw_interface_curvature [analyze_aw_interface_curvature]; vorpy.src.analyze.power_interface_curvature [analyze_power_interface_curvature]; vorpy.src.analyze.power_interface_cli [run_power_interface_pdb]; vorpy.src.analyze.power_vs_aw [_power_area]; vorpy.src.analyze.comparative_study_figures [make_figures]; vorpy.src.analyze.power_interface_cli [prepare_cazals_pdb]; vorpy.src.analyze.finite_solvent_boundary [solvent_hull, hull_clearances]; vorpy.src.analyze.finite_solvent_boundary [audit_feature_samples] | __future__, argparse, collections, csv, dataclasses, hashlib, json, math, numpy, os, pathlib, platform, re, shutil, statistics, sys, time |
| `comparative_study_figures.py` | — | __future__, csv, math, matplotlib, numpy, pathlib |
| `condition_b_geometry.py` | — | __future__, math, numpy |
| `curvature_showcase_2kai.py` | vorpy.src.calculations.edge_geometry [oriented_aw_face_edge_geodesic_curvature]; vorpy.src.calculations.vertex_geometry [aw_vertex_face_angle]; vorpy.src.network.edge_geometry_diagnostics [resolve_aw_network_edge]; vorpy.src.geometry.visualization.export [_edge_mesh]; vorpy.src.io [load_network]; vorpy.src.geometry.visualization.export [_surface_mesh]; vorpy.src.output.mesh [write_mesh] | __future__, argparse, collections, csv, difflib, json, math, matplotlib, numpy, pathlib, pickle |
| `dual_restriction_audit.py` | vorpy.src.analyze.power_interface_curvature [analyze_power_interface_curvature]; vorpy.src.analyze.power_interface_cli [run_power_interface_pdb]; vorpy.src.analyze.cazals_validation [power_facet_polygon_area]; vorpy.src.geometry.visualization.export [_edge_mesh]; vorpy.src.output.mesh [write_mesh]; vorpy.src.geometry.aw_alpha.births [surface_birth]; vorpy.src.io [load_network]; vorpy.src.geometry.visualization.export [_edge_mesh, _surface_mesh]; vorpy.src.analyze.apollonius [PrimalFeatureRef] | __future__, argparse, ast, csv, json, math, numpy, pathlib, types |
| `finite_solvent_boundary.py` | — | __future__, numpy, scipy |
| `matched_probe_power_aw.py` | vorpy.src.analyze.comparative_geometry_study [canonical_identity, heavy_atoms, read_pdb]; vorpy.src.analyze.apollonius [PrimalFeatureRef]; vorpy.src.analyze.cazals_validation [power_facet_polygon_area]; vorpy.src.geometry.visualization.export [_edge_mesh, _surface_mesh]; vorpy.src.output.mesh [MeshData, combine_mesh_parts, write_mesh]; vorpy.src.io [load_network]; vorpy.src.analyze.power_interface_curvature [analyze_power_interface_curvature]; vorpy.src.analyze.aw_interface_curvature [analyze_aw_interface_curvature]; vorpy.src.analyze.aw_interface_curvature [_oriented_surface_curvatures]; vorpy.src.geometry.aw_alpha.births [surface_birth, _pair_contact_candidate, _clearances] | __future__, argparse, collections, csv, hashlib, itertools, json, math, numpy, pathlib, pickle, types |
| `nonpolar_interface/__init__.py` | .classification [classify_atom, classify_surface_pair]; .export [export_geometry_tables]; .explore [summarize_surface_features]; .table [build_geometry_tables] | — |
| `nonpolar_interface/classification.py` | vorpy.src.chemistry.chemistry_interpreter [ion_residue_names, water_residue_names] | — |
| `nonpolar_interface/explore.py` | — | __future__, pandas |
| `nonpolar_interface/export.py` | — | __future__, json, numpy, pathlib |
| `nonpolar_interface/table.py` | .classification [NONPOLAR_RULE_VERSION, classify_atom, classify_surface_pair] | __future__, collections, numpy, pandas |
| `open_interface_curvature.py` | vorpy.src.calculations.curvature [integrate_quadratic_surface_gaussian_curvature]; vorpy.src.analyze.aw_interface_curvature [_validated_edge_geometry]; vorpy.src.calculations.edge_geometry [oriented_aw_face_edge_geodesic_curvature]; vorpy.src.calculations.vertex_geometry [aw_vertex_face_angle]; vorpy.src.network.edge_geometry_diagnostics [resolve_aw_network_edge]; vorpy.src.io [load_network] | __future__, argparse, collections, csv, dataclasses, json, math, numpy, pathlib, pickle, scipy, typing |
| `Part_I_Statistical_Ensembles/Fig1/AWVd_vs_Power/Surfaces.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs]; vorpy.src.calculations.vert [calc_vert]; vorpy.src.network.build_edge [build_edge]; vorpy.src.network.build_surf [build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/AWVd_vs_Power1/Surfaces.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs]; vorpy.src.calculations.vert [calc_vert]; vorpy.src.network.build_edge [build_edge]; vorpy.src.network.build_surf [build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/Edges.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges]; vorpy.src.calculations.vert [calc_vert]; vorpy.src.network.build_edge [build_edge_old]; vorpy.src.calculations.edge [calc_circ] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/ellipse_edge.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_circles, plot_edges, plot_surfs]; vorpy.src.calculations.edge [calc_circ]; vorpy.src.network.build_edge [build_edge] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/Inscribed_Circle.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_circles]; vorpy.src.calculations.edge [calc_circ] | matplotlib, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Edges/loop_edge.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_edges, plot_verts]; vorpy.src.network.build_edge [build_edge] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Encapsulation.py` | vorpy.src.visualize.mpl_visualize [plot_balls] | matplotlib, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Fill_Surface_points.py` | — | matplotlib, numpy |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Plane_Between_Centers.py` | vorpy.src.visualize.mpl_visualize [plot_surfs, plot_balls] | matplotlib, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Surfaces/Surfaces.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs]; vorpy.src.calculations.vert [calc_vert, calc_flat_vert]; vorpy.src.network.build_edge [build_edge_old]; vorpy.src.network.build_surf [build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Vertices/doublets.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs]; vorpy.src.calculations.vert [calc_vert]; vorpy.src.network.build_edge [build_edge]; vorpy.src.network.build_surf [build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Vertices/TangentialCircle.py` | — | matplotlib, numpy |
| `Part_I_Statistical_Ensembles/Fig1/Building_Blocks/Vertices/Vertices.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts]; vorpy.src.calculations.vert [calc_vert] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Cells/cells.py` | vorpy.src.system.system [System]; vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs] | matplotlib, numpy, os, pandas, sys |
| `Part_I_Statistical_Ensembles/Fig1/Edges/Edges.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges]; vorpy.src.calculations.vert [calc_vert]; vorpy.src.network.build_edge [build_edge]; vorpy.src.calculations.edge [calc_circ] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Edges/ellipse_edge.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_circles, plot_edges, plot_surfs]; vorpy.src.calculations.edge [calc_circ]; vorpy.src.network.build_edge [build_edge] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Edges/Inscribed_Circle.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_circles]; vorpy.src.calculations.edge [calc_circ] | matplotlib, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Edges/loop_edge.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_edges, plot_verts]; vorpy.src.network.build_edge [build_edge] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/proximity.py` | — | matplotlib, numpy, random, tkinter |
| `Part_I_Statistical_Ensembles/Fig1/SortBalls.py` | vorpy.src.calculations.calcs [calc_dist] | matplotlib, numpy, os, random, sys |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Encapsulation.py` | vorpy.src.visualize.mpl_visualize [plot_balls] | matplotlib, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Fill_Surface_points.py` | — | matplotlib, numpy |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Plane_Between_Centers.py` | vorpy.src.visualize.mpl_visualize [plot_surfs, plot_balls] | matplotlib, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Surfaces/Surfaces.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs]; vorpy.src.calculations.vert [calc_vert, calc_flat_vert]; vorpy.src.network.build_edge [build_edge]; vorpy.src.network.build_surf [build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Vertices/doublets.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts, plot_edges, plot_surfs]; vorpy.src.calculations.vert [calc_vert]; vorpy.src.network.build_edge [build_edge]; vorpy.src.network.build_surf [build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig1/Vertices/Vertices.py` | vorpy.src.visualize.mpl_visualize [plot_balls, plot_verts]; vorpy.src.calculations.vert [calc_vert] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig3/DistributionOfXYZ.py` | vorpy.src.analyze.tools.batch.compile_logs [get_logs_and_pdbs]; vorpy.src.system.system [System] | csv, os, sys |
| `Part_I_Statistical_Ensembles/Fig3/histofradiiCompare.py` | vorpy.src.analyze.tools.batch.compile_new_logs [get_logs_and_pdbs]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.system.system [System] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig3/HistofRadiiRaw.py` | vorpy.src.analyze.tools.batch.compile_new_logs [get_logs_and_pdbs]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.system.system [System] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig3/NumOlaps_csv.py` | — | csv, matplotlib, numpy, tkinter |
| `Part_I_Statistical_Ensembles/Fig3/NumOverlaps.py` | vorpy.src.system.system [System]; vorpy.src.calculations.calcs [calc_dist]; vorpy.src.analyze.tools.batch.get_files [get_files]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | csv, matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig3/OneOlapPlot.py` | — | csv, matplotlib, numpy, tkinter |
| `Part_I_Statistical_Ensembles/Fig3/OneXYZPlot.py` | — | csv, matplotlib, numpy |
| `Part_I_Statistical_Ensembles/Fig3/ranging_pdfs.py` | — | matplotlib, numpy, pandas, plotly, scipy |
| `Part_I_Statistical_Ensembles/Fig4/color_by_diff.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.inputs [read_pdb_line]; vorpy.src.output.pdb [make_pdb_line]; vorpy.src.analyze.tools.batch.get_files [get_files] | csv, matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig4/NetNeighbors.py` | vorpy.src.analyze.tools.batch.get_files [get_files]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.system.system [System] | matplotlib, numpy, os, random, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig4/overlaid_line_plots.py` | — | csv, matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/Fig4/plot_distributions_complete_cells.py` | vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.inputs [read_pdb_line]; vorpy.src.system.system [System] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig4/RadiusVsPercentDifference.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.inputs [read_pdb_line]; vorpy.src.output.pdb [make_pdb_line]; vorpy.src.analyze.tools.batch.get_files [get_files] | csv, matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig5/overlaid_line_plots_physical.py` | — | csv, matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/Fig5/physical_pdfs.py` | — | matplotlib, numpy, scipy |
| `Part_I_Statistical_Ensembles/Fig5/plot_hist_of_real_dists.py` | — | matplotlib, numpy, plotly, scipy |
| `Part_I_Statistical_Ensembles/Fig5/plot_physical_mpl.py` | vorpy.src.analyze.tools.plot_templates.line [line_plot] | csv, numpy, os, pandas, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig6/DivergenceVsSphericity.py` | vorpy.src.inputs [read_pdb_line]; vorpy.src.calculations [calc_sphericity, calc_isoperimetric_quotient]; vorpy.src.calculations [sort_lists]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | csv, matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig6/overlaid_line_plots_extra_data.py` | — | csv, matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/1_SphericityDifference_NerighborRadCVDiff.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/2_NeighborDistanceSphericity.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/3_Radius_UnsharedNeighbors.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/4_MissingNeighborDistance.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/5_pizza_pie.py` | vorpy.src.calculations [calc_dist, calc_aw_center, calc_pw_center, calc_tetra_vol, calc_surf_func]; vorpy.src.network [calc_surf_point, build_surf]; vorpy.src.visualize [plot_surfs, plot_edges, plot_balls] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig7/6_Overlap_Volume.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/7_RadialDistance_VolDiff.py` | vorpy.src.calculations [calc_dist, calc_aw_center, calc_pw_center, calc_tetra_vol, calc_surf_func]; vorpy.src.network [calc_surf_point, build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig7/B_Radius_UnsharedNeighbors.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/C_MissingNeighborDistance.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/D_ExclusiveNeighbors_SphericityDifference.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/ExclusiveNeighbors_SphericityDifference.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/F_RadialDistance_VolDiff.py` | vorpy.src.calculations [calc_dist, calc_aw_center, calc_pw_center, calc_tetra_vol, calc_surf_func]; vorpy.src.network [calc_surf_point, build_surf] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig7/G_pizza_pie.py` | vorpy.src.calculations [calc_dist, calc_aw_center, calc_pw_center, calc_tetra_vol, calc_surf_func]; vorpy.src.network [calc_surf_point, build_surf]; vorpy.src.visualize [plot_surfs, plot_edges, plot_balls] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig7/H_Overlap_Volume.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/Other/BallRadius_AbsoluteVolumeDifference.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/Other/LocalDensityPolydispersity.py` | vorpy.src.analyze.Part_I_Statistical_Ensembles.Review1.Other.NeighborsAverageRadii [get_data, select_file] | matplotlib, numpy, os, sys |
| `Part_I_Statistical_Ensembles/Fig7/Other/NeighborDistanceSphericity.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple]; vorpy.src.calculations [calc_aw_center, calc_pw_center, calc_dist] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/Other/NeighborsAverageRadii.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/Other/Radius_NeighborCount.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/Fig7/Other/SphericityDifference_NerighborRadCVDiff.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple] | matplotlib, numpy, os, sys, tkinter |
| `Part_I_Statistical_Ensembles/old/Check_Max_Vertex_Parameters.py` | Data.Analyze.tools.batch.compile_logs [get_logs_and_pdbs] | csv |
| `Part_I_Statistical_Ensembles/old/check_num_foams.py` | — | os, tkinter |
| `Part_I_Statistical_Ensembles/old/density_deviation.py` | — | matplotlib, numpy, os, scipy |
| `Part_I_Statistical_Ensembles/old/Open_density_factor.py` | Data.Analyze.tools.plot_templates.line [line_plot] | csv, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/old/plot_physical_plotly.py` | — | csv, numpy, pandas, plotly, tkinter |
| `Part_I_Statistical_Ensembles/old/plot_surf_graphy_plotly.py` | — | csv, numpy, os, pandas, plotly, tkinter |
| `Part_I_Statistical_Ensembles/old/StandardDeviationFoamGens.py` | — | matplotlib, scipy |
| `Part_I_Statistical_Ensembles/Other/check_deviation.py` | — | csv, itertools, matplotlib, numpy, tkinter |
| `Part_I_Statistical_Ensembles/Other/DensityVoxelization.py` | vorpy.src.system.system [System]; System.sys_funcs.calcs.calcs [calc_dist]; System.sys_funcs.calcs.sorting [get_balls, box_search]; Data.Analyze.tools.plot_templates.histogram [histogram] | csv, math, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/Other/Overlaid_Line_Plots_Density_Voxelization.py` | — | csv, matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/MissingBallsHeatMap.py` | Data.Analyze.tools.batch.get_files [get_files]; Data.Analyze.tools.compare.read_logs2 [read_logs2]; Data.Analyze.tools.compare.read_logs [read_logs]; vorpy.src.system.system [System] | matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/MonteCarloPlacement.py` | — | matplotlib, numpy |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/NeighborsMatrix.py` | Data.Analyze.tools.batch.get_files [get_files]; Data.Analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/OverlapsMatrix.py` | vorpy.src.system.system [System]; Data.Analyze.tools.batch.get_files [get_files]; Data.Analyze.tools.compare.read_logs2 [read_logs2]; System.sys_funcs.calcs.calcs [calc_dist] | csv, matplotlib, numpy, os, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/RadiusDistributionMatrix.py` | vorpy.src.system.system [System]; Data.Analyze.tools.batch.get_files [get_files]; Data.Analyze.tools.compare.read_logs2 [read_logs2]; System.sys_funcs.calcs.calcs [calc_dist]; System.sys_funcs.calcs.calcs [get_time] | matplotlib, numpy, os, scipy, time, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/SurfaceOverlapsMatrix.py` | Data.Analyze.tools.batch.compile_logs [get_logs_and_pdbs]; Data.Analyze.tools.compare.read_logs2 [read_logs2]; System.sys_funcs.calcs.calcs [calc_dist]; vorpy.src.system.system [System] | csv, matplotlib, numpy, os, pandas, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/XYZMatrix.py` | vorpy.src.system.system [System]; Data.Analyze.tools.batch.get_files [get_files]; Data.Analyze.tools.compare.read_logs2 [read_logs2]; System.sys_funcs.calcs.calcs [calc_dist] | matplotlib, numpy, os, scipy, tkinter |
| `Part_I_Statistical_Ensembles/SupplementaryMaterial/XYZMatrixFake.py` | vorpy.src.system.system [System]; Data.Analyze.tools.batch.get_files [get_files]; Data.Analyze.tools.compare.read_logs2 [read_logs2]; System.sys_funcs.calcs.calcs [calc_dist] | matplotlib, numpy, os, random, scipy, tkinter |
| `Part_I_Statistical_Ensembles/VerifyVorPy/VerifyVerts.py` | — | csv, matplotlib, matplotlib_venn, numpy, os, pandas, tkinter |
| `Part_I_Statistical_Ensembles/VerifyVorPy/VerifyVerts1.py` | — | csv, matplotlib, matplotlib_venn, numpy, os, pandas, tkinter |
| `Part_I_Statistical_Ensembles/VerifyVorPy/VerifyVolumes.py` | System.sys_funcs.calcs.calcs [calc_dist]; Data.Analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, matplotlib_venn, numpy, os, pandas, tkinter |
| `Part_II_Molecular_Analysis/Curvature/52ca_alanine_curvature_vs_area.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | dataclasses, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/52ca_alanine_dna_facing_violin.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | ast, dataclasses, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/5b_color_pdbs_by_curvature.py` | vorpy.src.output.pdb [make_pdb_line]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.inputs.pdb [read_pdb_line] | os, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/5c1_curvature_localization_along_sequence.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | collections, dataclasses, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/5c2_curvature_concentration.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | dataclasses, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/5c2a_alanine_curvature_vs_area_by_chain.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | dataclasses, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/5c2a_alanine_dna_facing.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/5c3_violin.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | dataclasses, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/curvature_assignments.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/curvature_name_bins.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs [read_logs] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/CurvatureByInterface.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.get_res_data [nucleics, proteins, ions, sols]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.box_whisker [box_whisker] | os, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_curvature_assignments/CurvatureDistributions.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.box_whisker [box_whisker]; vorpy.src.analyze.tools.plot_templates.line [line_plot] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4A_Peaks.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.plot_templates.line [line_plot] | collections, matplotlib, numpy, os, pandas, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4B_Peak_Assignment.py` | vorpy.src.analyze.tools.plot_templates.line [line_plot] | collections, matplotlib, numpy, os, pandas, re, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/4C_Combine_Curvature_Peak_Files.py` | — | __future__, argparse, dataclasses, numpy, os, pandas, tkinter, typing |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/Extra/surface_energy.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | ast, matplotlib, numpy, os, pandas, sys, tkinter |
| `Part_II_Molecular_Analysis/Curvature/F4_atomic_curvature/SurfaceEnergyCurvature.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | ast, matplotlib, numpy, os, pandas, sys, tkinter |
| `Part_II_Molecular_Analysis/extra/AtomAssignmentClustering.py` | — | matplotlib, numpy, pandas, sklearn, tkinter |
| `Part_II_Molecular_Analysis/extra/AtomAssignmentFileMaker.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.chemistry.chemistry_interpreter [amino_names, nucleo_names, ion_names, sol_names]; vorpy.src.objects.atom [get_element]; vorpy.src.calculations.calcs [calc_com] | csv, matplotlib, numpy, os, scipy, sklearn, sys, tkinter |
| `Part_II_Molecular_Analysis/extra/AvgAtomDiff.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.analyze.tools.compare.read_logs [read_logs] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_II_Molecular_Analysis/extra/ElementBinningDiffs.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.calculations.sorting [sort_lists] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_II_Molecular_Analysis/extra/ElementBinningTots.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.calculations.sorting [sort_lists] | matplotlib, numpy, os, scipy, sys, tkinter |
| `Part_II_Molecular_Analysis/extra/KMeansClustering.py` | — | matplotlib, numpy, os, pandas, seaborn, sklearn, tkinter |
| `Part_II_Molecular_Analysis/extra/TrainAWVolModel.py` | — | matplotlib, os, pandas, sklearn, tkinter |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/curvature_range_cartoons.py` | vorpy.src.network.fill [calc_surf_point]; vorpy.src.calculations.surf [calc_surf_func, calc_surf_func_reg]; vorpy.src.calculations.calcs [calc_dist]; vorpy.src.calculations.curvature [mean_curvature] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/curvature_range_plot.py` | vorpy.src.calculations.surf [calc_surf_func]; vorpy.src.calculations.curvature [mean_curvature]; vorpy.src.analyze.tools.plot_templates.line [line_plot] | matplotlib, numpy, os, sys |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/Curvature_Shapes_mpl.py` | — | matplotlib, numpy |
| `Part_II_Molecular_Analysis/F1_schemes_curvature_models/1C_curvature_representations/Show_Ratio_Dist.py` | — | matplotlib, numpy |
| `Part_II_Molecular_Analysis/F2_frame_comparison/AtomicCorrelation_multi_plots.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | __future__, matplotlib, numpy, pandas, pathlib, re, sys, tkinter |
| `Part_II_Molecular_Analysis/F2_frame_comparison/AtomicCorrelation_updated.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | __future__, matplotlib, numpy, pandas, pathlib, re, sys, tkinter |
| `Part_II_Molecular_Analysis/F2_frame_comparison/AtomicCorrelation.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | __future__, matplotlib, numpy, pandas, pathlib, re, sys, tkinter |
| `Part_II_Molecular_Analysis/F2_frame_comparison/CurvatureComparison_updated.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.read_logs2 [read_logs2] | __future__, matplotlib, numpy, os, pandas, pathlib, re, sys, tkinter |
| `Part_II_Molecular_Analysis/F2_frame_comparison/CurvatureComparison.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.read_logs2 [read_logs2] | __future__, matplotlib, numpy, os, pandas, pathlib, re, sys, tkinter |
| `Part_II_Molecular_Analysis/F2_frame_comparison/FrameAnalysis.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.plot_templates.line [line_plot] | matplotlib, os, sys, tkinter, warnings |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5A_terminal_effects.py` | .F5_config [KEY, UNITS, SCHEME_LABELS]; .F5_common [frame_summary, temporal_summary, write_csv, save_figure, empty_panel]; .F5_run_all [main] | matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5B_solvation_environments.py` | .F5_config [KEY]; .F5_common [frame_summary, write_csv, distribution_panels]; .F5_run_all [main] | pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5C_polymer_length.py` | .F5_config [KEY, UNITS, SCHEME_LABELS]; .F5_common [exposure, temporal_summary, write_csv, save_figure, empty_panel]; .F5_run_all [main] | matplotlib, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5D_internal_water.py` | .F5_config [KEY, SCHEME_LABELS]; .F5_common [write_csv, save_figure, empty_panel]; .F5_run_all [main] | matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5E_scheme_differences.py` | .F5_config [KEY, METRICS, UNITS, SCHEME_LABELS]; .F5_common [write_csv, save_figure, empty_panel, distribution_panels]; .F5_run_all [main]; .F5_config [UNITS] | matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/archive/F5F_water_coordination.py` | .F5_config [KEY, SCHEME_LABELS]; .F5_common [frame_summary, write_csv, save_figure, empty_panel]; .F5_run_all [main] | matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_common.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; .F5_config [Config, KEY, ATOM_KEY, UNIT_KEY, METRICS]; .F5_config [UNITS, SCHEME_LABELS] | collections, csv, dataclasses, hashlib, json, matplotlib, numpy, pandas, pathlib, pickle, re, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_config.py` | — | dataclasses, pathlib |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_peg_dataset.py` | .F5_common [read_pdb]; .F5C_water_contacts [water_contacts]; .F5B_water_burial [load_water_geometry] | numpy, pandas, pathlib |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5_run_all.py` | .F5_config [Config, KEY]; .F5_common [read_pdb, load_scheme, source_signature, independent_area_check, write_csv, attach_aw_exposure, match_schemes, setup_plots, exposure] | argparse, dataclasses, importlib, json, numpy, pandas, pathlib, pickle, sys, time, traceback |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5A_exposed_vs_buried.py` | .F5_config [Config, ATOM_KEY, KEY, METRICS]; .F5_common [match_schemes, setup_plots, save_figure, write_csv] | json, matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5A_PEG_chemistry_vs_hydration.py` | .F5A_exposed_vs_buried [assign_repeat_types]; .F5_common [load_scheme, match_schemes, read_pdb, save_figure, setup_plots]; .F5_config [ATOM_KEY, Config]; .F5_peg_dataset [screened_runs, paths, layer1_atoms] | matplotlib, numpy, pandas, pathlib, shutil, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5B_heparin_chemistry.py` | .F5_common [load_scheme, match_schemes, read_pdb]; .F5_config [Config] | matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5B_water_burial.py` | .F5_config [Config, KEY]; .F5_common [normalize_log, read_pdb, reconstruct_interfaces, save_figure, setup_plots, write_csv]; .F5_run_all [discover]; .F5_peg_dataset [screened_runs, matched_water_layers] | csv, dataclasses, matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5C_rerun_missing.py` | scripts.run_checkpoint_layers [build_command] | argparse, pathlib, subprocess, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/F5C_water_contacts.py` | .F5_common [normalize_log, read_pdb]; .F5_config [Config]; .F5B_water_burial [has_water_atom_rows]; .F5_peg_dataset [screened_runs, validated_water_layer] | matplotlib, numpy, pandas, pathlib, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_classification.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [exposure, attach_aw_exposure, reconstruct_interfaces, aggregate_residues]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config [Config]; .test_F5_loading [case] | numpy, pandas |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_direct_execution.py` | — | pathlib, pytest, subprocess, sys |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_loading.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config [Config]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [classify_atoms, normalize_log, audit_raw_log, topology_record, end_distances, annotate_units, read_pdb]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all [run]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [source_signature]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [load_scheme, independent_area_check]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all [FIGURES, run] | csv, dataclasses, numpy, pandas, pytest |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_scheme_comparison.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [reconstruct_interfaces, match_schemes]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config [Config]; .test_F5_loading [case]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [aggregate_residues]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all [add_residue_pairs] | pandas, pytest |
| `Part_II_Molecular_Analysis/F3_SimpleModels/tests/test_F5_surface_areas.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common [reconstruct_interfaces, aggregate_residues]; vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config [Config]; .test_F5_loading [case] | dataclasses, numpy, pandas, pytest |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/AllinaRow.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group] | matplotlib, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/compare_frames.py` | vorpy.src.analyze.tools.compare.read_logs [read_logs] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/Distribution_of_element_radii.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.plot_templates.histogram [histogram] | matplotlib, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/Distribution_Types.py` | — | matplotlib, numpy |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/frame_analysis.py` | vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.line [line_plot]; vorpy.src.system.system [System] | matplotlib, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/num_atoms_log.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.calculations.sorting [get_sys_type]; vorpy.src.analyze.tools.plot_templates.inline_plot [inline_plot] | os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/Plot_by_system.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs [read_logs] | matplotlib, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/ResidueBinningDiffs.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.get_res_data [residue_data]; vorpy.src.calculations.sorting [sort_lists]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.analyze.tools.compare.read_logs [read_logs] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/ResidueBinningTots.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.compare.get_res_data [residue_data]; vorpy.src.calculations.sorting [sort_lists] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Extra/vis_table_all_models.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.plot_templates.table [table]; vorpy.src.calculations.sorting [get_sys_type] | os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/config.py` | — | typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/pipeline.py` | .config [SMALL_MOL_ATOM_ALIASES, PROTEIN_ATOM_ALIASES, RNA_ATOM_ALIASES, DNA_ATOM_ALIASES, DIRECT_LABEL_GROUPS]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, pandas, sklearn, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/run_batch.py` | .config [DATA_ROOT, MOLECULE_GROUPS, OUTPUT_ROOT, SETTINGS_COMBINATIONS]; .pipeline [make_settings_name, main, resolve_group_folders] | os |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F2C_All_Data_Sets_Clustered_full/run_single.py` | .pipeline [main] | — |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_atomic_umap.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common [discover_systems, read_pair] | __future__, argparse, importlib, json, matplotlib, numpy, pandas, pathlib, sklearn, sys, umap |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_common.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | __future__, dataclasses, numpy, os, pandas, re, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_curvature_pca.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common [discover_systems, read_pair] | __future__, argparse, json, matplotlib, numpy, os, pandas, pathlib, re, seaborn, sklearn, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3A_atom_diffs.py` | vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3A_residue_diffs.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3B_molecule_totals.py` | vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.batch.get_files [get_all_files]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, sys, warnings |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3B_residue_diffs_updated.py` | vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3C_molecule_diffs_updated.py` | vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.batch.get_files [get_all_files]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, sys, warnings |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3C_residue_sphericity_vs_deviation.py` | vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.analyze.tools.batch.get_files [get_all_files]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.calculations.calcs [calc_sphericity] | numpy, os, sys, tkinter, warnings |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3C_Residue_type_diffs.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | math, matplotlib, numpy, os, sys, tkinter, warnings |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/3D_molecule_sphericity_vs_deviation.py` | vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.batch.get_files [get_all_files]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.calculations.calcs [calc_sphericity] | os, sys, warnings |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2_AnalysisTable.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2A_AtomicVolumeDifferences.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.analyze.tools.batch.get_files [get_all_files] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2Ac_atomic_deviations.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.plot_templates.scatter [scatter] | dataclasses, json, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_All_Data_Sets_Clustered.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.compare.sol_facing_percentage [compute_sol_facing_percent]; src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_Helper_Auto_Cluster [compare_manual_vs_ml]; src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_Helper_Auto_Cluster [make_sol_facing_binary]; src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_Helper_Auto_Cluster [print_cluster_summary]; src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_Helper_Auto_Cluster [run_clustering]; src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_Helper_Auto_Cluster [summarize_clusters] | collections, matplotlib, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper_Auto_Cluster.py` | — | hdbscan, kmodes, numpy, pandas, sklearn |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper_Protein_Assignment.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | json, matplotlib, os, pandas, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | math, numpy, os, pandas, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/F2C_Helper2.py` | — | math, matplotlib, numpy, os, pandas, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Ligand_Interface_Residues.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Mean_Vol_SA.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/NCP_Histone_DNA_Interfaces.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/NCP_Stacked_DNA_Contact_Profile.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/NCP_System_Analysis.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/p53tet_Chain_Equivalency_deviation.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Res_Vol_SA.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3_old/Tester_for_power_residues.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | numpy, os, pandas, sys, tkinter, warnings |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3A_multiscale_summary.py` | F3_common [POWER_COLOR, add_deviation_columns, build_atomic_metrics, build_molecule_metrics, build_residue_metrics, discover_systems, read_pair, summarize_scale] | __future__, matplotlib, numpy, os, pandas, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3B_atomic_groups.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common [add_deviation_columns, build_atomic_metrics, discover_systems, read_pair]; vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_All_Data_Sets_Clustered [canonicalize_atom_name]; vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_atomic_umap [analyze] | __future__, argparse, matplotlib, numpy, os, pandas, pathlib, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3C_residue_groups.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common [add_deviation_columns, build_atomic_metrics, build_residue_metrics, discover_systems, read_pair] | __future__, matplotlib, numpy, os, pandas, pathlib, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/F3E_solvent_exposure.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_common [add_deviation_columns, build_atomic_metrics, build_residue_metrics, discover_systems, read_pair, surface_pairs]; vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Multiscale_Scheme_Differences.F3_old.F2C_All_Data_Sets_Clustered [canonicalize_atom_name] | __future__, matplotlib, numpy, os, pandas, pathlib, sys, tkinter, typing |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/AvgSurfDiff.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.analyze.tools.plot_templates.bar [bar]; vorpy.src.analyze.tools.batch.get_files [get_all_files] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/InsideOutside.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.bar [bar] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/OtherSADiffs.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.compare.get_res_data [nucleics, proteins, ions, sols]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.plot_templates.bar [bar] | numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F4_atomic_deviations/Other/SolventAccessible.py` | vorpy.src.calculations.calcs [calc_com, calc_dist] | matplotlib, numpy, os, sys |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/color_residues.py` | — | csv, os, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/diff_by_res_vor_pow.py` | — | csv, numpy, os, sys, tkinter, traceback |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/coarse_grained_residue_volume_difference.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | matplotlib, numpy, os, pandas, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/f8_build_surfaces.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | numpy, os, pandas, pathlib, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/f8b_raincloud_interfaces.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.F8_CG_Curvature.f8_build_surfaces [build_all_surfaces] | matplotlib, numpy, os, pandas, pathlib, scipy, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/color_code_residues.py` | — | csv, matplotlib, numpy, os |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/curvature_assignments.py` | vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.chemistry.chemistry_interpreter [nucleo_names, amino_names]; vorpy.src.system.system [System] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/curvature_by_radius.py` | vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.system.system [System]; vorpy.src.analyze.tools.plot_templates.scatter [scatter]; vorpy.src.analyze.tools.calcs [sort_lists] | os, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/curvature_distribution.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.group.group [Group]; vorpy.src.analyze.tools.plot_templates.box_whisker [box_whisker] | matplotlib, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/F7_Coarse_grained_residues/F8_CG_Curvature/Old/max_avg_curv.py` | vorpy.src.analyze.tools.compare.compare_files [compare_files] | matplotlib, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme.py` | vorpy.src.analyze.tools.compare.compare_files [compare_files]; vorpy.src.analyze.tools.plot_templates.bar [bar] | matplotlib, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme2.py` | vorpy.src.analyze.tools.compare.compare_files [compare_files]; vorpy.src.analyze.tools.plot_templates.bar [bar] | csv, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/F5_Coarse_Grained_Models_and_Totals/totals_by_scheme3.py` | vorpy.src.analyze.tools.compare.compare_files [compare_files2]; vorpy.src.analyze.tools.plot_templates.bar [bar] | csv, numpy, os, sys, tkinter |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_2kai_curvature_forensics.py` | vorpy.src.analyze.power_interface_curvature [cazals_signed_beta]; vorpy.src.analyze.power_interface_cli [run_power_interface_pdb] | __future__, argparse, collections, copy, csv, itertools, math, numpy, pathlib |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_2kai_independent_audit.py` | vorpy.src.chemistry.Cazals [OTHER_RADIUS, PROBE_RADIUS, element_radii, special_radii]; vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb [prepare_cazals_pdb]; vorpy.src.analyze.power_interface_curvature [analyze_power_interface_curvature] | __future__, argparse, collections, csv, gudhi, json, math, numpy, pathlib |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_alpha0_audit.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation [cazals_power_interface]; vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb [prepare_cazals_pdb] | __future__, argparse, collections, csv, dataclasses, itertools, numpy, pathlib |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_pdb.py` | vorpy.src.chemistry.Cazals [OTHER_RADIUS, PROBE_RADIUS, base_radius, special_radii] | __future__, csv, dataclasses, numpy, pathlib, typing |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_power_interface_diagnostic.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation [cazals_power_interface]; vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb [prepare_cazals_pdb] | __future__, argparse, csv, dataclasses, numpy, pathlib |
| `Part_II_Molecular_Analysis/SuplimentalInformation/cazals_power_interface.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb [export_audit_csv, export_npz, prepare_cazals_pdb] | __future__, argparse, collections, csv, dataclasses, gudhi, numpy, pathlib, typing |
| `Part_II_Molecular_Analysis/SuplimentalInformation/comparative_figures/F5_alpha_selection.py` | .plot_common [configure_matplotlib, load_static_comparison, method_label, ordered_systems, require_methods, save_figure, select_output_directory, static_protein_only] | __future__, matplotlib, numpy |
| `Part_II_Molecular_Analysis/SuplimentalInformation/comparative_figures/plot_common.py` | — | __future__, matplotlib, pandas, pathlib, tkinter |
| `Part_II_Molecular_Analysis/SuplimentalInformation/Curvature/CubeVertexGaussianCurvature.py` | — | matplotlib, mpl_toolkits, numpy |
| `Part_II_Molecular_Analysis/SuplimentalInformation/Curvature/DulaunayDual.py` | — | itertools, matplotlib, numpy, pathlib |
| `Part_II_Molecular_Analysis/SuplimentalInformation/Curvature/visualize_interface_criteria.py` | — | matplotlib, numpy, pathlib |
| `Part_III_ML_Processes/Surfaces_Training.py` | vorpy.src.analyze.tools [read_logs2]; vorpy.src.inputs [read_pdb_simple] | csv, numpy, os, pandas, sklearn, sys, tkinter |
| `power_interface_cli.py` | vorpy.src.analyze.cazals_validation [facet_area_statistics, interface_atom_asymmetry, interface_atom_ids, interface_components, pdb_inventory, power_facet_polygon_area]; vorpy.src.analyze.power_interface_curvature [analyze_power_interface_curvature]; vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_pdb [prepare_cazals_pdb] | __future__, csv, math, pathlib, tempfile |
| `power_interface_curvature.py` | vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation [cazals_power_interface] | __future__, collections, csv, dataclasses, itertools, math, numpy, pathlib |
| `power_interface_selection.py` | vorpy.src.analyze.cazals_validation [power_facet_polygon_area] | collections, dataclasses, itertools |
| `power_vs_aw.py` | vorpy.src.analyze.cazals_validation [interface_components] | __future__, collections, csv, numpy, pathlib |
| `tools/__init__.py` | vorpy.src.analyze.tools.compare [*]; vorpy.src.analyze.tools.plot_templates [*]; vorpy.src.analyze.tools.CleanData [*]; vorpy.src.analyze.tools.batch [*]; vorpy.src.analyze.tools.calcs [*] | — |
| `tools/batch/__init__.py` | vorpy.src.analyze.tools.batch.compile_logs [*]; vorpy.src.analyze.tools.batch.compile_new_logs [*]; vorpy.src.analyze.tools.batch.generate_batch [*]; vorpy.src.analyze.tools.batch.get_files [*]; vorpy.src.analyze.tools.batch.InDepthLogAnalysis [*]; vorpy.src.analyze.tools.batch.IsoMetricQuotient [*]; vorpy.src.analyze.tools.batch.wipe_user_data [*] | — |
| `tools/batch/compile_logs.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.calculations [calc_sphericity, calc_dist, get_time] | csv, numpy, os, time, tkinter |
| `tools/batch/compile_new_logs.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.calculations [calc_sphericity, calc_dist, get_time] | csv, numpy, os, time, tkinter |
| `tools/batch/generate_batch.py` | — | os |
| `tools/batch/get_files.py` | — | os, tkinter |
| `tools/batch/InDepthLogAnalysis.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.calculations [calc_sphericity, calc_dist, get_time] | csv, numpy, os, time, tkinter |
| `tools/batch/IsoMetricQuotient.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.calculations [calc_sphericity, calc_dist, get_time] | csv, numpy, os, time, tkinter |
| `tools/batch/make_batch_fixed_missing_atoms.py` | vorpy.src.analyze.tools.CleanData.check_mol_data [check_mol_data] | os, sys, tkinter |
| `tools/batch/make_foam_runs.py` | — | datetime, os, platform, tkinter |
| `tools/batch/wipe_user_data.py` | — | os, shutil, tkinter |
| `tools/calcs.py` | — | numpy |
| `tools/CleanData/__init__.py` | vorpy.src.analyze.tools.CleanData.AdjustAfterRuns [*]; vorpy.src.analyze.tools.CleanData.CleanFoamData [*]; vorpy.src.analyze.tools.CleanData.CheckFoamFolders [*]; vorpy.src.analyze.tools.CleanData.CombineLogs [*]; vorpy.src.analyze.tools.CleanData.CreateFinalFoamData [*]; vorpy.src.analyze.tools.CleanData.CreateVertsFromLogs [*]; vorpy.src.analyze.tools.CleanData.FinalCleaning [*]; vorpy.src.analyze.tools.CleanData.GetBasicFiles [*]; vorpy.src.analyze.tools.CleanData.MakeFoamRuns [*]; vorpy.src.analyze.tools.CleanData.RenameFiles [*]; vorpy.src.analyze.tools.CleanData.RenameFolders [*]; vorpy.src.analyze.tools.CleanData.VerifyFiles [*]; vorpy.src.analyze.tools.CleanData.VerifyFoamData [*] | — |
| `tools/CleanData/AdjustAfterRuns.py` | — | os, shutil, sys, tkinter |
| `tools/CleanData/check_mol_data.py` | vorpy.src.inputs.pdb [read_pdb]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.system.system [System] | os, pathlib, sys, tkinter, typing |
| `tools/CleanData/check_neighbors.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | ast, dataclasses, numpy, os, pandas, re, sys, tkinter, typing |
| `tools/CleanData/CheckFoamFolders.py` | vorpy.src.analyze.tools.batch.get_files [get_files]; vorpy.src.analyze.tools.CleanData.VerifyFiles [verify_pdb, verify_logs] | os, sys, tkinter |
| `tools/CleanData/CleanFoamData.py` | — | csv, os, sys, tkinter |
| `tools/CleanData/combine_surfs.py` | — | os, tkinter |
| `tools/CleanData/CombineLogs.py` | vorpy.src.calculations [calc_com, round_func, calc_total_inertia_tensor, combine_inertia_tensors] | csv, datetime, numpy, os, sys, tkinter |
| `tools/CleanData/CombineLogs2.py` | vorpy.src.calculations [calc_com, round_func, calc_total_inertia_tensor, combine_inertia_tensors] | collections, csv, datetime, numpy, os, re, sys, tkinter |
| `tools/CleanData/CreateFinalFoamData.py` | — | os, sys, tkinter |
| `tools/CleanData/CreateVertsFromLogs.py` | vorpy.src.output.net [write_verts]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, shutil, sys, tkinter |
| `tools/CleanData/ExtractLogs.py` | vorpy.src.analyze.tools.CleanData.CombineLogs2 [read_logs_sections, write_logs_sections, parse_string_lists] | datetime, os, sys, tkinter |
| `tools/CleanData/FinalCleaning.py` | vorpy.src.analyze.tools.CleanData.CreateVertsFromLogs [write_log_to_vert] | numpy, os, sys, tkinter |
| `tools/CleanData/Fix_Logs_SA.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2]; vorpy.src.system.system [System] | csv, numpy, os, pathlib, sys, tkinter |
| `tools/CleanData/MakeFoamRuns.py` | — | datetime, os, platform, sys, tkinter |
| `tools/CleanData/recolor_surfs.py` | — | tkinter |
| `tools/CleanData/RenameFiles.py` | — | os, shutil, sys, tkinter |
| `tools/CleanData/RenameFolders.py` | vorpy.src.analyze.tools.batch.get_files [get_files]; vorpy.src.analyze.tools.CleanData.VerifyFiles [verify_logs, verify_pdb] | datetime, os, shutil, tkinter |
| `tools/CleanData/verify_frames.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, sys, tkinter |
| `tools/CleanData/VerifyFiles.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os, sys, tkinter |
| `tools/CleanData/VerifyFoamData.py` | — | os, sys, tkinter |
| `tools/compare/__init__.py` | vorpy.src.analyze.tools.compare.color_by_residue [*]; vorpy.src.analyze.tools.compare.compare_files [*]; vorpy.src.analyze.tools.compare.compare_verts [*]; vorpy.src.analyze.tools.compare.get_res_data [*]; vorpy.src.analyze.tools.compare.in_out [*]; vorpy.src.analyze.tools.compare.pdb_names [*]; vorpy.src.analyze.tools.compare.read_logs [*]; vorpy.src.analyze.tools.compare.read_logs2 [*] | — |
| `tools/compare/color_by_residue.py` | — | os |
| `tools/compare/compare_files.py` | vorpy.src.system.system [System]; vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.compare.get_res_data [residue_data]; vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | os |
| `tools/compare/compare_verts.py` | — | csv, pandas |
| `tools/compare/get_res_data.py` | vorpy.src.analyze.tools.compare.read_logs [read_logs]; vorpy.src.analyze.tools.compare.pdb_names [proteins, nucleics, ions, other, sols] | csv |
| `tools/compare/read_logs.py` | vorpy.src.analyze.tools.compare.read_logs2 [read_logs2] | csv, os, pandas |
| `tools/compare/read_logs2.py` | vorpy.src.interface.geometry_analysis [parse_geometry_record] | ast, csv, os, pandas, pathlib, re |
| `tools/compare/sol_facing_percentage.py` | — | ast, numpy, pandas, typing |
| `tools/plot_templates/__init__.py` | vorpy.src.analyze.tools.plot_templates.bar [*]; vorpy.src.analyze.tools.plot_templates.line [*]; vorpy.src.analyze.tools.plot_templates.scatter [*]; vorpy.src.analyze.tools.plot_templates.blank_axes [*]; vorpy.src.analyze.tools.plot_templates.blank_legend [*]; vorpy.src.analyze.tools.plot_templates.box_whisker [*]; vorpy.src.analyze.tools.plot_templates.histogram [*]; vorpy.src.analyze.tools.plot_templates.inline_plot [*]; vorpy.src.analyze.tools.plot_templates.starting_point [*]; vorpy.src.analyze.tools.plot_templates.table [*] | — |
| `tools/plot_templates/bar.py` | — | matplotlib, numpy |
| `tools/plot_templates/blank_axes.py` | — | matplotlib, numpy |
| `tools/plot_templates/blank_legend.py` | — | matplotlib |
| `tools/plot_templates/box_whisker.py` | — | matplotlib |
| `tools/plot_templates/histogram.py` | — | matplotlib, numpy |
| `tools/plot_templates/inline_plot.py` | — | matplotlib, numpy |
| `tools/plot_templates/line.py` | — | matplotlib, numpy |
| `tools/plot_templates/scatter.py` | — | matplotlib, numpy |
| `tools/plot_templates/starting_point.py` | vorpy.src.system.system [System]; vorpy.src.group.group [Group] | os, tkinter |
| `tools/plot_templates/table.py` | — | matplotlib |
