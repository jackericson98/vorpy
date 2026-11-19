# Figure Organization

Status: F3/F4 source consolidation complete; cleanup validated; no commit made.  
Last updated: 2026-11-16  
Authority: this document supersedes the figure-reorganization master plan, code map, Figure 3 planning reports, and F3/F4 dry-run manifest.

## Approved registry

| Collection | Scope |
|---|---|
| F1 | Methodology and molecular models |
| F2 | Frame persistence |
| F3 | Pure water and PEG; F3A remains reserved because no pure-water generator or result exists |
| F4 | Heparin |
| F5 | DNA and Hairpin only; NCP is not included automatically |
| F6 | Proteins, atomistic residue analyses, interfaces, and multiscale comparisons |
| F7 | Explicit coarse-grained representations and interfaces; tentative |
| Curvature | Separate research collection |
| Supplemental | Supporting and exploratory analyses, including Atomic UMAP |
| Archive | Historical or retired work |
| Shared_Utilities | One canonical implementation for reusable infrastructure |

AtomicCorrelation remains REVIEW. The complete monolithic legacy clustering implementation is the canonical executable pipeline; the incomplete split refactor remains preserved for historical comparison and shared helpers. Part I statistical-ensemble analyses are unchanged.

## Canonical F3/F4 locations

- F3 generators: `vorpy/src/analyze/Part_II_Molecular_Analysis/F3_Pure_Water_and_PEG/`
- F3 screened-data helper: `F3_Pure_Water_and_PEG/support/F5_peg_dataset.py`
- F3 preparation tools: `scripts/figure_reorganization/F3/`
- F3 tests: `F3_Pure_Water_and_PEG/tests/` and `vorpy/tests/analyze/F3/`
- F4 generator: `vorpy/src/analyze/Part_II_Molecular_Analysis/F4_Heparin/`
- Shared simple-model code: `vorpy/src/analyze/Part_II_Molecular_Analysis/Shared_Utilities/simple_models/`

## F3/F4 migration map

The SHA-256 column is the repaired pre-migration source/destination content hash and rollback fingerprint.

| Former source | Canonical destination | Figure | SHA-256 |
|---|---|---|---|
| `P/F3_SimpleModels/F3_PEG_prototypes.py` | `P/F3_Pure_Water_and_PEG/F3_PEG_prototypes.py` | F3 | `b680d37b05cc46d8e4d62fa2acb5f9037023a51f6f0db6d829ea78e15fc4c683` |
| `P/F3_SimpleModels/F3_round2_validation.py` | `P/F3_Pure_Water_and_PEG/F3_round2_validation.py` | F3 | `2839e4b324e02c1fb3fb2a7179498ac4257c253090e0656e21e64532b56ccaca` |
| `P/F3_SimpleModels/F3F_water_contact_prototype.py` | `P/F3_Pure_Water_and_PEG/F3F_water_contact_prototype.py` | F3 | `d634dfa8c2cbe9d18f80e32ecba4ddcb120b5f90e6adb5ebb0cc2734cfdb4a12` |
| `P/F3_SimpleModels/F5_peg_dataset.py` | `P/F3_Pure_Water_and_PEG/support/F5_peg_dataset.py` | F3 | `9ed71613cc8686c8c30a0901581258e477988b4be13fd8e67bdc7aa7bf73f28f` |
| `P/F3_SimpleModels/F5A_exposed_vs_buried.py` | `P/F3_Pure_Water_and_PEG/F5A_exposed_vs_buried.py` | F3 | `ffc1fe95f5337e6bd2ff7976989e4fc10c13e3177b7245617f1e9f7eb065c1ba` |
| `P/F3_SimpleModels/F5A_PEG_chemistry_vs_hydration.py` | `P/F3_Pure_Water_and_PEG/F5A_PEG_chemistry_vs_hydration.py` | F3 | `be083b907300042095c015a976667db43fa750b56bfc0d3e9f0ec0b3b5ba238b` |
| `P/F3_SimpleModels/F5B_heparin_chemistry.py` | `P/F4_Heparin/F5B_heparin_chemistry.py` | F4 | `c3be6122bdea0c29ab2a21da5ec99d498c814310e2f03411e70fbaac858db275` |
| `P/F3_SimpleModels/F5C_rerun_missing.py` | `P/F3_Pure_Water_and_PEG/F5C_rerun_missing.py` | F3 | `fcb492bb32baeb1767bd3a2aa7325d6dfc7bab1338aa2d6dd00ff4986ac4b23e` |
| `P/F3_SimpleModels/F5C_water_contacts.py` | `P/F3_Pure_Water_and_PEG/F5C_water_contacts.py` | F3 | `5d9ac685e923e55a047f4ac477bb085c0bbc421708e331d1220528dc1901ac0c` |
| `P/F3_SimpleModels/tests/test_F3_PEG_prototypes.py` | `P/F3_Pure_Water_and_PEG/tests/test_F3_PEG_prototypes.py` | F3 | `16b5294579c97eeef1171f8bbab2b03de7d143dee6206babe8f47478024cc4a4` |
| `P/F3_SimpleModels/tests/test_F3_round2.py` | `P/F3_Pure_Water_and_PEG/tests/test_F3_round2.py` | F3 | `5b5a5d3d40acca1824b036e8439b9f7d73f16168196f23cb36b0d6da0afd6cc4` |
| `scripts/audit_peg_boundary_patches.py` | `scripts/figure_reorganization/F3/audit_peg_boundary_patches.py` | F3 | `fba53fa9da59a79fc0ce835e9c10069d67809c09639e02bbae97200c440c8109` |
| `scripts/investigate_peg_first_layer.py` | `scripts/figure_reorganization/F3/investigate_peg_first_layer.py` | F3 | `3fb4a1234aafb1bfe37fec96de525e6e6bf4411070de6ad3984e35ab80b9a140` |
| `scripts/investigate_solvent_layers.py` | `scripts/figure_reorganization/F3/investigate_solvent_layers.py` | F3 | `8b67bf699efe2d585055362569bdb23a62d75e18f443b37767812e60eae11caa` |
| `scripts/probe_peg_missing_water.py` | `scripts/figure_reorganization/F3/probe_peg_missing_water.py` | F3 | `4a178d2fc23e7afad3e997660163b68501866f72283e89025c358fd70a438c12` |
| `scripts/reexport_peg_layer1.py` | `scripts/figure_reorganization/F3/reexport_peg_layer1.py` | F3 | `d17238d4b7de28882f52d7f42930810fd60d54e25ed8ac17d601fda062fb2509` |
| `scripts/summarize_peg_layer1_investigation.py` | `scripts/figure_reorganization/F3/summarize_peg_layer1_investigation.py` | F3 | `a1641cd298c7661b2db68fb87d5a160830216f012b542350af119e726f2a72c4` |
| `scripts/validate_peg_layer1_bundle.py` | `scripts/figure_reorganization/F3/validate_peg_layer1_bundle.py` | F3 | `0764a65e179f038527bb0c1411c02c501c56f0e46f46eb5583eabd2aa28031b1` |
| `vorpy/tests/command/test_peg_recovery_diagnostics.py` | `vorpy/tests/analyze/F3/test_peg_recovery_diagnostics.py` | F3 | `36656e22cfe90b6b53d360ec06db505a0ca2e16ee24a3b7bad7eea1fd8080236` |

Former analysis and preparation-script paths are compatibility aliases to these destinations. They contain no scientific implementation. The three duplicated relocated tests were removed from their former locations.

## Shared prerequisites and gates

| Module | Canonical status | Rollback SHA-256 |
|---|---|---|
| `F5_common.py` | Canonical implementation copied to Shared_Utilities; former path is an alias | `7932a20b16437b9c558cf08ba8d7a371705d26f731b6a73918e8a7af090baeed` |
| `F5_config.py` | Canonical implementation copied to Shared_Utilities; former path is an alias | `9867f4deb9b84a600021870ea42e534b36ac7203e110be8376acb7064461dbfd` |
| `F5_run_all.py` | Canonical implementation under Shared_Utilities; former path is an executable compatibility alias; canonical SHA-256 `79abb4834e0fd1a20be85ef1034da02bf7510fa6654f8073c729e6a73b632d8e` | `dd72d3ecc143f06243c355fde5ea5ca3a6363672ba358f741b5947bda69cf500` |
| `F5B_water_burial.py` | Canonical implementation under Shared_Utilities; former path is an executable compatibility alias; canonical SHA-256 `0f37d485c7e94e51974d764769a9b717ba2708cd2c63dc7ae6b790121d126db1` | `9e5c077532d8d332d564e5f3f00c26ee05215827a461c9fa2d5218ede00fbb22` |

No F3/F4-local copy of shared water geometry exists. `F5B_water_burial.py` remains the sole implementation. Its only relocation repair is a lazy import of the canonical F3 screened-dataset helper; calculations and defaults are unchanged.

## Compatibility aliases

All 20 former analysis and preparation entry points remain as thin aliases. They are retained because they protect historical module paths, direct-file/PyCharm commands, archive-script imports, and external workflows that cannot be exhaustively inventoried. Each alias bootstraps the repository path, imports the canonical module, and invokes canonical `main()` only when explicitly executed. Static inspection confirms that none defines an independent function or class implementation.

| Alias family | Verified consumer or compatibility need | Decision |
|---|---|---|
| Former `F5_common`, `F5_config`, `F5_run_all`, and `F5B_water_burial` paths | Historical F5 package imports, archive panel entry points, and documented direct-file commands | Retain |
| Former F3/F4 analysis paths | Fallback imports in relocated modules/tests plus historical module and PyCharm paths | Retain |
| Former root preparation-script paths | Fallback imports in the canonical preparation tools and external/importlib command paths | Retain |

No alias without a verified repository or historical/external compatibility purpose was found.

## Water-coverage discrepancy

The expectation was not stale. A custom PEG dataset was being ignored: the PEG branch always loaded the global screened registry, so a synthetic polymer-only dataset incorrectly reported 27 valid waters. The implementation now uses the screened registry only when `config.dataset_root` is the canonical organized dataset. Custom roots use their own discovered logs and correctly report `water_geometry_unavailable` when no water rows exist. Default screened-cohort calculations, tolerances, classifications, and normalization are unchanged.

## Validation evidence

- Pre-fix focused test reproduced two `ok` rows with 27 unrelated waters; post-fix test passed.
- Final focused suite: **50 passed**, one existing pandas `FutureWarning`.
- Import smoke: all 9 relocated analysis/support modules and 7 preparation modules imported without execution.
- Shared identity: former-path `read_pdb`, `Config`, runner `run`, and `water_geometry` are the same Python objects as the canonical shared exports.
- Copy integrity: all 19 migrated destinations matched their repaired source SHA-256 before aliases replaced former implementations; both shared canonical copies matched their recorded source hashes.
- Screened cohort: 207 scheme rows, 69 distinct model/frame sets, with recovered `PEG3000/frame_0005` AW bundle retained.
- Screened registry SHA-256: `5aeeba35121060c588e5102178023b2fdf581c2e64ec812db6b727491f81f770`.
- Tests wrote only to pytest temporary directories. No scientific generator, simulation, or production geometry ran.

The legacy direct-execution test now uses 13 isolated `runpy` import/symbol checks plus a runner `--help` check. All 14 pass without invoking `main()`, loading scientific datasets, or writing figures. This replaces the six generator-launching cases that previously timed out and the support-module execution case.

## Preserved output evidence

No migration command targeted an output directory. Representative validated files retained their existing modification times and completion hashes:

| Artifact | SHA-256 |
|---|---|
| `output/figure3_prototypes/F3_prototype_coverage.csv` | `dbffa7c3c16f550a281acbe445690229e276ee033c63fb0bbfedcde03987faf3` |
| `output/figure3_round2_validation/F3_round2_coverage.csv` | `dbffa7c3c16f550a281acbe445690229e276ee033c63fb0bbfedcde03987faf3` |
| `output/figure3_round2_water_contacts/F3F_coverage.csv` | `7f3953139215a35d6f645116d5ad72a2215bc03b548196193c33810f7b9f49e5` |
| legacy PEG F5A coverage | `3c19bc580ea4981a83a5a192d5de3d5c2bb8d820febb9326464c9c2d84236a58` |
| legacy PEG F5C coverage | `6b5bc78d66d41203ef83db9eb35d8fe881f75445c790d9f6860d0f726c4556bb` |

The dry run did not record a complete pre-copy output-hash baseline, so byte-for-byte before/after proof for every historical output is unavailable. The hashes above are the authoritative completion checkpoint; no output file was written by this migration.

## Cleanup and rollback

Removed as superseded:

- former copies of the two relocated F3 unit-test modules;
- former copy of `vorpy/tests/command/test_peg_recovery_diagnostics.py`;
- temporary master plan, Figure 3 planning/validation reports, code-map files, dry-run report, and migration manifest.

Retained intentionally:

- all scientific outputs, datasets, caches, screened-status records, and unresolved research material;
- all former analysis and preparation paths as compatibility aliases, including executable aliases for `F5_run_all.py` and `F5B_water_burial.py`;
- all unrelated dirty working-tree content.

Rollback evidence consists of the mapping hashes above, compatibility aliases, unchanged data/output locations, and Git history for tracked files. Alias removal requires an explicit external-compatibility decision; passing repository tests alone is insufficient.

## Remaining decisions

- F3A pure-water baseline remains planned, not implemented.
- Final F3–F7 panel assignments remain unresolved.
- External consumers of historical aliases are not fully enumerable, so those aliases remain intentionally.

## F6 — proteins, interfaces, and multiscale analysis

Status: canonical source migration complete for validated modules; no figures, simulations, or scientific results were generated.

### Canonical locations

- Core multiscale collection: `vorpy/src/analyze/Part_II_Molecular_Analysis/F6_Proteins_Interfaces_Multiscale/`
- Historical protein/interface/residue collection: `F6_Proteins_Interfaces_Multiscale/historical/`
- Figure3B recovery/export preparation: `scripts/figure_reorganization/F6/`
- Focused tests: `vorpy/tests/analyze/F6/`

Core canonical modules are `F3_common.py`, `F3A_multiscale_summary.py`, `F3B_atomic_groups.py`, `F3C_residue_groups.py`, `F3D_geometry_examples.py`, `F3E_solvent_exposure.py`, and `F3_build_figure.py`. The F3D and build modules remain explicit empty/draft slots; migration did not add analysis logic.

The historical collection contains 19 validated-import modules covering ligand interfaces, NCP protein–DNA interfaces, p53 chain comparisons, molecular totals, atomistic residue/type comparisons, and atomic deviations. NCP remains F6 historical work rather than F5 because these scripts combine protein and DNA.

The three canonical preparation tools are `plan_figure3b_recovery.py`, `validate_figure3b_recovery.py`, and `smoke_figure3b_pymol.py`. Former source and script locations are retained as 29 thin compatibility aliases with no function or class implementations.

Canonical source fingerprints at this checkpoint:

| Tree | Files | Bytes | SHA-256 |
|---|---:|---:|---|
| F6 analysis collection | 28 | 299,159 | `68d6a094ff0d5331db7df610887a7e5774bc32b6f94edd47ed7ab6dd0ba607bc` |
| F6 preparation tools | 4 | 19,288 | `869e441e0a640607ecdb2b69adb755f4360d8d527b2ba64d2a4b2d6e2d1fa2a2` |

### Dependencies and retained code

- `F3B_atomic_groups.py` and `F3E_solvent_exposure.py` still use the established `F2C_All_Data_Sets_Clustered.py` atom-name canonicalization. That legacy implementation and its helper remain in `F4_atomic_deviations/F3_old/`; only their deleted-package import was repaired.
- The complete `F4_atomic_deviations/F3_old/F2C_All_Data_Sets_Clustered.py` implementation is the canonical executable legacy clustering pipeline. The split `F2C_All_Data_Sets_Clustered_full` refactor remains intact as historical code and supplies its configuration, folder resolution, naming, PCA diagnostic, and plotting helpers. `run_single.py` and `run_batch.py` now import the complete monolith's existing `main` instead of the nonexistent `pipeline.main`; neither launcher was executed.
- Atomic UMAP remains Supplemental/exploratory and curvature PCA remains Curvature research. Their imports now use canonical F6 shared metrics, but their files were not moved in this F6 migration.
- Production `nonpolar_interface` and interface geometry/export APIs remain stable shared infrastructure outside a figure collection.
- `Extra/` and `Other/` were retained for review because import probing reached interactive top-level behavior. They were not copied or executed further.
- Historical `F2A_AtomicVolumeDifferences.py` retains explicit `E:/...` user data/output defaults. They are documented legacy settings, not safe paths to infer or rewrite.

### Legacy clustering implementation comparison

The implementations are not scientifically equivalent. They share the same configured data root, molecular-system groups, 20 batch settings, AW/Power volume-range filter, and many plotting/statistics helpers, but their executable scope and atom classification differ:

| Area | Complete monolith | Split refactor |
|---|---|---|
| Inputs and features | Reads AW/Power atom logs and surface logs; joins PDB identity and chemistry/bond context; derives AW/Power solvent-facing values, `DeltaV = Pow - AW`, coordinates, atom/residue names, and optional solvent binary | Joins AW and Power atom logs by index and carries a pre-existing `SolFacingPct`; has no complete orchestration |
| Classification/filtering | Expanded small-molecule aliases and environmental refinement; no forced protein groups | Smaller/different small-molecule aliases, different RNA O2/O4/O6 aliases, and 14 forced protein groups |
| Clustering | Calls the established helper for KMeans, agglomerative, DBSCAN, HDBSCAN, and K-Prototypes | Contains no clustering call or `main`; adds only standalone standardized PCA and solvent-facing summaries |
| Normalization | Numeric features use `StandardScaler`; categoricals are one-hot encoded except native categorical columns for K-Prototypes; missing selected features are dropped | PCA standardizes its selected numeric columns, but no end-to-end clustering normalization exists |
| Outputs | Writes ML assignments, cluster summaries, manual-versus-ML tables, group statistics/keys/points, and optional PNG/SVG | Plot helpers can save requested images/tables, but there is no runnable end-to-end output contract |

Consequently, routing the launchers to the monolith repairs importability without claiming that the split classifications reproduce monolithic results. Both historical source trees are retained, and no scientific output was regenerated.

The split batch configuration still embeds `E:\\Molecular` as its data root and an `E:\\OneDrive - Georgia State University\\...\\2C_Full_Plots_Clustered` output root. Direct execution therefore remains workstation-dependent and, by design, can create directories and write plots/tables. Validation imported the launchers only; it did not call either entry point or access those roots.

### F6 validation

- Canonical imports: 7 core and 19 historical modules passed.
- Focused pytest: 13 passed, 1 skipped. The skip is the existing optional `umap` dependency, which was not installed or changed.
- Figure3B nondestructive recovery regression: the bounded missing-vertex calculation passed. Geometry-building recovery tests were not run because this migration forbids simulations and regenerated outputs.
- All 32 canonical F6 Python files parsed. The two non-PyMOL preparation modules imported; the PyMOL tool was syntax-checked only because it requires the external PyMOL runtime.
- No scoped reference to the deleted `F3_Multiscale_Scheme_Differences` package remains.
- Legacy clustering resolution: 10 focused tests passed. The two launchers and both implementations imported without executing `main`; assertions confirmed both launchers resolve to the monolithic orchestrator, the 20 batch settings remain identical, and the known classification differences remain explicit.

### Preserved F6 scientific assets

The before/after fingerprints remained identical:

| Asset | Files | Bytes | SHA-256 |
|---|---:|---:|---|
| `output/Figure3B_revalidation` | 1,418 | 611,838,077 | `f354b9254ba597a4ba622f145bf316a48cf375a54699e5f7d37f5791b9a9dfb8` |
| Legacy non-Python `F4_atomic_deviations` asset tree | 1 | 1,573 | `712ec32f9bef982e2cde39053f6c3888240d7479ea3e1ec91a9294193dbf4032` |
| `vorpy/data/hammerhead.pdb` | 1 | 4,336,171 | `e9cec891c51aa240d3c3a9b66d48cba3d3bf83794c461ccabd97b709ef0e2cd1` |
| `vorpy/data/EDTA.pdb` | 1 | 50,189 | `3b233f7e98a1253f1506d562c54e192dcbfdf7a99a2efcbebe6aeccd305f888f` |
| `vorpy/data/T4LP.pdb` | 1 | 2,680,939 | `26dca4fe186d336139f08ca607ac075c0a2637bc6eec4b4a68943f445ca54e74` |

### F6 ownership decisions still required

- Decide whether a future migration should move the canonical monolith into F6 and leave a compatibility alias. It remains at its historical path now to avoid breaking external consumers and scientific provenance.
- Decide whether the partial split refactor should eventually be archived or reconciled. Its differing atom aliases and forced groups require scientific approval before any merge; it is not safe to treat as a drop-in equivalent.
- Assign `Extra/` and `Other/` scripts individually after removing or isolating interactive top-level execution.
- Relocate Atomic UMAP to Supplemental and curvature PCA to Curvature only in their respective approved migrations.
- Decide whether historical hardcoded external roots should become explicit command-line inputs; changing them was outside this organizational migration.

