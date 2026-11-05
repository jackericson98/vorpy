# VorPy Workbench Functionality Blueprint

**Status:** Product and implementation blueprint  
**Date:** 2026-09-04  
**Target:** Production-oriented PySide6/PyVistaQt application developed separately from VorPy until the integration contract is stable

## 1. Product definition

VorPy Workbench should be a project-based scientific application for building,
running, inspecting, comparing, and presenting molecular geometry analyses. It
should connect molecular structure, user-defined groups, VorPy results,
quantitative plots, and publication figures without attempting to replace the
entire feature set of PyMOL or VMD.

The central concept is:

> A project contains structures, reusable scientific selections, individual
> VorPy runs, comparisons between runs, and reproducible figure recipes.

The 3D viewer is the shared visual context for all of those objects. It should
not become the storage location for scientific state. Selections, groups,
interfaces, run settings, properties, and figures must exist as explicit domain
objects that can be saved, restored, tested, and reused.

### How to use this blueprint with an evolving GUI

This document defines **product goals, architectural boundaries, and desired
capabilities**. It is not a frozen description of the current screen layout and
must not be used to replace newer working code with an older mockup or earlier
implementation.

At the beginning of every implementation phase, Codex must:

1. inspect the current worktree, recent commits, `AGENTS.md`, tests, and current
   UI modules before proposing changes;
2. inventory improvements already present, including behavior, layout,
   appearance, performance, and tests;
3. map those improvements to the goals in this blueprint;
4. treat the current working implementation as the visual and behavioral
   baseline unless the user explicitly requests a redesign;
5. propose the smallest incremental change that advances the highest-priority
   unmet goal;
6. preserve unrelated improvements and user changes; and
7. add or update regression coverage before refactoring working behavior.

When current code and this document differ, preserve the current improvement
and adapt the blueprint's architecture around it. Ask the user before removing,
replacing, or substantially redesigning a working feature. Architectural
cleanup is not sufficient justification for a visible regression.

The implementation priority is therefore:

```text
preserve current improvements
→ identify the highest-value unmet goal
→ extend incrementally
→ verify behavior and appearance
→ update the checkpoint documentation
```

## 2. Confirmed requirements

### Workflows

- Load and inspect molecular structures.
- Configure and run VorPy.
- Load and explore completed VorPy results.
- Compare schemes or separate runs.
- Analyze trajectories across frames.
- Create publication figures.
- Create persistent groups and interfaces.
- Run one VorPy scheme at a time.
- Register previously completed runs in a project.

### First-class selectable objects

- Atoms
- Residues
- Chains
- Molecules
- Waters
- Ions
- Custom groups
- Interfaces between groups
- Surfaces
- Edges
- Vertices

### Group construction

- Manual atom and residue selection
- Rule-based selection
- Set union, intersection, subtraction, and inversion
- Reuse across trajectory frames

### Visualization and analysis

- Molecular representations
- VorPy geometry layers
- Property-based coloring
- Linked plots and data tables
- Side-by-side comparison
- Geometry: volume and surface area
- Topology: contacts and neighbors
- Mean and Gaussian curvature
- Integrated curvature and representative surface energy
- Interfaces and solvent exposure
- Build quality and completeness

## 3. Information architecture

The application should have five top-level workspaces. These are workflow
modes, not separate applications.

| Workspace | Purpose | Primary content |
|---|---|---|
| Structure | Load, inspect, and organize molecular inputs | 3D molecule, hierarchy, representations, frame controls |
| Build | Define groups/interfaces and configure a VorPy run | selection builder, group editor, scheme/settings, job queue |
| Results | Explore one completed run | geometry layers, property coloring, plots, tables, diagnostics |
| Compare | Compare two or more registered runs | linked viewers, difference maps, matched plots and tables |
| Figures | Assemble reproducible scientific graphics | viewport captures, plots, layout, labels, export settings |

The application should open in **Structure**. A loaded or completed run should
remain available when the user moves between workspaces.

## 4. Persistent application shell

The Analysis Studio shell should remain stable across workspaces:

### Left tool rail

- Project/home
- Structure
- Build
- Results
- Compare
- Figures
- Jobs
- Settings

The rail changes workspace. Atom/residue selection and camera tools belong in a
viewer toolbar, not in the global rail.

### Center stage

- One main viewer in Structure, Build, and Results
- Two synchronized viewers or a split-overlay mode in Compare
- Figure canvas in Figures
- Viewer-local toolbar along the top
- Visibility/representation strip below the viewport
- Optional frame strip when a trajectory is loaded

### Right inspector

Context-sensitive tabs:

- **Inspect:** metadata for the active atom, residue, chain, group, interface,
  surface, edge, or vertex
- **Display:** representation, visibility, color, opacity, clipping, and quality
- **Data:** property values, units, validity, provenance, and related records
- **Selection:** current selection expression and actions to create/update groups

The inspector should react to the active object. It should not duplicate an
entire settings page for every workspace.

### Bottom tray

Resizable and collapsible tabs:

- Summary
- Plots
- Data table
- Contacts
- Jobs/log
- Diagnostics

The bottom tray should never permanently consume most of the viewport. Its
height and collapsed state should be saved per project or user session.

## 5. Project tree

The project tree is the durable map of scientific objects:

```text
Project
├── Structures
│   ├── EDTA.pdb
│   └── trajectory.xtc + topology.pdb
├── Selections
│   ├── Active site
│   └── Solvent shell
├── Groups
│   ├── Protein
│   ├── Ligand
│   └── Waters within 4 Å
├── Interfaces
│   └── Protein ↔ Ligand
├── Runs
│   ├── EDTA — Atomic Voronoi
│   ├── EDTA — Power
│   └── Imported Primitive result
├── Comparisons
│   └── Atomic Voronoi vs Power
└── Figures
    └── Interface surface comparison
```

Tree rules:

- Clicking selects an object and updates viewer plus inspector.
- Checkboxes control visibility only; they do not change selection.
- Double-click opens the appropriate workspace/editor.
- Context menus expose rename, duplicate, recompute, reveal source, export, and
  remove-from-project actions.
- Deleting a referenced object requires dependency-aware confirmation.
- Imported files remain external by default; the project stores canonical paths
  and optionally supports copying them into a managed project directory.

## 6. Core domain model

Qt widgets must consume these models through controllers/services. Scientific
state must not live only in widgets or VTK actors.

### Project

- `id`
- `name`
- `project_file`
- `created_at`, `modified_at`
- collections of structures, selections, groups, interfaces, runs,
  comparisons, and figures
- UI state references, never raw widget state

### StructureSource

- stable ID and display name
- source paths
- format and topology source
- atoms, residues, chains, molecules
- optional trajectory metadata
- coordinate units
- structure fingerprint

### FrameSet

- topology/structure reference
- frame source
- frame count and time values
- current frame index
- lazy coordinate access
- atom-identity mapping status

### SelectionExpression

A serializable expression tree rather than a Python callback. Supported nodes:

- identity: atom IDs, residue IDs, chain IDs
- attributes: element, atom name, residue name, record type
- categories: protein, nucleic acid, water, ion, ligand, solvent
- spatial: within distance of another expression
- set operators: union, intersection, subtraction, inversion

### GroupDefinition

- ID, name, color, description
- `SelectionExpression`
- structure reference
- frame policy
- cached membership per relevant frame
- manually included/excluded identities

### InterfaceDefinition

- ID and name
- group A reference
- group B reference
- optional surrounding-solvent policy
- contact/surface inclusion rules
- validation preventing accidental self-interface unless explicitly allowed

### RunConfiguration

- source structure/frame set
- exactly one scheme: Atomic Voronoi, Power, or Primitive
- selected groups and interfaces
- solver/build settings
- analysis/property settings
- requested export artifacts
- output directory policy
- reproducible serialized command/configuration

### RunRecord

- ID, name, state, timestamps
- configuration snapshot
- VorPy version and environment metadata
- input and output paths
- result manifest reference
- progress/events
- warnings and errors
- imported-versus-executed origin

### ResultDataset

- owning run
- entity tables: atoms, residues, groups, interfaces, surfaces, edges, vertices
- property registry
- geometry assets
- provenance and units
- completeness/validation report

### ComparisonDefinition

- two or more compatible run IDs
- entity matching strategy
- selected properties
- reference run
- display layout and synchronized camera state

### FigureRecipe

- source run/comparison IDs
- selected entities and properties
- representation and color configuration
- camera state
- plot definition
- annotations
- canvas dimensions, DPI, and export format

## 7. Stable molecular identity

Trajectory reuse and comparisons require stable identity independent of VTK
point indices. Define a canonical `AtomKey` using the best available tuple:

```text
model / chain / residue sequence / insertion code / residue name /
atom name / alternate location
```

Use serial number only as supporting provenance because serials may change
between converted files. Every load should produce an identity report:

- unique atom keys
- duplicate/ambiguous keys
- missing identity fields
- trajectory-to-topology match percentage
- run-to-run entity match percentage

Groups reused across frames resolve through `AtomKey`, not coordinate or array
position.

## 8. Selection system

Selection is a service with one active mode and one active selection set.

### Modes

- Navigate
- Atom
- Residue
- Chain
- Molecule
- Surface
- Edge
- Vertex

### Required behavior

- Clicking replaces the current temporary selection by default.
- Ctrl-click toggles membership.
- Shift-click adds a range or connected selection where meaningful.
- Camera drags never change selection.
- Hidden categories cannot be picked.
- Selection highlighting never resets the camera.
- Selection survives workspace changes and property/color changes.
- Clear selection is explicit.
- “Create group from selection” converts temporary selection into a persistent
  `GroupDefinition`.

### Selection builder

The GUI should allow both structured controls and readable query text:

```text
protein and chain A and resid 20:45
water within 4.0 of group "Ligand"
(resname ASP or resname GLU) and exposed
```

Internally, parse this into the expression tree. Do not store raw executable
Python or evaluate arbitrary expressions.

## 9. Group and interface workflow

### Create group

1. Select atoms/residues manually or open the rule builder.
2. Preview membership and highlight it in the viewer.
3. Display atom/residue counts and unresolved identities.
4. Name and color the group.
5. Save the expression plus manual overrides.
6. For a trajectory, validate membership on representative frames.

### Combine groups

- Union: A ∪ B
- Intersection: A ∩ B
- Subtraction: A − B
- Inversion relative to the containing structure
- Derived groups retain references to parents, allowing optional live updates
- “Freeze membership” creates an identity snapshot

### Create interface

1. Choose two groups.
2. Preview spatial/contact relationship.
3. Choose whether interface results include touching solvent, interface atoms,
   interface surfaces, or whole contacting residues.
4. Validate overlap and empty-group conditions.
5. Save as a persistent `InterfaceDefinition`.

## 10. VorPy run workflow

Each run uses one scheme. Comparison is performed between separate run records.

### Run editor sections

1. **Input:** structure and optional frame(s)
2. **Scheme:** Atomic Voronoi, Power, or Primitive
3. **Targets:** whole system, groups, interfaces
4. **Build settings:** maximum vertices and scheme-specific options
5. **Analysis:** requested property families
6. **Outputs:** geometry format and artifact preset
7. **Destination:** project-managed or chosen output directory
8. **Review:** exact reproducible command/configuration

### Job states

```text
Draft → Queued → Preparing → Running → Exporting → Importing → Complete
                                      ↘ Failed
                 ↘ Cancelled
```

Jobs must be durable enough that reopening the project can distinguish:

- completed result
- interrupted local process
- missing output directory
- partially exported result
- result created by another VorPy version

## 11. VorPy integration contract

Start with a subprocess adapter rather than importing the VorPy GUI or solver
objects into the Qt process.

### Workbench-to-VorPy contract

Prefer a structured configuration file over constructing a long ambiguous CLI
command:

```bash
python -m vorpy --config /path/to/run_config.json --events-jsonl
```

Until VorPy supports this, the adapter may translate `RunConfiguration` into
the current CLI, but command construction must be isolated and unit tested.

### Progress protocol

VorPy should emit newline-delimited JSON events to stdout or a dedicated event
stream:

```json
{"event":"progress","stage":"build","label":"Constructing surfaces","current":42,"total":100}
{"event":"warning","code":"INCOMPLETE_CELL","message":"3 cells are incomplete"}
{"event":"artifact","kind":"surfaces","path":"EDTA/surfs.ply"}
{"event":"complete","manifest":"result_manifest.json"}
```

Human-readable logs can still be written separately. The GUI must not parse
terminal progress bars or depend on fixed prose.

### Result manifest

Every new run should produce `result_manifest.json` containing:

- schema version
- VorPy version and scheme
- input fingerprint
- complete run configuration
- timestamps and duration
- entity counts
- property table files and column metadata
- geometry artifacts with object scope and format
- group/interface definitions used by the run
- units
- warnings, failures, and completeness state

The Workbench should include a legacy importer for existing result directories
that infers artifacts from filenames. New runs should use the manifest as the
authoritative source.

## 12. Visualization architecture

Separate scientific data from rendering through these services:

- `SceneController`: viewer state and camera coordination
- `RepresentationRegistry`: available molecule/geometry representations
- `LayerController`: actor creation, visibility, opacity, and lifecycle
- `SelectionController`: pick modes, selection state, and highlights
- `PropertyRegistry`: property definitions, scopes, units, ranges, and palettes
- `ColorMappingService`: converts a property into per-entity display scalars
- `PlotController`: creates linked plots from property tables
- `ComparisonController`: entity matching, deltas, synchronized cameras
- `FigureController`: reproducible captures and layouts

Do not let each panel create its own VTK actors directly.

## 13. Representation registry

### Molecular representations

Initial production set:

- Ball and stick
- Space filling
- Sticks
- Points for very large systems
- Water and ion categories
- Optional backbone trace

Future/optional:

- Ribbon/cartoon representation, only if a reliable library-backed
  implementation is selected
- Molecular surface independent of VorPy geometry

### VorPy representations

- Surfaces
- Shell surfaces
- Separate per-entity surfaces
- Edges
- Shell edges
- Separate edges
- Vertices
- Shell vertices
- Separate vertices
- Atom-associated geometry
- Group geometry
- Interface geometry

Each registered representation declares:

- compatible entity scope
- required artifact type
- default style
- default opacity and color
- picking support
- large-data strategy

## 14. Property registry

Properties must be metadata-driven, not hard-coded into every widget.

Each property definition includes:

- canonical key
- display name
- entity scope
- unit
- description
- value type
- valid/missing-value rules
- recommended palette type
- default range policy
- whether values are local, integrated, normalized, or categorical

Initial property families:

| Family | Example properties | Typical scope |
|---|---|---|
| Geometry | volume, surface area | atom, residue, group, molecule |
| Topology | contact count, neighbor count, coordination | atom, residue, interface |
| Curvature | mean curvature, Gaussian curvature | triangle, surface |
| Integrated geometry | `int_mean_curv`, `int_mean_curv_sq`, `int_gauss_curv` | surface |
| Representative energy | `surf_energy` | surface |
| Interface | interface area, contact area, contacting entities | interface, residue |
| Solvent | solvent exposure, solvent-contact area | atom, residue, group |
| Quality | complete/incomplete, degeneracy, validation flags | cell, surface, run |

Preserve the canonical VorPy keys `int_mean_curv`, `int_mean_curv_sq`,
`int_gauss_curv`, and `surf_energy`. The GUI may show friendly names.
`surf_energy` must be described as representative/reference
curvature-dependent bending energy, not a calibrated molecular free energy.

### Coloring rules

- Sequential palettes for nonnegative magnitude.
- Diverging palettes only when zero has scientific meaning.
- Categorical palettes for status or entity class.
- Integrated surface quantities color the whole surface uniformly.
- Local curvature colors triangles or sampled surface locations.
- Missing values use a distinct neutral color and appear in the legend.
- Legends always show property name, unit, range, and scale mode.
- Manual, percentile, robust, logarithmic, and shared-comparison ranges.

## 15. Plots and linked data tables

Initial plot types:

- Histogram/distribution
- Scatter plot with selectable x/y properties
- Ranked bar plot
- Box/violin plot by residue, chain, group, or scheme
- Contact matrix/heatmap
- Neighbor/contact network summary
- Property correlation plot
- Trajectory time series with mean and uncertainty
- Scheme/run difference plot

Linked interaction:

- Selecting a plot point selects the corresponding 3D entity and table row.
- Selecting in the viewer filters or highlights plot/table records.
- Table filters update the active plot subset without hiding the full molecule
  unless the user requests it.
- Plot selections can be converted to persistent groups.

Use a single data-selection model shared by viewer, plots, and tables.

## 16. Comparison workspace

Comparison must operate on separate `RunRecord` objects.

### Layouts

- Side-by-side viewers with synchronized cameras
- Single viewer with switchable layers
- Difference coloring on a matched reference geometry
- Side-by-side plots/tables

### Compatibility gate

Before comparing, report:

- matching input structure/fingerprint
- atom/residue match percentage
- scheme and VorPy version
- property availability and unit compatibility
- group/interface definition equivalence
- unmatched entities

Never silently compare by row number.

### Difference operations

- absolute difference
- signed difference relative to a chosen reference
- percent difference where the denominator is valid
- ratio
- rank change
- contact gained/lost

Shared color ranges and synchronized cameras should be the default.

## 17. Trajectory analysis

The confirmed need is analysis across frames; smooth cinematic playback is not
required for the first production milestone.

Initial trajectory features:

- Register topology plus trajectory or a sequence of PDB frames.
- Scrub or step through frames.
- Apply stable groups to each frame.
- Queue one VorPy run per selected frame.
- Track per-frame run status.
- Aggregate scalar properties by atom, residue, group, or interface.
- Plot time series, mean, standard deviation, and confidence/quantile bands.
- Detect and report identity mismatches before running.
- Export frame-level and aggregate tables.

Do not load all frame geometry into VTK simultaneously. Cache a small working
set and unload old frame actors.

## 18. Publication figure workspace

Figures should be reproducible recipes, not screenshots with undocumented
manual state.

### Figure sources

- 3D viewport capture
- property legend
- plot
- data table excerpt
- text/annotation panel
- scale bar and orientation axes

### Controls

- preset sizes for journal column widths and slides
- physical dimensions and DPI
- light, dark, and transparent backgrounds
- camera lock and named views
- font and line scaling
- panel labels
- shared legends and color ranges
- PNG, TIFF, SVG/PDF where supported
- recipe JSON plus rendered output

For vector export, plots and text should remain vector. The VTK viewport can be
embedded as a high-resolution raster when true vector geometry is impractical.

## 19. Controller and service layout

Recommended package evolution:

```text
src/vorpy_workbench/
├── domain/
│   ├── project.py
│   ├── structure.py
│   ├── selection.py
│   ├── group.py
│   ├── run.py
│   ├── result.py
│   ├── comparison.py
│   └── figure.py
├── application/
│   ├── project_controller.py
│   ├── selection_controller.py
│   ├── run_controller.py
│   ├── comparison_controller.py
│   └── commands.py
├── services/
│   ├── structure_loader.py
│   ├── trajectory_loader.py
│   ├── project_store.py
│   ├── vorpy_subprocess.py
│   ├── result_importer.py
│   └── property_registry.py
├── visualization/
│   ├── scene_controller.py
│   ├── layer_controller.py
│   ├── representations.py
│   ├── color_mapping.py
│   └── picking.py
├── ui/
│   ├── main_window.py
│   ├── workspaces/
│   ├── inspectors/
│   ├── panels/
│   └── dialogs/
└── workers/
```

Avoid a single enormous `main_window.py`. The main window should compose
workspaces and route actions; it should not implement parsing, selection
semantics, subprocess management, or scientific calculations.

## 20. State and command architecture

Use explicit commands for meaningful mutations:

- add/remove structure
- create/update/delete group
- create/update/delete interface
- register/import run
- queue/cancel run
- create comparison
- update figure recipe

This enables:

- undo/redo for project edits
- dirty-state tracking
- autosave
- reproducible tests
- future scripting

Viewer-only actions such as camera rotation do not need project commands unless
the camera is saved into a figure recipe or named view.

## 21. Performance rules

- Group atoms by representation/material; never create one actor per atom.
- Combine bond and edge meshes where practical.
- Load large geometry lazily.
- Decimate only for display; preserve original data for analysis/export.
- Separate render quality from scientific data quality.
- Use background workers or processes for parsing, conversion, and aggregation.
- Never update Qt widgets from a worker thread directly.
- Cancel cooperatively and clean up subprocess trees.
- Cache property arrays and identity mappings by immutable dataset ID.
- Record load/render timings for large-system regression testing.

## 22. Persistence

Use a versioned project file, preferably JSON initially:

```text
project.vpyworkbench.json
```

Store:

- object IDs and relationships
- relative or canonical external paths
- selection/group/interface definitions
- registered run metadata
- comparison and figure recipes
- selected UI layout state

Do not embed large geometry or trajectories in the JSON. Reference artifacts
and optionally maintain a project-managed asset directory.

The loader must support schema migration and refuse destructive downgrades.

## 23. Error and diagnostics model

Separate three concepts:

- **User error:** invalid group, missing input, incompatible comparison
- **Run warning:** incomplete cells, missing optional artifact, ambiguous identity
- **Internal failure:** exception, corrupt manifest, renderer failure

The GUI should show concise messages with a details expander and copyable
technical information. Jobs and imported runs retain their warnings instead of
showing them only once in a modal dialog.

Diagnostics should be accessible but not dominate the normal interface.

## 24. Implementation roadmap

Each phase should end with a working application and focused tests.

### Phase 0 — Stabilize the current workbench

- Begin by documenting the current GUI as it exists now; do not assume the
  0.4.0 mockup or any earlier screenshot is authoritative.
- Create a preservation inventory of current visual and functional
  improvements before moving or rewriting code.
- Identify which blueprint goals are already partially or fully implemented so
  later phases do not duplicate or regress them.
- Finish visual QA of the Analysis Studio shell.
- Split current `main_window.py` into composable panels without changing behavior.
- Introduce IDs and controllers around the current `AnalysisResult`.
- Add GUI smoke tests for loading, visibility, repeated selection, and drag safety.
- Establish performance baselines for a small and medium molecule.

**Exit criterion:** current functionality is covered by tests and UI modules are
small enough to extend safely.

### Phase 1 — Project model and persistence

- Implement `Project`, `StructureSource`, IDs, and project store.
- Add New/Open/Save/Save As/Recent Projects.
- Add dirty-state and safe-close prompts.
- Move structure and geometry registration into the project tree.
- Persist window layout and relevant viewer state.

**Exit criterion:** close and reopen a project without losing structures,
registered results, or UI organization.

### Phase 2 — Unified selection, groups, and interfaces

- Introduce `SelectionController` and canonical atom identity.
- Add all selection modes and modifier behavior.
- Implement selection expression tree and rule builder.
- Create persistent groups, group algebra, and frozen/live membership.
- Implement interface definitions and validation.
- Link viewer, tree, inspector, and data table selection.

**Exit criterion:** create groups manually or by rule, combine them, define an
interface, save, reopen, and recover identical membership.

### Phase 3 — Result manifest and legacy import

- Define the versioned result manifest with VorPy.
- Implement manifest-backed result loading.
- Preserve a legacy directory importer for existing outputs.
- Build property registry and normalized entity tables.
- Add compatibility and completeness reports.

**Exit criterion:** imported and newly manifested results expose the same
Workbench-facing `ResultDataset` API.

### Phase 4 — VorPy subprocess runner

- Implement command/config translation.
- Add JSONL event parsing and durable job records.
- Add progress, cancellation, failure recovery, and output import.
- Configure one scheme per run.
- Add run duplication for easy scheme comparison.

**Exit criterion:** configure, execute, cancel, reopen, and inspect a VorPy run
without blocking or crashing the GUI.

### Phase 5 — Full result visualization

- Add representation and property registries.
- Render surface/edge/vertex layers by scope.
- Add property legends, range controls, missing values, and provenance.
- Add linked plots and tables.
- Add contact and neighbor views.
- Add quality/completeness visualization.

**Exit criterion:** every initial property family can be inspected numerically,
colored in 3D at its correct scope, plotted where meaningful, and exported.

### Phase 6 — Comparison workspace

- Implement entity matching and compatibility gate.
- Add synchronized side-by-side viewers.
- Add differences, ratios, percent changes, and contact gained/lost.
- Add matched plots/tables and shared range controls.

**Exit criterion:** compare Atomic Voronoi, Power, and Primitive run records
without relying on array row order.

### Phase 7 — Trajectory analysis

- Register frame sets and validate identities.
- Queue frame runs and track frame-level status.
- Add frame stepping and bounded geometry cache.
- Aggregate and plot properties across frames.
- Export frame and summary datasets.

**Exit criterion:** apply a group across frames, run selected frames, and plot a
group/interface property over time with mismatch reporting.

### Phase 8 — Publication figures

- Implement figure recipes and named cameras.
- Add plot/viewport/table panels.
- Add journal and slide presets.
- Export high-resolution figures and recipe metadata.

**Exit criterion:** reopen a recipe and regenerate a visually equivalent figure
from the same registered results.

### Phase 9 — Packaging and integration hardening

- Package for Linux and Windows.
- Validate clean installations.
- Add crash reporting/log bundles with opt-in handling.
- Document supported VorPy/result schema versions.
- Define migration point from separate prototype to official VorPy distribution.

## 25. Testing strategy

### Unit tests

- identity generation and ambiguity
- selection expression parsing/evaluation
- group algebra
- interface validation
- run command translation
- JSONL event parsing
- manifest schema and migrations
- property metadata and color range calculations
- run comparison matching and differences
- project serialization

### GUI tests

- action-to-workspace navigation
- project-tree selection and visibility independence
- repeated atom/residue selection
- drag-safe selection
- inspector updates
- group creation workflow
- run editor validation
- job progress/cancellation
- linked viewer/plot/table selection
- save/reopen layout

### Integration tests

- fake VorPy subprocess with deterministic events
- real small EDTA run
- import of representative legacy output directories
- scheme comparison with known matched entities
- multi-frame fixture with stable and intentionally broken identity

### Performance fixtures

- small: EDTA
- medium: cambrin or equivalent
- large: representative protein/solvated system

Track structure load, actor construction, geometry load, selection latency,
frame switch, and peak memory.

## 26. First Codex implementation prompt

Use this prompt from the `vorpy-workbench` repository with the sibling `vorpy`
repository added to the Codex workspace:

```text
Read AGENTS.md and VorPy_Workbench_Functionality_Blueprint.md completely.
Inspect the current worktree, recent commits, tests, and all current UI modules.
Report any uncommitted changes. The GUI has improved since this blueprint was
written: current working behavior and visual improvements are the baseline,
not the older version described here. Do not begin VorPy integration yet.

Implement Phase 0 only:
1. inventory all existing functionality, layout improvements, styling,
   interactions, performance work, and tests;
2. map the current implementation to blueprint goals and identify what is
   already complete, partial, or missing;
3. rank the unmet goals by dependency, scientific value, and regression risk;
4. propose the smallest safe module or architecture change needed for the
   highest-priority goal;
5. preserve every unrelated current improvement, including repeated drag-safe
   selection, water/ion visibility, result loading, layer controls,
   screenshots, background execution, and any newer functionality discovered;
6. add focused regression tests before refactoring behavior;
7. run compileall, pytest, ruff, and an interactive smoke test if the display is
   available; and
8. update AGENTS.md with the new checkpoint and preservation inventory.

Before editing, show the preservation inventory, goal-priority table, proposed
file changes, and test plan. Stop for confirmation if the proposal removes or
substantially changes an existing working feature. Do not modify the sibling
VorPy repository during Phase 0.
```

## 27. Decisions intentionally deferred

These should not block Phase 0–2:

- exact trajectory formats beyond an adapter interface
- whether ribbon/cartoon rendering is implemented internally or through a
  specialized dependency
- final binary packaging technology
- direct Python-library integration with VorPy
- remote/HPC job execution
- collaboration or shared project storage
- plugin architecture for third-party analyses

Defer them until the core project, selection, and result schemas are stable.

## 28. Product success criteria

The Workbench is conceptually successful when a scientist can:

1. create or open a project;
2. load a molecule or trajectory;
3. create reusable groups and an interface;
4. configure and run one VorPy scheme;
5. register existing results;
6. inspect surfaces, edges, vertices, and calculated properties;
7. move seamlessly between 3D selection, plots, and tables;
8. compare compatible runs without ambiguous entity matching;
9. analyze properties across frames; and
10. reproduce a publication figure from a saved recipe.

The GUI should make this sequence obvious without exposing the internal class
structure of VorPy to the user.
