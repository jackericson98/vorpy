# Workbench layout A checkpoint

Implemented 2026-09-12 against the integrated `vorpy/workbench` application.
The user's large layout A reference governs this change; B/C/D are deferred.

## Layout and ownership

- `ui/main_window.py`: authoritative application state, existing actions and
  panel controls, selection/group/interface operations, solver orchestration,
  summary population, empty/busy state presentation, project view-state hooks.
- `ui/panels.py`: reusable `WorkflowSidebar`, `ViewerPanel`, `ViewInspector`,
  and `ResultsInspector`, plus shared scroll containers and action buttons.
  These containers do not own scientific state.
- `ui/theme.py`: dark workbench styling with neutral cards and purple active
  controls. Existing assets and rendering remain unchanged.
- `ui/molecular_view.py`, worker, backend and project data model are unchanged.

The header links to existing Structure, Selection, Analysis, Results and
Settings controls. The horizontal splitter contains workflow, viewer and right
inspector. Solve is below System/Network display tabs. The vertical splitter
resizes Results, which has a separate collapse/expand control. Project view
state stores both splitter sizes and Results expansion, accepting older files
that omit those keys.

System remains the default display tab. Selection modes no longer switch to
Network; they preserve the chosen display tab. Selection actions are shared
between menus, the sidebar and viewer toolbar. Camera drag, delayed picking,
replacement/Shift selection, clipping, water/ion behavior and screenshot export
continue through the original viewer handlers.

## State and results

Empty Structure offers Open structure and Load example (packaged EDTA).
Loaded Structure retains the chemical information and multi-structure list.
Solve requires a loaded structure and is disabled during a job. File/session
switching is disabled during a job; cancellation and progress remain available.
Group and Interface creation explain their selection/group prerequisites.
Reset clears stale result/display presentation.

Overview prioritizes atoms, residues, chains, complete cells, surfaces,
volume, surface area and **defined interfaces**. Defined interfaces count saved
interface definitions, not computed interface measurements. Existing vertices,
edges, geometry, curvature, energy, classification and composition data remain
in detailed tables. The new Interfaces table lists definitions and their groups.

## Intentional compatibility boundaries

- Desktop structure loading and integrated solving still accept PDB. The CLI
  supports broader formats; the UI does not advertise them as desktop imports.
- Volume/area cards use `NetworkSummary` measurements with units when present.
  The current real solver adapter does not attach that metadata; unavailable
  values remain dashes. This redesign does not invent measurements or alter
  solver analysis behavior.
- Display controls scroll when the inspector cannot fit their full contents.
  Results scroll internally at small heights and can be expanded or collapsed.
- Backend scientific settings and defaults are unchanged. The maximum vertex
  radius help text now describes the existing radius parameter correctly.

## Validation

- Workbench tests: 34 passed, including project save/open, selections,
  Group/Interface behavior, representations, depth clipping and screenshots.
- Added regression checks for empty/loaded/busy/reset state, System default,
  maintaining the inspector tab when selecting, splitter round-trip, Results
  collapse/expand, and rejecting a pending pick after a drag.
- Updated a stale backend scalar-key assertion for the already implemented
  integrated mean/Gaussian curvature outputs.
- Actual X11 launch and rendered empty/loaded EDTA inspection.
- Actual GUI worker Power solve of EDTA's first residue (32 selected atoms):
  32 complete cells, 504 vertices, 425 surfaces. Verified Results population,
  enabled layer controls, geometry rendering and PNG screenshot export.
- Na5 Power fixture returned no displayable geometry; EDTA was used for the
  successful end-to-end check. No backend changes were made for this condition.
- Full repository pytest collection is blocked by missing `hypothesis` in
  `vorpy/tests/calculations/test_calcs.py`; this is not a UI test failure.
- UI compile and diff whitespace checks passed.

Validation captures and temporary scripts are under `/tmp/vorpy-layout-a-*`
and `/tmp/vorpy_*check.py`; they are not application dependencies.
