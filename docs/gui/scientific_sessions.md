# Scientific sessions (.vpy)

Workbench Save Session As writes a self-contained `.vpy` session. File → Open
Session / Project accepts sessions, older network `.vpy` archives, and legacy
`.vpyworkbench.json` projects. Saving an already-open legacy project migrates to
a sibling `.vpy`, leaving the original JSON untouched. The legacy JSON API remains
available for compatibility.

The original molecular file is not required to inspect a session. Loading does
not solve a network, construct a dual, calculate filtration births, select a
physical interface, triangulate, or integrate curvature. Missing legacy caches
remain missing. Existing export names, presets, and directory conventions are
unchanged.

## Archive contract

The existing ZIP/JSON/NumPy archive has format version **3** and
`document_kind: network | session`. Both kinds share the same field registry,
restricted codec, and scientific object graph. Version-1 and version-2 network
archives remain readable; readers predating version 3 cannot read archives that
use the new representation-cache fields. Each registered object record has
`schema_version: 1`. Unknown classes, fields, and record versions are rejected.
No pickle or archive-specified imports are used.

In addition to the existing molecular context and scientific DataFrames, both
document kinds persist attached interface analysis, DualComplex simplices and
primal references, AlphaFiltration records/blocked births, its AW source and
incidence, physical interface surfaces/edges/vertices, CurvatureSummary,
selection mappings, alpha incidence audit, cached coverage/certification and
unresolved records. Boundary configuration/generators, shared interface geometry,
and dry comparison caches are preserved. The session kind also stores multiple
systems/networks and Workbench results/project state.

Objects receive persistent UUID `archive_id` values. Archive-local integer
references preserve shared object identity; they are not Python addresses.
Network-local generator IDs, original atom mappings/labels and explicit generator
provenance remain distinct. Cached aggregate records carry separate
`archive_provenance` (network/system/interface IDs, scope, partition, alpha),
without rewriting their scientific metadata. Dedicated interface-network and
system-network scopes are explicit. Ownership mismatch is an error; inconsistent
scientific counts/metadata remain preserved and can be flagged by inventory.

Disposable table-identity cache keys are excluded and rebound after loading.
Stale cache state is preserved, rather than validating stale analysis by rebinding
it. Spatial acceleration structures, drawing caches, Qt/VTK actors, progress
callbacks and water-discovery closures remain runtime-only.

Workbench sessions preserve loaded results and cached trajectory frames (including export-group references),
active source/network, per-source groups/interfaces/selections, atom/bond/layer
arrays, visibility/styles, display and layout settings, backend/export settings,
active solve target and camera position. Trajectory file offsets may be stored,
but accessing unsaved frames still requires the trajectory source. Workers,
preload queues and live renderer objects are not stored.

## Programmatic use

```python
from vorpy.src.io import Session, save_session, load_session

save_session(Session.from_network(network), "analysis.vpy")
session = load_session("analysis.vpy")
network = session.active_network
interface = network.sys.ifaces[0]
analysis = interface.geometry_analysis

# Existing exporters consume the restored cached objects directly.
from vorpy.src.interface.selected_export import export_selected_interface
from vorpy.src.geometry.visualization.interface import export_interface_dual_visualization
interface.dir = "existing_export_destination"
export_selected_interface(interface)
export_interface_dual_visualization(interface)

capabilities = session.capabilities(
    interface=interface.interface_id,
    representation="Voronoi",
    result="mean_curvature",
)
```

The inventory reports AVAILABLE, PARTIAL, UNRESOLVED, NOT_CALCULATED,
NOT_SUPPORTED, or STALE using cached data only. Filters include network UUID,
interface ID, environment (`dry`, `solvent_competing`, or unknown), network scope
(`system` or `interface`), partition (`AW`, `Power`, `Primitive`), representation,
alpha query and result. Historical `wet` and `interface/dedicated` labels are
normalized when loaded. Existing incidence duals
use `dual`; future `Dual A` / `Dual B` are currently NOT_SUPPORTED. Requested
missing combinations are NOT_CALCULATED. Inventory does not repair discrepancies
or create missing analyses.

## Verification

`vorpy/tests/io/test_scientific_session.py` covers guarded cached round trips,
legacy version-1 loading, network/session document kinds, shared identity,
boundary state, partial/unresolved/stale preservation, inconsistent provenance,
Workbench reopen and both existing exporters. Its slow 2KAI test uses the local
validated network, matched-probe birth checkpoint and saved curvature component
tables. It materializes the reference before saving, then starts a fresh Python
process with scientific constructors/builders/integrators forbidden. Numerical
tables, mappings, metadata and geometry files must match; only absolute launcher
paths are normalized. The large local scientific artifacts are not test fixtures
distributed with the repository, so that acceptance case skips when absent.
