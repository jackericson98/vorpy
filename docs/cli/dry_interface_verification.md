# Dry interface regression verification

The exact command completed with exit code 0:

```console
python vorpy 2kai -g c 0 and c 1 -g c 2 -i
```

Output: `output/2KAI_7/`. The suffix reflects pre-existing output directories.
The complete generated directory tree is in `docs/cli/2kai-output-tree.txt` (98 files).

The first failing record was HETATM 2241, oxygen of HOH A 109. It is genuine
crystallographic water. Protein chain A already existed, so the old chain branch
did not initialize `sys.sol`, but the residue branch appended to it. PDB parsing
now initializes the established empty solvent container at the start and assigns
recognized water/ions to it independently of author chain labels. Protein and
generic hetero residues remain molecular topology. Shared residue metadata is
the sole classifier. Water analysis filters water names separately from ions.

`System('2kai')` now uses the existing bundled-input resolver for bare aliases.
The default surface export adds cached selected physical geometry separately
from the full physical `surfs.off`. Mapping IDs, selected surface IDs, and each
surface's generator pair are checked before export. No selection is recalculated
and no additional network is solved. Default archive writing omits only the
runtime water lookup closure in topology snapshots and per-water snapshots.
Legacy archives without Phase-1 caches retain their existing full geometry
export without inventing a new selection.

| Metric | Exact bundled input | Water-free regression input |
|---|---:|---:|
| Input atoms | 2246 | 2236 |
| Group 1 atoms (chains A+B) | 1799 | 1799 |
| Group 2 atoms (chain I) | 437 | 437 |
| Explicit solvent | yes, 10 HOH oxygens | no |
| Scheme / alpha | AW / additive 0 A | AW / additive 0 A |
| Eligible bicolor pairs | 475 | 480 |
| Alpha-selected pairs | 27 | 27 |
| Mapped selected physical surfaces | 27 | 27 |
| Selected area (A^2) | 127.46399209140874 | 127.67974816259849 |
| Topology V/E/F | 142/159/27 | 143/160/27 |
| Components | 14 | 14 |
| Boundary loops | unresolved | unresolved |
| Smooth mean partial contribution (A) | -0.130186 (info precision) | -0.12824789018852367 |
| Mean partial sum / total | unresolved / unresolved | unresolved / unresolved |
| Mean certification | no | no |
| Smooth Gaussian partial contribution | 0.111464 (info precision) | 0.11147673281714146 |
| Intrinsic Gaussian total | unresolved | unresolved |
| Gauss-Bonnet LHS / residual | unresolved | unresolved |
| Gaussian certification | no | no |

Selection and physical mapping have 100% coverage. Existing AW cell completeness
and selected-incidence requirements prevent certified curvature totals for this
unexpanded alpha-zero patch. Solvent absence does not gate these computations.
The existing mathematics and certification rules were preserved.

Open the selected physical interface in PyMOL:

```text
output/2KAI_7/c_0_and_c_1_c_2_interface/alpha_interface/selected_interface.pml
```

Open full/restricted Apollonius dual visualization:

```text
output/2KAI_7/c_0_and_c_1_c_2_interface/dual/apollonius_interface.pml
```

The selected dual layer is `apollonius_selected`; its mesh is
`dual/apollonius_selected_edges.off`. Both visualizations use the same solved
interface network. No PyMOL rendering was required.

Validation completed:

- Real water-free 2KAI command pipeline: passed, exactly one Network.build call.
- Broad parser/water/Phase-1/open-curvature/Cazals suite: 128 passed.
- Boundary and scheme-specific alpha suite: 43 passed.
- Frozen independent 2KAI Cazals audit and validation: 13 passed.
- Final focused geometry/water/new regression suite: 23 passed, slow test run separately.
- Actual solved archive callback round-trip: passed, runtime source context preserved.
- Final archive and pair/surface/mesh identity checks: passed, zero solves;
  27 surfaces, 4577 points and 7149 triangles for the bundled input.
- Dry mesh identity check: passed, 27 surfaces, 4584 points and 7157 triangles.
  Coordinates agree within the established OFF four-decimal serialization precision.
- Archive suite: 18 passed, one unrelated existing failure in
  `test_cli_archive_input_without_build`: the command pipeline has an existing
  `_run_archive` method but does not dispatch `.vpy` inputs to it, and consequently
  attempts molecular grouping on an unloaded System. This dispatch was left outside
  the parser/selected-interface scope.

Files changed for this task:

- `vorpy/src/inputs/pdb.py`
- `vorpy/src/system/system.py`
- `vorpy/src/interface/water.py`
- `vorpy/src/interface/export.py`
- `vorpy/src/interface/geometry_analysis.py` (export guidance only)
- `vorpy/src/interface/selected_export.py` (new)
- `vorpy/src/io/network_archive.py`
- `vorpy/tests/interface/test_dry_pdb_interface.py` (new)
- `docs/cli/dry_interfaces.md` (new)

Pre-existing and concurrent changes to dual visualization and other files were preserved.
