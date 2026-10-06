# Output Files

VorPy can export atom logs, surface logs, edges, vertices, shells, surrounding atoms, PDB structures, system information, group information, interface information, and curvature-derived quantities.

Exports can be requested through presets or individual components.

**TODO:** Add an authoritative file-by-file reference with filenames, schemas, units, and ownership level (system/group/interface).

## Sectioned run logs

`aw_logs.csv`, `pow_logs.csv`, and `prm_logs.csv` contain build information,
group information, atoms, surfaces, edges, and vertices. Full-network and
interface exports share column ordering rules in `vorpy/src/log_columns.py`.

Columns are grouped as follows:

- Build: source and version, network settings, vertex limits, timings.
- Group: volume/area/mass/density, mean curvature, Gaussian curvature,
  boundary and topology diagnostics, centers of mass and inertia tensors.
- Atoms: identity, input coordinates/radius/mass, cell geometry and shape,
  cell status, neighbors and contacts, distance statistics, mean curvature,
  Gaussian curvature, energy, spatial detail.
- Surfaces: identity, area/volume contributions/contact/overlap, mean curvature,
  Gaussian curvature, energy.
- Edges and vertices retain their compact identity/geometry/curvature order.

Column names and values are retained, including curvature total aliases.
Ordering does not change measurement definitions. Existing logs are not rewritten.

`vorpy.src.inputs.logs.read_logs` and the analysis readers `read_logs2` and
`read_logs` interpret modern columns by header name. The older analysis reader
preserves its historical output keys, and positional row parsers remain available
for legacy layouts. New consumers should use these readers instead of numeric
CSV column offsets. Archived analysis scripts that read raw CSV offsets directly
are not a supported schema interface.

## Cached alpha-selected interface geometry (Phase 1)

After a normal `Interface.build()`, `iface.geometry_analysis` contains an
`InterfaceGeometryAnalysis`. `iface.analyze_geometry(alpha=None, refresh=False)`
returns the cached result; `to_dict()` exposes its stable structured sections.
Exporters serialize this result and do not calculate it. The default query is
the existing native alpha-zero convention; `settings['interface_alpha']` or an
explicit method argument can supply another native value. AW uses additive
angstroms; Primitive and Power use their respective squared-radius and power
distance conventions in square angstroms. These values are not interchangeable
physical probe offsets.

The result separates alpha selection from selected physical Voronoi geometry.
It retains both orientations, mean-curvature partials, feature coverage and
certification, and Gaussian contributions. Missing totals are `None`.
Schema version 2 adapts the independently validated open-interface kernel:
intrinsic integrated Gaussian curvature is separate from exterior boundary
terms and the Gauss-Bonnet completed left-hand side. Certification requires
complete numerical accounting and manifold incidence, plus the kernel residual
tolerance for Gaussian accounting; it is not a formal interval-arithmetic proof.
Molecular union-of-balls contact patches remain `not_implemented`.

Standard interface logs append a three-row `Interface Geometry` section:
section label, stable metric headers, and one data row. Values are JSON literals
inside CSV cells, preserving booleans, nulls, lists, and strings. Both standard
readers expose the typed record as `result['interface geometry']`; existing
build/group/atom/surface/edge/vertex sections retain their schemas. The older
analysis reader also carries this extension through its modern-log adapter.

The cached side dictionaries now additionally contain:

- Topology: `topology_vertices`, `topology_edges`, `topology_faces`,
  `euler_characteristic`, `connected_components`, `boundary_loops`,
  `genus_by_component`, `manifold_incidence_certified`.
- Mean curvature (angstroms): `smooth_integrated_mean_curvature`,
  `internal_seam_integrated_mean_curvature`, `boundary_integrated_mean_curvature`,
  `partial_integrated_mean_curvature`, `total_integrated_mean_curvature`,
  `mean_curvature_certified`, `mean_supported_faces`, `mean_unsupported_faces`,
  `mean_supported_internal_seams`, `mean_unresolved_internal_seams`, `mean_coverage`.
  The boundary contribution is zero under the validated scalar convention.
- Intrinsic Gaussian curvature: `smooth_integrated_gaussian_curvature`,
  `intrinsic_seam_gaussian_curvature`, `interior_vertex_gaussian_curvature`,
  `intrinsic_integrated_gaussian_curvature`.
- Exterior boundary and completion: `boundary_geodesic_curvature`,
  `boundary_corner_turning`, `gauss_bonnet_lhs`, `gauss_bonnet_expected`,
  `gauss_bonnet_residual`, `gaussian_curvature_certified`.
- Gaussian support: `gaussian_supported_faces`, `gaussian_unsupported_faces`,
  `gaussian_coverage`, `gaussian_accounting_scope`,
  `supported_subcomplex_gauss_bonnet`. A supported-subcomplex diagnostic is
  explicitly distinct from completion of the full selected interface; full
  `gauss_bonnet_lhs` and residual remain null when only a subset was evaluated.

The Gaussian and angle-integral values are dimensionless. Existing keys remain
present: `total_K` and `total_gaussian_curvature` are compatibility aliases for
the certified **intrinsic** total, never for `gauss_bonnet_lhs` or `2*pi*chi`.
The old `smooth_K_partial` continues to expose the measured surface-only partial.
`boundary_geodesic_gaussian` now denotes only the exterior boundary contribution;
intrinsic seam terms have their own explicit field. Each side field is serialized
with a `vor_` prefix; orientation-specific explicit fields additionally have
`vor_g1_to_g2_` and `vor_g2_to_g1_` prefixes. Existing readers need no schema changes.

Reproduce the frozen matched-probe cache-adapter examples without tessellation
or alpha selection changes:

```powershell
python -m scripts.export_open_interface_geometry_examples
python -m pytest -o addopts='' vorpy/tests/interface/test_open_geometry_integration.py -q
```

Examples are written to `output/open_interface_integration/2KAI/{pow,aw}/` as
`info.txt`, `interface_geometry.csv`, and `geometry_analysis.json`. The example
adapter honors saved analyzer support flags, not merely finite diagnostic values.

Recognized physical solvent outside the two side definitions participates only
as competing generators. Virtual boundary sites do not count as chemical solvent.
No additional dry network is solved. Dry/solvent comparison is available only
when the same System holds an earlier dry analysis with matching non-solvent
generator geometry, radii, group memberships, scheme, alpha, and construction
settings. Otherwise its status explicitly reports unavailable. Gaussian deltas
remain null; mean deltas require two certified totals.

Reproduce the real solvated EDTA/Mg example from the saved source archive:

```powershell
python -m scripts.export_interface_geometry_example --network output/interface_geometry_phase1/EDTA/source_aw.vpy --group1-atoms 0-31 --group2-atoms 32 --build-interface --output-dir output/interface_geometry_phase1/EDTA_interface
```

The explicit `--build-interface` option exercises the normal Interface solve;
Phase 1 itself performs no additional solve. This input currently exposes an
existing later buried-water solver error. The reproduction helper reports that
error and exports the geometry already cached before the water stage. It does
not suppress errors in the normal application pipeline. The example's alpha
selection and area are available, but incomplete AW support prevents certified
mean totals. No unsupported mean value is replaced with zero.

For larger saved full networks, the default reuse mode can spend substantial
time in the existing complete AW pair-birth query. Phase 1 does not change its
mathematics or introduce a new restricted-query implementation.

To regenerate the source archive using existing solver and shell defaults:

```powershell
python -c "from pathlib import Path; from vorpy.workbench.services.vorpy_backend import VorPyBackend,VorPySolveSettings; from vorpy.src.io import save_network; destination=Path('output/interface_geometry_phase1/EDTA/source_aw.vpy').resolve(); result=VorPyBackend(VorPySolveSettings(boundary_mode='shell')).solve(Path('vorpy/data/EDTA.pdb').resolve(),lambda *_:None,lambda:False); save_network(result.export_group.net,destination)"
```

Focused verification (timeout addopts are overridden because this environment
lacks pytest-timeout):

```powershell
python -m pytest -o addopts='' vorpy/tests/interface/test_geometry_analysis.py vorpy/tests/geometry/test_scheme_alpha_filtrations.py vorpy/tests/geometry/test_alpha_interfaces.py vorpy/tests/geometry/test_aw_alpha_experimental.py vorpy/tests/analyze/test_aw_interface_curvature.py vorpy/tests/analyze/test_power_interface_curvature.py vorpy/tests/analyze/test_cazals_validation.py vorpy/tests/analyze/test_cazals_2kai_independent_audit.py vorpy/tests/analyze/test_cazals_2kai_curvature_forensics.py vorpy/tests/calculations/test_aw_interface_orientation.py vorpy/tests/interface/test_interface_water.py vorpy/tests/inputs/test_read_logs2.py vorpy/tests/inputs/test_log_column_order.py -q
```
