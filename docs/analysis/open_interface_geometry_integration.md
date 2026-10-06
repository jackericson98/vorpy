# Validated open-interface cache integration

The production integration lives in `vorpy/src/interface/geometry_analysis.py`.
It calls the frozen `open_interface_curvature` kernel on the alpha-selected
physical face incidence already owned by the Interface. `Interface.build()`
and `iface.analyze_geometry()` retain their cache API; `to_dict()`, normal
info export and the existing interface-geometry log record carry the extension.
Reverse orientation is derived from the forward cache. No additional network
solve, alpha selection change, radius change or mathematical implementation is
introduced. Schema version 2 invalidates an older analysis cache once, while
reusing its solved network and filtration source.

## Structured quantities

Each side now carries physical V/E/F, Euler characteristic, components,
boundary loops, genus by component and manifold-incidence certification.
Mean fields distinguish smooth, internal seam, zero free-boundary, partial
and certified total contributions, with supported/unsupported face and seam
counts and coverage. Intrinsic Gaussian fields distinguish smooth, seam and
interior-vertex contributions from exterior boundary geodesic and corner terms.
Gauss-Bonnet LHS, expected value, residual, support scope and certification are
separate. Precise field names and compatibility aliases are documented in
`docs/formats/output_files.md`.

The existing generic `total_K`/`total_gaussian_curvature` keys remain aliases
for certified intrinsic K, never Gauss-Bonnet completion. Unsupported totals
are null. Finite cached diagnostics do not override validated analyzer support
flags. The surface-only legacy `smooth_K_partial` remains available, including
the EDTA-Mg regression's 0.8904897507356364 value.

Certification denotes complete supported numerical accounting, not a formal
interval-arithmetic proof. It requires the validated manifold and contribution
checks, and the kernel's Gauss-Bonnet residual tolerance. No scheme name alone
is a certificate. The full selected AW interface cannot be certified by an
audit of a smaller supported subset. Its subset diagnostics explicitly retain
their own topology, quadrature order, contributions, target and residual.

## Frozen examples

Power: 350 selected pairs, area 898.1185793538328 square angstroms;
V=605, E=955, F=350, chi=0, two components, four boundary loops and genus [0,0].
All 835 seams are supported. H=-64.45049624200274 angstroms is certified.
Intrinsic K=5.188337486733097; boundary corner turning=-5.1883374867332295.
Completed Gauss-Bonnet is approximately zero, with residual -1.3233858453531866e-13.
Intrinsic K is not replaced with 2*pi*chi.

AW: 339 selected pairs, area 787.0925777403216 square angstroms;
V=584, E=924, F=339, chi=-1, two components, five boundary loops and genus [0,0].
Supported smooth H=1.6219633715872277, seam H=-55.77212428288777,
partial H=-54.150160911300546 angstroms. Five faces and 14 seams remain
unresolved; complete H and intrinsic K are null. The separate 334-face subset
has residual -0.0006654952459843599 at order 256. The full-interface completed
LHS and residual are null, and Gaussian certification is false.

Real EDTA-Mg: six eligible pairs, four selected pairs, five interface atoms,
one component, area 5.868263011906469 square angstroms, smooth K partial
0.8904897507356364. Both totals remain uncertified. The area agrees with the
old snapshot to one floating-point unit. Both standard log readers parse the
new normal export, and old logs remain readable.

## Reproduction

Run from the repository root, with Python's existing scientific dependencies:

```powershell
python -m scripts.export_open_interface_geometry_examples
python -m scripts.export_interface_geometry_example --network output/interface_geometry_phase1/EDTA/source_aw.vpy --group1-atoms 0-31 --group2-atoms 32 --output-dir output/open_interface_integration/EDTA
python -m pytest -o addopts='' vorpy/tests/interface/test_open_geometry_integration.py -q
python -m pytest -o addopts='' vorpy/tests/interface/test_geometry_analysis.py vorpy/tests/analyze/test_open_interface_curvature.py vorpy/tests/analyze/test_power_interface_curvature.py vorpy/tests/analyze/test_aw_interface_curvature.py vorpy/tests/analyze/test_cazals_2kai_curvature_forensics.py vorpy/tests/analyze/test_cazals_validation.py vorpy/tests/interface/test_interface_water.py -q
```

Power/AW examples appear under `output/open_interface_integration/2KAI/pow`
and `aw` as `info.txt`, typed `interface_geometry.csv` extension records, and
`geometry_analysis.json`. EDTA writes the normal detailed info and interface
log export. These commands reuse local frozen fixture artifacts; the static
and real EDTA fixture tests skip explicitly if their local files are absent.

The mathematical kernel SHA-256 remains
`472F9D3A6989668A3818AB077B774271A257959DC514738ED679968DD9CD0A53`.
The kernel and its mathematical tests were not edited for this integration.

Validation completed: 11 integration tests passed, including both real static
schemes, reader compatibility, both orientations and EDTA artifact regression.
The combined existing mathematical, Phase-1, Power/AW, water and frozen Cazals
suite passed 69 tests. The final cache/serializer/log-column subset passed
43 tests after the last cache-version and failure-status adjustments.
