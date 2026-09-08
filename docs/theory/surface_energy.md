# Representative Surface Energy

VorPy provides `surf_energy`, a representative curvature-dependent bending
quantity rather than a calibrated molecular free energy.

The current reference model is

```text
kappa_b = 1 kBT
C0 = 0
E_rep / kBT = 2 integral(H^2 dA)
             = 2 int_mean_curv_sq.
```

`int_mean_curv_sq` is evaluated over smooth triangulated surface patches. It is
orientation invariant and dimensionless when coordinates are in angstroms.
The current `surf_energy` therefore reports a reference bending score in units
of `kBT` under the stated unit-modulus convention.

This field does not include the planned edge correction to classical
integrated mean curvature. It also does not include electrostatics, dispersion,
solvent chemical potential, hydrophobicity, conformational entropy,
atom-specific surface tensions, fitted elastic constants, or calibrated
spontaneous curvature.

The later morphometric solvation model is a separate calibration layer:

```text
Delta G_np = p V + gamma A + kappa M + kappa_bar G.
```

Here `M` and `G` must eventually denote validated complete-boundary measures.
They must not be confused with the current reference bending energy.

Recommended reporting language is **representative surface energy**,
**reference bending energy**, or **curvature-dependent reference energy**.

Implementation: `calculations/surface_energy.py` and
`network/build_surfs.py`. See
[Morphometric Geometry Implementation](../development/morphometric_geometry.md)
for status and field ownership.
