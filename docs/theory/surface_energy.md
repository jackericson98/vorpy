# Representative Surface Energy

VorPy provides `surf_energy`, a representative/reference curvature-dependent bending quantity rather than a calibrated molecular free energy.

Current reference model:

```text
kappa_b = 1 kBT
C0 = 0
E_rep / kBT = 2 * integral(H^2 dA)
E_rep / kBT = 2 * int_mean_curv_sq
```

It can be used as a scalar measure of surface bending, but it does not by itself include electrostatics, dispersion, solvent chemical potential, hydrophobicity, conformational entropy, atom-specific surface tensions, fitted elastic constants, or system-specific spontaneous curvature.

Recommended reporting language: **representative surface energy**, **reference bending energy**, or **curvature-dependent reference energy**.
