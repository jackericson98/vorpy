# Apollonius dual and additive-weighted alpha complexes

These options are opt-in and run after the ordinary AW network has been
constructed. They use only the solved VorPy incidence tables; they do not build
an ordinary Delaunay triangulation or change the network solver.

Export the complete straight-line dual and per-dimension CSV tables:

```bash
python vorpy molecule.pdb -s mv 5 --apollonius
```

Calculate supported additive-weighted birth values and extract the subcomplex
at alpha 1.4 Å:

```bash
python vorpy molecule.pdb -s mv 5 --alpha-complex 1.4
```

The alpha value is an outward additive expansion of atomic radii,
`effective_radius = r_i + alpha`, for the distance
`||x - p_i|| - r_i`. Alpha implies construction/export of the Apollonius
complex. Files are written under the run output's `apollonius/<network>/`
directory, including per-dimension CSV files and a text summary.

The current exact calculations support cell births when the generator center
lies in that AW cell, bounded analytic AW edge minima, and AW vertex clearances
verified against all defining generators. Exact minima over trimmed pairwise
surface patches and cells whose generator center is excluded are reported as
unresolved. Alpha extraction omits cofaces whose faces are unavailable, so the
result remains closed under taking faces. No sampled-mesh approximation is
used to fill unresolved birth values.

The current network vertex table records four defining generators. Higher
order degeneracies are therefore reported as unsupported if the complete
incident set is not present in the network data.
