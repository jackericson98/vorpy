# Validation shapes

Reproduce the verified two-handle boundary from the repository root:

```sh
python -m vorpy.src.geometry.topology_diagnostic --torus-count 2 --major-radius 10 --minor-radius 3 --interior-count 16 --cross-section-resolution 12 --expect-genus 2 --output data/topology_audit/g2
```

The diagnostic constructs the weighted Voronoi network and derives selected
cell connectivity, boundary (V,E,F), Euler characteristic, genus, and
integrated Gaussian curvature from the assembled boundary. `--expect-genus`
checks those measurements and exits nonzero on failure; it does not change
the geometry or curvature calculations. It also checks that neighboring
sections share positive-area selected-selected Voronoi faces, and checks
boundary manifoldness and face winding independently.

The verified genus-2 junction uses one shared selected generator displaced by
`0.1 * minor_radius` perpendicular to both loop planes. Constraint-ring
phases vary between loops and samples to remove high-order co-spherical
junction vertices. Constraint spheres buried inside another loop's shell
are clipped. `--junction-offset` changes the displacement; start with the
default. The verified settings measure 31 complete selected cells, one
selected component, one closed orientable manifold boundary, Euler
characteristic -2, genus 2, and integrated Gaussian curvature -4*pi.

To write a geometry package without running the boundary diagnostic:

```sh
python -m vorpy.src.geometry.validation_shapes --shape multitorus --torus-count 2 --major-radius 10 --minor-radius 3 --interior-count 16 --cross-section-resolution 12
```

This writes XYZR, PDB, and PyMOL files to
`data/validation_shapes/multitorus_n2_R10_r3_i16_x12/`. Set `--torus-count`
to any positive integer. Counts above two are accepted by the generator but
remain unvalidated; the earlier five-loop geometry is connected through
positive-area faces, yet has an incomplete, nonmanifold boundary.

`--interior-count` controls centerline samples per loop. Three or more loops
round odd counts up to even so each middle loop has a sample at both joins;
the generated shape records requested and actual resolution. Interior PDB
generators are labeled `INT`; pass `--no-center` to export only surrounding
constraint spheres.

Python callers can use `multitorus(10.0, 3.0, torus_count=2,
interior_count=16)`. For a chain of `n` handles, the analytic expectation is
genus `n`, Euler characteristic `2 - 2*n`, and integrated Gaussian curvature
`4*pi*(1 - n)`. These expectations do not replace validation of the
constructed boundary. Junctions alter the mean curvature, so its reported
reference value is not an exact junction result.
