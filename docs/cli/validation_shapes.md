# Validation shapes

## Enclosed cavities

Generate and validate any nonnegative cavity count (subject to runtime and
memory), with an automatically sized exterior:

```sh
.venv/bin/python -m vorpy.src.geometry.topology_diagnostic --shape sphere_with_cavities --n-cavities 5 --cavity-radius 3 --resolution 32 --output data/topology_audit/cavity5
```

Python callers can use `sphere_with_cavities(n_cavities=5)`. To export generators
without building a network:

```sh
.venv/bin/python -m vorpy.src.geometry.validation_shapes --shape sphere_with_cavities --n-cavities 5
```

Cavity geometry defaults to resolution 32; the other shapes retain their existing
resolution defaults. Output basenames include the count, outer radius, cavity
radius and resolution. Set `--outer-radius` to request a fixed body size; an
undersized body raises an error rather than producing touching cavities.

For N>1 the cavity centers form a deterministic regular polygon, with nearest
spacing `6*cavity_radius`. Its circumradius is `3*r/sin(pi/N)`. Each excluded
center has a selected shell of `resolution` generators at distance `2*r`.
Different shell phases avoid aligned junctions. The default outer radius is
at least `placement_radius + 6*r`. Another `N*resolution` selected generators
lie on radius `outer_radius-r`, with `2*N*resolution+1` excluded exterior
supports on radius `outer_radius+r`. All generator radii are equal, so the
faces are planar. Selected support shells are separated from one another and
from the outer shell; the actual Voronoi cells fill the intervening material.
The regular polygon is symmetric and reproducible, though not an optimal
three-dimensional packing for large N. Radius and computational cost grow
with N. The diagnostic scales its maximum vertex radius with the body size;
if VorPy's settings limit is exceeded, reduce the physical cavity radius.

N=0 is a single selected central cell with spherical exterior supports. N=1
preserves the reviewed geometry: 32 selected generators at radius 6, one
excluded center, and 65 exterior generators at radius 14, for defaults r=3
and R=10. These are faceted approximations to spheres. Mean curvature
metadata is only the smooth reference `4*pi*(R-N*r)`.

The diagnostic measures the actual assembled boundary. It checks selected
connectivity and completeness, partitions retained faces, and sums measured
face, edge and vertex Gaussian contributions per component. Analytic targets
are used only after measurement. Winding checks use the solid-outward direction
from selected to excluded generators. Radial normal checks are relative to
each component's vertex centroid, so translated cavities are handled correctly.
Interior components have negative signed volume and positive integrated Gaussian
curvature. No curvature mathematics or signs are changed.

Verified results at resolution 32:

| Cavities | Boundary components | χ | Expected G | Measured G | Error |
|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 2 | 4π | 12.566370614359126 | -4.62e-14 |
| 1 | 2 | 4 | 8π | 25.13274122871821 | -1.35e-13 |
| 2 | 3 | 6 | 12π | 37.69911184307729 | -2.27e-13 |
| 3 | 4 | 8 | 16π | 50.265482457436306 | -3.84e-13 |
| 4 | 5 | 10 | 20π | 62.83185307179524 | -6.25e-13 |
| 5 | 6 | 12 | 24π | 75.39822368615434 | -6.96e-13 |

Every case has one connected selected solid, all selected cells complete, and
closed, manifold, orientable boundaries. Each component has χ=2 and G≈4π.
Automated topology regressions cover N=0,1,2,3; placement tests also cover larger
counts. A further N=8 build at resolution 16 also passes (9 boundary
components, χ=18, G≈36π). An enclosed cavity adds a boundary component (Δχ=+2, ΔG=+4π); a handle
instead changes genus (Δχ=-2, ΔG=-4π).

Reproduce the sequence:

```sh
for n in 0 1 2 3 4 5; do
  .venv/bin/python -m vorpy.src.geometry.topology_diagnostic \
    --shape sphere_with_cavities --n-cavities "$n" --resolution 32 \
    --output "data/topology_audit/cavity$n"
done
```

Each output directory contains `geometry.pdb`, `geometry.xyzr`, `report.json`,
`build.log`, and the assembled boundary as `boundary.obj` and PyMOL
`boundary.py`/`boundary.pml`. For example:

```sh
pymol data/topology_audit/cavity5/boundary.pml
```

The exterior is transparent blue and cavities orange. Toggle `boundary_0` to
hide the exterior; each cavity is its own object. The scripts contain actual
boundary triangles, not just generator spheres. A measured comparison table
is saved at `data/topology_audit/cavity_sequence.csv`.

## Handles

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
