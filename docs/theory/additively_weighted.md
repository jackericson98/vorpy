# Additively Weighted Voronoi

For a generator centered at `c_i` with radius `r_i`:

```text
d_AW(x, i) = ||x - c_i|| - r_i
```

A point is assigned to the generator whose surface is closest. This makes the construction directly sensitive to atomic radii. Additively weighted boundaries are not restricted to planes.

VorPy name: `aw`.
