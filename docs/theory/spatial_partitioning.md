# Spatial Partitioning

VorPy represents atoms or other spherical objects as generators of a three-dimensional spatial decomposition. Changing the distance rule changes the resulting cells, boundaries, neighbors, volumes, and interfaces.

VorPy emphasizes three constructions: additively weighted Voronoi, power/Laguerre, and primitive/unweighted Voronoi.

For molecular systems containing atoms of different effective sizes, these schemes can assign different local volumes, surface areas, neighbor identities, and interfaces. VorPy is designed to make those assumptions explicit and directly comparable.
