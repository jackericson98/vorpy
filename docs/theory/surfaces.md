# Surface Representation

VorPy represents solved cell boundaries using surfaces that can be triangulated for analysis and visualization. These surfaces support surface area, enclosed volume, curvature, interface area, and visualization meshes.

A solved network is organized around cells, surfaces, edges, and vertices. A shared surface indicates that the corresponding cells are Voronoi neighbors.

**TODO:** Verify the exact interpretation and units of the current `sr` surface-resolution setting.
