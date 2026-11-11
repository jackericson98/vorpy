# Changelog

## 3.5.3 (provisional release candidate)

- Added explicit export dependency validation for geometry-only execution
  profiles; incompatible analysis bundles fail cleanly after checkpoint save.
- Preserved source-independent checkpoint reload and cache-only export paths.
- Preserved bounded AW, POW, and PRM synthetic geometry, ordered-interface,
  Results, and Workbench regressions.
- Reduced redistributable package data to the synthetic `HO_cube.pdb` example;
  research structures, trajectories, benchmark datasets, and generated output
  remain external files.
- Updated the Workbench example action and CLI documentation for the reduced
  package resource set.

Known limitations remain documented in the release notes; this entry does not
claim molecular-scale scientific certification or complete open-interface
curvature validation.
