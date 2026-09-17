# Install the curvature foundation

Copy the files from this package into the root of the current VorPy repository,
preserving their paths.

New files:

- `vorpy/src/calculations/surface_geometry.py`
- `vorpy/src/calculations/edge_geometry.py`
- `vorpy/src/network/boundary_geometry.py`
- `vorpy/tests/calculations/test_surface_geometry.py`
- `vorpy/tests/calculations/test_edge_geometry.py`
- `vorpy/tests/network/test_boundary_geometry.py`
- `VORPY_CURVATURE_AUDIT.md`

Updated file:

- `vorpy/src/calculations/__init__.py`

Run from the repository root:

```bash
pytest -q vorpy/tests/calculations/test_surface_geometry.py
pytest -q vorpy/tests/calculations/test_edge_geometry.py
pytest -q vorpy/tests/network/test_boundary_geometry.py
pytest -q
```

This stage is additive and should not change any existing VorPy result or log.
