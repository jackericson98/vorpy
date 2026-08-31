# Testing and Python Compatibility

VorPy is scientific geometry software, so compatibility should be validated numerically and functionally rather than only by import success.

## Testing Layers

### Unit tests
Individual functions and geometric primitives.

### Numerical regression tests
Known molecular systems with expected counts and numerical outputs such as cells, vertices, edges, surfaces, selected volumes, surface areas, neighbors, interface areas, and curvature values. Floating-point comparisons should use tolerances.

### Integration tests
Complete subsystem workflows such as:

```text
parser -> system -> group -> network -> analysis -> export
```

All five input parsers should have integration coverage.

### End-to-end tests
Representative VorPy commands should complete, create expected outputs, preserve important numerical values within tolerance, and exercise AW, POW, and PRM workflows.

## Initial Python Matrix
- Python 3.10
- Python 3.11
- Python 3.12
- Python 3.13

Python 3.14 can initially be experimental.

## Platforms
- Linux
- Windows
- macOS

## CI Strategy
Every commit should run unit tests, small integration tests, and a small end-to-end smoke test across supported Python versions. Main/release validation should additionally run the full regression suite, all partitioning schemes, all parsers, interface workflows, representative exports, and larger reference systems where runtime permits.

A Python version should only be added to the README compatibility badge after the required suite passes consistently.

**TODO:** Review the existing `tests/` directory and convert this strategy into the actual GitHub Actions matrix.
