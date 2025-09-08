# tests/conftest.py

import pytest


@pytest.fixture(scope="session")
def tol():
    # Tight defaults you can loosen if needed for bigger scales
    return {"rel": 1e-12, "abs": 1e-12}
