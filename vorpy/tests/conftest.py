# vorpy/tests/conftest.py

import pytest
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning, module=r"vorpy\.src\.calculations\.calcs")
warnings.filterwarnings("ignore", category=RuntimeWarning, module=r"vorpy\.src\.calculations\.plane")


@pytest.fixture(scope="session")
def tol():
    # Tight defaults you can loosen if needed for bigger scales
    return {"rel": 1e-12, "abs": 1e-12}
