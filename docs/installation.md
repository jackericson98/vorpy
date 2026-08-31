# Installation

## PyPI
```bash
pip install vorpy3
```

## Conda
```bash
conda install vorpy3
```

## Python Version
VorPy currently targets Python 3.10+. Explicit compatibility claims for additional Python versions should be based on automated tests and end-to-end molecular regression calculations rather than import success alone.

## Development Installation
```bash
git clone https://github.com/jackericson98/vorpy.git
cd vorpy
pip install -e .
```

## Verification
```bash
python -c "import vorpy; print(vorpy)"
pytest -q
```

See [Testing and Compatibility](development/testing.md).
