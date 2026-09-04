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
python -m pip install -e ".[gui]"
```

This editable installation binds the `vorpy` command and the launchers in the
repository root to the current checkout. The included launchers are
`vorpy-linux.desktop` for Linux, `vorpy-mac.command` for macOS, and
`vorpy-windows.bat` for Windows.

## Verification
```bash
python -c "import vorpy; print(vorpy)"
pytest -q
```

See [Testing and Compatibility](development/testing.md).
