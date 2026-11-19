"""Compatibility entry point; implementation moved to scripts.figure_reorganization.F3.validate_peg_layer1_bundle."""
from pathlib import Path
import sys

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(root))
from scripts.figure_reorganization.F3.validate_peg_layer1_bundle import *  # noqa: F401,F403

if __name__ == '__main__':
    from scripts.figure_reorganization.F3.validate_peg_layer1_bundle import main
    main()
