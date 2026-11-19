"""Compatibility entry point; implementation moved to scripts.figure_reorganization.F3.probe_peg_missing_water."""
from pathlib import Path
import sys

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(root))
from scripts.figure_reorganization.F3.probe_peg_missing_water import *  # noqa: F401,F403

if __name__ == '__main__':
    from scripts.figure_reorganization.F3.probe_peg_missing_water import main
    main()
