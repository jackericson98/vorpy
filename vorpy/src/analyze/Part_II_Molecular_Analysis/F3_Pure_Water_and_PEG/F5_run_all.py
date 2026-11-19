"""Compatibility entry point; implementation moved to vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_run_all."""
from pathlib import Path
import sys

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(root))
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_run_all import *  # noqa: F401,F403

if __name__ == '__main__':
    from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_run_all import main
    main()
