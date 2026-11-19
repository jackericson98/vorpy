"""Compatibility entry point; implementation moved to vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F5A_exposed_vs_buried."""
from pathlib import Path
import sys

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(root))
from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F5A_exposed_vs_buried import *  # noqa: F401,F403

from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F5A_exposed_vs_buried import _self_test

if __name__ == '__main__':
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F5A_exposed_vs_buried import main
    main()
