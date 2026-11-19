"""Compatibility alias for vorpy.src.analyze.Part_II_Molecular_Analysis.F6_Proteins_Interfaces_Multiscale.historical.Ligand_Interface_Residues."""
from pathlib import Path
import importlib
import runpy
import sys

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(root))

if __name__ == '__main__':
    runpy.run_module('vorpy.src.analyze.Part_II_Molecular_Analysis.F6_Proteins_Interfaces_Multiscale.historical.Ligand_Interface_Residues', run_name='__main__')
else:
    _module = importlib.import_module('vorpy.src.analyze.Part_II_Molecular_Analysis.F6_Proteins_Interfaces_Multiscale.historical.Ligand_Interface_Residues')
    globals().update({name: value for name, value in vars(_module).items()
                      if not name.startswith('__')})
