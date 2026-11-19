"""Compatibility alias for scripts.figure_reorganization.F6.plan_figure3b_recovery."""
from pathlib import Path
import importlib
import runpy
import sys

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(root))

if __name__ == '__main__':
    runpy.run_module('scripts.figure_reorganization.F6.plan_figure3b_recovery', run_name='__main__')
else:
    _module = importlib.import_module('scripts.figure_reorganization.F6.plan_figure3b_recovery')
    globals().update({name: value for name, value in vars(_module).items()
                      if not name.startswith('__')})
