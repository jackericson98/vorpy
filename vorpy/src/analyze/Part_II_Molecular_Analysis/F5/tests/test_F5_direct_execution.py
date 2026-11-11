"""PyCharm-style file execution must not depend on the working directory."""
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPTS = sorted(Path(__file__).resolve().parents[1].glob('*.py'))


@pytest.mark.parametrize('script', SCRIPTS, ids=lambda path: path.name)
def test_direct_execution_outside_repository(script, tmp_path):
    # Isolated Python ignores PYTHONPATH and PyCharm's injected project roots.
    command = [sys.executable, '-I', str(script)]
    checks = {'F5A_exposed_vs_buried.py': '_self_test', 'F5B_water_burial.py': 'self_check'}
    if script.name in checks:
        # Check imports and calculations without loading datasets or overwriting plots.
        command = [sys.executable, '-I', '-c',
            f'import runpy, sys; runpy.run_path(sys.argv[1], run_name="f5_direct_test")["{checks[script.name]}"]()',
            str(script)]
    elif script.name == 'F5_run_all.py':
        command.append('--help')
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr
