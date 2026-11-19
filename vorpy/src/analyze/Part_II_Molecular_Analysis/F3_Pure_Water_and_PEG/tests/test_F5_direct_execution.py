"""Legacy entry points must import safely without launching generators."""
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1]
SMOKE_SYMBOLS = {
    'F3_PEG_prototypes.py': 'collect_observations',
    'F3_round2_validation.py': 'descriptive_trend',
    'F3F_water_contact_prototype.py': 'paired_differences',
    'F5_common.py': 'read_pdb',
    'F5_config.py': 'Config',
    'F5_peg_dataset.py': 'screened_runs',
    'F5_run_all.py': 'discover',
    'F5A_exposed_vs_buried.py': '_self_test',
    'F5A_PEG_chemistry_vs_hydration.py': 'load_frame',
    'F5B_heparin_chemistry.py': 'classify_chemistry',
    'F5B_water_burial.py': 'self_check',
    'F5C_rerun_missing.py': 'CASES',
    'F5C_water_contacts.py': 'water_contacts',
}


@pytest.mark.parametrize('name,symbol', sorted(SMOKE_SYMBOLS.items()))
def test_legacy_entry_point_imports_without_running_generator(name, symbol, tmp_path):
    script = SCRIPT_DIR / name
    command = [
        sys.executable, '-I', '-c',
        'import runpy, sys; namespace = runpy.run_path(sys.argv[1], run_name="f5_smoke"); '
        'assert sys.argv[2] in namespace',
        str(script), symbol,
    ]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr


def test_legacy_runner_help_is_nondestructive(tmp_path):
    result = subprocess.run(
        [sys.executable, '-I', str(SCRIPT_DIR / 'F5_run_all.py'), '--help'],
        cwd=tmp_path, capture_output=True, text=True, timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert '--output-root' in result.stdout
