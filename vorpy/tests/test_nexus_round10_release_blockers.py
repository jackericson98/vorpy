"""Round-10 regressions for release-blocking cache and project contracts."""

import json
import os
import subprocess
import sys
from pathlib import Path

from vorpy.src.io import save_session
from vorpy.tests.io.test_scientific_session import cached_session


ROOT = Path(__file__).resolve().parents[2]


def test_fresh_process_preserves_complete_cached_provenance(
    cached_session, tmp_path
):
    archive = save_session(cached_session, tmp_path / "provenance.vpy")
    code = r'''
import json
import sys
from vorpy.src.io import load_session
from vorpy.src.output.visualization import _cached_interface_result
from vorpy.src.network import Network

def fail_build(*args, **kwargs):
    raise AssertionError("cache replay called Network.build()")

Network.build = fail_build
session = load_session(sys.argv[1])
iface = next(item for item in session.active_network.sys.ifaces
             if item.net is session.active_network)
expected = dict(iface.geometry_analysis.archive_provenance)
actual = dict(_cached_interface_result(iface).provenance.details)
assert actual == expected, json.dumps({"expected": expected, "actual": actual}, default=str)
print(json.dumps(actual, sort_keys=True, default=str))
'''
    completed = subprocess.run(
        [sys.executable, "-c", code, str(archive)],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
