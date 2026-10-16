import base64
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
WATCHDOG = ROOT / "scripts" / "nexus_pow_watchdog.ps1"
POWERSHELL = shutil.which("powershell.exe")

pytestmark = pytest.mark.skipif(
    os.name != "nt" or POWERSHELL is None,
    reason="the integration watchdog is Windows PowerShell specific",
)


def _python_arguments(source):
    encoded = base64.b64encode(source.encode("utf-8")).decode("ascii")
    return (
        f"-c \"import base64;exec(base64.b64decode('{encoded}'))\" "
        "--fail-if-exists"
    )


def _run_watchdog(
    tmp_path,
    source,
    *,
    max_wall=5,
    max_output=1_073_741_824,
    require_archive=False,
    create_output=False,
):
    output = tmp_path / "output"
    archive = output / "result.vpy"
    if create_output:
        output.mkdir()
    command = [
        POWERSHELL,
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-File",
        str(WATCHDOG),
        "-FilePath",
        sys.executable,
        "-Arguments",
        _python_arguments(source),
        "-OutputRoot",
        str(output),
        "-ArchivePath",
        str(archive),
        "-MaxWallSeconds",
        str(max_wall),
        "-MaxCpuSeconds",
        "1800",
        "-MaxPeakRssBytes",
        str(4 * 1024**3),
        "-MaxTotalOutputBytes",
        str(max_output),
        "-MaxArchiveMembers",
        "100000",
        "-PollMilliseconds",
        "50",
    ]
    if require_archive:
        command.append("-RequireArchive")
    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=20,
    )
    payload = None
    for line in reversed(completed.stdout.splitlines()):
        try:
            payload = json.loads(line)
            break
        except json.JSONDecodeError:
            continue
    return completed, payload, output, archive


def test_watchdog_nonexistent_output_reports_zero(tmp_path):
    completed, payload, output, _ = _run_watchdog(
        tmp_path, "import time; time.sleep(0.6)"
    )

    assert completed.returncode == 0, completed.stderr
    assert payload["status"] == "PROCESS_EXITED"
    assert payload["total_output_bytes"] == 0
    assert not output.exists()


def test_watchdog_empty_output_reports_zero(tmp_path):
    source = f"import pathlib, time; pathlib.Path({str(tmp_path / 'output')!r}).mkdir(); time.sleep(0.6)"
    completed, payload, output, _ = _run_watchdog(tmp_path, source)

    assert completed.returncode == 0, completed.stderr
    assert payload["status"] == "PROCESS_EXITED"
    assert payload["total_output_bytes"] == 0
    assert output.is_dir()


@pytest.mark.parametrize(
    "files, expected",
    [
        ({"one.bin": 123}, 123),
        ({"nested/one.bin": 123, "nested/deeper/two.bin": 456}, 579),
    ],
)
def test_watchdog_totals_one_and_nested_files(tmp_path, files, expected):
    output = tmp_path / "output"
    writes = "; ".join(
        f"p=pathlib.Path({str(output / name)!r}); p.parent.mkdir(parents=True, exist_ok=True); p.write_bytes(b'x'*{size})"
        for name, size in files.items()
    )
    source = f"import pathlib, time; {writes}; time.sleep(0.8)"
    completed, payload, _, _ = _run_watchdog(tmp_path, source)

    assert completed.returncode == 0, completed.stderr
    assert payload["status"] == "PROCESS_EXITED"
    assert payload["total_output_bytes"] == expected


def test_watchdog_output_limit_violation_aborts(tmp_path):
    output = tmp_path / "output"
    source = (
        f"import pathlib, time; p=pathlib.Path({str(output / 'large.bin')!r}); "
        f"p.parent.mkdir(); p.write_bytes(b'x'*200000); time.sleep(5)"
    )
    completed, payload, _, _ = _run_watchdog(
        tmp_path, source, max_output=100_000
    )

    assert completed.returncode == 20, completed.stderr
    assert payload["status"] == "ABORTED_LIMIT"
    assert payload["trigger"] == "total_output_bytes"
    assert payload["total_output_bytes"] >= 200_000
    assert payload["remaining_process_ids"] == []


def test_watchdog_filesystem_scan_error_is_explicit_and_terminates(tmp_path):
    output = tmp_path / "output"
    source = (
        f"import pathlib, shutil, time; root=pathlib.Path({str(output)!r}); "
        f"root.mkdir(); shutil.rmtree(root); root.write_text('not a directory'); time.sleep(5)"
    )
    completed, payload, _, _ = _run_watchdog(tmp_path, source)

    assert completed.returncode == 22, completed.stderr
    assert payload["status"] == "FAILED_OUTPUT_SCAN"
    assert "not a directory" in payload["output_scan_error"]
    assert payload["remaining_process_ids"] == []


def test_watchdog_wall_time_terminates_parent_and_child_tree(tmp_path):
    source = (
        "import subprocess, sys, time; "
        "subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)']); "
        "time.sleep(30)"
    )
    completed, payload, _, _ = _run_watchdog(
        tmp_path, source, max_wall=1
    )

    assert completed.returncode == 20, completed.stderr
    assert payload["trigger"] == "wall_seconds"
    assert payload["child_process_count"] >= 1
    assert payload["termination_taskkill_exit_code"] == 0
    assert payload["remaining_process_ids"] == []


def test_watchdog_failed_process_and_required_archive_checks(tmp_path):
    failed, failed_payload, _, _ = _run_watchdog(
        tmp_path / "failed", "raise SystemExit(7)"
    )
    assert failed.returncode == 10, failed.stderr
    assert failed_payload["status"] == "PROCESS_FAILED"
    assert failed_payload["failure"] == "process_exit:7"

    missing, missing_payload, _, _ = _run_watchdog(
        tmp_path / "missing", "pass", require_archive=True
    )
    assert missing.returncode == 23, missing.stderr
    assert missing_payload["status"] == "FAILED_ARCHIVE_MISSING_OR_UNREADABLE"


def test_watchdog_captures_uncaught_exception_streams_and_summary(tmp_path):
    source = (
        "print('stdout-marker', flush=True); "
        "raise RuntimeError('uncaught dummy failure')"
    )
    completed, payload, output, _ = _run_watchdog(tmp_path, source)

    assert completed.returncode == 10, completed.stderr
    assert payload["status"] == "PROCESS_FAILED"
    assert payload["exit_code"] != 0
    assert payload["failure"].startswith("process_exit:")
    assert payload["started_at_utc"] < payload["ended_at_utc"]
    assert Path(payload["stdout_path"]).read_text(encoding="utf-8").strip() == "stdout-marker"
    stderr = Path(payload["stderr_path"]).read_text(encoding="utf-8")
    assert "RuntimeError: uncaught dummy failure" in stderr
    assert payload["stdout_bytes"] > 0
    assert payload["stderr_bytes"] > 0
    summary = Path(payload["summary_path"])
    assert summary.is_file()
    assert json.loads(summary.read_text(encoding="utf-8"))["status"] == "PROCESS_FAILED"
    assert output.with_name(output.name + ".watchdog").is_dir()


def test_watchdog_refuses_existing_output_root(tmp_path):
    completed, payload, output, _ = _run_watchdog(
        tmp_path, "raise SystemExit(99)", create_output=True
    )

    assert completed.returncode == 1
    assert payload is None
    assert "Refusing to overwrite an existing output root" in completed.stderr
    assert output.is_dir()
