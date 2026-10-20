import os
import os
import subprocess
import sys
from pathlib import Path

from vorpy.src.output.checkpoint import _load_checkpoint


def _terminal_lines(output):
    """Interpret carriage-return progress as a single terminal line."""
    lines = []
    current = []
    for char in output:
        if char == '\r':
            current = []
        elif char == '\n':
            lines.append(''.join(current))
            current = []
        else:
            current.append(char)
    if current:
        lines.append(''.join(current))
    return lines


def test_edta_default_cli_exports_clean_aw_bundle(tmp_path):
    output_parent = tmp_path / "output"
    root = Path(__file__).resolve().parents[3]
    completed = subprocess.run(
        [sys.executable, "vorpy", "edta", "-e", "dir", str(output_parent)],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=120,
        env={"PYTHONPATH": str(root), **os.environ},
    )

    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    output = completed.stdout
    assert 'network built - ' in output
    assert '32 complete cells' in output
    assert "Analysis - 100.00 %" not in output
    assert "Vertex Gaussian - 100%" not in output
    assert "Export - 100%" not in output
    assert "Solving AW network..." in output
    assert "Exporting results..." in output
    export_line = f"Files exported to: {output_parent / 'EDTA'}"
    assert output.rstrip().endswith(export_line)

    terminal_lines = _terminal_lines(output)
    solve_line = terminal_lines.index("Solving AW network...")
    assert terminal_lines[solve_line - 4:solve_line] == [
        "",
        "EDTA loaded",
        "  612 atoms, 195 residues, 1 chains, 193 waters",
        "",
    ]
    assert terminal_lines[solve_line + 1].startswith('"EDTA" network built - ')
    assert terminal_lines[solve_line + 2] == ""
    assert terminal_lines[solve_line + 3] == "Exporting results..."
    assert terminal_lines[solve_line + 4] == export_line
    assert terminal_lines[solve_line + 5] == ""

    directory = output_parent / "EDTA"
    system_dir = directory / "system"
    network_dir = directory / "groups" / "EDTA"
    assert (system_dir / "system_info.txt").is_file()
    system_info = (system_dir / "system_info.txt").read_text(encoding="utf-8")
    assert "Molecular Atoms:      32" in system_info
    assert "Waters:                193" in system_info
    assert "Total Atoms:            612" in system_info
    assert "Active Generators:      612" in system_info
    assert "Excluded Water:         0" in system_info
    assert (network_dir / "surfs.off").stat().st_size > 0
    # The public default is the small tier; full edge/vertex meshes remain in
    # the verified checkpoint or an explicit medium/large request.
    assert not (network_dir / "edges.off").exists()
    assert not (network_dir / "verts.off").exists()
    assert not list(directory.rglob("nonpolar_interface_geometry"))
    assert not (directory / "EDTA.vpy.timing.json").exists()

    archive = directory / "EDTA_checkpoint.vpy"
    checkpoint = _load_checkpoint(archive)
    assert checkpoint.document.source.filename == "EDTA.pdb"
    assert checkpoint.payloads
    system_states = [
        state for state in checkpoint.document.networks
        if state.scope == "group" or state.scope == "system"
    ]
    assert system_states
    state = max(system_states, key=lambda item: sum(item.completeness.values()))
    assert state.completeness["balls"] == 612
    assert state.completeness["verts"] > 0
    assert state.completeness["edges"] > 0
    assert state.completeness["surfs"] > 0
