"""Round 6 bounded multi-target construction and checkpoint acceptance."""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from vorpy.src.boundary import BoundaryConfig, BoundaryMode, SolventMode
from vorpy.src.group import Group
from vorpy.src.interface import Interface
from vorpy.src.io import Session, load_session, save_session
from vorpy.src.output.checkpoint import _load_checkpoint
from vorpy.src.system import System


def _four_group_system(tmp_path, *, order=("AB", "CD"), scheme="aw"):
    # A/B and C/D are two disjoint target pairs.  All four physical balls are
    # retained in each solve universe; the shell generators are boundary-only.
    balls = pd.DataFrame({
        "loc": [(0.0, 0.0, 0.0), (2.2, 0.0, 0.0),
                (0.0, 10.0, 0.0), (2.2, 10.0, 0.0)],
        "rad": [1.5] * 4,
        "mass": [12.0] * 4,
        "num": list(range(4)),
        "res_name": list("ABCD"),
        "name": ["CA"] * 4,
        "res": [None] * 4,
        "chn": [None] * 4,
    })
    system = System(
        spheres=balls,
        name="round6",
        output_directory=str(tmp_path / "system"),
        boundary_config=BoundaryConfig(mode=BoundaryMode.SHELL),
        solvent_mode=SolventMode.MOLECULE_ONLY,
    )
    groups = {}
    for name, index in zip("ABCD", range(4)):
        group = Group(
            system,
            name=name,
            group_id=name,
            atoms=[index],
            make_net=False,
            build_net=False,
            net_type=scheme,
            solvent_mode=SolventMode.MOLECULE_ONLY,
        )
        groups[name] = group
        system.groups.append(group)

    pairs = {"AB": (groups["A"], groups["B"]),
             "CD": (groups["C"], groups["D"])}
    system.make_interfaces(
        [pairs[key] for key in order],
        analyze_waters=False,
        calculate_curvature=False,
    )
    return system


def _network_signature(interface):
    network = interface.net
    physical_ids = tuple(
        int(value) for value in network.generator_system_indices
        if value is not None
    )
    return {
        "interface_id": interface.interface_id,
        "network_id": network.archive_id,
        "scheme": (network.settings or {}).get("net_type"),
        "source_ids": physical_ids,
        "counts": tuple(len(getattr(network, key)) for key in
                         ("balls", "verts", "edges", "surfs")),
        "area": float(network.surfs["sa"].sum()) if len(network.surfs) else 0.0,
    }


@pytest.mark.timeout(60)
@pytest.mark.parametrize("scheme", ("aw", "pow", "prm"))
def test_two_disjoint_targets_have_full_source_and_order_independent_geometry(tmp_path, scheme):
    first = _four_group_system(tmp_path / "first", order=("AB", "CD"), scheme=scheme)
    second = _four_group_system(tmp_path / "second", order=("CD", "AB"), scheme=scheme)

    first_by_id = {item["interface_id"]: item for item in
                   map(_network_signature, first.ifaces)}
    second_by_id = {item["interface_id"]: item for item in
                    map(_network_signature, second.ifaces)}
    assert set(first_by_id) == {"A__B", "C__D"}
    assert set(first_by_id) == set(second_by_id)
    assert len({item["network_id"] for item in first_by_id.values()}) == 2
    for interface_id in first_by_id:
        left = first_by_id[interface_id]
        right = second_by_id[interface_id]
        assert left["scheme"] == right["scheme"] == scheme
        assert left["source_ids"] == (0, 1, 2, 3)
        assert right["source_ids"] == (0, 1, 2, 3)
        assert left["counts"] == right["counts"]
        assert left["area"] == pytest.approx(right["area"])


@pytest.mark.timeout(60)
def test_multi_target_session_reloads_and_exports_without_build(tmp_path):
    system = _four_group_system(tmp_path / "build")
    session = Session.from_network(system.ifaces[0].net)
    archive = tmp_path / "multi-session.vpy"
    save_session(session, archive, fail_if_exists=True)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()

    root = Path(__file__).resolve().parents[2]
    child = r'''
import hashlib
import json
import sys
from pathlib import Path

from vorpy.src.network import Network
from vorpy.src.io import load_session
from vorpy.src.interface.export import interface_exports

archive = Path(sys.argv[1])
export_root = Path(sys.argv[2])
def forbidden(*args, **kwargs):
    raise AssertionError("Network.build() called during multi-target replay")
Network.build = forbidden
session = load_session(archive)
observed = []
for interface in session.systems[0].ifaces:
    directory = export_root / interface.interface_id
    directory.mkdir(parents=True, exist_ok=False)
    interface.dir = str(directory)
    interface_exports(interface, surfs=True, edges=True, verts=True, info=True)
    observed.append({
        "interface_id": interface.interface_id,
        "network_id": interface.net.archive_id,
        "counts": [len(getattr(interface.net, key)) for key in ("balls", "verts", "edges", "surfs")],
        "files": sorted(path.name for path in directory.iterdir()),
    })
print("NEXUS_ROUND6_REPLAY=" + json.dumps({
    "network_count": len(session.networks),
    "interfaces": observed,
}, sort_keys=True))
'''
    result = subprocess.run(
        [sys.executable, "-c", child, str(archive), str(tmp_path / "exports")],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=45,
        env={**os.environ, "PYTHONPATH": str(root)},
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    marker = next(line for line in result.stdout.splitlines()
                  if line.startswith("NEXUS_ROUND6_REPLAY="))
    observed = json.loads(marker.split("=", 1)[1])
    assert observed["network_count"] == 2
    assert {item["interface_id"] for item in observed["interfaces"]} == {"A__B", "C__D"}
    assert all("interface_surfs.off" in item["files"] for item in observed["interfaces"])
    assert all("interface_edges.off" in item["files"] for item in observed["interfaces"])
    assert all("interface_verts.off" in item["files"] for item in observed["interfaces"])
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == digest


def test_multi_target_rejects_duplicate_overlap_and_empty_targets(tmp_path):
    system = _four_group_system(tmp_path / "invalid")
    a, b = system.groups[0], system.groups[1]
    with pytest.raises(ValueError, match="Duplicate interface target identity"):
        system.make_interfaces([(a, b)])

    c = Group(system, name="overlap", group_id="overlap", atoms=[0],
              make_net=False, build_net=False,
              solvent_mode=SolventMode.MOLECULE_ONLY)
    d = Group(system, name="other", group_id="other", atoms=[1],
              make_net=False, build_net=False,
              solvent_mode=SolventMode.MOLECULE_ONLY)
    with pytest.raises(ValueError, match="must be disjoint"):
        system.make_interfaces([(c, c)])
    empty = Group.__new__(Group)
    empty.name = "empty"
    empty.group_id = "empty"
    empty.ball_ndxs = []
    with pytest.raises(ValueError, match="must contain at least one"):
        system.make_interfaces([(empty, d)])


def test_public_cli_plan_accepts_two_ordered_interface_targets():
    from vorpy.src.command.calculation_plan import build_calculation_plan

    plan = build_calculation_plan([
        "input.pdb",
        "-g", "c", "A", "-g", "c", "B",
        "-g", "c", "C", "-g", "c", "D",
        "-i", "g", "c", "A", "g", "c", "B",
        "-i", "g", "c", "C", "g", "c", "D",
        "-c", "system",
    ])
    assert plan.full_system_request is True
    assert [request.operands for request in plan.interface_requests] == [
        (("c", "A"), ("c", "B")),
        (("c", "C"), ("c", "D")),
    ]


@pytest.mark.timeout(90)
def test_public_cli_builds_two_disjoint_interface_targets(tmp_path):
    source = tmp_path / "four_groups.pdb"
    atoms = (
        (1, "A", 1, 0.0, 0.0, 0.0),
        (2, "B", 1, 2.0, 0.0, 0.0),
        (3, "C", 1, 0.0, 2.0, 0.0),
        (4, "D", 1, 0.0, 0.0, 2.0),
    )
    source.write_text(
        "".join(
            f"ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            for serial, chain, residue, x, y, z in atoms
        ) + "END\n",
        encoding="ascii",
    )
    output_root = tmp_path / "cli-output"
    archive = tmp_path / "cli-session.vpy"
    command = [
        sys.executable, "-m", "vorpy", str(source),
        "-g", "c", "A", "-g", "c", "B", "-g", "c", "C", "-g", "c", "D",
        "-i", "g", "c", "A", "g", "c", "B",
        "-i", "g", "c", "C", "g", "c", "D",
        "-c", "system", "-s", "nt", "aw", "and", "mv", "5",
        "--profile", "geometry",
        "--boundary", "shell", "-e", "none", "-e", "dir", str(output_root),
        "--save-network", str(archive), "--fail-if-exists",
    ]
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        command,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=80,
        env={**os.environ, "PYTHONPATH": str(root)},
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert archive.is_file()
    assert "c_A_c_B_interface" in result.stdout
    assert "c_C_c_D_interface" in result.stdout
    checkpoint = _load_checkpoint(archive)
    assert len(checkpoint.document.interfaces) == 2
    assert {
        tuple(item.display_orientation) for item in checkpoint.document.interfaces
    } == {("c_A", "c_B"), ("c_C", "c_D")}
    # The canonical checkpoint has one shared source geometry payload and two
    # independently identified materialized interface targets.  Their source
    # views are replayable from this payload; no second source solve is needed.
    assert len(checkpoint.document.networks) == 1
    assert (output_root / source.stem).is_dir()
