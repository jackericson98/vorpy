from pathlib import Path

from vorpy.src.command.calculation_plan import build_export_plan
from vorpy.src.command.command_export import _run_component_plan
from vorpy.src.output.checkpoint import _load_checkpoint
from vorpy.src.output.sys import export_sys
from vorpy.tests.test_nexus_interface_output_contract import _contract_case


def _run_tier(root, preset):
    root.mkdir(parents=True, exist_ok=True)
    _, interface = _contract_case(root)
    system = interface.sys
    output = root / "out"
    system.files["dir"] = str(output)
    system.update_progress = lambda **kwargs: None
    system.exports = lambda **kwargs: export_sys(system, **kwargs)
    plan = build_export_plan([(preset,)], scope="interface")
    system._export_plan = plan
    _run_component_plan(system, plan)
    return system, output


def test_canonical_interface_tiny_is_source_and_checkpoint_complete(tmp_path, monkeypatch):
    from vorpy.src.network import Network

    def forbidden_build(*args, **kwargs):
        raise AssertionError("cache-only export rebuilt Network")

    monkeypatch.setattr(Network, "build", forbidden_build)
    system, output = _run_tier(tmp_path, "tiny")

    assert (output / "synthetic_interface.pdb").read_bytes() == (
        tmp_path / "synthetic_interface.pdb"
    ).read_bytes()
    assert (output / "system" / "system_info.txt").is_file()
    assert (output / "groups" / "A" / "A.pdb").is_file()
    assert (output / "groups" / "B" / "B.pdb").is_file()
    interface_dir = output / "interfaces" / "dna__ligand"
    assert (interface_dir / "interface_atoms.pdb").is_file()
    assert (interface_dir / "interface_edges.off").is_file()
    assert (interface_dir / "interface_info.txt").is_file()
    assert not (interface_dir / "logs.csv").exists()
    assert not (interface_dir / "logs.provenance.json").exists()
    checkpoint = _load_checkpoint(output / f"{system.name}_checkpoint.vpy")
    assert checkpoint.document.source.filename == "synthetic_interface.pdb"
    assert checkpoint.payloads


def test_canonical_interface_tiers_are_cumulative(tmp_path):
    inventories = {}
    checkpoint_sizes = {}
    for preset in ("tiny", "small", "medium", "large", "all"):
        system, output = _run_tier(tmp_path / preset, preset)
        inventories[preset] = {
            str(path.relative_to(output))
            for path in output.rglob("*")
            if path.is_file()
        }
        checkpoint_sizes[preset] = (
            output / f"{system.name}_checkpoint.vpy"
        ).stat().st_size

    assert inventories["tiny"] <= inventories["small"]
    assert inventories["small"] <= inventories["medium"]
    assert inventories["medium"] <= inventories["large"]
    assert inventories["large"] <= inventories["all"]
    assert all(size > 0 for size in checkpoint_sizes.values())
    assert not any("logs.csv" in item for item in inventories["all"])
    assert not any("logs.provenance.json" in item for item in inventories["all"])
