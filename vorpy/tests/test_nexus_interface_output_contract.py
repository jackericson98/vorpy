from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from vorpy.src.interface import Interface
from vorpy.src.interface.export import interface_exports
from vorpy.src.network.views import NetworkView
from vorpy.src.output.export_scope import export_scope
from vorpy.src.output.pdb import make_pdb_line
from vorpy.tests.network.test_views import _network


def _register(group, interface_id, interface=None, other_group=None, side=None,
              interface_name=None):
    group.interfaces[interface_id] = interface
    group.interface_metadata[interface_id] = {
        "interface_id": interface_id,
        "interface_name": interface_name,
        "group_name": group.name,
        "other_group_name": getattr(other_group, "name", None),
        "side": side,
    }


def _contract_case(tmp_path):
    source = _network(scheme="aw")
    source.solvent_mode = "FULL_SOLVENT"
    source.settings["solvent_mode"] = "FULL_SOLVENT"
    source.balls.loc[2, "stable_id"] = "water-1"
    source.balls.loc[3, "stable_id"] = "water-2"
    source.surfs.at[0, "balls"] = [0, 1, 2]

    water_surface = source.surfs.iloc[1].copy()
    water_surface["balls"] = [2, 3]
    dna_water_surface = source.surfs.iloc[1].copy()
    dna_water_surface["balls"] = [0, 2]
    source.surfs = pd.concat(
        [source.surfs, pd.DataFrame([water_surface, dna_water_surface], index=[2, 3])]
    )

    rows = []
    for index in range(14):
        rows.append({
            "num": index,
            "loc": (float(index), 0.0, 0.0),
            "name": f"C{index}",
            "res_name": "DNA" if index == 10 else "MOL",
            "res_seq": 1,
            "chain_name": "A" if index == 10 else "C",
            "element": "C",
        })
    system_balls = pd.DataFrame(rows)
    pdb = tmp_path / "synthetic_interface.pdb"
    pdb.write_text(
        "".join(
            make_pdb_line(
                ser_num=index + 1,
                name=f"C{index}",
                res_name="DNA" if index == 10 else "MOL",
                chain="A" if index == 10 else "C",
                res_seq=1,
                x=float(index),
                y=0.0,
                z=0.0,
                elem="C",
            )
            for index in range(14)
        ) + "END\n",
        encoding="ascii",
    )

    group_a = SimpleNamespace(
        name="dna", group_id="dna", ball_ndxs=[10],
        settings=dict(source.settings), net=None, interfaces={},
        interface_metadata={},
    )
    group_b = SimpleNamespace(
        name="ligand", group_id="ligand", ball_ndxs=[11],
        settings=dict(source.settings), net=None, interfaces={},
        interface_metadata={},
    )
    group_a.register_interface = lambda **kwargs: _register(group_a, **kwargs)
    group_b.register_interface = lambda **kwargs: _register(group_b, **kwargs)
    system = SimpleNamespace(
        name="synthetic-interface",
        frame_index=1,
        frame_count=1,
        type="pdb",
        groups=[group_a, group_b],
        ifaces=[],
        balls=system_balls,
        round_to=3,
        files={"dir": str(tmp_path), "base_file": str(pdb)},
    )
    group_a.sys = system
    group_b.sys = system

    interface = Interface(system, group_a, group_b)
    interface.net = source
    interface.source_network = source
    interface.source_view = NetworkView.for_interface(source, [10], [11])
    interface.group1_indices = {10}
    interface.group2_indices = {11}
    interface.group2_name = "ligand"
    interface.dir = str(tmp_path / "interface")
    system.ifaces.append(interface)
    return source, interface


def _export_case(tmp_path):
    source, interface = _contract_case(tmp_path)
    interface_exports(
        interface,
        atoms=True,
        surfs=True,
        edges=True,
        verts=True,
        logs=True,
        info=True,
    )
    return source, interface, Path(interface.dir)


def test_nexus_interface_public_artifact_contract(tmp_path):
    source, interface, directory = _export_case(tmp_path)
    expected = {
        "interface_surfs.off",
        "interface_edges.off",
        "interface_verts.off",
        "interface_info.txt",
        "group_a_atoms.pdb",
        "group_b_atoms.pdb",
        "surr_atoms.pdb",
        "surr_residues.pdb",
    }
    missing = sorted(name for name in expected if not (directory / name).is_file())
    assert not missing, (
        "NEXUSA interface output contract is incomplete; missing: "
        + ", ".join(missing)
    )

    assert not (directory / "visualization_manifest.json").exists()
    assert not list(directory.glob("*.geometry.vpy"))
    assert len(source.balls) == 5
    assert export_scope(interface).source_ids("surface") == (0,)


def test_nexus_interface_current_physical_scope_and_pdb_lineage(tmp_path):
    source, interface, directory = _export_case(tmp_path)
    scope = export_scope(interface)

    assert scope.source_ids("surface") == (0,)
    assert scope.source_ids("edge") == (0, 1, 2, 3)
    assert scope.source_ids("vertex") == (0, 1, 2, 3)
    assert len(source.surfs) == 4

    # Until the public aliases are integrated, inspect the existing physical
    # files to prove the scoped geometry and source PDB identity independently.
    assert (directory / "interface_surfs.off").read_text(encoding="utf-8").splitlines()[1] == "4 1 0"
    assert " C10 " in (directory / "group_a_atoms.pdb").read_text(encoding="ascii")
    assert " C11 " in (directory / "group_b_atoms.pdb").read_text(encoding="ascii")
    assert "COMPLETE" in (directory / "interface_surfaces.csv").read_text(encoding="utf-8")
    assert "Surface Area (Å²)" in (directory / "interface_surfaces.csv").read_text(encoding="utf-8")
    assert "12" in (directory / "surr_atoms.pdb").read_text(encoding="ascii")
    assert "12" in (directory / "surr_residues.pdb").read_text(encoding="ascii")
