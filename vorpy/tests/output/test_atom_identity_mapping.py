"""Identity contracts; topology IDs need not equal source PDB serials."""
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.command.group import get_group_spheres
from vorpy.src.group import Group
from vorpy.src.group.export import _surface_neighbor_system_indices
from vorpy.src.inputs.pdb import read_pdb
from vorpy.src.output.pdb import make_pdb_line, write_pdb
from vorpy.src.system.system import System


@pytest.mark.parametrize("scheme", ["aw", "pow", "prm"])
@pytest.mark.parametrize("mode,local", [("SOLVENT_CONTEXT", 2), ("MOLECULE_ONLY", 1)])
def test_filtered_and_regrouped_identity(tmp_path, monkeypatch, scheme, mode, local):
    source = tmp_path / "identity.pdb"
    records = [make_pdb_line(ser_num=serial, name=name, res_name=res,
                             res_seq=i + 1, x=i * 2, elem="O")
               for i, (serial, name, res) in enumerate(
                   [(17, "O", "ALA"), (99, "O", "HOH"), (407, "OXT", "GLY")])]
    source.write_text("HEADER identity\n" + records[0] + "TER\n" + "".join(records[1:]))
    loaded = SimpleNamespace(files={"base_file": str(source), "vpy_dir": None, "dir": None},
                             type="mol", print_actions=False)
    read_pdb(loaded)
    assert get_group_spheres(loaded.balls, ["2"]) == [2]
    assert loaded.balls.iloc[2]["name"] == "OXT"
    system = System(spheres=loaded.balls, name="identity", make_dir=False, solvent_mode=mode)
    group = Group(system, atoms=[2], make_net=True, net_type=scheme, print_metrics=False)
    assert group.ball_ndxs == [2]
    assert int(group.net.balls.iloc[local]["system_num"]) == 2
    assert tuple(group.net.balls.iloc[local]["loc"]) == tuple(loaded.balls.iloc[2]["loc"])
    # A known face tests endpoint attribution independently of geometry solving.
    group.net.surfs = pd.DataFrame({"balls": [[0, local]], "sa": [1.0], "contact_area": [0.0]})
    assert _surface_neighbor_system_indices(group) == [0]
    monkeypatch.chdir(tmp_path)
    write_pdb(group.ball_ndxs, "group_atoms", loaded)
    write_pdb(_surface_neighbor_system_indices(group), "surr_atoms", loaded)
    assert records[2] in (tmp_path / "group_atoms.pdb").read_text()
    assert records[0] in (tmp_path / "surr_atoms.pdb").read_text()
    regrouped = Group(system, atoms=[0, 2], make_net=True, net_type=scheme, print_metrics=False)
    assert int(regrouped.net.balls.iloc[local]["system_num"]) == 2


def test_hammerhead_cli_selection_matches_original_record(tmp_path, monkeypatch):
    source = Path(__file__).resolve().parents[2] / "data" / "hammerhead.pdb"
    system = SimpleNamespace(files={"base_file": str(source), "vpy_dir": None, "dir": None},
                             type="mol", print_actions=False)
    read_pdb(system)
    assert get_group_spheres(system.balls, ["51505"]) == [51505]
    atom = system.balls.iloc[51505]
    assert (atom["name"], atom["res_name"], atom["res_seq"], atom["chain_name"]) == ("O3'", "A", 21, "A")
    assert tuple(atom["loc"]) == (73.85, 71.43, 21.21)
    monkeypatch.chdir(tmp_path)
    write_pdb([51505], "group_atoms", system)
    record = next(line for line in (tmp_path / "group_atoms.pdb").read_text().splitlines() if line.startswith("ATOM"))
    assert int(record[6:11]) == 51506
    assert record[12:16].strip() == atom["name"]
