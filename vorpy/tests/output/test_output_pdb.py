from pathlib import Path
import pandas as pd
import pytest
from types import SimpleNamespace

from vorpy.src.output.pdb import write_pdb


class DummySystem:
    def __init__(self, base_file, out_dir):
        self.name = "test"
        self.files = {
            "base_file": str(base_file),
            "dir": str(out_dir),
        }

        # Deliberately use one-based PDB serial values in ``num`` while
        # VorPy's internal/system indices remain zero-based.
        self.balls = pd.DataFrame(
            {
                "num": [1, 2, 3, 4],
                "name": ["A1", "A2", "B1", "B2"],
                "res_name": ["AAA", "AAA", "BBB", "BBB"],
                "res_seq": [1, 1, 2, 2],
                "chain": ["A", "A", "B", "B"],
                "loc": [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
                "element": ["C", "N", "O", "S"],
            }
        )


def test_source_backed_full_export_does_not_change_cwd(tmp_path):
    source = tmp_path / 'source.pdb'
    source.write_text('ATOM      7  OW  SOL    2       1.000   2.000   3.000  1.00  0.00           O\n')
    output = tmp_path / 'new'
    output.mkdir()
    system = DummySystem(source, output)
    system.balls = system.balls.iloc[:1]
    original = Path.cwd()
    write_pdb([0], 'full', system, directory=str(output))
    assert Path.cwd() == original


def test_source_less_molecular_export_uses_attested_serials_and_raw_chain(tmp_path):
    system = SimpleNamespace(type='mol', files={'base_file': None}, balls=pd.DataFrame({
        'num': [0], 'pdb_serial': [202], 'name': ['OW'], 'loc': [[1., 2., 3.]], 'rad': [1.5],
        'res': [SimpleNamespace(name='SOL')], 'chn': [SimpleNamespace(name='H2O')],
        'chain_name': [' '], 'res_seq': [2], 'element': ['O']}))
    write_pdb([0], 'attested', system, directory=str(tmp_path))
    line = (tmp_path / 'attested.pdb').read_text().splitlines()[0]
    assert int(line[6:11]) == 202 and line[21] == ' '
    system.balls = system.balls.drop(columns='pdb_serial')
    original = Path.cwd()
    with pytest.raises(ValueError, match='Original PDB serial identities are unavailable'):
        write_pdb([0], 'unverified', system, directory=str(tmp_path))
    assert not (tmp_path / 'unverified.pdb').exists()
    assert Path.cwd() == original


def test_pdb_reader_preserves_nonconsecutive_original_serials(tmp_path):
    from vorpy.src.output.pdb import make_pdb_line
    from vorpy.src.system import System
    source = tmp_path / 'original.pdb'
    source.write_text(make_pdb_line(ser_num=202, name='OW', res_name='SOL', chain=' ', res_seq=2, elem='O')
                      + make_pdb_line(ser_num=909, name='HW1', res_name='SOL', chain=' ', res_seq=2, x=1., elem='H'))
    system = System(file=str(source), make_dir=False)
    assert system.balls['pdb_serial'].tolist() == [202, 909]


def test_subset_indexing_ignores_ter_and_other_non_atom_records(tmp_path):
    base = tmp_path / "input.pdb"
    base.write_text(
        "HEADER test\n"
        "ATOM      1  A1  AAA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  A2  AAA A   1       1.000   0.000   0.000  1.00  0.00           N\n"
        "TER\n"
        "REMARK this line must not shift atom indexing\n"
        "ATOM      3  B1  BBB B   2       2.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      4  B2  BBB B   2       3.000   0.000   0.000  1.00  0.00           S\n"
        "END\n",
        encoding="utf-8",
    )

    sys = DummySystem(base, tmp_path)

    # Zero-based system index 2 is B1, despite TER/REMARK between chains.
    write_pdb([2], "subset", sys, directory=str(tmp_path))

    text = (tmp_path / "subset.pdb").read_text(encoding="utf-8")

    assert " B1  BBB B" in text
    assert " A2  AAA A" not in text
    assert " B2  BBB B" not in text


def test_subset_indexing_is_zero_based_not_pdb_serial(tmp_path):
    base = tmp_path / "input.pdb"
    base.write_text(
        "HEADER test\n"
        "ATOM      1  A1  AAA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  A2  AAA A   1       1.000   0.000   0.000  1.00  0.00           N\n"
        "TER\n"
        "ATOM      3  B1  BBB B   2       2.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      4  B2  BBB B   2       3.000   0.000   0.000  1.00  0.00           S\n",
        encoding="utf-8",
    )

    sys = DummySystem(base, tmp_path)

    # System index 1 must select the second atom, whose PDB serial is 2.
    write_pdb([1], "subset", sys, directory=str(tmp_path))

    text = (tmp_path / "subset.pdb").read_text(encoding="utf-8")

    assert " A2  AAA A" in text
    assert " A1  AAA A" not in text
    assert " B1  BBB B" not in text


def test_subset_indexing_handles_hetatm_records(tmp_path):
    base = tmp_path / "input.pdb"
    base.write_text(
        "HEADER test\n"
        "ATOM      1  A1  AAA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  ZN  ZN  A   2       1.000   0.000   0.000  1.00  0.00          ZN\n"
        "TER\n"
        "ATOM      3  B1  BBB B   3       2.000   0.000   0.000  1.00  0.00           O\n"
        "END\n",
        encoding="utf-8",
    )

    sys = DummySystem(base, tmp_path)
    sys.balls = sys.balls.iloc[:3].copy()

    # HETATM participates in the same positional atom-record index space.
    write_pdb([1], "subset", sys, directory=str(tmp_path))

    text = (tmp_path / "subset.pdb").read_text(encoding="utf-8")

    assert "HETATM" in text
    assert " ZN " in text
    assert " A1  AAA A" not in text
    assert " B1  BBB B" not in text
