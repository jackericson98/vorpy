from pathlib import Path
import pandas as pd

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
