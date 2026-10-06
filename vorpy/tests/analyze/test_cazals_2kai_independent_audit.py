import csv
import json
from pathlib import Path

import pytest

from vorpy.src.analyze.Part_II_Molecular_Analysis.SuplimentalInformation.cazals_2kai_independent_audit import (
    run_independent_2kai_audit,
)


@pytest.mark.slow
def test_independent_2kai_audit_reconstructs_frozen_s_h(tmp_path):
    pdb = Path(__file__).resolve().parents[2] / "data" / "2KAI.pdb"
    output = tmp_path / "audit"
    result = run_independent_2kai_audit(pdb, output)

    assert result["stop_stage"] is None
    assert result["pdb_census"]["total_coordinate_records"] == 2246
    assert result["pdb_census"]["ATOM_records"] == 2236
    assert result["pdb_census"]["HETATM_records"] == 10
    stages = {row["stage"]: row for row in result["stages"]}
    assert stages["01_pdb_atoms_radii"]["match"]
    assert stages["02_gudhi_regular_alpha"]["independent_counts"] == {
        0: 2236, 1: 17153, 2: 29762, 3: 14844,
    }
    assert stages["04_alpha0_bicolor_m5_selection"]["independent_final_pairs"] == 362
    assert stages["05_independent_edge_lengths"]["edge_count"] == 864
    assert stages["07_cazals_sign_and_final_sH"]["s_H_signed_deg"] == pytest.approx(
        -7.63423724994, abs=1e-9
    )
    assert (output / "stage_01_atoms_and_radii.csv").is_file()
    assert (output / "stage_02_gudhi_regular_and_alpha_simplices.csv").is_file()
    assert (output / "stage_04_bicolor_m5_facets.csv").is_file()
    with (output / "stage_05_independent_interface_edges.csv").open(newline="", encoding="utf-8") as stream:
        assert sum(1 for _ in csv.DictReader(stream)) == 864
    final = json.loads((output / "stage_07_final_sH.json").read_text(encoding="utf-8"))
    assert final["signed_numerator_A_rad"] == pytest.approx(-127.0500153782, abs=1e-9)
