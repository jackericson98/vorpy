from pathlib import Path

from vorpy.src.analyze.curvature_showcase_2kai import _parse_off, audit_solvated


ROOT = Path(__file__).resolve().parents[3]


def test_centered_2kai_solvent_audit_counts_real_records():
    result = audit_solvated(ROOT / "vorpy/data/2KAI_MD_centered.pdb")
    assert result["all_atom_records"] == 44520
    assert result["protein_atoms"] == 4340
    assert result["water_atoms"] == 40086
    assert result["water_residues"] == 10000
    assert result["ion_or_other_atoms"] == 94
    assert result["available_solvated_geometry_logs_atom_count"] == 4340
    assert result["available_logs_include_water_or_ions"] is False


def test_aw_off_mesh_loader_ignores_blank_lines_and_face_colors():
    mesh = ROOT / "output/2KAI_curvature_showcase/static_aw_alpha_physical.off"
    points, faces = _parse_off(mesh)
    assert points.ndim == 2 and points.shape[1] == 3
    assert len(points) > 0 and len(faces) > 0
    assert all(3 <= len(face) for face in faces)
    assert min(min(face) for face in faces) >= 0
    assert max(max(face) for face in faces) < len(points)


def test_gaussian_totals_are_explicitly_unresolved_not_zero_filled():
    import csv
    table = ROOT / "output/2KAI_curvature_showcase/curvature_comparison.csv"
    rows = {r["representation"]: r for r in csv.DictReader(table.open(encoding="utf-8-sig"))}
    for key in ("static_power_alpha_interface", "static_aw_alpha_interface",
                "solvated_power_alpha_interface", "solvated_aw_alpha_interface"):
        assert rows[key]["integrated_gaussian"] == ""
        assert rows[key]["certified_gaussian_total"].lower() == "false"
        assert "UNRESOLVED" in rows[key]["status"] or "unresolved" in rows[key]["status"].lower()


def test_alpha_contact_export_does_not_claim_whole_spheres_are_patches():
    status = (ROOT / "output/2KAI_curvature_showcase/alpha_contact_patch_status.txt").read_text(encoding="utf-8")
    assert "unresolved" in status.lower()
    assert "Whole atom spheres" not in status
