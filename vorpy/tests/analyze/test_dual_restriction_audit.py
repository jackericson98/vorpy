from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest

from vorpy.src.analyze.dual_restriction_audit import (
    classify_power_edges, export_power_visualization, summarize_distances,
)
from vorpy.src.analyze.power_interface_curvature import analyze_power_interface_curvature
from vorpy.src.analyze.power_interface_cli import run_power_interface_pdb


def _atoms(points, groups):
    result = []
    for i, (point, group) in enumerate(zip(points, groups)):
        result.append(SimpleNamespace(
            index=i, chain="A" if group == "A" else "I", residue_name="GLY",
            residue_number=i + 1, insertion_code="", atom_name="CA", pdb_serial=i + 1,
            group=group, base_radius=1.5, expanded_radius=2.9, weight=2.9 ** 2,
            x=float(point[0]), y=float(point[1]), z=float(point[2]),
        ))
    return result


def _run(points, group_a, group_b):
    pytest.importorskip("gudhi")
    radii = np.full(len(points), 2.9)
    production = analyze_power_interface_curvature(points, radii, group_a, group_b)
    groups = ["A" if i in group_a else "B" for i in range(len(points))]
    return production, classify_power_edges(production, _atoms(points, groups))


def test_power_stage_audit_is_exhaustive_and_obeys_stage_order():
    points = np.random.default_rng(224).uniform(-4.0, 4.0, (24, 3))
    result, rows = _run(points, set(range(12)), set(range(12, 24)))
    assert len(rows) == len(result.full_simplices[1])
    assert {tuple(sorted((row["generator_i"], row["generator_j"]))) for row in rows} == result.full_simplices[1]
    assert all(row["alpha0_membership"] for row in rows if row["final_cazals_membership"])
    assert all(row["alpha0_membership"] and row["bicolor"]
               for row in rows if row["stage"] == "rejected_M5")
    assert all(row["condition_b_applicable"] for row in rows if row["final_cazals_membership"])
    assert len(summarize_distances(rows)[0:1]) == 1


def test_power_stage_classification_invariant_to_rigid_transform():
    points = np.random.default_rng(91).normal(size=(18, 3))
    group_a, group_b = set(range(9)), set(range(9, 18))
    result1, rows1 = _run(points, group_a, group_b)
    q, _ = np.linalg.qr(np.random.default_rng(8).normal(size=(3, 3)))
    transformed = points @ q + np.array([40.0, -13.0, 7.0])
    result2, rows2 = _run(transformed, group_a, group_b)
    classification1 = {tuple(sorted((r["generator_i"], r["generator_j"]))): r["stage"] for r in rows1}
    classification2 = {tuple(sorted((r["generator_i"], r["generator_j"]))): r["stage"] for r in rows2}
    assert classification1 == classification2
    assert {f.generator_ids for f in result1.final_facets} == {f.generator_ids for f in result2.final_facets}
    assert result1.edge_totals["signed"] == pytest.approx(result2.edge_totals["signed"], abs=1e-8)


def test_power_stage_classification_invariant_to_generator_permutation():
    points = np.random.default_rng(7).normal(size=(17, 3))
    a, b = set(range(8)), set(range(8, 17))
    _, rows = _run(points, a, b)
    permutation = np.random.default_rng(3).permutation(len(points))
    inverse = np.empty_like(permutation)
    inverse[permutation] = np.arange(len(permutation))
    perm_a = {int(inverse[index]) for index in a}
    perm_b = {int(inverse[index]) for index in b}
    _, perm_rows = _run(points[permutation], perm_a, perm_b)

    def physical_stages(values, xyz):
        output = {}
        for row in values:
            i, j = int(row["generator_i"]), int(row["generator_j"])
            key = tuple(sorted((tuple(np.round(xyz[i], 10)), tuple(np.round(xyz[j], 10)))))
            output[key] = row["stage"]
        return output

    assert physical_stages(rows, points) == physical_stages(perm_rows, points[permutation])


def test_stage_visualization_export_is_read_only(tmp_path):
    points = np.random.default_rng(13).normal(size=(12, 3))
    result, rows = _run(points, set(range(6)), set(range(6, 12)))
    atoms = _atoms(points, ["A"] * 6 + ["B"] * 6)
    before = (tuple(f.generator_ids for f in result.final_facets), result.edge_totals.copy())
    metadata = export_power_visualization(rows, atoms, tmp_path, "synthetic.pdb")
    assert metadata["distance_cutoff_used"] is False
    assert metadata["counts"]["2KAI_full_regular"] == len(rows)
    assert tuple(f.generator_ids for f in result.final_facets) == before[0]
    assert result.edge_totals == before[1]
    assert (tmp_path / "2KAI_dual_restriction_stages.pml").is_file()


@pytest.fixture(scope="module")
def frozen_2kai_audit(tmp_path_factory):
    pdb = Path(__file__).resolve().parents[2] / "data" / "2KAI.pdb"
    out = tmp_path_factory.mktemp("dual_restriction_2kai")
    run = run_power_interface_pdb(pdb, {"A", "B"}, {"I"}, out / "production",
                                  probe_radius=1.4, alpha=0.0, condition_beta_m=5.0)
    result, prepared = run["result"], run["prepared"]
    snapshot = (tuple(f.generator_ids for f in result.final_facets),
                tuple((e.regular_triangle, e.length, e.length_beta) for e in result.included_edges),
                result.edge_totals.copy())
    rows = classify_power_edges(result, prepared.atoms)
    yield result, rows, snapshot


def test_frozen_2kai_selection_and_curvature_unchanged_by_stage_audit(frozen_2kai_audit):
    result, rows, snapshot = frozen_2kai_audit
    final = {tuple(sorted((row["generator_i"], row["generator_j"]))) for row in rows
             if row["final_cazals_membership"]}
    assert len(rows) == len(result.full_simplices[1])
    assert sum(row["bicolor"] for row in rows) == 690
    assert sum(row["bicolor"] and row["alpha0_membership"] for row in rows) == 362
    assert final == set(snapshot[0])
    assert len(result.included_edges) == 864
    assert result.edge_totals["length"] == pytest.approx(953.5241609499, abs=2e-6)
    assert result.edge_totals["signed"] == pytest.approx(-127.0500153782, abs=2e-6)
    assert result.edge_totals["unsigned"] == pytest.approx(607.3649734659, abs=2e-6)
    assert tuple((e.regular_triangle, e.length, e.length_beta) for e in result.included_edges) == snapshot[1]
