import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.Part_II_Molecular_Analysis.F4_atomic_deviations import F3_atomic_umap as analysis


def atom_table():
    return pd.DataFrame({
        'Index': [1, 2, 3, 4, 5, 6], 'Name': ['CA'] * 6, 'Residue': ['ALA'] * 6,
        'Residue Sequence': [1, 2, 3, 4, 5, 6], 'Chain': ['A'] * 6,
        'Volume': [10., 12, 15, 20, 25, 30], 'Surface Area': [3., 5, 4, 9, 11, 10],
        'Number of Neighbors': [2, 3, 2, 5, 6, 4], 'Maximum Mean Curvature': [1., 2, 3, 4, 5, 6],
        'Radius': [1., 1, 1.5, 1.5, 2, 2], 'Mass': [12.] * 6,
        'Complete Cell?': [True, True, True, True, True, False],
        'X': [100, 200, 300, 400, 500, 600], 'Closest Neighbor': [9, 8, 7, 6, 5, 4],
        'neighbors': [[1, 2]] * 6,
        'Moment of Inertia Tensor': [np.diag([1., 2., 3.]).tolist()] * 6,
    })


def test_all_physical_scalar_features_and_inertia_are_available_without_identity_leakage():
    features, report = analysis.prepare_atoms(atom_table(), 'AW')
    assert len(features) == 5
    assert report['incomplete_atoms_excluded'] == 1
    assert {'AW Volume', 'AW Surface Area', 'AW Number of Neighbors',
            'AW Maximum Mean Curvature', 'AW Radius', 'AW Mass'}.issubset(features)
    assert 'AW X' not in features and 'AW Closest Neighbor' not in features
    assert 'AW Residue Sequence' not in features and 'AW neighbors' not in features
    assert features['AW Inertia Eigenvalue 3'].eq(3).all()
    duplicated = pd.concat([atom_table(), atom_table().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match='Duplicate'):
        analysis.prepare_atoms(duplicated, 'AW')


def test_matching_uses_compound_identity_and_reports_unmatched_atoms(monkeypatch):
    from types import SimpleNamespace
    aw = atom_table()
    power = aw.iloc[::-1].copy()
    power.loc[power['Index'] == 2, 'Residue'] = 'GLY'
    power['Maximum Mean Curvature'] += 0.5
    monkeypatch.setattr(analysis, 'discover_systems', lambda *a, **k: [SimpleNamespace(name='D_test')])
    monkeypatch.setattr(analysis, 'read_pair', lambda *a, **k: ({'atoms': aw}, {'atoms': power}))
    atoms, report = analysis.collect_descriptors('.', [])
    assert atoms['Index'].tolist() == [1, 3, 4, 5]
    assert 'Power Maximum Mean Curvature' in atoms
    assert report[0]['aw_unmatched'] == 1
    assert atoms['System'].eq('D_test').all()


def test_scaling_imputation_and_exact_duplicate_removal_preserve_row_mapping():
    features, _ = analysis.prepare_atoms(atom_table(), 'AW')
    atoms = features.reset_index()
    atoms['Power Volume'] = atoms['AW Volume']
    atoms['AW Unavailable'] = np.nan
    atoms.loc[1, 'AW Surface Area'] = np.inf
    metadata, scaled, report = analysis.prepare_features(atoms)
    assert metadata['Index'].tolist() == [1, 2, 3, 4, 5]
    assert scaled.shape[0] == 5 and np.isfinite(scaled).all()
    np.testing.assert_allclose(scaled.mean(axis=0), 0, atol=1e-14)
    np.testing.assert_allclose(scaled.std(axis=0), 1)
    assert report['duplicate_features']['Power Volume'] == 'AW Volume'
    assert 'AW Unavailable' in report['all_missing_features']
    assert 'AW Mass' in report['constant_features']
    complete_metadata, _, _ = analysis.prepare_features(atoms, missing='complete')
    assert complete_metadata['Index'].tolist() == [1, 3, 4, 5]


def test_actual_umap_and_hdbscan_export_reproducible_assignments(tmp_path):
    pytest.importorskip('umap')
    rng = np.random.default_rng(17)
    raw = np.vstack([rng.normal(-3, 0.15, (20, 4)), rng.normal(3, 0.15, (20, 4))])
    atoms = pd.DataFrame(raw, columns=['AW Volume', 'AW Surface Area',
                                     'AW Number of Neighbors', 'Power Volume'])
    atoms.insert(0, 'Index', np.arange(40))
    atoms.insert(0, 'System', 'synthetic')
    metadata, scaled, report = analysis.prepare_features(atoms)
    results = analysis.export_results(metadata, scaled, report, [], tmp_path,
                                      neighbors=5, min_cluster_size=5)
    second = analysis.fit_clusters(scaled, neighbors=5, min_cluster_size=5)
    np.testing.assert_array_equal(results['Cluster'], second[2])
    np.testing.assert_allclose(results[['UMAP1', 'UMAP2']], second[1])
    assert results['Index'].tolist() == list(range(40))
    assert len(set(results['Cluster']) - {-1}) >= 2
    assert (tmp_path / 'F3_atomic_umap_report.json').is_file()
    assert (tmp_path / 'F3_atomic_umap_metrics.png').is_file()


def test_f3b_umap_option_passes_the_selected_original_logs_and_output(tmp_path, monkeypatch):
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F6_Proteins_Interfaces_Multiscale import F3B_atomic_groups as f3b
    calls = []
    monkeypatch.setattr(f3b, 'print_atomic_cache_status', lambda _: None)
    monkeypatch.setattr(f3b, 'load_or_collect_atomic_data', lambda **kwargs: atom_table())
    monkeypatch.setattr(f3b, 'export_atomic_outliers', lambda *args, **kwargs: None)
    monkeypatch.setattr(f3b, 'make_figure', lambda *args, **kwargs: None)
    monkeypatch.setattr(analysis, 'analyze', lambda *args, **kwargs: calls.append((args, kwargs)))
    f3b.main(data_root='original-logs', figure_dir=str(tmp_path), umap=True)
    assert calls == [(('original-logs', tmp_path / 'F3_atomic_umap'), {'exclude_keys': f3b.EXCLUDE_KEYS})]

