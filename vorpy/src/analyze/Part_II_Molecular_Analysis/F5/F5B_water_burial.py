"""Compare AW and Power water-molecule geometry from water-inclusive F5 logs.

Run this file directly in PyCharm or with Python. Set DATASET_ROOT to a dataset
with the same F5 layout and complete water atom cells and incident faces.
Polymer-scoped logs cannot supply the external-area denominator or water volume.
Polymer-surrounded waters are internal-hydration candidates. Enclosure and
accessibility require an exterior-connectivity or solvent-accessibility test.
"""
import csv
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

if not __package__:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from .F5_config import Config, KEY
from .F5_common import normalize_log, read_pdb, reconstruct_interfaces, save_figure, setup_plots, write_csv
from .F5_run_all import discover


DATASET_ROOT = Config.dataset_root
OUTPUT_DIR = Path(__file__).parent / 'output/F5B'
EXTERIOR_MAX = 0.10
SURROUNDED_MIN = 0.90
CATEGORIES = ['Exterior-like', 'Intermediate', 'Polymer-surrounded']
MOLECULE_KEY = KEY + ['residue_id']
METRICS = {
    'volume': 'Water-molecule volume',
    'external_area': 'External water-molecule surface area',
    'neighbor_count': 'Unique neighboring molecules',
}


def has_water_atom_rows(path, config):
    """Check log scope before loading geometry. F5 validates the full log on load."""
    water_names = {name.upper()[:3] for name in config.water_residues}
    with path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.reader(stream)
        for row in reader:
            if row and row[0].strip().casefold() == 'atoms':
                header = [column.strip().casefold() for column in next(reader, [])]
                residue_column = header.index('residue')
                for atom in reader:
                    if atom and atom[0].strip().casefold() == 'surfaces':
                        return False
                    if len(atom) > residue_column and atom[residue_column].strip().upper() in water_names:
                        return True
                return False
    raise ValueError('Missing Atoms section')


def water_geometry(cells, faces, pdb, polymer, config):
    """Pool complete water cells, dropping both sides of intramolecular faces."""
    # An occurrence ID includes chain and residue identity and distinguishes wrapped numbers.
    expected_members = pdb[pdb.residue_id.isin(cells.residue_id)].groupby('residue_id').size()
    observed_members = cells.groupby('residue_id').size()
    complete_cells = cells.primary_valid & cells.neighbor_list_consistent
    complete_molecules = complete_cells.groupby(cells.residue_id).all()
    complete_molecules &= observed_members.eq(expected_members)
    valid_ids = complete_molecules[complete_molecules].index

    molecule_id = pdb.residue_id.copy()
    if polymer == 'Heparin':
        polymer_atoms = pdb.kind.eq('polymer')
        if not pdb.loc[polymer_atoms, 'order_valid'].all():
            raise ValueError('Heparin polymer molecule identity could not be validated')
        # F5 validates these Heparin models as one connected polysaccharide.
        molecule_id.loc[polymer_atoms] = 'Heparin polymer'

    incident_faces = []
    for side, other in [('i', 'j'), ('j', 'i')]:
        water = faces[f'residue_id_{side}']
        neighbor = faces[f'atom_{other}'].map(molecule_id)
        external = faces[f'kind_{side}'].eq('water') & water.isin(valid_ids) & water.ne(neighbor)
        incident_faces.append(pd.DataFrame({
            'residue_id': water[external],
            'neighbor_id': neighbor[external],
            'area': faces.loc[external, 'surface_area'],
            'polymer_neighbor': faces.loc[external, f'kind_{other}'].eq('polymer'),
        }))
    external_faces = pd.concat(incident_faces, ignore_index=True)
    external_faces['polymer_area'] = external_faces.area.where(external_faces.polymer_neighbor, 0.)
    result = external_faces.groupby('residue_id').agg(
        external_area=('area', 'sum'), polymer_area=('polymer_area', 'sum'),
    )
    active = external_faces[external_faces.area.gt(config.minimum_interface_area)]
    result['neighbor_count'] = active.groupby('residue_id').neighbor_id.nunique().reindex(result.index, fill_value=0)
    volume = cells.volume.where(np.isfinite(cells.volume) & cells.volume.ge(0))
    volume_by_water = volume.groupby(cells.residue_id)
    total_volume = volume_by_water.sum(min_count=1)
    result['volume'] = total_volume.where(volume_by_water.count().eq(observed_members))
    result = result[result.external_area.gt(0)].copy()
    result['polymer_contact_fraction'] = result.polymer_area / result.external_area
    return result.reset_index()


def load_water_geometry(path, pdb, polymer, config):
    atoms, faces, _, _ = normalize_log(path, pdb, config)
    water_atoms = atoms[atoms.kind.eq('water')]
    if water_atoms.empty:
        raise ValueError('No recognized water atom rows after identity validation')
    cells, faces, _ = reconstruct_interfaces(water_atoms, faces, pdb, config)
    return water_geometry(cells, faces, pdb, polymer, config)


def burial_category(fraction):
    # At the boundaries, 10% is exterior-like and 90% is polymer-surrounded.
    return np.select(
        [fraction.le(EXTERIOR_MAX), fraction.ge(SURROUNDED_MIN)],
        [CATEGORIES[0], CATEGORIES[2]], default=CATEGORIES[1],
    )


def summarize(pairs):
    rows = []
    for metric, label in METRICS.items():
        for category in CATEGORIES:
            group = pairs.loc[pairs.category.eq(category)]
            aw = group[f'{metric}_aw'].astype(float)
            power = group[f'{metric}_pow'].astype(float)
            valid = np.isfinite(aw) & np.isfinite(power) & aw.ne(0)
            differences = 100 * ((power[valid] - aw[valid]) / aw[valid]).abs()
            rows.append(dict(
                metric=label, category=category, mean_absolute_percent_difference=differences.mean(),
                n_matched=len(group), n_valid=int(valid.sum()), n_excluded=int((~valid).sum()),
            ))
    return pd.DataFrame(rows)


def plot_results(aw_waters, summary, config):
    setup_plots(config)
    import matplotlib.pyplot as plt

    for metric, label in METRICS.items():
        rows = summary[summary.metric.eq(label)]
        if not rows.n_valid.any():
            continue
        fig, ax = plt.subplots(figsize=config.figure_size)
        ax.bar(rows.category, rows.mean_absolute_percent_difference, color='#2166ac')
        ax.set_ylabel('Mean absolute difference from AW (%)')
        ax.set_title(f'{label} · Power vs AW')
        save_figure(fig, config.output_root / 'F5B' / f'F5B_{metric}', config)
    if aw_waters.empty:
        return
    fig, ax = plt.subplots(figsize=config.figure_size)
    ax.hist(100 * aw_waters.polymer_contact_fraction, bins=np.linspace(0, 100, 21), color='#2166ac')
    for threshold in [EXTERIOR_MAX, SURROUNDED_MIN]:
        ax.axvline(100 * threshold, color='#444444', linestyle='--')
    ax.set(xlabel='AW polymer-contact fraction (%)', ylabel='Water/frame observations',
           title='Water-molecule polymer contact in AW')
    save_figure(fig, config.output_root / 'F5B' / 'F5B_polymer_contact_histogram', config)

    aw_waters = aw_waters.assign(polymer_surrounded=aw_waters.category.eq(CATEGORIES[2]))
    counts = aw_waters.groupby(['polymer', 'model'], as_index=False).agg(
        polymer_surrounded_observations=('polymer_surrounded', 'sum'),
        valid_water_observations=('residue_id', 'size'),
        frames_with_valid_water=('frame', 'nunique'),
    )
    write_csv(counts, config.output_root / 'F5B' / 'F5B_surrounded_counts.csv')
    fig, ax = plt.subplots(figsize=config.figure_size)
    ax.bar(counts.model, counts.polymer_surrounded_observations, color='#2166ac')
    ax.tick_params(axis='x', labelrotation=40)
    ax.set(ylabel='Polymer-surrounded water/frame observations',
           title='Internal-hydration candidates by model in AW')
    save_figure(fig, config.output_root / 'F5B' / 'F5B_surrounded_counts', config)


def self_check():
    # Two waters and one polymer. The 100 Å² internal face must never enter the denominator.
    pdb = pd.DataFrame({'residue_id': ['a', 'a', 'a', 'b', 'polymer']})
    cells = pd.DataFrame({'residue_id': ['a', 'a', 'a', 'b'], 'volume': [1., 2., 3., 4.],
                          'primary_valid': True, 'neighbor_list_consistent': True})
    faces = pd.DataFrame({
        'atom_i': [0, 0, 1, 2], 'atom_j': [1, 4, 4, 3],
        'residue_id_i': ['a'] * 4, 'residue_id_j': ['a', 'polymer', 'polymer', 'b'],
        'kind_i': ['water'] * 4, 'kind_j': ['water', 'polymer', 'polymer', 'water'],
        'surface_area': [100., 6., 3., 1.],
    })
    config = Config()
    result = water_geometry(cells, faces, pdb, 'PEG', config).set_index('residue_id')
    assert result.loc['a', ['volume', 'external_area', 'polymer_area', 'neighbor_count']].tolist() == [6., 10., 9., 2]
    assert result.loc['b', 'polymer_contact_fraction'] == 0
    assert 'a' not in water_geometry(cells.iloc[1:], faces, pdb, 'PEG', config).residue_id.tolist()
    assert burial_category(pd.Series([EXTERIOR_MAX, (EXTERIOR_MAX + SURROUNDED_MIN) / 2, SURROUNDED_MIN])).tolist() == CATEGORIES


def main(config=None):
    if not 0 <= EXTERIOR_MAX < SURROUNDED_MIN <= 1:
        raise ValueError('Thresholds must satisfy 0 <= EXTERIOR_MAX < SURROUNDED_MIN <= 1')
    self_check()
    config = replace(config, schemes=('aw', 'pow')) if config else Config(
        dataset_root=Path(DATASET_ROOT), output_root=OUTPUT_DIR.parent,
        schemes=('aw', 'pow'), figure_size=(8, 5))
    results = {'aw': [], 'pow': []}
    coverage = []
    for frame_key, logs in discover(config).groupby(KEY, sort=False):
        metadata = dict(zip(KEY, frame_key))
        pdb = None
        for row in logs.itertuples(index=False):
            record = dict(metadata, scheme=row.scheme, source=str(row.log), status=row.status,
                          reason=row.missing_files, n_valid_water_observations=0)
            if row.status == 'available':
                try:
                    if not has_water_atom_rows(row.log, config):
                        record.update(status='water_geometry_unavailable', reason='No water atom rows in log')
                    else:
                        if pdb is None:
                            pdb = read_pdb(row.pdb, row.polymer, config)
                        waters = load_water_geometry(row.log, pdb, row.polymer, config)
                        if waters.empty:
                            record.update(status='water_geometry_unavailable', reason='No complete water molecules')
                        else:
                            results[row.scheme].append(waters.assign(**metadata))
                            record.update(status='ok', n_valid_water_observations=len(waters))
                except (OSError, csv.Error, ValueError, KeyError, IndexError, TypeError, StopIteration) as error:
                    record.update(status='invalid', reason=str(error))
            coverage.append(record)
    coverage = pd.DataFrame(coverage)
    write_csv(coverage, config.output_root / 'F5B' / 'F5B_coverage.csv')

    columns = MOLECULE_KEY + list(METRICS) + ['polymer_area', 'polymer_contact_fraction']
    aw, power = [pd.concat(results[scheme], ignore_index=True) if results[scheme]
                 else pd.DataFrame(columns=columns) for scheme in config.schemes]
    aw['category'] = burial_category(aw.polymer_contact_fraction)
    # Both schemes were checked against the same PDB. Each row contains every atom of that water.
    pairs = aw.merge(power, on=MOLECULE_KEY, suffixes=('_aw', '_pow'), validate='one_to_one')
    summary = summarize(pairs)
    write_csv(summary, config.output_root / 'F5B' / 'F5B_summary.csv')
    if not aw.empty:
        plot_results(aw, summary, config)
    print(coverage.status.value_counts().to_string())
    print(f'Valid AW water/frame observations: {len(aw):,}. Matched AW–Power observations: {len(pairs):,}.')
    if aw.empty:
        print('Water burial plots require water-inclusive AW logs with complete water cells and incident faces.')
    print(f"Output: {(config.output_root / 'F5B').resolve()}")


if __name__ == '__main__':
    main()
