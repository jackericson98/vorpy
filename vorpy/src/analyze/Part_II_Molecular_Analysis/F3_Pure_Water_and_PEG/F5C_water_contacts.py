"""Figure 5C from AW solvent logs; one observation per complete water/frame.

Areas are in Angstrom squared. Use the existing F5 coordinate and area tolerances
to validate PDB identity and reconcile each cell with all its incident surfaces.
No geometry is rebuilt and no source data are changed.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

if not __package__:
    import sys
    root = next(parent for parent in Path(__file__).resolve().parents
                if (parent / 'pyproject.toml').is_file())
    sys.path.insert(0, str(root))
    __package__ = ('vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG'
                   if 'F3_Pure_Water_and_PEG' in Path(__file__).parts else
                   'vorpy.src.analyze.Part_II_Molecular_Analysis.F3_SimpleModels')

from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import normalize_log, read_pdb
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import LEGACY_OUTPUT_ROOT, Config
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5B_water_burial import has_water_atom_rows
try:
    from .support.F5_peg_dataset import screened_runs, validated_water_layer
except ModuleNotFoundError:
    from .F5_peg_dataset import screened_runs, validated_water_layer


# Input paths follow the rerun layout: MODEL/frame_NNNN/MODEL/solvent/aw_logs.csv.
INPUT_ROOT = Path(
    r'C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack'
    r'\Ericson Voronoi DNA\P2\NewModels\PEG\Reruns\layers_4\aw'
)
# Some rerun folders omit the PDB; the F5 reader validates this fallback against
# every logged atom's identity and coordinates before it can be used.
PDB_ROOT = Config.dataset_root / 'PEG'
OUTPUT_DIR = LEGACY_OUTPUT_ROOT / 'PEG' / 'F5C'
RECOVERY_ROOT = OUTPUT_DIR / 'reruns'
MODELS = ('EG', 'PEG300', 'PEG1000', 'PEG2000', 'PEG3000', 'PEG4000', 'PEG6000')
# Include this model and all larger models in the ordered list above.
MIN_MODEL = 'EG'
FIGURE_SIZE = (8, 5)
AXIS_LABEL_SIZE = 18
TICK_LABEL_SIZE = 14
TITLE_SIZE = 20
LEGEND_SIZE = 11
TITLE = f'Water contacts ({MIN_MODEL} and larger)'
DPI = 300
MARKER_SIZE = 7
ALPHA = 0.25
TREND_BIN_COUNT = 20  # Equal-width Water SA bins; 20 gives 5-percentage-point bins.
WATER_COLOR = 'blue'
POLYMER_COLOR = 'red'


def water_contacts(log_path, pdb_path, config, layer_atoms=None):
    pdb = read_pdb(pdb_path, 'PEG', config)
    cells, faces, _, _ = normalize_log(log_path, pdb, config)
    if layer_atoms is not None:
        if not set(layer_atoms).issubset(cells.index):
            raise ValueError(f'Layer-1 cells absent from native log: {log_path}')
        cells = cells.loc[sorted(layer_atoms)]

    # Deduplicate repeated records by surface ID, never by neighboring atom pair:
    # distinct surfaces between the same pair must each count as a contact.
    faces[['atom_i', 'atom_j']] = np.sort(faces[['atom_i', 'atom_j']], axis=1)
    face_columns = ['surface_index', 'atom_i', 'atom_j', 'surface_area']
    faces = faces.drop_duplicates(face_columns)
    if faces.surface_index.duplicated().any():
        raise ValueError(f'Conflicting records for a surface ID: {log_path}')
    if faces.atom_i.eq(faces.atom_j).any() or faces.surface_area.lt(0).any():
        raise ValueError(f'Self-interface or negative surface area: {log_path}')
    if not faces.atom_i.isin(pdb.index).all() or not faces.atom_j.isin(pdb.index).all():
        raise ValueError(f'Surface endpoint missing from PDB: {log_path}')

    # Each face has one incidence on each atomic cell. A water-water face will
    # correctly contribute once to each of its two different water molecules.
    incident = pd.concat([
        faces.rename(columns={'atom_i': 'atom_index', 'atom_j': 'neighbor'}),
        faces.rename(columns={'atom_j': 'atom_index', 'atom_i': 'neighbor'}),
    ], ignore_index=True)
    incident['residue_id'] = incident.atom_index.map(pdb.residue_id)
    incident['neighbor_residue'] = incident.neighbor.map(pdb.residue_id)
    incident['neighbor_kind'] = incident.neighbor.map(pdb.kind)
    incident['valid_area'] = np.isfinite(incident.surface_area)

    # Require every cell to be complete, its entire incident-area sum to match
    # the logged cell area, and its observed neighbors to match the log list.
    summed = incident.groupby('atom_index').agg(
        area=('surface_area', 'sum'), valid_area=('valid_area', 'all'),
        neighbors=('neighbor', set),
    ).reindex(cells.index)
    valid_cells = (
        cells.complete.fillna(False) & cells.identity_valid
        & summed.valid_area.fillna(False)
        & np.isfinite(cells.surface_area) & cells.surface_area.gt(0)
        & np.isclose(summed.area, cells.surface_area,
                     rtol=config.area_rtol, atol=config.area_atol)
    )
    valid_cells &= pd.Series([
        isinstance(logged, list) and observed == set(logged)
        for logged, observed in zip(cells.logged_neighbors, summed.neighbors)
    ], index=cells.index)

    # PDB residue occurrence identity includes segment, chain, and insertion
    # code, and distinguishes wrapped residue numbers. Require OW/HW1/HW2.
    waters = pdb[pdb.kind.eq('water')]
    names = waters.groupby('residue_id').atom_name.agg(list)
    three_sites = names.map(lambda values: sorted(values) == ['HW1', 'HW2', 'OW'])
    observed = cells.groupby('residue_id').size().reindex(names.index, fill_value=0)
    complete = valid_cells.groupby(cells.residue_id).all().reindex(names.index, fill_value=False)
    valid_ids = names.index[three_sites & observed.eq(3) & complete]

    # Exclude O-H and H-H faces inside the same water before summing external
    # area. Other species, if present, remain in the external-area denominator.
    external = incident[
        incident.residue_id.isin(valid_ids)
        & incident.residue_id.ne(incident.neighbor_residue)
    ].copy()
    external['water_contact'] = external.neighbor_kind.eq('water')
    external['polymer_contact'] = external.neighbor_kind.eq('polymer')
    external['water_area'] = external.surface_area.where(external.water_contact, 0.)
    result = external.groupby('residue_id').agg(
        external_area=('surface_area', 'sum'), water_area=('water_area', 'sum'),
        water_contacts=('water_contact', 'sum'), polymer_contacts=('polymer_contact', 'sum'),
    )
    # Never turn an absent cell/surface set or a zero denominator into contacts.
    result = result[np.isfinite(result.external_area) & result.external_area.gt(0)].copy()
    internal = faces[faces.atom_i.map(pdb.residue_id).eq(faces.atom_j.map(pdb.residue_id))]
    internal_area = internal.groupby(internal.atom_i.map(pdb.residue_id)).surface_area.sum()
    expected = cells.groupby('residue_id').surface_area.sum() - 2 * internal_area.reindex(result.index, fill_value=0)
    if not np.isclose(result.external_area, expected.reindex(result.index),
                      rtol=config.area_rtol, atol=config.area_atol).all():
        raise ValueError(f'Inconsistent external water area: {log_path}')
    result['water_sa_percent'] = 100 * (result.water_area / result.external_area)
    if not result.water_sa_percent.between(0, 100).all():
        raise ValueError(f'Water SA percentage outside [0, 100]: {log_path}')
    return result


def main():
    if TREND_BIN_COUNT < 1:
        raise ValueError('TREND_BIN_COUNT must be positive')
    config = Config()
    observations = []
    runs = screened_runs()
    coverage = []
    for model in MODELS[MODELS.index(MIN_MODEL):]:
        model_results = []
        selected = runs[runs.molecule.eq(model) & runs.scheme.eq('aw')]
        for row in selected.itertuples(index=False):
            result, _, atoms = validated_water_layer(row, config)
            result = result.assign(model=model, frame=row.frame)
            model_results.append(result)
            coverage.append(dict(model=model, frame=row.frame, atoms=len(atoms), water_molecules=len(result)))
        count = sum(len(result) for result in model_results)
        print(f'{model}: {len(model_results)} frames, {3 * count} atoms, {count} layer-1 water molecules', flush=True)
        if not count:
            raise ValueError(f'No valid waters for {model}; check INPUT_ROOT and complete solvent layers')
        observations.extend(model_results)

    data = pd.concat(observations, ignore_index=True)
    # Join bin means with straight segments: each water/frame has equal weight
    # within its bin. Empty bins are omitted; no trend is extrapolated beyond
    # the mean SA values of the populated bins. Nonnegative counts give
    # nonnegative means and connecting segments, without clipping fitted values.
    bins = pd.cut(data.water_sa_percent, np.linspace(0, 100, TREND_BIN_COUNT + 1),
                  include_lowest=True)
    trend = data.groupby(bins, observed=True)[
        ['water_sa_percent', 'water_contacts', 'polymer_contacts']
    ].mean()
    assert np.isfinite(trend.to_numpy()).all() and trend.ge(0).all().all()
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=FIGURE_SIZE, layout='constrained')
    for column, color, label in [
        ('water_contacts', WATER_COLOR, 'Water–water contacts'),
        ('polymer_contacts', POLYMER_COLOR, 'Water–polymer contacts'),
    ]:
        ax.scatter(data.water_sa_percent, data[column], s=MARKER_SIZE,
                   alpha=ALPHA, c=color, edgecolors='none', label=label)
        ax.plot(trend.water_sa_percent, trend[column], color=color, linewidth=2.5,
                linestyle='--', label=label.replace('contacts', 'polyline (bin means)'))
    ax.set(xlabel='Water SA (%)', ylabel='Number of contacts (surfaces)',
           xlim=(0, 100), ylim=(0, None))
    ax.xaxis.label.set_size(AXIS_LABEL_SIZE)
    ax.yaxis.label.set_size(AXIS_LABEL_SIZE)
    ax.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)
    ax.set_title(TITLE, fontsize=TITLE_SIZE)
    ax.legend(fontsize=LEGEND_SIZE)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUTPUT_DIR / 'F5C_water_contacts.csv', index_label='residue_id')
    trend.to_csv(OUTPUT_DIR / 'F5C_water_contacts_bins.csv', index_label='water_sa_bin')
    pd.DataFrame(coverage).to_csv(OUTPUT_DIR / 'F5C_coverage.csv', index=False)
    for extension in ('png', 'svg'):
        path = OUTPUT_DIR / f'F5C_water_contacts.{extension}'
        fig.savefig(path, dpi=DPI)
        print(f'Saved {path}', flush=True)
    plt.close(fig)


if __name__ == '__main__':
    main()
