"""Generate Figure 3F from matched, validated first-layer PEG waters.

Contacts are counts of external Voronoi interface faces, matching the validated
F5C definition.  They are not hydrogen bonds or counts of unique neighboring
molecules.  Waters are intersected by full PDB identity across AW, Power, and
Primitive before AW/Power comparison.
"""
from pathlib import Path

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

from .F3_PEG_prototypes import MODELS, REPO_ROOT
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import save_figure, setup_plots
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import Config
try:
    from .support.F5_peg_dataset import SCHEMES, screened_runs, validated_water_layer, water_identity
except ModuleNotFoundError:
    from .F5_peg_dataset import SCHEMES, screened_runs, validated_water_layer, water_identity


OUTPUT_DIR = REPO_ROOT / 'output' / 'figure3_round2_water_contacts'
SCHEME_COLORS = {'aw': '#2166ac', 'pow': '#d73027'}
SCHEME_LABELS = {'aw': 'AW', 'pow': 'Power'}
CONTACT_METRICS = {
    'water_contacts': 'Water–water interface faces',
    'polymer_contacts': 'Polymer–water interface faces',
    'water_sa_percent': 'Water-facing external area (%)',
}


def match_frame_layers(rows, config):
    """Validate each scheme and retain the three-scheme identity intersection."""
    layers = {}
    identities = {}
    for row in rows.itertuples(index=False):
        waters, pdb, _ = validated_water_layer(row, config)
        identity = {water_identity(pdb, residue): residue for residue in waters.index}
        if len(identity) != len(waters):
            raise ValueError(f'Ambiguous water identity in {row.molecule}/{row.frame}/{row.scheme}')
        layers[row.scheme] = waters
        identities[row.scheme] = identity
    if set(layers) != set(SCHEMES):
        raise ValueError(f'Expected {SCHEMES}, found {tuple(layers)}')
    common = sorted(set.intersection(*(set(identities[scheme]) for scheme in SCHEMES)))
    matched = []
    for scheme in ('aw', 'pow'):
        residues = [identities[scheme][key] for key in common]
        table = layers[scheme].loc[residues].reset_index(drop=True).copy()
        table['water_id'] = np.arange(len(common))
        table['scheme'] = scheme
        matched.append(table)
    coverage = {
        f'{scheme}_validated_waters': len(layers[scheme]) for scheme in SCHEMES
    }
    coverage['matched_waters'] = len(common)
    return pd.concat(matched, ignore_index=True), coverage


def collect_matched_waters(config=None):
    config = config or Config()
    runs = screened_runs()
    observations, coverage = [], []
    for (model, frame), rows in runs.groupby(['molecule', 'frame'], sort=True):
        matched, counts = match_frame_layers(rows, config)
        matched['model'] = model
        matched['frame'] = frame
        observations.append(matched)
        coverage.append(dict(model=model, frame=frame, **counts))
    coverage = pd.DataFrame(coverage)
    if len(coverage) != 69:
        raise ValueError(f'Expected 69 screened frames, found {len(coverage)}.')
    return pd.concat(observations, ignore_index=True), coverage


def summarize_contacts(observations):
    frame = (observations.groupby(['model', 'frame', 'scheme'], observed=True)
             .agg(mean_water_contacts=('water_contacts', 'mean'),
                  mean_polymer_contacts=('polymer_contacts', 'mean'),
                  mean_water_sa_percent=('water_sa_percent', 'mean'),
                  n_matched_waters=('water_id', 'size')).reset_index())
    model = (frame.groupby(['model', 'scheme'], observed=True)
             .agg(mean_water_contacts=('mean_water_contacts', 'mean'),
                  sd_frame_mean_water_contacts=('mean_water_contacts', 'std'),
                  mean_polymer_contacts=('mean_polymer_contacts', 'mean'),
                  sd_frame_mean_polymer_contacts=('mean_polymer_contacts', 'std'),
                  mean_water_sa_percent=('mean_water_sa_percent', 'mean'),
                  sd_frame_mean_water_sa_percent=('mean_water_sa_percent', 'std'),
                  n_frames=('frame', 'nunique'),
                  n_water_frame_observations=('n_matched_waters', 'sum')).reset_index())
    return frame, model


def paired_differences(observations):
    keys = ['model', 'frame', 'water_id']
    aw = observations[observations.scheme.eq('aw')]
    power = observations[observations.scheme.eq('pow')]
    pairs = aw.merge(power, on=keys, suffixes=('_aw', '_pow'), validate='one_to_one')
    pairs['delta_water_contacts_faces'] = pairs.water_contacts_pow - pairs.water_contacts_aw
    pairs['delta_polymer_contacts_faces'] = pairs.polymer_contacts_pow - pairs.polymer_contacts_aw
    pairs['delta_water_sa_percentage_points'] = pairs.water_sa_percent_pow - pairs.water_sa_percent_aw
    pairs['delta_external_area_pct'] = 100 * (
        pairs.external_area_pow - pairs.external_area_aw) / pairs.external_area_aw
    frame = (pairs.groupby(['model', 'frame'], observed=True)
             .agg(mean_delta_water_contacts_faces=('delta_water_contacts_faces', 'mean'),
                  mean_delta_polymer_contacts_faces=('delta_polymer_contacts_faces', 'mean'),
                  mean_delta_water_sa_percentage_points=('delta_water_sa_percentage_points', 'mean'),
                  mean_delta_external_area_pct=('delta_external_area_pct', 'mean'),
                  n_matched_waters=('water_id', 'size')).reset_index())
    model = (frame.groupby('model', observed=True)
             .agg(mean_delta_water_contacts_faces=('mean_delta_water_contacts_faces', 'mean'),
                  sd_frame_mean_delta_water_contacts_faces=('mean_delta_water_contacts_faces', 'std'),
                  mean_delta_polymer_contacts_faces=('mean_delta_polymer_contacts_faces', 'mean'),
                  sd_frame_mean_delta_polymer_contacts_faces=('mean_delta_polymer_contacts_faces', 'std'),
                  mean_delta_water_sa_percentage_points=('mean_delta_water_sa_percentage_points', 'mean'),
                  sd_frame_mean_delta_water_sa_percentage_points=('mean_delta_water_sa_percentage_points', 'std'),
                  mean_delta_external_area_pct=('mean_delta_external_area_pct', 'mean'),
                  sd_frame_mean_delta_external_area_pct=('mean_delta_external_area_pct', 'std'),
                  n_frames=('frame', 'nunique'),
                  n_water_frame_observations=('n_matched_waters', 'sum')).reset_index())
    return pairs, frame, model


def interior_support_audit(observations):
    diagnostic = observations.assign(
        zero_water_faces=observations.water_contacts.eq(0),
        low_water_facing_area=observations.water_sa_percent.le(10),
        polymer_face_dominated=observations.polymer_contacts.gt(observations.water_contacts),
    )
    summary = (diagnostic.groupby(['model', 'scheme'], observed=True)
               .agg(n_water_frame_observations=('water_id', 'size'),
                    zero_water_face_candidates=('zero_water_faces', 'sum'),
                    low_water_facing_area_candidates=('low_water_facing_area', 'sum'),
                    polymer_face_dominated_candidates=('polymer_face_dominated', 'sum'),
                    n_frames=('frame', 'nunique')).reset_index())
    summary['complete_atomic_geometry_required'] = True
    summary['exterior_connectivity_test_available'] = False
    summary['interior_or_buried_classification_supported'] = False
    summary['reason'] = ('Complete local cells do not establish enclosure or accessibility; '
                         'candidate diagnostics are not interior-water classifications.')
    return summary


def plot_contacts(model, output_dir):
    config = Config(figure_size=(12, 4.2), dpi=240, font_size=10, export_formats=('png', 'svg'))
    setup_plots(config)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=config.figure_size)
    x = np.arange(len(MODELS))
    for ax, (metric, label) in zip(axes, CONTACT_METRICS.items()):
        for offset, scheme in [(-.08, 'aw'), (.08, 'pow')]:
            rows = model[model.scheme.eq(scheme)].set_index('model').reindex(MODELS)
            mean = rows[f'mean_{metric}']
            sd = rows[f'sd_frame_mean_{metric}']
            ax.errorbar(x + offset, mean, yerr=sd, color=SCHEME_COLORS[scheme],
                        marker='o', capsize=2, label=SCHEME_LABELS[scheme])
        ax.set_xticks(x, MODELS, rotation=40, ha='right')
        ax.set(title=label, xlabel='Molecular system', ylabel='Frame-mean value')
    axes[0].legend(frameon=False)
    fig.suptitle('Matched first-layer water contact environments')
    save_figure(fig, output_dir / 'F3F_matched_water_contacts', config)


def plot_differences(model, output_dir):
    config = Config(figure_size=(12, 4.2), dpi=240, font_size=10, export_formats=('png', 'svg'))
    setup_plots(config)
    import matplotlib.pyplot as plt
    columns = [
        ('delta_water_contacts_faces', 'Water–water faces (Power − AW)'),
        ('delta_polymer_contacts_faces', 'Polymer–water faces (Power − AW)'),
        ('delta_water_sa_percentage_points', 'Water-facing area (Power − AW, pp)'),
    ]
    fig, axes = plt.subplots(1, 3, figsize=config.figure_size)
    rows = model.set_index('model').reindex(MODELS)
    x = np.arange(len(MODELS))
    for ax, (metric, label) in zip(axes, columns):
        ax.errorbar(x, rows[f'mean_{metric}'], yerr=rows[f'sd_frame_mean_{metric}'],
                    color='#59636e', marker='o', capsize=2)
        ax.axhline(0, color='#777777', linewidth=.8)
        ax.set_xticks(x, MODELS, rotation=40, ha='right')
        ax.set(title=label, xlabel='Molecular system', ylabel='Frame-mean difference')
    fig.suptitle('AW–Power differences for identity-matched first-layer waters')
    save_figure(fig, output_dir / 'F3F_matched_water_scheme_differences', config)


def main(output_dir=OUTPUT_DIR):
    output_dir = Path(output_dir)
    observations, coverage = collect_matched_waters()
    frame, model = summarize_contacts(observations)
    pairs, difference_frame, difference_model = paired_differences(observations)
    interior = interior_support_audit(observations)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        'F3F_coverage.csv': coverage,
        'F3F_matched_water_observations.csv': observations,
        'F3F_contact_frame_means.csv': frame,
        'F3F_contact_model_summary.csv': model,
        'F3F_paired_water_differences.csv': pairs,
        'F3F_difference_frame_means.csv': difference_frame,
        'F3F_difference_model_summary.csv': difference_model,
        'F3F_interior_support_audit.csv': interior,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / name, index=False)
    plot_contacts(model, output_dir)
    plot_differences(difference_model, output_dir)
    print(f'Generated Figure 3F from {len(coverage)} screened frames and '
          f'{len(pairs)} matched water/frame observations.')


if __name__ == '__main__':
    main()
