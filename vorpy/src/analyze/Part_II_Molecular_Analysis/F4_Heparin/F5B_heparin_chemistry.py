"""Figure 5B: Heparin chemistry, AW water coordination, and chain position.

Use existing polymer-cell logs only: volume is in angstrom cubed and atomic
Surface Area (not Contact Area) is in angstrom squared. Each signed percentage
is 100 * (Power - AW) / AW, calculated before averaging on the same valid atoms.
SD means sample SD (ddof=1) of pooled atom/frame observations, not uncertainty
from independent replicates. A repeat here is one monosaccharide PDB residue
(dpN has N repeats), consistent with the existing F5 chain-order reader.
"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd


# Inputs and all presentation settings. No solver, layer export, or cache writes.
INPUT_ROOT = Path(
    r'C:\Users\Ryzen-3600\Versiti Dropbox\John Ericson\Versiti lab\Jack'
    r'\Ericson Voronoi DNA\OrganizedModels\Heparin'
)
MODELS = ('dp2', 'dp12', 'dp18', 'dp24', 'dp30')
DPI = 300
# Typography defaults for all six plots (points); per-plot overrides stay below.
FONT_SIZE = 13
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 13
LEGEND_SIZE = 13
TITLE_SIZE = 17
COUNT_SIZE = 12  # Counts, annotations, and figure notes.
POINT_SIZE = 6
POINT_ALPHA = 0.10
POINT_JITTER = 0.10
RANDOM_SEED = 0
MEAN_SIZE = 48
MEAN_MIN_N = 10
MEAN_LINE_WIDTH = 0.9
ERROR_CAP_SIZE = 2
BAR_WIDTH = 0.75
POSITION_WIDTH = 0.12
ZERO_COLOR = '#555555'
# Each entry controls ONE figure independently. Limits of None use autoscaling.
# point_size/mean_size are scatter areas (points squared); marker_size is a
# diameter in points for terminal/legend markers. Error bars show sample SD.
# legend_loc uses Matplotlib figure-legend locations, e.g. 'upper center'.
PLOT_SETTINGS = {
    'volume_water_coordination': {
        'figsize': (8, 6), 'xlim': None, 'ylim': None, 'show_title': True,
        'xlabel_size': AXIS_LABEL_SIZE, 'ylabel_size': AXIS_LABEL_SIZE,
        'xtick_size': TICK_LABEL_SIZE, 'ytick_size': TICK_LABEL_SIZE,
        'show_legend': True, 'legend_loc': 'upper center', 'marker_size': 5,
        'point_size': POINT_SIZE, 'mean_size': MEAN_SIZE, 'show_error_bars': False,
    },
    'volume_functional_groups': {
        'figsize': (8, 7), 'xlim': None, 'ylim': None, 'show_title': True,
        'xlabel_size': AXIS_LABEL_SIZE, 'ylabel_size': AXIS_LABEL_SIZE,
        'xtick_size': TICK_LABEL_SIZE, 'ytick_size': TICK_LABEL_SIZE,
        'show_legend': True, 'legend_loc': 'upper center', 'marker_size': 5,
        'show_error_bars': True,
    },
    'volume_terminal_effects': {
        'figsize': (9, 6), 'xlim': None, 'ylim': None, 'show_title': True,
        'xlabel_size': AXIS_LABEL_SIZE, 'ylabel_size': AXIS_LABEL_SIZE,
        'xtick_size': TICK_LABEL_SIZE, 'ytick_size': TICK_LABEL_SIZE,
        'show_legend': True, 'legend_loc': 'upper center', 'marker_size': 5,
        'show_error_bars': True,
    },
    'surface_area_water_coordination': {
        'figsize': (8, 6), 'xlim': None, 'ylim': None, 'show_title': True,
        'xlabel_size': AXIS_LABEL_SIZE, 'ylabel_size': AXIS_LABEL_SIZE,
        'xtick_size': TICK_LABEL_SIZE, 'ytick_size': TICK_LABEL_SIZE,
        'show_legend': True, 'legend_loc': 'upper center', 'marker_size': 5,
        'point_size': POINT_SIZE, 'mean_size': MEAN_SIZE, 'show_error_bars': False,
    },
    'surface_area_functional_groups': {
        'figsize': (8, 7), 'xlim': None, 'ylim': None, 'show_title': True,
        'xlabel_size': AXIS_LABEL_SIZE, 'ylabel_size': AXIS_LABEL_SIZE,
        'xtick_size': TICK_LABEL_SIZE, 'ytick_size': TICK_LABEL_SIZE,
        'show_legend': True, 'legend_loc': 'upper center', 'marker_size': 5,
        'show_error_bars': True,
    },
    'surface_area_terminal_effects': {
        'figsize': (9, 6), 'xlim': None, 'ylim': None, 'show_title': True,
        'xlabel_size': AXIS_LABEL_SIZE, 'ylabel_size': AXIS_LABEL_SIZE,
        'xtick_size': TICK_LABEL_SIZE, 'ytick_size': TICK_LABEL_SIZE,
        'show_legend': True, 'legend_loc': 'upper center', 'marker_size': 5,
        'show_error_bars': True,
    },
}
METRIC_LABELS = {'volume': 'volume', 'surface_area': 'surface area'}
ANALYSES = {'A': 'water_coordination', 'B': 'functional_groups', 'C': 'terminal_effects'}
COLORS = {
    'Carbon': '#484848',
    'Hydroxyl oxygen': '#0072B2',
    'Glycosidic oxygen': '#009E73',
    'Sulfate oxygen': '#D55E00',
    'Carboxylate oxygen': '#CC79A7',
    'Ring oxygen': '#56B4E9',
    'Sulfur': '#A68A00',
}
POSITIONS = ('Terminal', 'Near-terminal', 'Interior')
# Named covalent bonds are validated against coordinates, not Voronoi neighbors.
# These bounds (angstrom) match the existing F5 heavy/O-H connectivity checks.
HEAVY_BOND_RANGE = (0.9, 1.9)
OH_BOND_RANGE = (0.6, 1.35)
BASE_BONDS = 'C1-C2 C2-C3 C3-C4 C4-C5 C5-O5 O5-C1 C5-C6 C3-O3 C4-O4'
RESIDUE_BONDS = {
    'AGLC': 'C2-N N-S2 S2-OS22 S2-OS23 S2-OS24 C6-O6 O6-S6 '
            'S6-OS62 S6-OS63 S6-OS64',
    'AIDO': 'C2-O2 O2-S2 S2-OS22 S2-OS23 S2-OS24 C6-O61 C6-O62',
}

if not __package__:
    import sys
    root = next(parent for parent in Path(__file__).resolve().parents
                if (parent / 'pyproject.toml').is_file())
    sys.path.insert(0, str(root))
    __package__ = ('vorpy.src.analyze.Part_II_Molecular_Analysis.F4_Heparin'
                   if 'F4_Heparin' in Path(__file__).parts else
                   'vorpy.src.analyze.Part_II_Molecular_Analysis.F3_SimpleModels')

from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import (
    load_scheme, match_schemes, read_pdb,
)
from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_config import (
    LEGACY_OUTPUT_ROOT, Config,
)

OUTPUT_DIR = LEGACY_OUTPUT_ROOT / 'Heparin' / 'F5B'


def classify_chemistry(pdb):
    """Require known residue/atom names AND the expected local bond graph.

    PDBs have explicit hydrogens but no required CONECT records. Recover bonds
    from unwrapped coordinates using the bounds above; ambiguous or broken
    motifs remain Unresolved. No periodic bond repair is guessed. Sulfate
    oxygen includes O-sulfate ester bridges and nonbridging O atoms on both
    O-sulfates and N-sulfamate groups; ring O5 is NOT a glycosidic oxygen.
    """
    p = pdb[pdb.kind.eq('polymer')].copy()
    heavy = p[~p.is_hydrogen]
    labels = list(zip(heavy.residue_id, heavy.atom_name))
    if len(labels) != len(set(labels)):
        raise ValueError('Duplicate polymer atom names within a residue')
    xyz = heavy[['x', 'y', 'z']].to_numpy()
    distances = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
    bonded = (distances >= HEAVY_BOND_RANGE[0]) & (distances <= HEAVY_BOND_RANGE[1])
    actual = {key: {labels[j] for j in np.flatnonzero(bonded[i])}
              for i, key in enumerate(labels)}
    expected = {key: set() for key in labels}

    def add_bond(left, right):
        # Retain missing named endpoints so a missing atom cannot pass validation.
        expected.setdefault(left, set()).add(right)
        expected.setdefault(right, set()).add(left)

    for residue, group in p.groupby('residue_id', sort=False):
        name = group.residue_name.iloc[0]
        if name not in RESIDUE_BONDS:
            raise ValueError(f'Unsupported Heparin residue {name}')
        bonds = BASE_BONDS + ' ' + RESIDUE_BONDS[name]
        if group.atom_name.eq('O1').any():
            bonds += ' C1-O1'
        for bond in bonds.split():
            left, right = bond.split('-')
            add_bond((residue, left), (residue, right))

    # Inter-residue glycosidic bonds require named C1--O4 endpoints. Valence
    # checks below reject missing or multiple connections at either endpoint.
    chains = dict(zip(labels, heavy.chain))
    for key, neighbors in actual.items():
        if key[1] == 'O4':
            for other in neighbors:
                if other[0] != key[0] and other[1] == 'C1' and chains[key] == chains[other]:
                    add_bond(key, other)
    elements = dict(zip(labels, heavy.element))
    valid = {key: actual[key] == expected[key] and bool(expected[key])
                  and elements[key] == key[1][0] for key in labels}
    for key in labels:
        if key[1] == 'C1':
            valid[key] &= len(expected[key]) == 3
        if key[1] == 'O4':
            has_ho4 = ((p.residue_id == key[0]) & (p.atom_name == 'HO4')).any()
            valid[key] &= len(expected[key]) == (1 if has_ho4 else 2)

    h = p[p.is_hydrogen]
    p['atom_class'] = 'Unresolved'
    p['classification_reason'] = 'named covalent motif not validated'
    for idx, atom in p.iterrows():
        key = (atom.residue_id, atom.atom_name)
        if atom.is_hydrogen:
            p.loc[idx, ['atom_class', 'classification_reason']] = ['Hydrogen', 'not plotted']
            continue
        # Check adjacent motifs too: e.g. O61 needs the complete C6(O61,O62) group.
        if not valid.get(key, False) or not all(valid.get(n, False) for n in expected[key]):
            continue
        name, residue = atom.atom_name, atom.residue_name
        role = None
        if name in {f'C{i}' for i in range(1, 7)} and atom.element == 'C':
            role = 'Carbon'
        elif name in {'S2', 'S6'} and atom.element == 'S':
            role = 'Sulfur'
        elif name == 'N' and atom.element == 'N':
            role = 'Nitrogen'
        elif atom.element == 'O':
            dh = np.linalg.norm(h[['x', 'y', 'z']].to_numpy() - atom[['x', 'y', 'z']].to_numpy(dtype=float), axis=1)
            attached = set(zip(h.loc[(dh >= OH_BOND_RANGE[0]) & (dh <= OH_BOND_RANGE[1]), 'residue_id'],
                               h.loc[(dh >= OH_BOND_RANGE[0]) & (dh <= OH_BOND_RANGE[1]), 'atom_name']))
            hydroxyl = name in {'O1', 'O3'} or (name == 'O4' and len(expected[key]) == 1)
            if attached != ({(key[0], 'H' + name)} if hydroxyl else set()):
                p.loc[idx, 'classification_reason'] = 'O-H attachment not validated'
                continue
            if hydroxyl:
                role = 'Hydroxyl oxygen'
            elif name == 'O4':
                role = 'Glycosidic oxygen'
            elif name == 'O5':
                role = 'Ring oxygen'
            elif (name in {'OS22', 'OS23', 'OS24', 'OS62', 'OS63', 'OS64'}
                  or (residue == 'AGLC' and name == 'O6')
                  or (residue == 'AIDO' and name == 'O2')):
                role = 'Sulfate oxygen'
            elif residue == 'AIDO' and name in {'O61', 'O62'}:
                role = 'Carboxylate oxygen'
        if role is not None:
            p.loc[idx, ['atom_class', 'classification_reason']] = [role, 'validated named bond graph']

    # The shared reader establishes chain order from a single C1--O4 path.
    # Terminal takes precedence for dp2, which has no near-terminal/interior.
    p['position'] = np.select(
        [p.order_valid & p.end_distance.eq(0), p.order_valid & p.end_distance.eq(1),
         p.order_valid & p.end_distance.ge(2)], POSITIONS, default='Unresolved')
    return p.reset_index(drop=True)


def load_frame(folder, model, config):
    """Match existing complete cells; AW interfaces supply unique water IDs."""
    pdb = read_pdb(folder / 'frame.pdb', 'Heparin', config)
    chemistry = classify_chemistry(pdb)
    if chemistry.residue_id.nunique() != int(model[2:]):
        raise ValueError(f'{model}: PDB residue count disagrees with dp')
    metadata = dict(polymer='Heparin', model=model, frame=folder.name)
    atoms = pd.concat([
        load_scheme(folder / f'{scheme}_logs.csv', pdb, config,
                    dict(metadata, scheme=scheme))['atoms'] for scheme in ('aw', 'pow')
    ], ignore_index=True)
    pairs, _ = match_schemes(atoms, config)
    valid = pairs.pair_valid.copy()
    excluded = {'invalid_or_unmatched_pairs': int((~valid).sum())}
    for metric in METRIC_LABELS:
        aw, power = pairs[f'{metric}_aw'], pairs[f'{metric}_s']
        denominator_ok = np.isfinite(aw) & aw.gt(0)
        power_ok = np.isfinite(power) & power.gt(0)
        # Reason counts can overlap cell-invalid counts. Require both metrics
        # to be valid so volume and surface-area plots use identical atoms.
        excluded[f'{metric}_invalid_aw_denominator'] = int((pairs.identity_match & ~denominator_ok).sum())
        excluded[f'{metric}_invalid_power_value'] = int((pairs.identity_match & ~power_ok).sum())
        valid &= denominator_ok & power_ok
    excluded['total_excluded'] = int((~valid).sum())
    points = pairs.loc[valid, ['model', 'frame', 'atom_index', 'volume_aw', 'volume_s',
                               'surface_area_aw', 'surface_area_s', 'water_coordination_aw']].copy()
    # The shared AW reader counts a SET of PDB water-molecule occurrence IDs
    # over positive-area polymer--water interfaces, including water H contacts.
    # Multiple contacted sites on the same water therefore count only once.
    points = points.rename(columns={'water_coordination_aw': 'water_contacts'})
    if not (points.water_contacts.ge(0) & points.water_contacts.eq(points.water_contacts.round())).all():
        raise ValueError('Invalid AW water coordination count')
    for metric in METRIC_LABELS:
        points[f'delta_{metric}_pct'] = 100 * (points[f'{metric}_s'] - points[f'{metric}_aw']) / points[f'{metric}_aw']
    points = points.merge(chemistry[['atom_index', 'residue_name', 'atom_name', 'atom_class',
                                     'classification_reason', 'position']],
                          on='atom_index', validate='one_to_one')
    if not np.isfinite(points[['delta_volume_pct', 'delta_surface_area_pct']]).all().all():
        raise ValueError('Nonfinite signed volume/surface-area differences')
    return points, excluded


def summarize(points):
    """One CSV: plotted pooled summaries, per-model summaries, and name audit."""
    points = points.assign(frame_id=points.model + '/' + points.frame)
    plotted = points[points.atom_class.isin(COLORS)]
    tables = []
    selections = [('A', plotted[~plotted.atom_class.eq('Sulfur')], ['atom_class', 'water_contacts']),
                  ('B', plotted, ['atom_class']),
                  ('C', plotted[~plotted.model.eq('dp2') & ~plotted.atom_class.eq('Sulfur')
                                & plotted.position.isin(POSITIONS)], ['atom_class', 'position']),
                  ('chemistry_audit', points, ['residue_name', 'atom_name', 'atom_class', 'classification_reason'])]
    for panel, data, keys in selections:
        # Pooled observations weight longer chains more heavily. Per-model rows
        # retain chain-length information without adding panels or output files.
        for by_model in (False, True):
            grouped = data.groupby((['model'] if by_model else []) + keys, observed=True)
            table = grouped.agg(mean_delta_volume_pct=('delta_volume_pct', 'mean'),
                                sd_delta_volume_pct=('delta_volume_pct', 'std'),
                                mean_delta_surface_area_pct=('delta_surface_area_pct', 'mean'),
                                sd_delta_surface_area_pct=('delta_surface_area_pct', 'std'),
                                n=('delta_volume_pct', 'size'), n_frames=('frame_id', 'nunique'),
                                n_models=('model', 'nunique')).reset_index()
            table['panel'] = panel
            if not by_model:
                table['model'] = 'pooled'
            table['draw_mean'] = table.n.ge(MEAN_MIN_N) if panel == 'A' else True
            tables.append(table)
    return pd.concat(tables, ignore_index=True)


def plot_figure(points, summary, metric, panel):
    """Draw one independent figure, keeping the original A/B/C visual styles."""
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    key = f'{metric}_{ANALYSES[panel]}'
    settings = PLOT_SETTINGS[key]
    mean_col, sd_col = f'mean_delta_{metric}_pct', f'sd_delta_{metric}_pct'
    plt.rcParams.update({'font.size': FONT_SIZE, 'axes.titlesize': TITLE_SIZE,
                         'legend.fontsize': LEGEND_SIZE, 'svg.fonttype': 'none'})
    fig, ax = plt.subplots(figsize=settings['figsize'])
    pooled = summary[summary.model.eq('pooled') & summary.panel.eq(panel)]
    roles = [role for role in COLORS if role in set(pooled.atom_class)]
    if panel == 'A':
        rng = np.random.default_rng(RANDOM_SEED)
        for role in roles:
            color = COLORS[role]
            group = points[points.atom_class.eq(role)]
            ax.scatter(group.water_contacts + rng.uniform(-POINT_JITTER, POINT_JITTER, len(group)),
                       group[f'delta_{metric}_pct'], s=settings['point_size'], alpha=POINT_ALPHA,
                       color=color, linewidths=0, rasterized=True)
            means = pooled[pooled.atom_class.eq(role) & pooled.draw_mean]
            ax.plot(means.water_contacts, means[mean_col], color=color, linewidth=MEAN_LINE_WIDTH)
            if settings['show_error_bars']:
                ax.errorbar(means.water_contacts, means[mean_col], yerr=means[sd_col],
                            fmt='none', color=color, capsize=ERROR_CAP_SIZE, elinewidth=.8)
            ax.scatter(means.water_contacts, means[mean_col], s=settings['mean_size'],
                       color=color, edgecolors='white', linewidths=.7, zorder=3)
        xlabel = 'Neighboring water molecules in AW'
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.text(.03, .97, f'Large means: n ≥ {MEAN_MIN_N}\nHorizontal jitter for visibility',
                transform=ax.transAxes, va='top', fontsize=COUNT_SIZE)
    elif panel == 'B':
        bars = pooled.set_index('atom_class').reindex(COLORS).dropna(subset=['n'])
        bars = bars[(bars.index != 'Sulfur') | bars.n.ge(MEAN_MIN_N)]
        roles = list(bars.index)
        x = np.arange(len(bars))
        errors = bars[sd_col] if settings['show_error_bars'] else None
        ax.bar(x, bars[mean_col], yerr=errors, color=[COLORS[c] for c in roles],
               width=BAR_WIDTH, capsize=ERROR_CAP_SIZE, error_kw={'elinewidth': .8})
        for i, (_, row) in enumerate(bars.iterrows()):
            height = row[mean_col] + (row[sd_col] if settings['show_error_bars'] else 0)
            # Keep enlarged black counts clear of the filled negative bars.
            ax.annotate(f'n={int(row.n):,}', (i, max(0, height)), xytext=(0, 5),
                        textcoords='offset points', ha='center', fontsize=COUNT_SIZE, rotation=90)
        ax.set_xticks(x, [c.replace(' oxygen', '\noxygen') for c in roles], rotation=45, ha='right')
        xlabel = 'Atomic chemical class'
    else:
        # Compare the SAME classes across unchanged repeat regions; dp2 remains
        # excluded from C. No pooled elemental oxygen average is introduced.
        for j, role in enumerate(roles):
            group = pooled[pooled.atom_class.eq(role)].set_index('position').reindex(POSITIONS)
            x = np.arange(len(POSITIONS)) + (j - (len(roles) - 1) / 2) * POSITION_WIDTH
            errors = group[sd_col] if settings['show_error_bars'] else None
            ax.errorbar(x, group[mean_col], yerr=errors, fmt='o', markersize=settings['marker_size'],
                        color=COLORS[role], capsize=ERROR_CAP_SIZE, elinewidth=.8)
            for i, (_, row) in enumerate(group.iterrows()):
                if pd.notna(row.n):
                    height = row[mean_col] + (row[sd_col] if settings['show_error_bars'] else 0)
                    ax.annotate(f'n={int(row.n):,}', (x[i], height), xytext=(0, 5),
                                textcoords='offset points', ha='center', fontsize=COUNT_SIZE, rotation=90)
        ax.set_xticks(range(3), ['Terminal\nfirst / last', 'Near-terminal\nsecond / penultimate',
                                'Interior\nremaining'])
        xlabel = 'Monosaccharide repeat position (dp12–dp30)'

    ax.axhline(0, color=ZERO_COLOR, linewidth=.8, zorder=0)
    if settings['show_title']:
        title = ANALYSES[panel].replace('_', ' ').capitalize()
        ax.set_title(f'{panel}  Heparin: {title}\nAtomic {METRIC_LABELS[metric]}',
                     loc='left', fontweight='bold', pad=12)
    ax.set_xlabel(xlabel, fontsize=settings['xlabel_size'])
    ax.set_ylabel(f'Signed {METRIC_LABELS[metric]} difference\nfrom AW (%)', fontsize=settings['ylabel_size'])
    ax.tick_params(axis='x', labelsize=settings['xtick_size'])
    ax.tick_params(axis='y', labelsize=settings['ytick_size'])
    ax.xaxis.get_offset_text().set_fontsize(settings['xtick_size'])
    ax.yaxis.get_offset_text().set_fontsize(settings['ytick_size'])
    ax.spines[['top', 'right']].set_visible(False)
    # Leave room above the SD endpoints for the vertical observation counts.
    if panel != 'A':
        low, high = ax.get_ylim()
        ax.set_ylim(low, high + .50 * (high - low))
    if settings['xlim'] is not None:
        ax.set_xlim(settings['xlim'])
    if settings['ylim'] is not None:
        ax.set_ylim(settings['ylim'])
    if settings['show_legend']:
        handles = [Line2D([0], [0], marker='o', linestyle='none', color=COLORS[role], label=role,
                          markersize=settings['marker_size']) for role in roles]
        fig.legend(handles=handles, loc=settings['legend_loc'], ncol=3,
                   fontsize=LEGEND_SIZE, frameon=False)
    systems = 'dp12, dp18, dp24, dp30 (dp2 excluded)' if panel == 'C' else 'dp2, dp12, dp18, dp24, dp30'
    error_note = 'Error bars: sample SD. ' if settings['show_error_bars'] else ''
    fig.text(.5, .025,
             f'{systems}\nn: atom/frame observations. {error_note}Frames are correlated.\n'
             'Pooled statistics weight longer chains more heavily.',
             ha='center', fontsize=COUNT_SIZE)
    # Reserve space for the larger two-line titles, rotated labels, and notes.
    fig.subplots_adjust(left=.17, right=.98, top=.72 if settings['show_legend'] else .85,
                        bottom=.35 if panel == 'B' else .30)
    for extension in ('png', 'svg'):
        path = OUTPUT_DIR / f'F5B_{key}.{extension}'
        fig.savefig(path, dpi=DPI, facecolor='white')
    print(f'F5B_{key} (PNG + SVG)', flush=True)
    plt.close(fig)


def main():
    config = Config()
    frames = []
    for model in MODELS:
        folders = sorted((INPUT_ROOT / model).glob('frame_*'))
        if not folders:
            raise FileNotFoundError(f'No frames for {model} under {INPUT_ROOT}')
        model_points, excluded = [], {}
        for folder in folders:
            # Fail visibly on missing or inconsistent inputs rather than silently
            # dropping a requested model/frame from the figure.
            points, counts = load_frame(folder, model, config)
            model_points.append(points)
            for reason, count in counts.items():
                excluded[reason] = excluded.get(reason, 0) + count
        combined = pd.concat(model_points, ignore_index=True)
        if not combined.atom_class.isin(COLORS).any():
            raise ValueError(f'{model}: no validated plotted observations')
        print(f'{model}: {len(folders)} frames, {len(combined)} matched valid atom/frame observations; '
              f'{excluded["total_excluded"]} pairs excluded from both metrics', flush=True)
        print(f'  Exclusion reasons (may overlap): {excluded}', flush=True)
        frames.append(combined)
    points = pd.concat(frames, ignore_index=True)
    print('Missing input files/frames: none detected in the discovered frame folders.', flush=True)
    summary = summarize(points)
    print('Matched valid observations by chemical class:', points.atom_class.value_counts().to_dict(), flush=True)
    unresolved = points[points.atom_class.eq('Unresolved')]
    print('Unresolved chemical classifications:', len(unresolved), flush=True)
    if not unresolved.empty:
        print(unresolved.groupby(['model', 'residue_name', 'atom_name', 'classification_reason']).size().to_string(), flush=True)
    print('Unresolved chain positions:', int(points.position.eq('Unresolved').sum()), flush=True)
    print('Hydrogen and nitrogen are calculated/audited but not plotted; sulfur is shown only in B; '
          'dp2 is excluded from C.', flush=True)
    for panel in ANALYSES:
        included = summary.loc[summary.model.eq('pooled') & summary.panel.eq(panel), 'n'].sum()
        print(f'{panel}, each metric: {int(included)} plotted observations; '
              f'{len(points) - int(included)} excluded by class/region selection', flush=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / 'F5B_heparin_chemistry_summary.csv'
    summary.to_csv(path, index=False)
    print(path, flush=True)
    for metric in METRIC_LABELS:
        for panel in ANALYSES:
            plot_figure(points, summary, metric, panel)
    print(f'Saved all 12 requested figure files to {OUTPUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
