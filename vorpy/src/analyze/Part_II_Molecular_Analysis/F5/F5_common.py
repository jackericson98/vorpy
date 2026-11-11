"""F5 scientific tables. Source logs and solver state are never modified."""
from collections import Counter
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import csv
import json
import pickle
import re

import numpy as np
import pandas as pd

if not __package__:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[5]))
    __package__ = 'vorpy.src.analyze.Part_II_Molecular_Analysis.F5'

from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from .F5_config import Config, KEY, ATOM_KEY, UNIT_KEY, METRICS


def classify_atoms(atoms, polymer, config):
    """Residue AND atom names identify water; unsupported water sites are flagged."""
    a = atoms.copy()
    res = a.residue_name.str.upper()
    name = a.atom_name.str.upper()
    water_res = res.isin(config.water_residues)
    water_name = name.isin(config.water_oxygens + config.water_hydrogens)
    polymer_res = config.peg_residues if polymer == 'PEG' else config.heparin_residues
    a['kind'] = np.select([res.isin(polymer_res), water_res & water_name,
                           res.isin(config.ion_residues)], ['polymer', 'water', 'ion'], default='other')
    a['identity_valid'] = ~(water_res & ~water_name)
    a['is_hydrogen'] = a.element.str.upper().isin(['H', 'D'])
    a['atom_class'] = np.select([a.is_hydrogen, a.element.eq('C'), a.element.eq('O'),
                                 a.element.eq('S'), a.element.eq('N')],
                                ['hydrogen', 'carbon', 'oxygen', 'sulfur', 'nitrogen'], default='other')
    return a


def read_pdb(path, polymer, config):
    """Read coordinate identity only; retain order, serial, full residue name and wraps.

    This is deliberately not a VorPy log parser. Molecule occurrence IDs avoid
    conflating wrapped PDB residue numbers. Only a single MODEL is permitted.
    """
    rows = []
    previous = None
    occurrence = -1
    models = 0
    segment = 0
    with Path(path).open(encoding='utf-8-sig') as stream:
        for line in stream:
            record = line[:6].strip()
            if record == 'MODEL':
                models += 1
                if models > 1:
                    raise ValueError('frame.pdb contains multiple MODELs')
            if record == 'TER':
                segment += 1
                previous = None
            if record not in {'ATOM', 'HETATM'}:
                continue
            if line[16:17].strip() not in {'', 'A'}:
                raise ValueError('Alternate PDB conformers require an explicit selection')
            name, residue, chain = line[12:16].strip(), line[17:21].strip(), line[21:22].strip()
            seq, insertion = int(line[22:26]), line[26:27].strip()
            identity = (segment, chain, residue, seq, insertion)
            if identity != previous:
                occurrence += 1
                previous = identity
            element = line[76:78].strip().upper()
            if not element:
                element = re.sub(r'^[0-9]+', '', name).upper()[:1]
            rows.append((len(rows), int(line[6:11]), name, residue, seq, chain, insertion,
                         f'{segment}:{occurrence}', element,
                         float(line[30:38]), float(line[38:46]), float(line[46:54])))
    a = pd.DataFrame(rows, columns=['atom_index', 'pdb_serial', 'atom_name', 'residue_name',
        'residue_sequence', 'chain', 'insertion', 'residue_id', 'element', 'x', 'y', 'z'])
    if a.empty or not np.isfinite(a[['x', 'y', 'z']]).all().all():
        raise ValueError('Empty or nonfinite PDB coordinates')
    a = classify_atoms(a, polymer, config).set_index('atom_index', drop=False)
    water = a[a.kind.eq('water')]
    counts = water.groupby('residue_id').atom_name.agg(list)
    invalid_molecules = {i for i, names in counts.items()
                         if sum(n.upper() in config.water_oxygens for n in names) != 1
                         or sum(n.upper() in config.water_hydrogens for n in names) not in {0, 2}}
    a.loc[a.residue_id.isin(invalid_molecules), 'identity_valid'] = False
    return annotate_units(a, polymer)


def end_distances(order):
    order = np.asarray(order)
    return np.minimum(order, len(order) - 1 - order)


def annotate_units(atoms, polymer):
    """Validate a PEG O(CCO)n path or Heparin inter-residue C1--O4 path.

    No periodic unwrapping is guessed: discontinuous paths stay unavailable.
    """
    a = atoms.copy()
    a['unit_id'] = a.residue_id
    a['repeat_index'] = np.nan
    a['repeat_count'] = np.nan
    a['end_distance'] = np.nan
    a['order_valid'] = False
    a['order_reason'] = 'not a polymer atom'
    p = a[a.kind.eq('polymer')]
    if p.empty:
        raise ValueError('No recognized polymer residues in PDB')
    a.loc[p.index, 'order_reason'] = 'backbone topology not validated'
    if polymer == 'PEG':
        # Each PEG molecule must be represented by one residue in these models.
        for residue, group in p.groupby('residue_id', sort=False):
            heavy = group[~group.is_hydrogen].sort_index()
            n = (len(heavy) - 1) // 3
            expected = ['O'] + ['C', 'C', 'O'] * n
            if n < 1 or heavy.element.tolist() != expected:
                continue
            names = [f'{el}{i+1}' for i, el in enumerate(expected)]
            if heavy.atom_name.tolist() != names:
                continue
            xyz = heavy[['x', 'y', 'z']].to_numpy()
            distances = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
            if not np.all((distances >= .9) & (distances <= 1.85)):
                continue
            units = np.maximum(0, (np.arange(len(heavy)) - 1) // 3)
            assignments = dict(zip(heavy.index, units))
            valid = True
            for idx, atom in group[group.is_hydrogen].iterrows():
                dist = np.linalg.norm(xyz - atom[['x', 'y', 'z']].to_numpy(dtype=float), axis=1)
                nearest = int(dist.argmin())
                if not .6 <= dist[nearest] <= 1.35:
                    valid = False
                    break
                assignments[idx] = units[nearest]
            if not valid:
                continue
            for idx, unit in assignments.items():
                a.loc[idx, ['unit_id', 'repeat_index', 'repeat_count', 'end_distance', 'order_valid', 'order_reason']] = [
                    f'{residue}:repeat{unit}', unit, n, min(unit, n-1-unit), True, 'validated O(CCO)n path']
            internal_oxygen = heavy.index[heavy.element.eq('O')][1:-1]
            a.loc[internal_oxygen, 'atom_class'] = 'PEG ether oxygen'
            a.loc[heavy.index[[0, -1]], 'atom_class'] = 'PEG terminal oxygen'
    else:
        # Discover connectivity across all residues, rather than sorting sequence numbers.
        residues = list(p.residue_id.drop_duplicates())
        graph = {r: set() for r in residues}
        for _, c in p[p.atom_name.eq('C1')].iterrows():
            for _, o in p[p.atom_name.eq('O4')].iterrows():
                if c.residue_id == o.residue_id or c.chain != o.chain:
                    continue
                dist = np.linalg.norm(c[['x', 'y', 'z']].to_numpy(dtype=float) - o[['x', 'y', 'z']].to_numpy(dtype=float))
                if .9 <= dist <= 1.9:
                    graph[c.residue_id].add(o.residue_id)
                    graph[o.residue_id].add(c.residue_id)
        ends = [r for r in residues if len(graph[r]) == 1]
        if len(ends) == 2 and all(1 <= len(g) <= 2 for g in graph.values()):
            order, previous, current = [], None, ends[0]
            while current is not None and current not in order:
                order.append(current)
                following = graph[current] - {previous}
                previous, current = current, next(iter(following), None)
            if len(order) == len(residues):
                for i, residue in enumerate(order):
                    a.loc[p.index[p.residue_id.eq(residue)], ['repeat_index', 'repeat_count', 'end_distance', 'order_valid', 'order_reason']] = [
                        i, len(order), min(i, len(order)-1-i), True, 'validated inter-residue C1--O4 path']
    return a


def audit_raw_log(path):
    """Narrow integrity pass, not another log reader. Reject silent numeric coercions.

    read_logs2 defaults missing surface areas to zero; record those IDs to restore
    NaN after reading. Missing atom numerics are already omitted by the reader.
    """
    section, headers = '', []
    counts = Counter()
    missing_surface_areas = set()
    sections = {'atoms', 'surfaces', 'edges', 'vertices', 'group information', 'build information', 'build informaiton'}
    numeric = {'index', 'ball 1', 'ball 2', 'volume', 'surface area', 'contact area', 'x', 'y', 'z',
               'residue sequence', 'euler characteristic', 'integrated gaussian curvature',
               'gauss-bonnet expected', 'gauss-bonnet error', 'boundary components'}
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.reader(stream)
        for row in reader:
            if not row or not any(v.strip() for v in row):
                continue
            first = row[0].strip().lower()
            if first in sections:
                section = first
                headers = [v.strip().lower() for v in next(reader)]
                if section in {'atoms', 'surfaces'}:
                    required = {'index', 'name', 'residue', 'residue sequence', 'x', 'y', 'z', 'surface area'} if section == 'atoms' else {'index', 'ball 1', 'ball 2', 'surface area'}
                    if not required.issubset(headers):
                        raise ValueError(f'{section}: missing headers {required-set(headers)}')
                continue
            if section not in {'atoms', 'surfaces', 'group information'}:
                continue
            counts[section] += 1
            data = dict(zip(headers, row))
            for key in numeric.intersection(headers):
                value = data.get(key, '').strip()
                if not value:
                    if key in {'index', 'ball 1', 'ball 2', 'residue sequence', 'x', 'y', 'z'} and section in {'atoms', 'surfaces'}:
                        raise ValueError(f'{section}: missing {key}')
                    continue
                number = float(value)
                if not np.isfinite(number):
                    raise ValueError(f'{section}: nonfinite {key}')
                if key in {'index', 'ball 1', 'ball 2', 'residue sequence'} and number != int(number):
                    raise ValueError(f'{section}: noninteger {key}')
            for key in ('complete cell?', 'boundary complete', 'boundary closed', 'boundary manifold', 'boundary orientable'):
                if data.get(key, '').strip().lower() not in {'', 'true', 'false'}:
                    raise ValueError(f'Invalid boolean {key}')
            if section == 'surfaces' and not data.get('surface area', '').strip():
                missing_surface_areas.add(int(data['index']))
    if not counts['atoms'] or not counts['surfaces']:
        raise ValueError('Empty atom or surface section')
    return counts, missing_surface_areas


def normalize_log(path, pdb, config):
    counts, missing = audit_raw_log(path)
    data = read_logs2(path, all_=False, balls=True, surfs=True)
    a = data['atoms'].rename(columns={'Index': 'atom_index', 'Name': 'atom_name', 'Residue': 'log_residue_name',
        'Residue Sequence': 'residue_sequence', 'Chain': 'chain', 'X': 'x', 'Y': 'y', 'Z': 'z',
        'Volume': 'volume', 'Surface Area': 'surface_area', 'Complete Cell?': 'complete', 'Neighbors': 'logged_neighbors'})
    if len(a) != counts['atoms'] or len(data['surfs']) != counts['surfaces']:
        raise ValueError('Reader/raw row count mismatch')
    if a.atom_index.duplicated().any() or not a.atom_index.isin(pdb.index).all():
        raise ValueError('Duplicate atom indices or unresolved PDB indices')
    p = pdb.loc[a.atom_index].reset_index(drop=True)
    a = a.reset_index(drop=True)
    if 'chain' not in a:
        a['chain'] = ''
    if not a.atom_name.eq(p.atom_name).all() or not a.residue_sequence.eq(p.residue_sequence).all():
        raise ValueError('Log index does not match PDB atom/residue identity')
    if not a.log_residue_name.str.strip().eq(p.residue_name.str[:3]).all():
        raise ValueError('Log/PDB residue names disagree (including 3-column truncation)')
    if not a.chain.fillna('').str.strip().eq(p.chain).all():
        raise ValueError('Log/PDB chains disagree')
    if not np.allclose(a[['x','y','z']], p[['x','y','z']], atol=config.coordinate_tolerance, rtol=0):
        raise ValueError('Log/PDB coordinates disagree; no index shift or wrapping guessed')
    for col in ('volume', 'surface_area', 'complete'):
        if col not in a:
            a[col] = np.nan
    for col in p.columns:
        if col not in a:
            a[col] = p[col].to_numpy()
    a['complete'] = a.complete.astype('boolean')
    a = a.set_index('atom_index', drop=False)
    s = data['surfs'].rename(columns={'Index':'surface_index', 'Surface Area':'surface_area', 'Contact Area':'contact_area'})
    s['atom_i'] = s.Balls.map(lambda pair: int(pair[0]))
    s['atom_j'] = s.Balls.map(lambda pair: int(pair[1]))
    s.loc[s.surface_index.isin(missing), 'surface_area'] = np.nan
    if 'contact_area' not in s:
        s['contact_area'] = np.nan
    return a, s[['surface_index','atom_i','atom_j','surface_area','contact_area']], data['group data'], counts


def reconstruct_interfaces(atoms, surfaces, context, config):
    """Unique undirected pairs and incident areas. Row order carries no identity."""
    s = surfaces.copy()
    if s.atom_i.eq(s.atom_j).any():
        raise ValueError('Self-interface in surface records')
    s[['atom_i','atom_j']] = np.sort(s[['atom_i','atom_j']].to_numpy(), axis=1)
    if not s.atom_i.isin(context.index).all() or not s.atom_j.isin(context.index).all():
        raise ValueError('Surface endpoint absent from PDB')
    duplicate_count = int(s.duplicated(['atom_i','atom_j']).sum())
    duplicates = s[s.duplicated(['atom_i','atom_j'], keep=False)]
    for _, group in duplicates.groupby(['atom_i','atom_j']):
        for col in ['surface_area','contact_area']:
            if not np.allclose(group[col], group[col].iloc[0], rtol=config.area_rtol, atol=config.area_atol, equal_nan=True):
                raise ValueError('Conflicting repeated surface pair; disconnected faces cannot be inferred')
    # One surface ID cannot denote different pairs, even when the numeric areas agree.
    if s.drop_duplicates(['surface_index','atom_i','atom_j']).surface_index.duplicated().any():
        raise ValueError('Surface index denotes multiple endpoint pairs')
    s = s.drop_duplicates(['atom_i','atom_j']).reset_index(drop=True)
    if s.surface_area.lt(0).any():
        raise ValueError('Negative interface area')
    for suffix in ['i','j']:
        endpoints = context.loc[s[f'atom_{suffix}']].reset_index(drop=True)
        for col in ['kind','residue_id','unit_id','repeat_index','chain','atom_name','residue_name','residue_sequence',
                    'identity_valid','x','y','z']:
            s[f'{col}_{suffix}'] = endpoints[col].to_numpy()
    pairkind = s.apply(lambda r: tuple(sorted((r.kind_i, r.kind_j))), axis=1)
    s['interface_type'] = pairkind.map({('polymer','polymer'):'polymer–polymer', ('polymer','water'):'polymer–water',
        ('ion','polymer'):'polymer–other', ('other','polymer'):'polymer–other', ('water','water'):'water–water'}).fillna('other')
    incidence = []
    for side, other in [('i','j'),('j','i')]:
        part = pd.DataFrame({'atom_index':s[f'atom_{side}'], 'neighbor':s[f'atom_{other}'],
            'neighbor_kind':s[f'kind_{other}'], 'neighbor_residue':s[f'residue_id_{other}'],
            'neighbor_unit':s[f'unit_id_{other}'], 'neighbor_repeat':s[f'repeat_index_{other}'],
            'neighbor_chain':s[f'chain_{other}'], 'area':s.surface_area,
            'identity_valid':s.identity_valid_i & s.identity_valid_j})
        incidence.append(part)
    inc = pd.concat(incidence, ignore_index=True)
    out = atoms.copy()
    by_atom = dict(tuple(inc.groupby('atom_index', sort=False)))
    derived = []
    for idx, atom in out.iterrows():
        values = {'atom_index':idx}
        edges = by_atom.get(idx, inc.iloc[:0])
        areas = {}
        for kind in ['water','polymer','other']:
            mask = edges.neighbor_kind.eq(kind) if kind != 'other' else ~edges.neighbor_kind.isin(['water','polymer'])
            vals = edges.loc[mask, 'area']
            areas[kind] = vals.sum() if vals.notna().all() else np.nan
            values[f'{kind}_area'] = areas[kind]
        observed = edges.area.sum() if edges.area.notna().all() else np.nan
        values['interface_area_sum'] = observed
        partition = bool(np.isfinite(atom.surface_area) and atom.surface_area > 0 and
                         np.isclose(observed, atom.surface_area, rtol=config.area_rtol, atol=config.area_atol) and
                         edges.identity_valid.all() and atom.identity_valid)
        complete = pd.notna(atom.complete) and bool(atom.complete)
        values['partition_valid'] = partition
        values['fraction_valid'] = partition and complete
        values['primary_valid'] = partition and (complete or config.include_incomplete)
        for kind, area in areas.items():
            values[f'{kind}_fraction'] = area / atom.surface_area if pd.notna(atom.surface_area) and atom.surface_area > 0 else np.nan
        active = edges[edges.area > config.minimum_interface_area]
        water = active[active.neighbor_kind.eq('water')]
        poly = active[active.neighbor_kind.eq('polymer')]
        waters = frozenset(water.neighbor_residue)
        neighbors = frozenset(map(int, active.neighbor))
        polymers = frozenset(map(int, poly.neighbor))
        for col, value in [('water_coordination',len(waters)), ('polymer_coordination',len(polymers)),
                           ('polymer_residue_coordination',poly.neighbor_residue.nunique()),
                           ('same_residue_contacts',int(poly.neighbor_residue.eq(atom.residue_id).sum())),
                           ('same_repeat_contacts',int(poly.neighbor_unit.eq(atom.unit_id).sum()))]:
            values[col] = value
        local = poly.neighbor_chain.eq(atom.chain) & poly.neighbor_repeat.sub(atom.repeat_index).abs().eq(1)
        nonlocal_mask = (~poly.neighbor_chain.eq(atom.chain) | poly.neighbor_repeat.sub(atom.repeat_index).abs().gt(1))
        values['adjacent_repeat_contacts'] = int(local.sum()) if atom.order_valid else np.nan
        values['nonlocal_contacts'] = int(nonlocal_mask.sum()) if atom.order_valid else np.nan
        values['area_per_water'] = areas['water']/len(waters) if waters else np.nan
        values['area_per_polymer'] = areas['polymer']/len(polymers) if polymers else np.nan
        # Set-valued columns are kept in pickle caches; CSVs use explicit JSON arrays.
        for col, value in [('water_neighbors',waters),('polymer_neighbors',polymers),('neighbors',neighbors),
                           ('polymer_residue_neighbors',frozenset(poly.neighbor_residue)),
                           ('nonlocal_unit_neighbors',frozenset(poly.loc[nonlocal_mask,'neighbor_unit']) if atom.order_valid else frozenset())]:
            values[col] = value
        values['logged_neighbor_count'] = len(atom.get('logged_neighbors', [])) if isinstance(atom.get('logged_neighbors'),list) else np.nan
        values['surface_neighbor_count'] = edges.neighbor.nunique()
        values['neighbor_list_consistent'] = (set(atom.logged_neighbors) == set(edges.neighbor)) if isinstance(atom.get('logged_neighbors'),list) else False
        derived.append(values)
    out = out.reset_index(drop=True).merge(pd.DataFrame(derived),on='atom_index',validate='one_to_one')
    if not config.include_hydrogens:
        out = out[~out.is_hydrogen].copy()
    for col in ['partition_valid','fraction_valid','primary_valid','neighbor_list_consistent']:
        out[col] = out[col].astype(bool)
    return out, s, duplicate_count


def union(values):
    return frozenset().union(*values)


def aggregate_residues(atoms, unit=False):
    group_col = 'unit_id' if unit else 'residue_id'
    rows = []
    for identity, group in atoms.groupby(group_col, sort=False):
        row = {key:group.iloc[0][key] for key in ['polymer','model','frame','scheme','residue_name',
            'residue_sequence','chain','residue_id','repeat_count','order_valid','order_reason'] if key in group}
        row[group_col] = identity
        row['unit_id'] = identity if unit else group.iloc[0].unit_id
        row['repeat_index'] = group.repeat_index.iloc[0] if unit else np.nan
        row['end_distance'] = group.end_distance.iloc[0] if unit else np.nan
        row['atom_count'] = len(group)
        row['member_identity'] = frozenset(tuple(v) for v in group[
            ['atom_index','atom_name','residue_name','residue_sequence','chain','x','y','z']].itertuples(index=False,name=None))
        for col in ['volume','surface_area','water_area','polymer_area','other_area']:
            row[col] = group[col].sum(min_count=len(group))
        row['mean_atomic_volume'] = group.volume.mean() if group.volume.notna().all() else np.nan
        for kind in ['water','polymer','other']:
            row[f'{kind}_fraction'] = row[f'{kind}_area']/row['surface_area'] if row['surface_area'] > 0 else np.nan
        row['water_neighbors'] = union(group.water_neighbors)
        row['polymer_neighbors'] = union(group.polymer_neighbors) - set(group.atom_index)
        row['polymer_residue_neighbors'] = union(group.polymer_residue_neighbors) - set(group.residue_id)
        row['neighbors'] = union(group.neighbors) - set(group.atom_index)
        row['water_coordination'] = len(row['water_neighbors'])
        row['polymer_coordination'] = len(row['polymer_neighbors'])
        row['polymer_residue_coordination'] = len(row['polymer_residue_neighbors'])
        row['nonlocal_contacts'] = len(union(group.nonlocal_unit_neighbors)) if group.order_valid.all() else np.nan
        for col in ['same_residue_contacts','same_repeat_contacts','adjacent_repeat_contacts']:
            row[col] = group[col].sum(min_count=len(group))
        row['area_per_water'] = row['water_area']/row['water_coordination'] if row['water_coordination'] else np.nan
        row['area_per_polymer'] = row['polymer_area']/row['polymer_coordination'] if row['polymer_coordination'] else np.nan
        row['primary_valid'] = bool(group.primary_valid.all())
        row['complete'] = bool(group.complete.fillna(False).all())
        row['zero_water'] = row['water_area'] == 0
        rows.append(row)
    return pd.DataFrame(rows)


def exposure(values, config):
    """Provisional bins: [0,.1], (.1,.4], (.4,1]. Invalid values stay missing."""
    return pd.cut(values, [-1e-12,config.buried_threshold,config.exposed_threshold,1.],
                  labels=['Buried','Intermediate','Exposed']).astype(object)


def attach_aw_exposure(atoms, units, config):
    ref = units[units.scheme.eq('aw')].copy()
    ref['aw_water_fraction'] = ref.water_fraction.where(ref.primary_valid)
    ref['environment'] = exposure(ref.aw_water_fraction, config)
    ref['aw_zero_water'] = ref.water_area.eq(0) & ref.primary_valid
    ref = ref[UNIT_KEY+['aw_water_fraction','environment','aw_zero_water']]
    return (atoms.merge(ref,on=UNIT_KEY,how='left',validate='many_to_one'),
            units.merge(ref,on=UNIT_KEY,how='left',validate='many_to_one'))


def match_schemes(table, config, level='atom'):
    keys = ATOM_KEY if level == 'atom' else UNIT_KEY
    identity = ['atom_name','residue_name','residue_sequence','chain','x','y','z'] if level == 'atom' else ['residue_name','chain','atom_count','repeat_index','member_identity']
    columns=list(dict.fromkeys(keys+identity+METRICS+['scheme','primary_valid','neighbors']+
        [c for c in ['environment','atom_class','aw_water_fraction'] if c in table]))
    table=table[columns]
    reference = table[table.scheme.eq('aw')]
    pairs, reports = [], []
    for scheme in ('pow','prm'):
        other = table[table.scheme.eq(scheme)]
        if other.empty:
            continue
        merged = reference.merge(other,on=keys,suffixes=('_aw','_s'),how='outer',indicator=True,validate='one_to_one')
        matched = merged['_merge'].eq('both')
        identical = matched.copy()
        for col in identity:
            left, right = merged[f'{col}_aw'], merged[f'{col}_s']
            if col in {'x','y','z'}:
                identical &= np.isclose(left,right,atol=config.coordinate_tolerance,rtol=0)
            else:
                identical &= left.eq(right) | (left.isna() & right.isna())
        valid = identical & merged.primary_valid_aw.eq(True) & merged.primary_valid_s.eq(True)
        merged['identity_match'] = identical
        merged['pair_valid'] = valid
        merged['scheme'] = scheme
        merged['level'] = level
        for metric in METRICS:
            scale = 100 if metric == 'water_fraction' else 1
            merged[f'delta_{metric}'] = scale*(merged[f'{metric}_s']-merged[f'{metric}_aw']).where(valid)
            merged[f'abs_delta_{metric}'] = merged[f'delta_{metric}'].abs()
        merged['neighbor_jaccard'] = [len(a & b)/len(a | b) if a | b else 1. for a,b in zip(
                merged.neighbors_aw.map(lambda v:v if isinstance(v,frozenset) else frozenset()),
                merged.neighbors_s.map(lambda v:v if isinstance(v,frozenset) else frozenset()))]
        merged.loc[~valid, 'neighbor_jaccard'] = np.nan
        for metadata, rows in merged.groupby(KEY, dropna=False):
            reports.append(dict(zip(KEY,metadata),scheme=scheme,level=level,
                n_aw=int(rows['_merge'].isin(['both','left_only']).sum()),
                n_other=int(rows['_merge'].isin(['both','right_only']).sum()),
                n_matched=int(rows['_merge'].eq('both').sum()), n_identity_mismatch=int((rows['_merge'].eq('both') & ~rows.identity_match).sum()),
                n_valid=int(rows.pair_valid.sum())))
        pairs.append(merged)
    return pd.concat(pairs,ignore_index=True) if pairs else pd.DataFrame(), pd.DataFrame(reports)


TOPOLOGY = ['Euler Characteristic','Integrated Gaussian Curvature','Gauss-Bonnet Expected','Gauss-Bonnet Error',
    'Gauss-Bonnet Relative Error','Boundary Complete','Boundary Closed','Boundary Manifold',
    'Boundary Orientable','Boundary Components','Missing Faces','Missing Edges','Missing Vertices']


def topology_record(group):
    result = {col:group.get(col,np.nan) for col in TOPOLOGY}
    result['topology_available'] = pd.notna(result['Euler Characteristic'])
    result['gauss_bonnet_applicable'] = all(result[col] is True or isinstance(result[col],np.bool_) and result[col]
        for col in ['Boundary Complete','Boundary Closed','Boundary Manifold','Boundary Orientable'])
    result['internal_occupancy_status'] = 'unavailable: no validated envelope or solvent-connectivity segmentation'
    return result


def load_scheme(path, pdb, config, metadata):
    a, s, group, counts = normalize_log(path,pdb,config)
    a, s, duplicates = reconstruct_interfaces(a,s,pdb,config)
    if not a.kind.eq('polymer').all():
        raise ValueError('Expected polymer-scoped atom rows; mixed/all-system scope requires an explicit selection')
    for key, value in metadata.items():
        a[key], s[key] = value, value
    validation = dict(metadata, status='ok', atom_rows=len(a), surface_rows=len(s), duplicate_surfaces=duplicates,
        incomplete_atoms=int((~a.complete.fillna(False)).sum()), invalid_partitions=int((~a.partition_valid).sum()),
        excluded_primary=int((~a.primary_valid).sum()), unknown_water_sites=int((~pdb.identity_valid).sum()),
        neighbor_list_mismatches=int((~a.neighbor_list_consistent).sum()),
        repeat_order_valid=bool(a.order_valid.all()), repeat_count=a.repeat_count.max(),
        raw_atom_rows=counts['atoms'], raw_surface_rows=counts['surfaces'], index_convention='zero-based PDB atom record order',
        source=str(path), contact_area_differs=int((~np.isclose(s.surface_area,s.contact_area,equal_nan=True)).sum()))
    inventory=[]
    for (residue,kind),rows in pdb.groupby(['residue_name','kind']):
        inventory.append(dict({k:v for k,v in metadata.items() if k!='scheme'},residue_name=residue,kind=kind,n_atoms=len(rows),
            n_molecules=rows.residue_id.nunique(),n_hydrogens=int(rows.is_hydrogen.sum()),identity_invalid=int((~rows.identity_valid).sum())))
    return {'atoms':a, 'residues':aggregate_residues(a), 'units':aggregate_residues(a,unit=True),
            'interfaces':s, 'topology':pd.DataFrame([dict(metadata,**topology_record(group))]),
            'validation':pd.DataFrame([validation]), 'pdb_inventory':pd.DataFrame(inventory)}


def source_signature(paths, config):
    source = [Path(__file__),Path(__file__).with_name('F5_config.py'),Path(read_logs2.__code__.co_filename)]
    payload = {'config':asdict(config), 'inputs':[(str(p.resolve()),p.stat().st_size,p.stat().st_mtime_ns) for p in paths],
               'code':[(str(p),sha256(p.read_bytes()).hexdigest()) for p in source]}
    return sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest()


def write_csv(table, path):
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    out = table.copy()
    for col in out.select_dtypes('object'):
        if out[col].map(lambda v:isinstance(v,(set,frozenset,list,tuple))).any():
            out[col] = out[col].map(lambda v:json.dumps(sorted(v) if isinstance(v,(set,frozenset)) else v) if isinstance(v,(set,frozenset,list,tuple)) else v)
    out.to_csv(path,index=False)


def independent_area_check(path, atoms, context, config, metadata):
    """Independent raw-surface sums for five selected polymer atoms, per real log."""
    selected = atoms.sort_values('atom_index').iloc[np.unique(np.linspace(0,len(atoms)-1,min(5,len(atoms))).astype(int))]
    sums = {int(i):0. for i in selected.atom_index}
    seen, section, headers = set(), '', []
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        reader = csv.reader(stream)
        for row in reader:
            if not row:
                continue
            first = row[0].strip().lower()
            if first == 'surfaces':
                section, headers = 'surfaces', next(reader)
                continue
            if first in {'edges','vertices'}:
                section = ''
            if section != 'surfaces' or not row[0].strip():
                continue
            raw = dict(zip(headers,row))
            i,j = int(raw['Ball 1']),int(raw['Ball 2'])
            pair = tuple(sorted((i,j)))
            if pair in seen:
                continue
            seen.add(pair)
            area = float(raw['Surface Area']) if raw.get('Surface Area','').strip() else np.nan
            for left,right in [(i,j),(j,i)]:
                if left in sums and context.at[right,'kind'] == 'water':
                    sums[left] += area
    checks = []
    for idx, value in sums.items():
        expected = atoms.set_index('atom_index').at[idx,'water_area']
        checks.append(dict(metadata,atom_index=idx,raw_water_area=value,pipeline_water_area=expected,
            passed=bool(np.isclose(value,expected,rtol=config.area_rtol,atol=config.area_atol,equal_nan=True))))
    return pd.DataFrame(checks)


def frame_summary(table, groups, metrics):
    """One descriptive observation per model/frame/group; atoms are not replicates."""
    valid = table[table.primary_valid].copy()
    rows = []
    for keys, group in valid.groupby(groups,observed=True,dropna=False):
        keys = keys if isinstance(keys,tuple) else (keys,)
        row = dict(zip(groups,keys))
        row.update(n_observations=len(group), n_atoms=int(group.atom_count.sum()) if 'atom_count' in group else len(group),
                   n_units=group.unit_id.nunique())
        for metric in metrics:
            row[metric] = group[metric].mean()
            row[f'n_valid_{metric}'] = group[metric].notna().sum()
        rows.append(row)
    return pd.DataFrame(rows,columns=groups+['n_observations','n_atoms','n_units']+list(metrics)+[f'n_valid_{m}' for m in metrics])


def temporal_summary(frames, groups, metrics):
    rows = []
    for keys, group in frames.groupby(groups,observed=True,dropna=False):
        keys = keys if isinstance(keys,tuple) else (keys,)
        row = dict(zip(groups,keys),n_frames=group.frame.nunique())
        for metric in metrics:
            row[f'{metric}_mean'] = group[metric].mean()
            row[f'{metric}_std'] = group[metric].std(ddof=1)
            row[f'{metric}_n_frames'] = group.loc[group[metric].notna(),'frame'].nunique()
        rows.append(row)
    return pd.DataFrame(rows,columns=groups+['n_frames']+[f'{metric}_{suffix}' for metric in metrics for suffix in ['mean','std','n_frames']])


def save_figure(fig, path, config):
    import matplotlib.pyplot as plt
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.tight_layout(rect=(0,.04,1,.95))
    for extension in config.export_formats:
        fig.savefig(path.with_suffix('.'+extension),dpi=config.dpi,bbox_inches='tight')
    if config.show:
        plt.show()
    plt.close(fig)


def setup_plots(config):
    import matplotlib
    if not config.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':config.font_size,
        'axes.labelsize':config.font_size,'axes.titlesize':config.font_size+1,
        'legend.fontsize':max(config.font_size-2,6),'lines.linewidth':1.4,
        'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})


def empty_panel(ax, text='No valid observations'):
    ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes,fontsize=9)


def distribution_panels(frames, metrics, title, path, config, category='environment'):
    """Boxes of model/frame means; no atom-level pseudo-replication."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from .F5_config import UNITS, SCHEME_LABELS
    fig,axes=plt.subplots(2,3,figsize=config.figure_size,squeeze=False)
    categories=['Buried','Intermediate','Exposed']
    colors=dict(config.colors)
    for ax,metric in zip(axes.flat,metrics):
        found=False
        for k,scheme in enumerate(config.schemes):
            for i,cat in enumerate(categories):
                vals=frames.loc[frames.scheme.eq(scheme) & frames[category].eq(cat),metric].dropna() if not frames.empty else pd.Series(dtype=float)
                if vals.empty:
                    continue
                found=True
                bp=ax.boxplot(vals,positions=[i+(k-1)*.23],widths=.20,patch_artist=True,showfliers=False,
                    medianprops={'color':'black','linewidth':1.1})
                bp['boxes'][0].set_facecolor(colors[scheme]);bp['boxes'][0].set_alpha(.55)
        ax.set_xticks(range(3),categories,rotation=15)
        ax.set_ylabel(UNITS.get(metric,metric))
        if not found: empty_panel(ax)
    for ax in list(axes.flat)[len(metrics):]: ax.set_visible(False)
    fig.legend([Line2D([0],[0],color=colors[s],lw=3) for s in config.schemes],
               [SCHEME_LABELS[s] for s in config.schemes],loc='lower center',ncol=3)
    fig.suptitle(title+'\nBoxes describe model/frame means; frames may be correlated',fontsize=config.font_size+1)
    save_figure(fig,path,config)
