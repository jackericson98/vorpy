"""Screened PEG paths and identity-checked first-layer membership; no solver calls."""
from pathlib import Path

import pandas as pd

from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5_common import read_pdb


REPOSITORY_ROOT = next(parent for parent in Path(__file__).resolve().parents
                       if (parent / 'pyproject.toml').is_file())
STATUS_PATH = REPOSITORY_ROOT / 'output/peg_one_layer_investigation_20261112/one_layer_run_status.csv'
SCHEMES = ('aw', 'pow', 'prm')


def screened_runs():
    rows = pd.read_csv(STATUS_PATH, dtype=str).fillna('')
    if rows.duplicated(['molecule', 'frame', 'scheme']).any():
        raise ValueError('Duplicate screened run identities')
    valid = rows.usable_for_screened_layer1_atom_analysis.str.lower().eq('true')
    selected = []
    for _, group in rows.groupby(['molecule', 'frame'], sort=True):
        if len(group) == 3 and set(group.scheme) == set(SCHEMES) and valid.loc[group.index].all():
            selected.append(group)
    result = pd.concat(selected, ignore_index=True)
    if len(result) != 69 * 3:
        raise ValueError('Expected the documented 69 screened three-scheme frame sets')
    recovered = result[(result.molecule == 'PEG3000') & (result.frame == 'frame_0005') & (result.scheme == 'aw')]
    expected = STATUS_PATH.parent / 'recovered_PEG3000_aw_frame0005_v2'
    if Path(recovered.iloc[0].analysis_directory).resolve() != expected.resolve():
        raise ValueError('PEG3000 frame 5 must use the verified recovered AW bundle')
    return result


def paths(row, *, solvent=False):
    root = Path(row.analysis_directory)
    if solvent:
        log = root / 'solvent' / f'{row.scheme}_logs.csv'
    else:
        candidates = [root / 'groups' / name / f'{row.scheme}_logs.csv'
                      for name in (row.molecule, f'{row.molecule}_full_system')]
        existing = [path for path in candidates if path.is_file()]
        if len(existing) != 1:
            raise ValueError(f'Expected one polymer group log: {root}')
        log = existing[0]
    pdb = root / f'{row.molecule}.pdb'
    if not log.is_file() or not pdb.is_file():
        raise ValueError(f'Screened log/source PDB missing: {root}')
    return log, pdb


def layer1_atoms(row, pdb):
    """Join layer PDB records to full source identity, never serial-minus-one."""
    columns = ['atom_name', 'residue_name', 'chain', 'residue_sequence', 'insertion', 'x', 'y', 'z']
    keys = {}
    for atom in pdb.itertuples():
        key = tuple(getattr(atom, column) for column in columns)
        if key in keys:
            raise ValueError('Ambiguous source atom identity')
        keys[key] = atom.atom_index
    root = Path(row.analysis_directory) / 'solvent'
    files = list(root.glob('layer_01.pdb'))
    if len(files) != 1:
        raise ValueError(f'Expected one layer-1 PDB: {root}')
    selected = []
    with files[0].open() as stream:
        for line in stream:
            if line[:6].strip() not in {'ATOM', 'HETATM'}:
                continue
            key = (line[12:16].strip(), line[17:21].strip(), line[21:22].strip(),
                   int(line[22:26]), line[26:27].strip(),
                   float(line[30:38]), float(line[38:46]), float(line[46:54]))
            index = keys.get(key)
            if index is None or int(line[6:11]) != pdb.loc[index, 'pdb_serial']:
                raise ValueError(f'Unverified layer atom identity: {files[0]}')
            selected.append(index)
    atoms = pdb.loc[selected]
    if (not selected or len(set(selected)) != len(selected) or not atoms.kind.eq('water').all()
            or len(selected) != int(row.layer1_atoms)):
        raise ValueError('Empty, duplicate or non-water layer membership')
    names = atoms.groupby('residue_id').atom_name.agg(lambda group: sorted(group))
    if not names.map(lambda group: group == ['HW1', 'HW2', 'OW']).all():
        raise ValueError('Layer 1 must contain whole three-site waters')
    return set(selected)


def validated_water_layer(row, config):
    # Retain every incident face while restricting cells to layer 1.
    try:
        from ..F5C_water_contacts import water_contacts
    except ImportError:
        from .F5C_water_contacts import water_contacts
    log, source = paths(row, solvent=True)
    pdb = read_pdb(source, 'PEG', config)
    atoms = layer1_atoms(row, pdb)
    waters = water_contacts(log, source, config, layer_atoms=atoms)
    if len(waters) * 3 != len(atoms):
        raise ValueError(f'Screened layer lost incomplete/inconsistent waters: {log}')
    return waters, pdb, atoms


def water_identity(pdb, residue):
    columns = ['atom_name', 'residue_name', 'chain', 'residue_sequence', 'insertion', 'x', 'y', 'z']
    return tuple(pdb.loc[pdb.residue_id.eq(residue)].sort_values('atom_name')[columns].itertuples(index=False, name=None))


def matched_water_layers(rows, config):
    """Use the three-scheme intersection of validated whole-water identities."""
    import numpy as np
    from vorpy.src.analyze.Part_II_Molecular_Analysis.Shared_Utilities.simple_models.F5B_water_burial import load_water_geometry
    layers, identities = {}, {}
    for row in rows.itertuples(index=False):
        checked, pdb, atoms = validated_water_layer(row, config)
        log, _ = paths(row, solvent=True)
        geometry = load_water_geometry(log, pdb, 'PEG', config, layer_atoms=atoms).set_index('residue_id')
        if set(geometry.index) != set(checked.index) or not np.isclose(
                geometry.loc[checked.index].external_area, checked.external_area,
                rtol=config.area_rtol, atol=config.area_atol).all():
            raise ValueError('Cross-scheme reader disagrees with validated water external area')
        identity = {water_identity(pdb, residue): residue for residue in geometry.index}
        if len(identity) != len(geometry):
            raise ValueError('Ambiguous water identity')
        layers[row.scheme], identities[row.scheme] = geometry, identity
    common = set.intersection(*(set(identities[scheme]) for scheme in SCHEMES))
    matched = {}
    for scheme in SCHEMES:
        keys = sorted(common)
        table = layers[scheme].loc[[identities[scheme][key] for key in keys]].copy()
        # Keep the AW occurrence label only after matching the full PDB identity.
        table.index = [identities['aw'][key] for key in keys]
        table.index.name = 'residue_id'
        matched[scheme] = table.reset_index()
    return matched
