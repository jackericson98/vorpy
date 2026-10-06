"""Export the persisted 2KAI AW example without solving or recalculating births.

Requires the existing full_aw.vpy, matched-probe birth checkpoint, identity
audit, and saved curvature component tables. No replacement data is invented.
"""
import csv
from dataclasses import replace
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vorpy.src.geometry.duals import build_dual
from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex
from vorpy.src.geometry.interfaces import build_alpha_interface
from vorpy.src.geometry.visualization.interface import export_interface_dual_visualization
from vorpy.src.interface.interface import Interface
from vorpy.src.io import load_network


def read(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def main():
    cache = ROOT / 'output/comparative_study/matched_probe_control/2KAI'
    net = load_network(ROOT / 'output/comparative_study/cases/2KAI/H/full_aw.vpy')
    atoms = read(cache / 'input_identity_audit.csv')
    for row in atoms:
        generator = int(row['generator_id'])
        np.testing.assert_allclose(net.balls.loc[generator, 'loc'],
            [float(row[key]) for key in ('x_A', 'y_A', 'z_A')], rtol=0., atol=1e-12)
        np.testing.assert_allclose(net.balls.loc[generator, 'rad'], float(row['base_radius_A']), rtol=0., atol=1e-12)
    a = frozenset(int(row['generator_id']) for row in atoms if row['chain'] in {'A', 'B'})
    b = frozenset(int(row['generator_id']) for row in atoms if row['chain'] == 'I')
    stable_metadata = {int(row['generator_id']): {
        'chain': row['chain'], 'residue_number': row['residue_number'],
        'insertion_code': row['insertion_code'], 'residue_name': row['residue_name'],
        'atom_name': row['atom_name'], 'stable_atom_id': row['stable_atom_id'], 'element': row['element']} for row in atoms}
    checkpoint = read(cache / 'aw_pair_birth_checkpoint.csv')
    thresholds = {float(row['threshold']) for row in checkpoint}
    if len(thresholds) != 1:
        raise ValueError('Persisted AW checkpoint has inconsistent alpha queries')
    alpha = thresholds.pop()
    records = {d: {} for d in range(4)}
    for row in checkpoint:
        pair = tuple(sorted((int(row['generator_i']), int(row['generator_j']))))
        birth = float(row['surface_birth_A']) if row['surface_birth_A'] else None
        records[1][pair] = AlphaSimplex(1, pair, birth, birth, row['birth_source'], birth is not None,
            notes=row['birth_reason'], diagnostics={
                'restriction_status': row['restriction_status'],
                'restriction_certificate': row['restriction_certificate'],
                'birth_lower_bound_A': float(row['birth_lower_bound_A']) if row['birth_lower_bound_A'] else None})
    filtration = AlphaFiltration(net, 'aw', records, {
        'native_alpha_units': 'A', 'filtration_kind': 'persisted experimental AW restricted-cell pair births',
        'birth_dimensions_calculated': 1, 'stable_atom_metadata': stable_metadata,
        'source': str(cache / 'aw_pair_birth_checkpoint.csv')}, tolerance=1e-6)
    dual = build_dual(net)
    rep = build_alpha_interface(net, dual, filtration, alpha, a, b, calculate_curvature=False)
    selected = {row.generator_tuple for row in filtration.interface_at(alpha, a, b)}
    persisted = {tuple(sorted((int(row['generator_i']), int(row['generator_j'])))) for row in checkpoint
                 if row['restriction_status'] == 'retained' and
                 ((int(row['generator_i']) in a and int(row['generator_j']) in b) or
                  (int(row['generator_j']) in a and int(row['generator_i']) in b))}
    if selected != persisted:
        raise ValueError('Cached query does not reproduce the persisted selected pairs')
    curvature = {int(row['surface_id']): row for row in read(
        ROOT / 'output/2KAI_curvature_showcase/aw_surface_curvature_components.csv')}
    rep.surfaces = [replace(surface, smooth_H_A=(float(curvature[surface.feature_id]['integral_H_dA_A'])
        if surface.feature_id in curvature and curvature[surface.feature_id]['support_status'] == 'included'
        and curvature[surface.feature_id]['integral_H_dA_A'] else None)) for surface in rep.surfaces]
    iface = Interface.__new__(Interface)
    iface.net, iface.name = net, '2KAI_cached_AB_I'
    iface.dir = str(ROOT / 'output/apollonius_interface_cached_2KAI')
    iface.group1, iface.group2 = SimpleNamespace(name='A+B', group_id='A+B'), SimpleNamespace(name='I', group_id='I')
    iface._geometry_dual, iface._geometry_filtration, iface._geometry_representation = dual, filtration, rep
    iface.geometry_analysis = SimpleNamespace(metadata={'alpha_value': alpha}, alpha_selection={
        'selected_generator_pairs': [list(pair) for pair in sorted(selected)],
        'eligible_bicolor_pairs': sum(surface.supported and surface.complete and surface.bounded
                                      and surface.area_A2 is not None for surface in rep.candidate_surfaces)})
    # Fail loudly if any visualization path tries to solve again.
    from vorpy.src.network import Network
    def forbidden(*args, **kwargs):
        raise RuntimeError('Cached visualization attempted a network solve')
    Network.build = forbidden
    result = export_interface_dual_visualization(iface)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
