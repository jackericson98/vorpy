"""Serialize the frozen matched-probe fixtures through the production cache adapter.

No tessellation or alpha selection is recomputed. Saved physical polygons,
incidence and validated analyzer contributions are the inputs.
"""
from collections import defaultdict
import csv
import json
from pathlib import Path
import pickle
from types import SimpleNamespace

import pandas as pd
from scipy.spatial import cKDTree

from vorpy.src.geometry.interfaces.model import InterfaceSurface, InterfaceEdge, InterfaceRepresentation
from vorpy.src.interface.geometry_analysis import (
    InterfaceGeometryAnalysis, finite, summarize_physical_interface,
    reverse_perspective, write_geometry_info, write_geometry_log,
)

ROOT = Path(__file__).resolve().parents[1]
MATCHED = ROOT / 'output/comparative_study/matched_probe_control/2KAI'
COMPONENTS = ROOT / 'output/2KAI_curvature_showcase'


def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))


def load_fixture(scheme):
    pair_rows = rows(MATCHED / 'bicolor_pair_comparison.csv')
    pairs = {tuple(sorted((int(r['generator_i']), int(r['generator_j']))))
             for r in pair_rows if r['in_Power' if scheme == 'pow' else 'in_AW'] == 'True'}
    atoms = rows(MATCHED / 'input_identity_audit.csv')
    group_a = frozenset(int(r['generator_id']) for r in atoms if r['chain'] in {'A', 'B'})
    group_b = frozenset(int(r['generator_id']) for r in atoms if r['chain'] == 'I')
    prefix = 'power' if scheme == 'pow' else 'aw'
    surfaces = rows(COMPONENTS / f'{prefix}_surface_curvature_components.csv')
    edge_rows = rows(COMPONENTS / f'{prefix}_edge_curvature_components.csv')
    if scheme == 'pow':
        _, power = pickle.loads((MATCHED / 'power_geometry_cache.pkl').read_bytes())
        keys = list(power.power_vertices)
        positions = [power.power_vertices[k].position for k in keys]
        tree = cKDTree(positions)
        cache = json.loads((MATCHED / 'power_facet_area_cache.json').read_text())['areas']
        face_edges, face_vertices, edge_faces = {}, {}, defaultdict(list)
        edge_keys = {}
        for sid, pair in enumerate(sorted(pairs)):
            cycle = []
            for point in cache[f'{pair[0]}:{pair[1]}']['polygon_vertices']:
                distances, indices = tree.query(point, k=2)
                if distances[0] > 2e-5 or distances[1] < 1e-8:
                    raise ValueError('Unresolved Power dual vertex identity')
                cycle.append(int(indices[0]))
            face_vertices[sid] = cycle
            face_edges[sid] = []
            for i, a in enumerate(cycle):
                b = cycle[(i+1) % len(cycle)]
                key = tuple(sorted((a,b)))
                eid = edge_keys.setdefault(key, len(edge_keys))
                face_edges[sid].append(eid); edge_faces[eid].append(sid)
        edge_table = []
        for (a,b), eid in edge_keys.items():
            generators = tuple(sorted(set(keys[a]) & set(keys[b])))
            if len(generators) != 3:
                raise ValueError('Power physical edge is not a solved dual triangle')
            edge_table.append({'id': eid, 'verts':[a,b], 'balls':generators, 'surfs':edge_faces[eid]})
        by_pair = {(int(r['generator_i']),int(r['generator_j'])):r for r in surfaces}
        surfaces = [dict(by_pair[pair], surface_id=sid) for sid,pair in enumerate(sorted(pairs))]
        network = SimpleNamespace(settings={'net_type':'pow'},
            surfs=pd.DataFrame([{'balls':pair, 'edges':face_edges[sid], 'verts':face_vertices[sid],
                                'int_gauss_curv':0.} for sid,pair in enumerate(sorted(pairs))]),
            edges=pd.DataFrame(edge_table).set_index('id'),
            verts=pd.DataFrame({'loc':positions}))
        source_edges = {tuple(map(int,r['generator_tuple'].split(';'))):r for r in edge_rows}
        edge_rows = [dict(source_edges[tuple(row['balls'])], edge_id=eid)
                     for eid,row in network.edges.iterrows()]
    else:
        from vorpy.src.io import load_network
        network = load_network(str(ROOT / 'output/comparative_study/cases/2KAI/H/full_aw.vpy'))
    rep = InterfaceRepresentation(scheme, 'frozen_matched_probe', 0. if scheme=='pow' else 1.4,
                                  group_a, group_b)
    for row in surfaces:
        sid = int(row['surface_id'])
        pair = tuple(sorted((int(row['generator_i']),int(row['generator_j']))))
        if pair not in pairs:
            raise ValueError('Saved surface selection differs from matched-probe selection')
        rep.surfaces.append(InterfaceSurface(sid,pair,None,finite(row['area_A2']),True,True,
            True,True,True,'selected', boundary_edge_ids=tuple(network.surfs.loc[sid,'edges']),
            boundary_vertex_ids=tuple(network.surfs.loc[sid,'verts']),
            smooth_H_A=finite(row['integral_H_dA_A']) if scheme=='pow' or row['support_status']=='included' else None))
    selected_ids = {s.feature_id for s in rep.surfaces}
    for row in edge_rows:
        eid = int(row['edge_id'])
        incidence = tuple(sid for sid in selected_ids if eid in network.surfs.loc[sid,'edges'])
        raw = finite(row['raw_signed_contribution_A_rad'])
        if scheme=='aw' and row['support_status']!='included':
            raw = None
        rep.edges.append(InterfaceEdge(eid,tuple(map(int,row['generator_tuple'].split(';'))),None,
            incidence,incidence,True,True,raw is not None,'saved',
            endpoint_vertex_ids=tuple(network.edges.loc[eid,'verts']),
            classification=row['classification'],length_A=finite(row['length_A']),
            signed_integral_A_rad=raw,unsigned_integral_A_rad=abs(raw) if scheme=='pow' and raw is not None else None,
            curvature_status='integrated' if raw is not None else 'unresolved'))
    # Signed H comes directly from the saved validated analyzer. Saved unsigned
    # integrals are not reconstructed from a signed net integral.
    metrics = next(r for r in rows(MATCHED/'total_interface_metrics.csv')
                   if r['comparison']=='total_restricted_interface' and r['scheme']==prefix)
    rep.curvature = SimpleNamespace(positive_edge_A_rad=float(metrics['positive_edge_A_rad']),
                                   negative_edge_A_rad=float(metrics['negative_edge_A_rad']),
                                   raw_unsigned_edge_A_rad=float(metrics['raw_unsigned_edge_A_rad']))
    return rep, network


def fixture_analysis(scheme):
    rep, network = load_fixture(scheme)
    side = summarize_physical_interface(rep, network)
    side['unsigned_edge'] = rep.curvature.raw_unsigned_edge_A_rad
    analysis = InterfaceGeometryAnalysis(
        metadata={'schema_version':2,'interface_id':'2KAI_matched_probe','scheme':scheme,
                  'alpha_convention':'power_distance' if scheme=='pow' else 'additive',
                  'alpha_value':rep.alpha,'alpha_units':'A^2' if scheme=='pow' else 'A',
                  'solvent_context':False,'additional_network_solves':0,
                  'molecular_contact_patch_status':'not_implemented'},
        alpha_selection={'selected_pairs':len(rep.selected_surfaces),'coverage':1.,'status':'frozen selection'},
        voronoi_side_1=side,voronoi_side_2=reverse_perspective(side),
        dry_solvent_comparison={'available':False,'status':'not applicable'},
        unresolved=side['unresolved_reasons'])
    return analysis


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default=str(ROOT/'output/open_interface_integration/2KAI'))
    args = parser.parse_args()
    out = Path(args.output_dir); out.mkdir(parents=True,exist_ok=True)
    for scheme in ('pow','aw'):
        analysis = fixture_analysis(scheme)
        folder = out/scheme; folder.mkdir(exist_ok=True)
        with (folder/'info.txt').open('w',encoding='utf-8') as stream:
            write_geometry_info(stream,analysis)
        with (folder/'interface_geometry.csv').open('w',encoding='utf-8',newline='') as stream:
            write_geometry_log(csv.writer(stream),analysis)
        (folder/'geometry_analysis.json').write_text(json.dumps(analysis.to_dict(),indent=2,allow_nan=False))
        print(folder)


if __name__ == '__main__':
    main()
