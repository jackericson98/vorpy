"""Validate explicitly named saved AW/Power runs; stage portable Figure 3B pairs.

Does not solve geometry. Originals and raw CLI outputs are read-only.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

repository = next(parent for parent in Path(__file__).resolve().parents
                  if (parent / 'pyproject.toml').is_file())
sys.path.insert(0, str(repository))
from vorpy.src.io.network_archive import load_network
from vorpy.src.network.analyze import diagnose_geometry_completeness


def check_network(network, scheme):
    central, = network.group
    report = diagnose_geometry_completeness(network)['cells'][0]
    assert report['complete'], report
    if scheme == 'aw':
        assert report['bounded'], report
    assert not any(report[k] for k in ('missing_faces','missing_edges','missing_vertices')), report
    assert not any(report['closure_relations'].values()), report
    locations = np.asarray(list(network.balls['loc']))
    radii = np.asarray(network.balls['rad'])
    actual = network.balls.loc[central]
    assert np.isfinite(actual['vol']) and actual['vol'] > 0
    assert np.isfinite(actual['sa']) and actual['sa'] > 0
    maximum_residual = maximum_intrusion = 0.
    for _, face in network.surfs.iterrows():
        left, right = face['balls']
        points = np.asarray(face['points'])
        assert len(points) and len(face['tris']) and face['sa'] > 0
        if scheme == 'aw':
            values = np.linalg.norm(points-locations[left],axis=1)-radii[left]
            other = np.linalg.norm(points-locations[right],axis=1)-radii[right]
            for start in range(0,len(points),128):
                chunk = points[start:start+128]
                distances = np.linalg.norm(chunk[:,None,:]-locations[None,:,:],axis=2)-radii
                maximum_intrusion = max(maximum_intrusion, float(np.max(values[start:start+128]-distances.min(axis=1))))
        else:
            values = np.sum((points-locations[left])**2,axis=1)-radii[left]**2
            other = np.sum((points-locations[right])**2,axis=1)-radii[right]**2
        maximum_residual = max(maximum_residual, float(np.max(np.abs(values-other))))
    # AW triangulation inputs are normalized at 10 decimal places in the existing mesher.
    assert maximum_residual < 1e-7, maximum_residual
    assert maximum_intrusion < 1e-7, maximum_intrusion
    result = dict(system_index=int(actual['system_num']), local_id=int(central),
                  generators=len(network.balls), vertices=len(network.verts), edges=len(network.edges),
                  faces=len(network.surfs), volume=float(actual['vol']), area=float(actual['sa']),
                  maximum_bisector_residual=maximum_residual, maximum_aw_competitor_intrusion=maximum_intrusion,
                  completeness=report)
    if scheme == 'pow':
        delta = locations-locations[central]
        offsets = -np.sum(delta*delta,axis=1)-radii[central]**2+radii*radii
        keep = np.arange(len(locations)) != central
        assert np.all(offsets[keep] < 0), 'Reference needs a strict interior point'
        halfspaces = np.column_stack((2*delta[keep],offsets[keep]))
        reference_vertices = HalfspaceIntersection(halfspaces,np.zeros(3)).intersections
        reference = ConvexHull(reference_vertices)
        np.testing.assert_allclose([actual['vol'],actual['sa']], [reference.volume,reference.area],rtol=1e-11)
        actual_vertices = np.asarray(list(network.verts['loc']))-locations[central]
        np.testing.assert_allclose(np.min(np.linalg.norm(reference_vertices[:,None,:]-actual_vertices[None,:,:],axis=2),axis=1),0,atol=1e-10)
        assert len(reference_vertices) == len(actual_vertices)
        active = np.sum(np.abs(reference_vertices @ halfspaces[:,:3].T + halfspaces[:,3]) < 1e-9,axis=0) >= 3
        expected_neighbors = set(np.flatnonzero(keep)[active])
        actual_neighbors = {int(i) for pair in network.surfs['balls'] for i in pair if i != central}
        assert actual_neighbors == expected_neighbors
        result['independent_halfspace_reference'] = dict(vertices=len(reference_vertices),volume=reference.volume,area=reference.area)
        result['independent_halfspace_reference']['face_neighbors'] = sorted(map(int,expected_neighbors))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples',type=Path,required=True,help='Original Figure3B_AtomExamples, for the generator only')
    parser.add_argument('--runs',type=Path,default=Path('output/Figure3B_revalidation'))
    parser.add_argument('--example',action='append',required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('figure3b',args.examples/'generate_figure3b_pymol.py')
    figure = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(figure)
    summary = {}
    for name in args.example:
        index = name.rsplit('_',1)[1]
        staged = args.runs/'validated'/name
        summary[name] = {}
        for scheme in ('aw','pow'):
            archives = list((args.runs/name/scheme).glob('*/*.vpy'))
            if len(archives) != 1:
                raise ValueError('Expected one CLI archive: '+str(args.runs/name/scheme))
            archive = archives[0]
            network = load_network(archive)
            summary[name][scheme] = check_network(network, scheme)
            destination = staged/scheme
            folder = destination/('a_'+index)
            folder.mkdir(parents=True,exist_ok=True)
            shutil.copy2(archive.with_suffix('.pdb'),destination/archive.with_suffix('.pdb').name)
            raw = archive.parent/('a_'+index)
            for source,target in [('surfs.off',scheme+'_surfs.off'),('edges.off',scheme+'_edges.off'),
                                  ('surr_atoms.pdb',scheme+'_surr_atoms.pdb'),('group_atoms.pdb','group_atoms.pdb'),
                                  ('surr_resids.pdb','surr_resids.pdb'),(scheme+'_logs.csv',scheme+'_logs.csv'),
                                  (scheme+'_atom_identity.csv',scheme+'_atom_identity.csv')]:
                shutil.copy2(raw/source,folder/target)
            row, data = figure.audit(staged,scheme)
            assert row['recommended_action'] == 'READY', row
            atomic, = figure.read_logs(folder/(scheme+'_logs.csv'))['Atoms']
            assert int(atomic['Index']) == summary[name][scheme]['local_id']
            # CSV metrics are serialized to six decimal places by the existing writer.
            for field, value in [('Volume',summary[name][scheme]['volume']),('Surface Area',summary[name][scheme]['area'])]:
                assert abs(float(atomic[field])-value) <= 0.500001e-6
            summary[name][scheme]['export_audit'] = row
        status = figure.main(['--root',str(args.runs/'validated'),'--example',name])
        assert status == 0
    path = args.runs/'real_pair_validation.json'
    path.write_text(json.dumps(summary,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
    print('Independent geometry and export checks passed:',path)


if __name__ == '__main__':
    main()
