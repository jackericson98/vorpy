from dataclasses import replace
import csv

import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config import Config
from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import (
    classify_atoms, normalize_log, audit_raw_log, topology_record, end_distances,
    annotate_units, read_pdb,
)


def case():
    context = pd.DataFrame({
        'atom_index':[10,30,40,41,42,50], 'atom_name':['O1','C2','OW','HW1','HW2','SOD'],
        'residue_name':['PEG','PEG','SOL','SOL','SOL','SOD'], 'residue_sequence':[1,1,2,2,2,3],
        'residue_id':['p','p','w','w','w','ion'], 'unit_id':['u0','u1','w','w','w','ion'],
        'chain':['A']*6, 'element':['O','C','O','H','H','NA'],
        'x':[0.,1.,2.,2.,2.,4.], 'y':[0.]*6, 'z':[0.]*6,
        'repeat_index':[0,1,np.nan,np.nan,np.nan,np.nan], 'repeat_count':[2]*6,
        'end_distance':[0,0,np.nan,np.nan,np.nan,np.nan], 'order_valid':[True]*6,
        'order_reason':['fixture']*6,
    }).set_index('atom_index',drop=False)
    context = classify_atoms(context,'PEG',Config())
    atoms = context.loc[[30,10]].copy()
    atoms['volume'] = [20.,10.]
    atoms['surface_area'] = [4.,10.]
    atoms['complete'] = True
    atoms['logged_neighbors'] = [[10],[30,40,41,42,50]]
    for col,value in [('polymer','PEG'),('model','fixture'),('frame','frame_0001'),('scheme','aw')]:
        atoms[col] = value
    surfaces = pd.DataFrame({'surface_index':[9,2,8,6,1], 'atom_i':[10]*5,
        'atom_j':[40,30,41,42,50], 'surface_area':[2.,4.,1.,1.,2.], 'contact_area':[.1]*5})
    return context,atoms,surfaces


def test_classification_requires_residue_and_water_name():
    c,_,_ = case()
    assert c.at[10,'kind'] == 'polymer'
    assert c.loc[[40,41,42],'kind'].eq('water').all()
    assert c.at[50,'kind'] == 'ion'
    c.loc[40,'atom_name'] = 'UNKNOWN'
    result = classify_atoms(c,'PEG',Config())
    assert result.at[40,'kind'] == 'other'
    assert not result.at[40,'identity_valid']


def test_missing_topology_remains_unavailable():
    row = topology_record({})
    assert np.isnan(row['Euler Characteristic'])
    assert not row['gauss_bonnet_applicable']
    assert row['internal_occupancy_status'].startswith('unavailable')
    row = topology_record({'Boundary Complete':True,'Boundary Closed':True,'Boundary Manifold':True,'Boundary Orientable':False})
    assert not row['gauss_bonnet_applicable']


def test_end_distances():
    assert end_distances(range(5)).tolist() == [0,1,2,1,0]
    assert end_distances(range(1)).tolist() == [0]


def test_peg_repeat_order_requires_named_connected_backbone():
    c,_,_ = case()
    a = pd.concat([c.loc[[10,30]],c.loc[[30,10]]],ignore_index=True)
    a['atom_index'] = range(4)
    a['atom_name'] = ['O1','C2','C3','O4']
    a['element'] = ['O','C','C','O']
    a['x'] = [0.,1.4,2.8,4.2]
    a['unit_id'] = a.residue_id = 'p'
    result = annotate_units(a,'PEG')
    assert result.order_valid.all()
    assert result.repeat_count.eq(1).all()
    a.loc[2,'x'] = 10.
    assert not annotate_units(a,'PEG').order_valid.any()


def test_reader_preserves_missing_surface_area_and_rejects_bad_numeric(tmp_path):
    path = tmp_path/'log.csv'
    with path.open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['Atoms']); w.writerow(['Index','Name','Residue','Residue Sequence','Chain','X','Y','Z','Surface Area','Volume','Complete Cell?'])
        w.writerow([0,'O1','PEG',1,'A',0,0,0,10,'',True])
        w.writerow(['Surfaces']); w.writerow(['Index','Ball 1','Ball 2','Surface Area'])
        w.writerow([7,0,1,''])
    c,_,_=case()
    c=c.iloc[:2].copy().reset_index(drop=True)
    c['atom_index']=[0,1]; c=c.set_index('atom_index',drop=False)
    a,s,_,_=normalize_log(path,c,Config())
    assert np.isnan(a.at[0,'volume']) and np.isnan(s.surface_area.iloc[0])
    path.write_text(path.read_text().replace('0,O1,PEG','oops,O1,PEG'))
    with pytest.raises(ValueError): audit_raw_log(path)


def test_pdb_molecule_ids_survive_repeated_residue_numbers(tmp_path):
    path=tmp_path/'frame.pdb'
    rows=[]
    for i,(name,res,seq) in enumerate([('O1','PEG',1),('OW','SOL',2),('HW1','SOL',2),('HW2','SOL',2),('OW','SOL',3),('OW','SOL',2)]):
        rows.append(f'ATOM  {i+1:5d} {name:>4s} {res:4s} {seq:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           O\n')
    path.write_text(''.join(rows))
    a=read_pdb(path,'PEG',Config())
    assert a.at[1,'residue_id'] == a.at[3,'residue_id']
    assert a.at[1,'residue_id'] != a.at[5,'residue_id']


def mini_dataset(tmp_path):
    folder=tmp_path/'source/PEG/EG/frame_0001'
    folder.mkdir(parents=True)
    records=[('O1','EG_',1,0.),('C2','EG_',1,1.4),('C3','EG_',1,2.8),('O4','EG_',1,4.2),
             ('OW','SOL',2,6.),('OW','SOL',3,8.)]
    pdb=[]
    for i,(name,res,seq,x) in enumerate(records):
        pdb.append(f'ATOM  {i+1:5d} {name:>4s} {res:4s} {seq:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00          {name[0]:>2s}\n')
    (folder/'frame.pdb').write_text(''.join(pdb))
    for scheme in ['aw','pow','prm']:
        with (folder/f'{scheme}_logs.csv').open('w',newline='') as stream:
            w=csv.writer(stream)
            w.writerow(['Atoms']);w.writerow(['Index','Name','Residue','Residue Sequence','X','Y','Z','Volume','Surface Area','Complete Cell?','Neighbors','Number of Neighbors'])
            for i in range(4):w.writerow([i,records[i][0],'EG_',1,records[i][3],0,0,1,2,True,[[1,4],[0,2],[1,3],[2,5]][i],2])
            w.writerow(['Surfaces']);w.writerow(['Index','Ball 1','Ball 2','Surface Area','Contact Area'])
            for i,(a,b) in enumerate([(0,1),(1,2),(2,3),(0,4),(3,5)]):w.writerow([i,a,b,1,.5])
    return folder


def test_runner_continues_failed_scheme_and_reuses_valid_cache(tmp_path):
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all import run
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import source_signature
    folder=mini_dataset(tmp_path)
    (folder/'pow_logs.csv').write_text('broken')
    config=replace(Config(),dataset_root=tmp_path/'source',output_root=tmp_path/'out',
        polymers=('PEG',),models=('EG',),frames=('frame_0001',))
    result=run(config,[])
    assert result['loaded_logs']==2 and result['failed_or_missing']==1
    before=pd.read_csv(config.output_root/'validation/pdb_residue_inventory.csv')
    rerun=run(config,[])
    after=pd.read_csv(config.output_root/'validation/pdb_residue_inventory.csv')
    pd.testing.assert_frame_equal(before,after)
    assert rerun['loaded_logs']==2
    paths=[folder/'aw_logs.csv',folder/'frame.pdb']
    old=source_signature(paths,config)
    paths[0].write_text(paths[0].read_text()+'\n')
    assert source_signature(paths,config)!=old


@pytest.mark.parametrize('family,model',[('PEG','PEG3000'),('Heparin','dp24')])
def test_representative_real_frame_reader_and_raw_area_validation(family,model):
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import load_scheme, independent_area_check
    config=Config()
    folder=config.dataset_root/family/model/'frame_0001'
    if not folder.is_dir():pytest.skip('External organized dataset is unavailable')
    pdb=read_pdb(folder/'frame.pdb',family,config)
    metadata=dict(polymer=family,model=model,frame='frame_0001',scheme='aw')
    result=load_scheme(folder/'aw_logs.csv',pdb,config,metadata)
    assert result['atoms'].partition_valid.all()
    assert result['atoms'].order_valid.all()
    assert independent_area_check(folder/'aw_logs.csv',result['atoms'],pdb,config,metadata).passed.all()



def test_runner_generates_only_retained_panels_in_custom_output(tmp_path):
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all import FIGURES, run
    mini_dataset(tmp_path)
    config = replace(Config(), dataset_root=tmp_path/'source', output_root=tmp_path/'out',
        polymers=('PEG',), models=('EG',), frames=('frame_0001',),
        export_formats=('png',), dpi=40)
    assert FIGURES == {'A': 'F5A_exposed_vs_buried', 'B': 'F5B_water_burial'}
    result = run(config, ['A', 'B'])
    assert result['figure_failures'] == []
    assert (config.output_root/'F5A/F5A_summary.csv').is_file()
    assert (config.output_root/'F5A/F5A_volume.png').is_file()
    coverage = pd.read_csv(config.output_root/'F5B/F5B_coverage.csv')
    assert coverage.status.eq('water_geometry_unavailable').all()
    assert (config.output_root/'F5B/F5B_summary.csv').is_file()
    assert sorted(p.name for p in config.output_root.glob('F5*')) == ['F5A', 'F5B']
