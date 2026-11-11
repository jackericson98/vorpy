from dataclasses import replace
import numpy as np
import pandas as pd
import pytest

from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import reconstruct_interfaces, aggregate_residues
from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config import Config
from .test_F5_loading import case


def test_index_matching_area_partition_and_unique_water_count():
    c,a,s=case()
    out,interfaces,_=reconstruct_interfaces(a,s,c,Config())
    o=out.set_index('atom_index').loc[10]
    assert o.water_area == 4. and o.water_fraction == .4
    assert o.polymer_area == 4. and o.other_area == 2.
    assert o.water_coordination == 1 and o.area_per_water == 4.
    assert o.same_residue_contacts == 1 and o.adjacent_repeat_contacts == 1 and o.nonlocal_contacts == 0
    assert o.primary_valid and o.neighbor_list_consistent
    assert len(interfaces)==5 and interfaces.surface_area.sum()==10.
    assert out.interface_area_sum.sum()==14.  # The polymer-polymer face contributes to both atoms.


def test_residue_fraction_uses_area_weights_and_deduplicates_molecules():
    c,a,s=case(); out,_,_=reconstruct_interfaces(a,s,c,Config())
    r=aggregate_residues(out).iloc[0]
    assert r.water_fraction == pytest.approx(4/14)
    assert r.volume == 30. and r.mean_atomic_volume == 15.
    assert r.water_coordination==1 and r.polymer_coordination==0


def test_zero_water_contacts_are_valid():
    c,a,s=case(); out,_,_=reconstruct_interfaces(a,s,c,Config())
    r=out.set_index('atom_index').loc[30]
    assert r.water_fraction==0. and r.water_coordination==0
    assert r.primary_valid and np.isnan(r.area_per_water)


def test_duplicate_records_count_once_but_conflicts_fail():
    c,a,s=case(); dup=pd.concat([s,s.iloc[[0]]],ignore_index=True)
    out,interfaces,count=reconstruct_interfaces(a,dup,c,Config())
    assert count==1 and len(interfaces)==5
    assert out.set_index('atom_index').at[10,'water_area']==4.
    dup.loc[len(dup)-1,'surface_area']=30.
    with pytest.raises(ValueError,match='Conflicting'): reconstruct_interfaces(a,dup,c,Config())


@pytest.mark.parametrize('defect',['incomplete','missing_area','bad_partition','missing_face'])
def test_invalid_cells_retained_and_flagged(defect):
    c,a,s=case()
    if defect=='incomplete': a.loc[10,'complete']=False
    if defect=='missing_area': s.loc[0,'surface_area']=np.nan
    if defect=='bad_partition': a.loc[10,'surface_area']=1.
    if defect=='missing_face': s=s.iloc[1:]
    out,_,_=reconstruct_interfaces(a,s,c,Config())
    r=out.set_index('atom_index').loc[10]
    assert len(out)==2 and not r.primary_valid
    if defect=='bad_partition': assert r.water_fraction==4.  # Never clip.
    if defect=='incomplete':
        out,_,_=reconstruct_interfaces(a,s,c,replace(Config(),include_incomplete=True))
        assert out.set_index('atom_index').at[10,'primary_valid']


def test_hydrogen_filter_keeps_water_hydrogen_interface_area():
    c,a,s=case(); out,_,_=reconstruct_interfaces(a,s,c,replace(Config(),include_hydrogens=False))
    assert out.set_index('atom_index').at[10,'water_area']==4.


def test_nonlocal_contacts_follow_repeat_identity():
    c,a,s=case()
    c.at[30,'repeat_index']=4
    a.at[30,'repeat_index']=4
    out,_,_=reconstruct_interfaces(a,s,c,Config())
    row=out.set_index('atom_index').loc[10]
    assert row.adjacent_repeat_contacts==0 and row.nonlocal_contacts==1
    assert row.nonlocal_unit_neighbors==frozenset({'u1'})
