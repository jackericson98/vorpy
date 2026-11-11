import pandas as pd
import pytest

from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import reconstruct_interfaces, match_schemes
from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config import Config
from .test_F5_loading import case


def paired():
    c,a,s=case(); a,_,_=reconstruct_interfaces(a,s,c,Config())
    p=a.copy(); p['scheme']='pow'; p['water_fraction']+=.1; p['volume']+=2
    return a,p


def test_paired_differences_percentage_points_and_zero_counts():
    a,p=paired(); out,report=match_schemes(pd.concat([a,p]),Config())
    assert out.pair_valid.all() and out.delta_water_fraction.tolist()==pytest.approx([10.,10.])
    assert out.delta_volume.eq(2).all() and out.delta_water_coordination.eq(0).all()
    assert out.neighbor_jaccard.eq(1).all() and report.n_valid.sum()==2


@pytest.mark.parametrize('column,value',[('x',100.),('atom_name','wrong'),('residue_name','wrong')])
def test_scheme_identity_mismatch_cannot_generate_differences(column,value):
    a,p=paired(); p.loc[p.index[0],column]=value
    out,report=match_schemes(pd.concat([a,p]),Config())
    assert out.delta_volume.notna().sum()==1
    assert report.n_identity_mismatch.sum()==1


def test_unmatched_atoms_reported_and_empty_neighbor_sets_have_jaccard_one():
    a,p=paired()
    a['neighbors']=[frozenset()]*2; p['neighbors']=[frozenset()]*2
    out,report=match_schemes(pd.concat([a,p.iloc[:1]]),Config())
    assert out.pair_valid.sum()==1 and report.n_aw.sum()==2 and report.n_other.sum()==1
    assert out.loc[out.pair_valid,'neighbor_jaccard'].iloc[0]==1.


def test_repeat_matches_verify_member_coordinates():
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import aggregate_residues
    a,p=paired()
    p.loc[p.index[0],'x']+=1
    units=pd.concat([aggregate_residues(a,unit=True),aggregate_residues(p,unit=True)])
    out,report=match_schemes(units,Config(),level='unit')
    assert out.pair_valid.sum()==1 and report.n_identity_mismatch.sum()==1


def test_original_residue_pairs_retain_residue_key_and_summed_volumes():
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import aggregate_residues
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_run_all import add_residue_pairs
    a,p=paired()
    data={'residues':pd.concat([aggregate_residues(a),aggregate_residues(p)])}
    report=add_residue_pairs(data,Config())
    pairs=data['residue_pairs']
    assert pairs.residue_id.tolist()==['p']
    assert pairs.delta_volume.iloc[0]==4.
    assert pairs.level.eq('residue').all() and report.n_valid.sum()==1
