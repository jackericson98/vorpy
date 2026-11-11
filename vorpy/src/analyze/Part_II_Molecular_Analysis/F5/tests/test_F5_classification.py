import numpy as np
import pandas as pd

from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_common import exposure, attach_aw_exposure, reconstruct_interfaces, aggregate_residues
from vorpy.src.analyze.Part_II_Molecular_Analysis.F5.F5_config import Config
from .test_F5_loading import case


def test_provisional_exposure_boundaries_and_invalid_values():
    result=exposure(pd.Series([0,.1,.1001,.4,.4001,1.,-1.,1.1,np.nan]),Config())
    assert result.iloc[:6].tolist()==['Buried','Buried','Intermediate','Intermediate','Exposed','Exposed']
    assert result.iloc[6:].isna().all()


def test_environment_is_aw_repeat_reference_for_other_schemes():
    c,a,s=case(); a,_,_=reconstruct_interfaces(a,s,c,Config())
    units=aggregate_residues(a,unit=True)
    other=units.copy(); other['scheme']='pow'; other['water_fraction']=.9
    allunits=pd.concat([units,other],ignore_index=True)
    _,out=attach_aw_exposure(a,allunits,Config())
    assert out.loc[out.unit_id.eq('u0'),'environment'].eq('Intermediate').all()
    assert out.loc[out.unit_id.eq('u1'),'environment'].eq('Buried').all()
