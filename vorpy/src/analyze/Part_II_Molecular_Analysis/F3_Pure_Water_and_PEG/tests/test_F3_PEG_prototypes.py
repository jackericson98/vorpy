import pandas as pd
import pytest

try:
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F3_PEG_prototypes import (
        chemistry_tables, infer_repeat_units, length_tables,
    )
except ModuleNotFoundError:
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_SimpleModels.F3_PEG_prototypes import (
        chemistry_tables, infer_repeat_units, length_tables,
    )


def test_infer_repeat_units_requires_complete_peg_formula():
    assert infer_repeat_units(10) == 1
    assert infer_repeat_units(45) == 6
    with pytest.raises(ValueError):
        infer_repeat_units(44)


def test_summaries_weight_frames_then_systems_instead_of_atoms():
    observations = pd.DataFrame([
        dict(model='short', frame='f1', repeat_units=1, metric='volume',
             atom_class_aw='Carbon', atom_index=0, delta_pct=10., abs_delta_pct=10.),
        dict(model='short', frame='f1', repeat_units=1, metric='volume',
             atom_class_aw='Carbon', atom_index=1, delta_pct=20., abs_delta_pct=20.),
        dict(model='long', frame='f1', repeat_units=2, metric='volume',
             atom_class_aw='Carbon', atom_index=0, delta_pct=30., abs_delta_pct=30.),
    ])
    frame, model = length_tables(observations)
    assert frame.loc[frame.model.eq('short'), 'mean_absolute_difference_pct'].item() == 15
    assert model.loc[model.model.eq('long'), 'mean_absolute_difference_pct'].item() == 30

    _, _, chemistry = chemistry_tables(observations)
    # Equal system weighting gives (15 + 30) / 2, not the pooled-atom mean of 20.
    assert chemistry.equal_system_mean_pct.item() == 22.5
    assert chemistry.n_systems.item() == 2
