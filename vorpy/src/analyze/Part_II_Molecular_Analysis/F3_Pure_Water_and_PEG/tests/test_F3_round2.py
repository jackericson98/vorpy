import pandas as pd
import pytest

try:
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F3_round2_validation import descriptive_trend
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_Pure_Water_and_PEG.F3F_water_contact_prototype import paired_differences
except ModuleNotFoundError:
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_SimpleModels.F3_round2_validation import descriptive_trend
    from vorpy.src.analyze.Part_II_Molecular_Analysis.F3_SimpleModels.F3F_water_contact_prototype import paired_differences


def test_descriptive_trend_reports_slope_without_replicate_inference():
    slope, rho, count = descriptive_trend([0, 1, 2], [1, 3, 5])
    assert slope == pytest.approx(2)
    assert rho == 1
    assert count == 3


def test_water_differences_use_raw_faces_percentage_points_and_aw_denominator():
    common = dict(model='EG', frame='f1', water_id=0)
    observations = pd.DataFrame([
        dict(common, scheme='aw', water_contacts=10, polymer_contacts=4,
             water_sa_percent=80., external_area=50.),
        dict(common, scheme='pow', water_contacts=12, polymer_contacts=3,
             water_sa_percent=75., external_area=55.),
    ])
    pairs, frame, model = paired_differences(observations)
    assert pairs.delta_water_contacts_faces.item() == 2
    assert pairs.delta_polymer_contacts_faces.item() == -1
    assert pairs.delta_water_sa_percentage_points.item() == -5
    assert pairs.delta_external_area_pct.item() == 10
    assert frame.n_matched_waters.item() == 1
    assert model.n_water_frame_observations.item() == 1
