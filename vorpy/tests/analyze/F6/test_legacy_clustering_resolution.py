import importlib


PACKAGE = "vorpy.src.analyze.Part_II_Molecular_Analysis.F4_atomic_deviations"


def test_split_launchers_use_complete_monolithic_orchestrator():
    monolith = importlib.import_module(
        f"{PACKAGE}.F3_old.F2C_All_Data_Sets_Clustered"
    )
    run_single = importlib.import_module(
        f"{PACKAGE}.F2C_All_Data_Sets_Clustered_full.run_single"
    )
    run_batch = importlib.import_module(
        f"{PACKAGE}.F2C_All_Data_Sets_Clustered_full.run_batch"
    )

    assert run_single.main is monolith.main
    assert run_batch.main is monolith.main


def test_split_refactor_remains_scientifically_distinct():
    monolith = importlib.import_module(
        f"{PACKAGE}.F3_old.F2C_All_Data_Sets_Clustered"
    )
    split_config = importlib.import_module(
        f"{PACKAGE}.F2C_All_Data_Sets_Clustered_full.config"
    )
    split_pipeline = importlib.import_module(
        f"{PACKAGE}.F2C_All_Data_Sets_Clustered_full.pipeline"
    )

    assert not hasattr(split_pipeline, "main")
    assert split_config.SETTINGS_COMBINATIONS == monolith.SETTINGS_COMBINATIONS
    assert split_config.SMALL_MOL_ATOM_ALIASES != monolith.SMALL_MOL_ATOM_ALIASES
    assert split_config.RNA_ATOM_ALIASES != monolith.RNA_ATOM_ALIASES
    assert split_config.FORCE_INCLUDE_GROUPS != monolith.FORCE_INCLUDE_GROUPS
