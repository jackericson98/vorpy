from pathlib import Path

from vorpy.src.boundary import BoundaryConfig, SolventMode, prepare_network_geometry
from vorpy.src.group import Group
from vorpy.src.io.network_archive import load_network, save_network
from vorpy.src.system.system import System
from vorpy.tests.network.test_views import _network


ROOT = Path(__file__).resolve().parents[2]
EDTA = ROOT / "data" / "EDTA.pdb"


def _edta_group(retained=(), boundary_config=None, max_vert=40):
    system = System(
        file=str(EDTA),
        make_dir=False,
        retained_system_indices=retained,
        boundary_config=boundary_config,
    )
    group = Group(
        system,
        name="edta",
        atoms=list(range(32)),
        max_vert=max_vert,
        make_net=True,
        build_net=False,
    )
    return system, group


def test_edta_mg_can_be_retained_without_retaining_water():
    system, group = _edta_group((32,))
    geometry = prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)

    assert len(system.balls) == 612
    assert geometry.system_indices == tuple(range(33))
    assert geometry.retained_system_indices == (32,)
    assert geometry.retained_generator_ids == ("system:32",)
    assert geometry.solvent_local_indices == ()
    assert len(group.net.balls) == 33
    assert group.net.group == list(range(32))
    assert group.net.retained_system_indices == (32,)
    assert group.net.retained_generator_ids == ("system:32",)
    assert group.net.retained_generator_provenance["boundary_independent"] is True


def test_retention_changes_solve_identity_and_keeps_default_exclusion():
    _, default_group = _edta_group()
    _, retained_group = _edta_group((32,))

    default_context = default_group.net.get_solve_context_key(frame_id="static")
    retained_context = retained_group.net.get_solve_context_key(frame_id="static")

    assert default_context.physical_generator_ids[-1] == "system:31"
    assert retained_context.physical_generator_ids[-1] == "system:32"
    assert default_context.digest != retained_context.digest
    assert retained_group.net.solvent_generator_ids == ()
    assert retained_group.net.retained_generator_digest


def test_retained_edta_mg_shell_network_marks_domain_boundary_separately():
    _, group = _edta_group(
        (32,),
        boundary_config=BoundaryConfig("shell", shell_spacing=4.0),
        max_vert=5,
    )
    group.build(calculate_curvature=False)
    network = group.net
    physical = network.balls.iloc[:33]
    boundary_indices = set(network.boundary_indices)

    assert len(physical) == 33
    assert boundary_indices
    assert network.solvent_local_indices == ()
    assert network.retained_generator_ids == ("system:32",)
    # The EDTA group cells are bounded and complete; Mg is a retained
    # competitor, not an implicitly selected/analysed cell.
    assert physical.iloc[:32]["complete"].astype(bool).all()
    assert bool(physical.iloc[32]["complete"]) is False
    assert all(
        entry[4] is True
        for entry in network.generator_membership_evidence
        if entry[0] in boundary_indices
    )
    assert not any(
        entry[3] is True
        for entry in network.generator_membership_evidence
        if entry[0] == 32
    )


def test_retention_metadata_survives_archive_round_trip(tmp_path):
    network = _network(scheme="aw")
    network.group = [0, 1, 2, 3]
    network.sys = System(
        simple=True,
        make_dir=False,
        retained_system_indices=(32,),
    )
    network.sys.balls = network.balls.copy(deep=True)
    network.sys.balls["res"] = None
    network.sys.balls["chn"] = None
    network.solvent_mode = SolventMode.MOLECULE_ONLY.value
    network.retained_system_indices = (32,)
    network.retained_generator_ids = ("system:32",)
    network.refresh_solvent_provenance()
    before_context = network.get_solve_context_key(frame_id="static")
    before_retention = network.retained_generator_provenance

    archive = save_network(network, tmp_path / "edta-mg-retained.vpy")
    loaded = load_network(archive)

    assert loaded.retained_system_indices == (32,)
    assert loaded.retained_generator_ids == ("system:32",)
    assert loaded.retained_generator_provenance == before_retention
    assert loaded.get_solve_context_key(frame_id="static") == before_context
