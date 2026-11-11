from types import SimpleNamespace

from vorpy.src.output.visualization import _cached_result_identity


def _interface(source_network_id=None, *, solve_universe_id=None):
    network = SimpleNamespace(
        archive_id="archive:network",
        source_network_id=source_network_id,
        network_scope="system",
        settings={"net_type": "power"},
    )
    if solve_universe_id is not None:
        network.get_solve_context_key = lambda: SimpleNamespace(
            solve_universe_id=solve_universe_id
        )
    system = SimpleNamespace(name="synthetic", archive_id="archive:system")
    return SimpleNamespace(
        archive_id="archive:interface",
        interface_id="A_B",
        name="A_B",
        net=network,
        sys=system,
        group1=SimpleNamespace(name="A"),
        group2=SimpleNamespace(name="B"),
        geometry_analysis=None,
        _geometry_filtration=None,
        _geometry_representation=None,
    )


def test_cached_identity_keeps_producer_and_archive_ids_distinct():
    iface = _interface("producer:network-7")

    identity = _cached_result_identity(iface)

    assert identity.network_id == "archive:network"
    assert identity.source_network_id == "producer:network-7"
    assert identity.source_identity_status == "AUTHORITATIVE"
    assert identity.source_interface_id == "A_B"
    assert identity.network_id != identity.source_network_id


def test_cached_identity_does_not_fabricate_producer_id_for_legacy_archive():
    iface = _interface()

    identity = _cached_result_identity(iface)

    assert identity.network_id == "archive:network"
    assert identity.source_network_id is None
    assert identity.source_identity_status == "UNRESOLVED"


def test_cached_identity_does_not_promote_provisional_namespace():
    iface = _interface(
        "network:context-1",
        solve_universe_id="solve:context-1",
    )

    identity = _cached_result_identity(iface)

    assert identity.source_network_id == "network:context-1"
    assert identity.source_identity_status == "PROVISIONAL"


def test_cached_identity_rejects_archive_id_as_producer_identity():
    iface = _interface("archive:network")

    identity = _cached_result_identity(iface)

    assert identity.source_network_id is None
    assert identity.source_identity_status == "UNRESOLVED"


def test_legacy_molecular_identity_can_remain_explicitly_unresolved():
    iface = _interface()
    surface = SimpleNamespace(
        contact_selection_source="legacy-cache",
        contact_selection_id="selection:1",
        convention="expanded_union_of_balls",
        definition_version="mcs-v1",
    )

    identity = _cached_result_identity(
        iface,
        result_kind="molecular_contact_surface",
        representation="molecular_contact_surface",
        side="A",
        surface=surface,
    )

    assert identity.source_network_id is None
    assert identity.source_identity_status == "UNRESOLVED"
