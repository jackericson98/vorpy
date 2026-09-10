from types import SimpleNamespace

from vorpy.src.command.vpy_cmnd import Command


class DummySystem:
    def __init__(self, groups=None, interfaces=None):
        self.groups = [] if groups is None else groups
        self.interfaces = [] if interfaces is None else interfaces


def test_diagnose_edges_flag_is_removed_before_grouped_parsing(monkeypatch):
    command = Command(sys=DummySystem())
    monkeypatch.setattr(
        "sys.argv",
        [
            "vorpy", "edta", "-s", "mv", "5",
            "--diagnose-edges", "-e", "small",
        ],
    )

    command.parse_commands()

    assert command.diagnose_edges is True
    assert command.settings_cmnds == [["mv", "5"]]
    assert command.exports == [["small"]]
    assert command.load_commands == []


def test_edge_diagnostic_prints_summary_and_problem_edges(capsys):
    records = (
        SimpleNamespace(
            edge_index=4,
            status="matched_curve",
            ball_indices=(0, 1, 2),
            vertex_indices=(3, 8),
            reason="matched",
        ),
        SimpleNamespace(
            edge_index=9,
            status="turning_point",
            ball_indices=(2, 3, 4),
            vertex_indices=(10, 11),
            reason="rho interval reaches a turning point",
        ),
    )
    report = SimpleNamespace(
        records=records,
        summary=lambda: {
            "total_edges": 2,
            "matched_edges": 1,
            "matched_fraction": 0.5,
            "status_counts": {
                "matched_curve": 1,
                "matched_nonsingular": 0,
                "turning_point": 1,
            },
            "maximum_sample_error": 1.25e-5,
            "rms_sample_error": 4.0e-6,
            "maximum_absolute_length_difference": 8.0e-5,
        },
    )
    net = SimpleNamespace(
        settings={"net_type": "aw"},
        diagnose_aw_edges=lambda: report,
    )
    command = Command(
        sys=DummySystem(groups=[SimpleNamespace(name="EDTA", net=net)])
    )

    command.run_edge_diagnostics()
    output = capsys.readouterr().out

    assert "ANALYTIC AW EDGE DIAGNOSTIC: EDTA" in output
    assert "Matched analytic curves:      1" in output
    assert "Matched turning conics:       0" in output
    assert "Match coverage:               50.00 %" in output
    assert "edge 9: turning_point" in output


def test_edge_diagnostic_skips_non_aw_network(capsys):
    net = SimpleNamespace(settings={"net_type": "pow"})
    command = Command(
        sys=DummySystem(groups=[SimpleNamespace(name="Power", net=net)])
    )

    command.run_edge_diagnostics()
    output = capsys.readouterr().out

    assert "Skipping analytic edge diagnostic for Power" in output
    assert "network type 'pow' is not AW" in output


def test_diagnostic_networks_deduplicates_shared_network_objects():
    net = SimpleNamespace(settings={"net_type": "aw"})
    system = DummySystem(
        groups=[SimpleNamespace(name="Group", net=net)],
        interfaces=[SimpleNamespace(name="Interface", net=net)],
    )
    command = Command(sys=system)

    assert list(command._diagnostic_networks()) == [("Group", net)]
