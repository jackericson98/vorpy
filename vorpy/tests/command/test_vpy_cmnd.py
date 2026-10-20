import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.boundary import SolventMode
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
from vorpy.src.system.system import System
from vorpy.workbench.services.vorpy_backend import (
    InterfaceSolveSelection,
    VorPyBackend,
    VorPySolveSettings,
)


@pytest.mark.parametrize('all_frames', [False, True])
@pytest.mark.parametrize('custom_directory', [False, True])
def test_cli_creates_only_selected_output_directory(tmp_path, monkeypatch, all_frames, custom_directory):
    source = tmp_path / 'sample.pdb'
    atom = 'ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n'
    source.write_text(atom + 'ENDMDL\n' + (atom + 'ENDMDL\n' if all_frames else ''))
    original_set_files = System.set_files

    def set_files(system, *args, **kwargs):
        original_set_files(system, *args, **kwargs)
        system.files['vpy_dir'] = str(tmp_path)

    monkeypatch.setattr(System, 'set_files', set_files)
    args = ['vorpy', str(source), '-e', 'only', 'pdb']
    if all_frames:
        args.append('--all-frames')
    if custom_directory:
        args.extend(['-e', 'dir', str(tmp_path / 'custom results')])
    monkeypatch.setattr('sys.argv', args)
    monkeypatch.setattr(Command, 'create_groups', lambda self: None)

    Command().run()

    root = tmp_path / ('custom results' if custom_directory else 'output') / 'sample'
    directories = [root / 'frames' / f'frame_{index:04d}' for index in (1, 2)] if all_frames else [root]
    for directory in directories:
        assert (directory / 'sample.pdb').is_file()
    assert (tmp_path / 'output').exists() is (not custom_directory)
    assert list(root.parent.iterdir()) == [root]


class DummySystem:
    def __init__(self, groups=None, interfaces=None):
        self.groups = [] if groups is None else groups
        self.interfaces = [] if interfaces is None else interfaces


def test_legacy_cli_interface_result_state_is_explicitly_unresolved(tmp_path, capsys):
    class InterfaceSystem(DummySystem):
        def __init__(self):
            super().__init__(
                groups=[SimpleNamespace(name='A'), SimpleNamespace(name='B')]
            )
            self.files = {'dir': str(tmp_path)}

        def make_interfaces(self, pairs, *, analyze_waters):
            self.interface_pairs = pairs
            self.analyze_waters = analyze_waters

    system = InterfaceSystem()
    command = Command(sys=system)
    command.base_file = tmp_path / 'synthetic.input'
    command.interface_mode = True
    command.run_exports = lambda: None
    command._save_archive = lambda: None
    command.load_files = lambda: None
    command.apply_settings = lambda: None
    command.create_groups = lambda: None
    command._build_full_system_target = lambda: None

    command._process_system()

    assert system.interface_pairs == [(system.groups[0], system.groups[1])]
    assert command.interface_result_state == {
        'status': 'UNRESOLVED',
        'typed_result': False,
        'provenance': False,
        'dependency': False,
        'reason': command.interface_result_state['reason'],
    }
    assert 'typed Workbench InterfaceResult' in command.interface_result_state['reason']
    assert 'Warning: interface Results unresolved; see interface_info.txt' in capsys.readouterr().out


def test_cli_direct_ordered_interface_handoff_keeps_partial_geometry(
    tmp_path, capsys, monkeypatch
):
    system = System(
        spheres=pd.DataFrame({
            'loc': [(0.0, 0.0, 0.0), (2.0, 0.2, 0.1), (0.2, 2.1, 0.3), (0.1, 0.4, 2.2)],
            'rad': [1.0, 1.1, 0.9, 1.05],
            'mass': [12.0, 16.0, 18.0, 14.0],
            'num': [100, 101, 102, 103],
            'res_name': ['ALA', 'GLY', 'SER', 'THR'],
            'name': ['CA', 'O', 'CA', 'CA'],
        }),
        name='cli-interface-handoff',
        make_dir=False,
        solvent_mode=SolventMode.MOLECULE_ONLY,
        retained_system_indices=(),
    )
    selection = InterfaceSolveSelection(
        interface_id='chains-A-B__chain-I',
        group_a_id='chains-A-B',
        group_a_system_indices=(0, 1),
        group_b_id='chain-I',
        group_b_system_indices=(2, 3),
    )
    interface = VorPyBackend(
        VorPySolveSettings(network_type='pow')
    ).make_interface_network(system, selection)
    interface.build(analyze_waters=False, calculate_curvature=False)
    from copy import copy as shallow_copy
    incomplete_source = shallow_copy(interface.net)
    incomplete_source.network_mode = 'complete'
    incomplete_source.group = list(range(len(incomplete_source.balls)))
    incomplete_source.source_network_id = 'producer:cli-incomplete-view'
    interface.net.source_network_id = incomplete_source.source_network_id
    system.ifaces.clear()
    interface.geometry_analysis = InterfaceGeometryAnalysis(
        metadata={'result_status': 'PARTIAL'},
        voronoi_side_1={},
    )
    system.files = {'dir': str(tmp_path)}

    def make_interfaces(pairs, *, analyze_waters):
        assert pairs == [(system.groups[0], system.groups[1])]
        system.ifaces.append(interface)

    system.make_interfaces = make_interfaces
    saved = {}

    def save_network(network, path, **kwargs):
        saved.update(network=network, path=path, kwargs=kwargs)

    monkeypatch.setattr('vorpy.src.io.save_network', save_network)
    command = Command(sys=system)
    command.base_file = tmp_path / 'synthetic_cli.pdb'
    command.interface_mode = True
    command.create_groups = lambda: None
    command.run_exports = lambda: None
    command.save_network_path = tmp_path / 'live-interface.vpy'
    command._full_system_group = SimpleNamespace(net=incomplete_source)

    command._process_system()

    assert len(command.interface_results) == 1
    assert command.interface_result_state['status'] == 'PARTIAL'
    assert command.interface_result_state['typed_result'] is True
    assert command.interface_results[0].provenance.geometry_lineage.origin.value == (
        'DIRECT_SOLVE'
    )
    assert system.cli_interface_results == command.interface_results
    assert saved['network'] is system
    assert saved['network'].cli_interface_results == command.interface_results
    assert saved['path'] == command.save_network_path
    output = capsys.readouterr().out
    assert 'Warning: interface extraction failed:' not in output
    assert 'CERTIFIED' not in output


def test_actual_cli_ordered_pow_interface_archive_reload_and_cache_export(tmp_path):
    """Exercise the validated live producer through two fresh CLI processes."""
    source = tmp_path / 'synthetic_ab_pow.pdb'
    atoms = (
        ('A', 1, 0.0, 0.0, 0.0), ('A', 2, 1.6, 0.0, 0.0),
        ('A', 3, 0.0, 1.6, 0.0), ('A', 4, 0.0, 0.0, 1.6),
        ('B', 5, 2.0, 1.0, 0.0), ('B', 6, 3.6, 1.0, 0.0),
        ('B', 7, 2.0, 2.6, 0.0), ('B', 8, 2.0, 1.0, 1.6),
        ('I', 9, 1.0, 0.5, 0.8), ('I', 10, 2.6, 0.5, 0.8),
        ('I', 11, 1.0, 2.1, 0.8), ('I', 12, 1.0, 0.5, 2.4),
    )
    source.write_text(
        ''.join(
            f'ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}'
            f'    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n'
            for serial, (chain, residue, x, y, z) in enumerate(atoms, 1)
        ) + 'END\n',
        encoding='ascii',
    )
    output = tmp_path / 'solve-output'
    archive = tmp_path / 'synthetic-ab-pow.vpy'
    timing_report = tmp_path / 'synthetic-ab-pow.timing.json'
    command = [
        sys.executable, '-m', 'vorpy', '--fail-if-exists', str(source),
        '-g', 'c', 'A', 'and', 'c', 'B', '-g', 'c', 'I', '-i',
        '-c', 's',
        '-s', 'nt', 'pow', 'and', 'sr', '0.2', 'and', 'mv', '30',
        '--boundary', 'shell', '--solvent-mode', 'MOLECULE_ONLY',
        '--profile', 'full', '-e', 'only', 'logs', '-e', 'dir',
        str(output), '--save-network', str(archive),
        '--timing-report', str(timing_report),
    ]
    root = Path(__file__).resolve().parents[3]
    solved = subprocess.run(
        command, cwd=root, capture_output=True, text=True, timeout=120,
        env={**os.environ, 'PYTHONPATH': str(root)},
    )
    assert solved.returncode == 0, solved.stdout + '\n' + solved.stderr
    assert solved.stdout.count('network built -') == 1, solved.stdout
    assert '"c_A_and_c_B_c_I_interface" built - ' in solved.stdout
    assert 'Interface Result status:' not in solved.stdout
    assert 'Extracting ' not in solved.stdout
    assert 'CERTIFIED' not in solved.stdout
    assert 'VORPY_BUILD_TIMING' not in solved.stdout
    assert 'EXPORT TIMING' not in solved.stdout
    assert archive.is_file()
    timing = json.loads(timing_report.read_text(encoding='utf-8'))
    assert timing['format'] == 'vorpy.cli_timing'
    assert set(timing['phases']) >= {
        'loading_initialization', 'geometry_network', 'analysis',
        'interface_extraction', 'export', 'archive_save', 'total_wall',
    }
    assert timing['phases']['geometry_network']['source'].startswith('FORGE')
    assert timing['phases']['export']['source'].startswith('VECTOR')
    assert timing['phases']['archive_save']['source'].startswith('ATLAS')

    with zipfile.ZipFile(archive) as bundle:
        state = json.loads(bundle.read('state.json'))
    interface_record = next(
        record for record in state['objects'] if record['kind'] == 'interface'
    )
    assert 'typed_results' in interface_record['fields']
    assert 'results_dependency' in interface_record['fields']
    typed_snapshot = interface_record['fields']['typed_results']['tuple'][0]
    snapshot_map = dict(typed_snapshot['map'])
    assert snapshot_map['format'] == 'vorpy.interface_result_snapshot'

    reload_code = (
        'import json, sys; '
        'from vorpy.src.io import load_network; '
        'network = load_network(sys.argv[1]); '
        'iface = network.sys.ifaces[0]; result = iface.typed_results[0]; '
        'print(json.dumps({"group_A": result.identity.group_A, '
        '"group_B": result.identity.group_B, '
        '"orientation": result.identity.orientation, '
        '"partition": result.identity.partition, '
        '"source_network_id": result.identity.source_network_id, '
        '"archive_network_id": result.identity.network_id, '
        '"status": result.status.value, '
        '"lineage_origin": result.provenance.geometry_lineage.origin.value, '
        '"lineage_source": result.provenance.geometry_lineage.source_network_id, '
        '"dependency_status": iface.results_dependency.status.value}))'
    )
    reloaded = subprocess.run(
        [sys.executable, '-c', reload_code, str(archive)],
        cwd=root, capture_output=True, text=True, timeout=30,
        env={**os.environ, 'PYTHONPATH': str(root)},
    )
    assert reloaded.returncode == 0, reloaded.stdout + '\n' + reloaded.stderr
    restored = json.loads(reloaded.stdout)
    assert restored['group_A'] == 'c_A_and_c_B'
    assert restored['group_B'] == 'c_I'
    assert restored['orientation'] == 'group1 -> group2'
    assert restored['partition'] == 'power'
    assert restored['source_network_id']
    assert restored['archive_network_id'] != restored['source_network_id']
    assert restored['status'] == 'PARTIAL'
    assert restored['lineage_origin'] == 'MATERIALIZED_VIEW'
    assert restored['lineage_source'] == restored['source_network_id']
    assert restored['dependency_status'] == 'CERTIFIED'

    cache_output = tmp_path / 'cache-output'
    cached = subprocess.run(
        [
            sys.executable, '-m', 'vorpy', '--fail-if-exists',
            '--profile-cache-export', str(archive), '-e', 'dir',
            str(cache_output), '-e', 'only', 'logs',
        ],
        cwd=root, capture_output=True, text=True, timeout=60,
        env={**os.environ, 'PYTHONPATH': str(root)},
    )
    assert cached.returncode == 0, cached.stdout + '\n' + cached.stderr
    profiles = list(cache_output.rglob('cache_export_profile.json'))
    assert len(profiles) == 1
    profile = json.loads(profiles[0].read_text(encoding='utf-8'))
    assert profile['cache_state']['interface_count'] == 1
    assert profile['export_phases']['operations']
    logs = list(cache_output.rglob('*logs.csv'))
    assert logs


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
