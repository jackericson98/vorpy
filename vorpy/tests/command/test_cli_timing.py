import json
from types import SimpleNamespace

from vorpy.src.command.vpy_cmnd import Command, timing_report_cli_options


def test_timing_report_option_is_explicit_and_parser_safe(tmp_path):
    report = tmp_path / 'timing.json'
    remaining, parsed = timing_report_cli_options([
        'input.pdb', '--timing-report', str(report), '-i',
    ])

    assert remaining == ['input.pdb', '-i']
    assert parsed == report.resolve()


def test_timing_report_uses_forge_vector_and_atlas_events(tmp_path):
    network = SimpleNamespace(
        metrics={'tot': 8.0, 'anal': 2.5},
        group_name='source',
    )
    system = SimpleNamespace(
        name='synthetic',
        groups=[SimpleNamespace(name='source', net=network)],
        ifaces=[],
        files={'dir': str(tmp_path)},
    )
    command = Command(sys=system)
    command.base_file = tmp_path / 'synthetic.pdb'
    command.settings_dict = {'net_type': 'pow'}
    command.solvent_mode = 'MOLECULE_ONLY'
    command.timing_report_path = tmp_path / 'timing.json'
    command._cli_timing = {
        'schema_version': 1,
        'unit': 'seconds',
        'phases': {},
    }
    command._cli_timing_started = (0.0, 0.0)

    command._record_network_timings()
    system.export_timing = {
        'system info': 3.0,
        'network archive': 1.0,
        'total': 4.0,
    }
    system.archive_timing = {'total': 1.0}
    command._record_export_archive_timings(None)
    command._write_timing_report(tmp_path)

    payload = json.loads(command.timing_report_path.read_text(encoding='utf-8'))
    assert payload['format'] == 'vorpy.cli_timing'
    assert payload['unit'] == 'seconds'
    assert payload['phases']['geometry_network']['wall_seconds'] == 5.5
    assert payload['phases']['analysis']['wall_seconds'] == 2.5
    assert payload['phases']['export']['wall_seconds'] == 3.0
    assert payload['phases']['archive_save']['wall_seconds'] == 1.0
    assert payload['phases']['geometry_network']['source'].startswith('FORGE')
    assert payload['phases']['export']['source'].startswith('VECTOR')
    assert payload['phases']['archive_save']['source'].startswith('ATLAS')
    assert payload['phases']['total_wall']['source'].startswith('PRISM')
