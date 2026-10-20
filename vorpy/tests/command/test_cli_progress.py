import io
import sys
import time
from types import SimpleNamespace

import pandas as pd
import pytest

from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.io import ArchiveError
from vorpy.src.output.output import ExportProgress
from vorpy.src.system.system import System


class _Output(io.StringIO):
    def __init__(self, *, tty):
        super().__init__()
        self._tty = tty

    def isatty(self):
        return self._tty


def _synthetic_system(tmp_path):
    return System(
        spheres=pd.DataFrame({
            'loc': [
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (0.0, 2.0, 0.0), (0.0, 0.0, 2.0),
                (0.8, 0.8, 0.8),
            ],
            'rad': [0.75, 0.5, 0.75, 0.75, 0.35],
            'mass': [1.0] * 5,
            'num': list(range(5)),
            'res_name': ['ALA', 'ALA', 'ALA', 'ALA', 'HOH'],
            'name': ['CA'] * 5,
        }),
        name='synthetic-cli-progress',
        output_directory=str(tmp_path),
        make_dir=False,
        solvent_mode='MOLECULE_ONLY',
    )


@pytest.mark.parametrize('scheme', ['aw', 'pow', 'prm'])
@pytest.mark.parametrize('tty', [False, True])
def test_interface_cli_progress_has_single_stage_trace(tmp_path, monkeypatch, scheme, tty):
    system = _synthetic_system(tmp_path)
    command = Command(sys=system)
    command.base_file = str(tmp_path / 'synthetic.pdb')
    command._arguments = [
        command.base_file,
        '-g', 'a', '0', '-g', 'a', '4', '-i',
        '-s', 'nt', scheme, 'and', 'mv', '40',
        '--profile', 'geometry', '--boundary', 'none',
    ]
    command.parse_commands()
    command.run_exports = lambda: None
    command._save_archive = lambda: None

    stream = _Output(tty=tty)
    monkeypatch.setattr(sys, 'stdout', stream)
    command._process_system()
    output = stream.getvalue()

    assert f'Solving {scheme.upper()} interface network...' in output
    assert 'Network construction' in output
    assert 'Analyzing geometry' in output
    assert 'Export preparation...' in output
    assert 'Final summary...' in output
    assert 'complete (elapsed ' in output
    if not tty:
        assert '\r' not in output


def test_export_progress_reports_categories_counts_and_timings():
    events = []
    system = SimpleNamespace(
        verbose=False,
        update_progress=lambda **kwargs: events.append(kwargs),
    )
    progress = ExportProgress(3, system)
    for name in ('interface: info', 'interface: surfaces', 'network archive'):
        category = progress.show(name)
        progress.record(category, 0.25)
        progress.step()
    progress.finish()

    labels = [event['process'] for event in events]
    assert labels == [
        'Exporting scientific tables: interface: info (1/3)',
        'Exporting physical geometry: interface: surfaces (2/3)',
        'Saving archive: network archive (3/3)',
        'Exports complete',
    ]
    assert system.export_stage_counts == {
        'Exporting scientific tables': 1,
        'Exporting physical geometry': 1,
        'Saving archive': 1,
    }
    assert system.export_stage_timing == {
        'Exporting scientific tables': 0.25,
        'Exporting physical geometry': 0.25,
        'Saving archive': 0.25,
    }


def test_export_failure_reports_named_stage_and_keeps_original_exception(capsys):
    from vorpy.src.output.output import _run_export

    command = Command(sys=SimpleNamespace(
        verbose=False,
        update_progress=lambda **kwargs: None,
    ))
    command._cli_stage_enabled = True
    command.sys._cli_export_stage_callback = command._cli_export_stage_callback
    progress = ExportProgress(1, command.sys)

    def fail_export(**kwargs):
        raise RuntimeError('synthetic canonical-log failure')

    with pytest.raises(RuntimeError, match='synthetic canonical-log failure'):
        _run_export(
            progress, 'interface: canonical logs', fail_export,
        )

    assert command.sys.export_failure['stage'] == 'interface: canonical logs'
    assert 'Export stage failed: Canonical-log preparation and serialization' in capsys.readouterr().err


def test_slow_export_stage_is_visible_without_fabricated_percentage(capsys):
    from vorpy.src.output.output import _run_export

    command = Command(sys=SimpleNamespace(
        verbose=False,
        update_progress=lambda **kwargs: None,
    ))
    command._cli_stage_enabled = True
    command.sys._cli_export_stage_callback = command._cli_export_stage_callback
    progress = ExportProgress(1, command.sys)

    def slow_export(**kwargs):
        time.sleep(0.02)

    _run_export(progress, 'interface: canonical logs', slow_export)
    output = capsys.readouterr().out
    assert 'Canonical-log preparation and serialization (interface: canonical logs)...' in output
    assert 'Canonical-log preparation and serialization: complete (elapsed ' in output
    assert '%' not in output


def test_explicit_archive_failure_reports_stage_and_preserves_exception(
    tmp_path, monkeypatch, capsys,
):
    system = SimpleNamespace(
        frame_count=1,
        files={'dir': str(tmp_path)},
        _last_export_archive_path=None,
    )
    command = Command(sys=system)
    command.save_network_path = tmp_path / 'failed.vpy'
    command._cli_stage_enabled = True

    def fail_save(*args, **kwargs):
        raise ArchiveError('synthetic archive failure')

    monkeypatch.setattr('vorpy.src.io.save_network', fail_save)
    with pytest.raises(ArchiveError, match='synthetic archive failure'):
        command._save_archive()

    assert system.cli_stage_failure['stage'] == 'Saving archive'
    assert system.cli_stage_failure['exception_type'] == 'ArchiveError'
    assert 'Stage failed: Saving archive' in capsys.readouterr().err


def test_interface_bundle_reports_substages_without_profile_dump(tmp_path, capsys):
    from vorpy.src.command.vpy_cmnd import Command
    from vorpy.src.interface.export import interface_exports
    from vorpy.tests.test_nexus_interface_output_contract import _contract_case

    _, interface = _contract_case(tmp_path)
    system = interface.sys
    command = Command(sys=system)
    command.interface_mode = True
    command._cli_stage_enabled = True
    system._cli_export_stage_callback = command._cli_export_stage_callback

    interface_exports(
        interface,
        atoms=True,
        surfs=True,
        edges=True,
        verts=True,
        info=True,
        contact_tables=True,
    )

    output = capsys.readouterr().out
    assert 'Export preparation...' in output
    assert 'Contact/summary exports...' in output
    assert 'Exporting physical geometry...' in output
    assert 'INTERFACE EXPORT PROFILE' not in output
    assert 'complete (elapsed ' in output


def test_final_summary_marks_unresolved_science_separately_from_directory(
    tmp_path, capsys,
):
    system = SimpleNamespace(files={'dir': str(tmp_path)}, ifaces=[])
    command = Command(sys=system)
    command.interface_mode = True
    command._cli_stage_enabled = True
    command.interface_result_errors = ('FULL_SOLVENT extraction failed',)

    command._cli_report_final_summary(str(tmp_path), status='UNRESOLVED')

    output = capsys.readouterr().out
    assert 'Final summary: UNRESOLVED' in output
    assert 'Export directory:' in output
    assert 'complete' not in output.split('Final summary:', 1)[1].splitlines()[0]


def test_missing_alpha_capability_is_partial_not_failed(capsys, tmp_path):
    interface = SimpleNamespace(
        name='dna_ligand_interface',
        group1=SimpleNamespace(name='dna'),
        group2=SimpleNamespace(name='ligand'),
        net=SimpleNamespace(surfs=None, balls=None),
        source_view=None,
        geometry_analysis=SimpleNamespace(
            metadata={
                'alpha_status': 'UNRESOLVED',
                'alpha_failure_kind': 'GUDHI_UNAVAILABLE',
                'alpha_failure_reason': 'GUDHI is unavailable',
                'curvature_requested': False,
            },
            alpha_selection={'status': 'unresolved'},
            voronoi_side_1={},
        ),
    )
    system = SimpleNamespace(files={'dir': str(tmp_path)}, ifaces=[interface])
    command = Command(sys=system)
    command.interface_mode = True
    command._cli_stage_enabled = True
    command.interface_result_state = {'statuses': ('PARTIAL',)}

    assert command._cli_scientific_status() == 'PARTIAL'
    command._cli_report_final_summary(str(tmp_path))
    output = capsys.readouterr().out
    assert 'Final summary: PARTIAL' in output
    assert 'Alpha status: UNRESOLVED (optional capability unavailable:' in output
