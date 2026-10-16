from types import SimpleNamespace

import pytest

from vorpy.src.io import ArchiveError
from vorpy.src.io.archive_codec import Writer
from vorpy.src.output import output


def _progress():
    system = SimpleNamespace(
        verbose=False,
        update_progress=lambda **kwargs: None,
    )
    return output.ExportProgress(1, system)


def test_export_failure_records_stage_and_reraises_original_traceback(capsys):
    progress = _progress()

    def failing_export():
        raise RuntimeError('synthetic canonical log failure')

    with pytest.raises(RuntimeError, match='synthetic canonical log failure') as caught:
        output._run_export(progress, 'interface: canonical logs', failing_export)

    diagnostic = progress.sys.export_failure
    assert diagnostic == {
        'stage': 'interface: canonical logs',
        'callable': 'test_export_failure_records_stage_and_reraises_original_traceback.<locals>.failing_export',
        'exception_type': 'RuntimeError',
        'message': 'synthetic canonical log failure',
    }
    assert any(frame.name == 'failing_export' for frame in __import__('traceback').extract_tb(caught.value.__traceback__))
    assert progress.current == 0
    assert "stage='interface: canonical logs'" in capsys.readouterr().err


def test_canonical_log_failure_is_stage_labeled_and_does_not_create_file(tmp_path, monkeypatch):
    progress = _progress()
    interface = SimpleNamespace(dir=str(tmp_path))

    def fail_rows(_iface):
        raise TypeError('unsupported interface-log value: function')

    monkeypatch.setattr(output, '_export_canonical_interface_log', output._export_canonical_interface_log)
    from vorpy.src.output import visualization
    monkeypatch.setattr(visualization, '_compact_result_rows', fail_rows)

    with pytest.raises(TypeError, match='unsupported interface-log value'):
        output._run_export(
            progress,
            'interface: canonical logs',
            output._export_canonical_interface_log,
            iface=interface,
        )

    assert not (tmp_path / 'logs.csv').exists()
    assert progress.sys.export_failure['stage'] == 'interface: canonical logs'
    assert progress.sys.export_failure['exception_type'] == 'TypeError'


def test_archive_failure_is_stage_labeled_and_does_not_mark_archive_complete(tmp_path, monkeypatch):
    progress = _progress()
    system = SimpleNamespace(
        name='synthetic',
        files={'dir': str(tmp_path)},
        _last_export_archive_path=None,
        update_progress=lambda **kwargs: None,
    )
    progress.sys = system

    def fail_save(*args, **kwargs):
        raise ArchiveError('unsupported scientific value type: function')

    import vorpy.src.io as io
    monkeypatch.setattr(io, 'save_network', fail_save)

    with pytest.raises(ArchiveError, match='unsupported scientific value type'):
        output._run_export(
            progress,
            'network archive',
            output._export_network_archive,
            sys=system,
        )

    assert not (tmp_path / 'synthetic.vpy').exists()
    assert system._last_export_archive_path is None
    assert progress.sys.export_failure['stage'] == 'network archive'
    assert progress.sys.export_failure['exception_type'] == 'ArchiveError'


def test_archive_encoder_rejects_unsupported_scientific_values():
    with pytest.raises(ArchiveError, match='Unsupported scientific value type: function'):
        Writer(None, {}).encode(lambda: None)
