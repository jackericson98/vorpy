from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from vorpy.src.command.vpy_cmnd import (
    _atomic_copy,
    _check_archive_target_available,
    _prepare_output_directory,
    _write_cache_export_profile,
    fail_if_exists_cli_options,
)


ROOT = Path(__file__).resolve().parents[2]
WATCHDOG = ROOT / 'scripts' / 'nexus_pow_watchdog.ps1'


def test_fail_if_exists_is_opt_in_and_removed_before_legacy_parsing():
    remaining, enabled = fail_if_exists_cli_options(
        ['input.pdb', '--fail-if-exists', '-e', 'small']
    )

    assert remaining == ['input.pdb', '-e', 'small']
    assert enabled is True
    remaining, enabled = fail_if_exists_cli_options(['input.pdb', '-e', 'small'])
    assert remaining == ['input.pdb', '-e', 'small']
    assert enabled is False


def test_existing_output_and_archive_targets_fail_without_mutation(tmp_path):
    output = tmp_path / 'pilot'
    output.mkdir()
    with pytest.raises(FileExistsError, match='existing output directory'):
        _prepare_output_directory(output, fail_if_exists=True)

    archive = tmp_path / 'pilot.vpy'
    archive.write_bytes(b'keep')
    with pytest.raises(SystemExit, match='existing archive'):
        _check_archive_target_available(archive)
    assert archive.read_bytes() == b'keep'

    # The reservation is opt-in; the legacy copy helper still replaces by default.
    source = tmp_path / 'replacement.vpy'
    source.write_bytes(b'new')
    _atomic_copy(source, archive)
    assert archive.read_bytes() == b'new'


def test_output_directory_reservation_allows_one_racing_creator(tmp_path):
    target = tmp_path / 'racing-output'

    def reserve():
        try:
            _prepare_output_directory(target, fail_if_exists=True)
        except FileExistsError:
            return 'exists'
        return 'created'

    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = list(pool.map(lambda _: reserve(), range(2)))

    assert sorted(outcomes) == ['created', 'exists']
    assert target.is_dir()


def test_archive_copy_reservation_allows_one_racing_writer(tmp_path):
    sources = [tmp_path / 'one.vpy', tmp_path / 'two.vpy']
    sources[0].write_bytes(b'one')
    sources[1].write_bytes(b'two')
    destination = tmp_path / 'reserved.vpy'

    def copy(source):
        try:
            _atomic_copy(source, destination, fail_if_exists=True)
        except FileExistsError:
            return 'exists'
        return 'created'

    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = list(pool.map(copy, sources))

    assert sorted(outcomes) == ['created', 'exists']
    assert destination.read_bytes() in {b'one', b'two'}


def test_cache_profile_reservation_preserves_existing_target(tmp_path):
    profile = tmp_path / 'cache_export_profile.json'
    profile.write_bytes(b'keep')
    system = type('SystemStub', (), {
        'files': {'dir': str(tmp_path)},
        '_fail_if_exists': True,
    })()

    with pytest.raises(FileExistsError, match='cache export profile'):
        _write_cache_export_profile(system, tmp_path / 'input.vpy', {'timing': {}})
    assert profile.read_bytes() == b'keep'


def test_approved_watchdog_requires_reservation_and_timing_sidecar_checks():
    if not WATCHDOG.is_file():
        pytest.skip('integration-owned watchdog is absent in this checkout')
    source = WATCHDOG.read_text(encoding='utf-8')
    assert '--fail-if-exists' in source
    assert 'Watchdog requires the approved --fail-if-exists reservation flag.' in source
    assert 'Refusing to overwrite an existing output root' in source
    assert 'Refusing to overwrite an existing archive timing sidecar' in source
