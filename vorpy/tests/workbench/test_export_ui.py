import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from pathlib import Path
import pytest
from PySide6.QtWidgets import QApplication, QDialog

from vorpy.workbench.domain import AnalysisResult
from vorpy.workbench.export_config import ExportConfig, OPTION_LABELS
from vorpy.workbench.services.export import PRESETS
from vorpy.workbench.ui.export_contents_dialog import ExportContentsDialog
from vorpy.workbench.ui import export_panel
from vorpy.workbench.ui.theme import STUDIO_STYLESHEET


@pytest.fixture(scope='module')
def app():
    app = QApplication.instance() or QApplication([])
    app.setStyleSheet(STUDIO_STYLESHEET)
    return app


@pytest.mark.parametrize('preset,expected', [
    ('Small', {'info', 'shell_surfs', 'logs'}),
    ('Medium', {'info', 'shell_surfs', 'surfs', 'shell_edges', 'edges', 'shell_verts', 'verts',
                'logs', 'atoms', 'surr_atoms', 'surr_resids'}),
    ('Large', {'shell_verts', 'shell_edges', 'shell_surfs', 'info', 'edges', 'verts',
               'atoms', 'surr_atoms', 'surr_resids', 'logs', 'atom_surfs', 'atom_edges', 'atom_verts'}),
])
def test_preset_cards_and_dialog_preserve_legacy_options(app, preset, expected):
    panel = export_panel.ExportPanel(lambda: None)
    panel.size_buttons[preset].click()
    assert set(panel.config.options) == expected == set(PRESETS[preset])
    dialog = ExportContentsDialog(ExportConfig())
    dialog.preset_buttons[preset].click()
    assert {key for key, check in dialog.checks.items() if check.isChecked()} == expected
    assert len(dialog.checks) == len(OPTION_LABELS) == 18
    panel.close()
    dialog.close()


def test_custom_recognition_and_reset(app):
    original = ExportConfig()
    dialog = ExportContentsDialog(original)
    dialog.checks['info'].click()
    assert dialog.config.preset == 'Custom'
    assert dialog.preset_buttons['Custom'].isChecked()
    assert original.options.get('info', False)
    dialog.checks['info'].click()
    assert dialog.config.preset == 'Medium'
    dialog.checks['sep_verts'].click()
    dialog.reset()
    assert dialog.config == original
    assert dialog.count_label.text() == '11 export components selected'
    dialog.close()


@pytest.mark.parametrize('accepted', [True, False])
def test_dialog_apply_cancel_and_summary(app, monkeypatch, accepted):
    panel = export_panel.ExportPanel(lambda: None)
    def edit(dialog):
        dialog.checks['sep_verts'].click()
        return QDialog.Accepted if accepted else QDialog.Rejected
    monkeypatch.setattr(ExportContentsDialog, 'exec', edit)
    panel.edit_custom()
    assert bool(panel.config.options.get('sep_verts')) == accepted
    assert ('Separate vertices' in panel.summary.text()) == accepted
    assert ('CUSTOM' in panel.summary_title.text()) == accepted
    assert panel.config.preset == ('Custom' if accepted else 'Medium')
    panel.close()


@pytest.mark.parametrize('custom', [False, True])
def test_export_uses_authoritative_config(app, monkeypatch, tmp_path, custom):
    result = AnalysisResult(tmp_path / 'sample.pdb', 'sample', export_group=object())
    panel = export_panel.ExportPanel(lambda: result)
    panel.update_state(result, False)
    assert panel.size_buttons['Medium'].text() == 'Standard\nRecommended'
    assert all('Medium' not in button.text() for button in panel.size_buttons.values())
    assert panel.config.preset == 'Medium'
    if custom:
        panel.config.options = {'sep_verts': True, 'atoms': True}
        panel.refresh_summary()
    calls = []
    monkeypatch.setattr(export_panel, 'export_result', lambda *args, **kwargs: calls.append((args, kwargs)) or tmp_path / 'written')
    panel.directory.setText(str(tmp_path / 'chosen'))
    panel.round_to.setValue(8)
    panel.file_type.setCurrentText('VTP')
    panel.export_button.click()
    assert calls == [((result, str(tmp_path / 'chosen'), 'Custom' if custom else 'Medium',
                      panel.config.options, 8, 'vtp'), {'atom_format': 'pdb', 'color_overrides': {}})]
    assert str(tmp_path / 'chosen') in panel.destination.text()
    assert 'Exported to' in panel.hint.text()
    panel.close()


def test_unsolved_configuration_and_compact_layout(app):
    panel = export_panel.ExportPanel(lambda: None)
    panel.resize(540, 280)
    panel.show()
    app.processEvents()
    assert not panel.export_button.isEnabled()
    assert panel.custom_button.isEnabled() and panel.directory.isEnabled()
    panel.size_buttons['Small'].click()
    assert panel.config.preset == 'Small'
    assert panel.summary.text() == 'Network info • Logs • Shell surfaces'
    assert 'Solve' in panel.hint.text()
    assert panel.scroll.horizontalScrollBar().maximum() == 0
    assert panel.rect().contains(panel.export_button.geometry())
    assert panel.export_button.isVisibleTo(panel)
    for width in (680, 1000):
        panel.resize(width, 300)
        app.processEvents()
        assert panel.scroll.horizontalScrollBar().maximum() == 0
        widths = [button.width() for button in panel.size_buttons.values()]
        assert max(widths) - min(widths) <= 1
    panel.close()


def test_new_export_controls_and_session_action(app, monkeypatch):
    from dataclasses import asdict
    from PySide6.QtGui import QAction, QColor
    calls = []
    action = QAction('Save session')
    action.triggered.connect(lambda: calls.append('save'))
    panel = export_panel.ExportPanel(lambda: None, action)
    panel.atom_format.setCurrentText('XYZ')
    assert panel.config.atom_format == 'xyz'
    monkeypatch.setattr(export_panel.QColorDialog, 'getColor', lambda *args: QColor('#123456'))
    panel.color_buttons['edges'].click()
    assert panel.config.color_overrides == {'edges': '#123456'}
    dialog = ExportContentsDialog(panel.config)
    dialog.checks['full_molecule'].click()
    assert dialog.config.options['full_molecule']
    assert dialog.config.preset == 'Custom'
    panel.config.options = dialog.config.options
    state = asdict(panel.config)
    panel.restore_config({})
    assert panel.config.atom_format == 'pdb'
    panel.restore_config(state)
    assert asdict(panel.config) == state
    assert 'Full molecule' in panel.summary.text()
    panel.save_session_button.click()
    assert calls == ['save']
    action.setEnabled(False)
    assert not panel.save_session_button.isEnabled()
    panel.viewer_colors_button.click()
    assert panel.config.color_overrides == {}
    assert 'Uses current viewer colors' in panel.color_hint.text()
    panel.close()
    dialog.close()
