"""Export workflow backed by a single configuration and the legacy exporters."""
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QApplication, QButtonGroup, QComboBox, QDialog, QFileDialog, QColorDialog,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton,
    QSizePolicy, QSpinBox, QVBoxLayout, QWidget,
)
from vorpy.workbench.export_config import ExportConfig, OPTION_LABELS, PRESET_LABELS
from vorpy.workbench.services.export import export_result
from vorpy.workbench.ui.export_contents_dialog import ExportContentsDialog, section_label
from vorpy.workbench.ui.panels import scroll_panel


class ExportPanel(QWidget):
    def __init__(self, result_provider, save_action=None):
        super().__init__()
        self.result_provider = result_provider
        self.config = ExportConfig()
        self._busy = False
        self._result = None
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(4)
        content = QWidget()
        body = QVBoxLayout(content)
        body.setContentsMargins(4, 2, 4, 2)
        body.setSpacing(6)
        body.addWidget(section_label('EXPORT PRESET'))
        cards = QHBoxLayout()
        self.size_group = QButtonGroup(self)
        self.size_buttons = {}
        for key, subtitle in (('Small', 'Essentials'), ('Medium', 'Recommended'), ('Large', 'Complete geometry')):
            button = QPushButton(f'{PRESET_LABELS[key]}\n{subtitle}')
            button.setObjectName('exportPreset')
            button.setCheckable(True)
            button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            button.clicked.connect(lambda checked=False, key=key: self.select_preset(key))
            self.size_group.addButton(button)
            self.size_buttons[key] = button
            cards.addWidget(button, 1)
        body.addLayout(cards)
        contents = QHBoxLayout()
        contents.addWidget(section_label('EXPORT CONTENTS'))
        self.custom_button = QPushButton('Customize contents…')
        self.custom_button.clicked.connect(self.edit_custom)
        contents.addWidget(self.custom_button)
        contents.addStretch()
        body.addLayout(contents)
        body.addWidget(section_label('OUTPUT'))
        output = QWidget()
        output.setMaximumWidth(820)
        output_layout = QVBoxLayout(output)
        output_layout.setContentsMargins(0, 0, 0, 0)
        output_layout.setSpacing(4)
        location = QHBoxLayout()
        location.addWidget(QLabel('Directory'))
        self.directory = QLineEdit()
        self.directory.setMinimumWidth(0)
        self.directory.setPlaceholderText('Default output directory (beside source)')
        self.directory.setAccessibleName('Output directory')
        self.directory.textChanged.connect(self.directory_changed)
        location.addWidget(self.directory, 1)
        self.browse_button = QPushButton('Browse…')
        self.browse_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.browse_button.clicked.connect(self.browse)
        location.addWidget(self.browse_button)
        output_layout.addLayout(location)
        options = QHBoxLayout()
        options.addWidget(QLabel('Decimal precision'))
        self.round_to = QSpinBox()
        self.round_to.setRange(0, 12)
        self.round_to.setValue(self.config.precision)
        self.round_to.setAccessibleName('Decimal precision for numeric logs')
        self.round_to.setToolTip('Decimal precision for numeric log output')
        self.round_to.valueChanged.connect(lambda value: setattr(self.config, 'precision', value))
        options.addWidget(self.round_to)
        options.addSpacing(12)
        options.addWidget(QLabel('Surface format'))
        self.file_type = QComboBox()
        self.file_type.addItems(['OFF', 'PLY', 'VTP'])
        self.file_type.setAccessibleName('Surface format')
        self.file_type.setToolTip('File format for surfaces, edges, and vertices')
        self.file_type.currentTextChanged.connect(lambda value: setattr(self.config, 'file_type', value.lower()))
        options.addWidget(self.file_type)
        options.addStretch()
        output_layout.addLayout(options)
        atom_options = QHBoxLayout()
        atom_options.addWidget(QLabel('Atom structure format'))
        self.atom_format = QComboBox()
        self.atom_format.addItems(['PDB', 'XYZ'])
        self.atom_format.setAccessibleName('Atom structure format')
        self.atom_format.currentTextChanged.connect(lambda value: setattr(self.config, 'atom_format', value.lower()))
        atom_options.addWidget(self.atom_format)
        atom_options.addStretch()
        output_layout.addLayout(atom_options)
        precision_hint = QLabel('Decimal places in numeric logs; atom formats retain their coordinate precision.')
        precision_hint.setWordWrap(True)
        precision_hint.setObjectName('sectionLabel')
        output_layout.addWidget(precision_hint)
        body.addWidget(output)
        body.addWidget(section_label('GEOMETRY COLORS'))
        colors = QHBoxLayout()
        self.color_buttons = {}
        for kind, label in (('surfaces', 'Surfaces'), ('edges', 'Edges'), ('vertices', 'Vertices')):
            button = QPushButton(label)
            button.clicked.connect(lambda checked=False, kind=kind: self.choose_color(kind))
            self.color_buttons[kind] = button
            colors.addWidget(button)
        self.viewer_colors_button = QPushButton('Use viewer colors')
        self.viewer_colors_button.clicked.connect(self.use_viewer_colors)
        colors.addWidget(self.viewer_colors_button)
        colors.addStretch()
        body.addLayout(colors)
        self.color_hint = QLabel()
        self.color_hint.setWordWrap(True)
        body.addWidget(self.color_hint)
        self.refresh_color_controls()
        self.summary_title = section_label('')
        body.addWidget(self.summary_title)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setMinimumWidth(0)
        self.summary.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        body.addWidget(self.summary)
        self.destination = QLabel()
        self.destination.setWordWrap(True)
        self.destination.setMinimumWidth(0)
        self.destination.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        self.destination.setTextFormat(Qt.PlainText)
        body.addWidget(self.destination)
        body.addStretch()
        self.scroll = scroll_panel(content)
        layout.addWidget(self.scroll, 1)
        footer = QHBoxLayout()
        self.hint = QLabel()
        self.hint.setTextFormat(Qt.PlainText)
        self.hint.setWordWrap(True)
        self.hint.setMinimumWidth(0)
        self.hint.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        footer.addWidget(self.hint, 1)
        self.save_session_button = QPushButton('Save session…')
        self.save_session_button.setToolTip('Save the Workbench project, including current export settings')
        self.save_session_button.setEnabled(save_action is not None and save_action.isEnabled())
        if save_action is not None:
            self.save_session_button.clicked.connect(save_action.trigger)
            save_action.changed.connect(lambda: self.save_session_button.setEnabled(save_action.isEnabled()))
        footer.addWidget(self.save_session_button, alignment=Qt.AlignBottom)
        self.export_button = QPushButton('Export Results')
        self.export_button.setObjectName('primaryAction')
        self.export_button.clicked.connect(self.export)
        footer.addWidget(self.export_button, alignment=Qt.AlignRight | Qt.AlignBottom)
        layout.addLayout(footer)
        self.refresh_summary()
        self.update_state(None, False)

    def choose_color(self, kind):
        layer = next((layer for layer in (self._result.layers if self._result else []) if layer.kind == kind), None)
        initial = self.config.color_overrides.get(kind, layer.color if layer else '#806df0')
        color = QColorDialog.getColor(QColor(initial), self, f'Export {kind} color')
        if color.isValid():
            self.config.color_overrides[kind] = color.name()
            self.refresh_color_controls()

    def use_viewer_colors(self):
        self.config.color_overrides.clear()
        self.refresh_color_controls()

    def refresh_color_controls(self):
        for kind, button in self.color_buttons.items():
            color = self.config.color_overrides.get(kind)
            button.setStyleSheet(f'border: 2px solid {color};' if color else '')
            button.setToolTip(f'Export override: {color}' if color else 'Uses current viewer colors; click to override')
        overrides = ', '.join(self.config.color_overrides)
        self.color_hint.setText(f'Solid overrides: {overrides}. Other geometry follows the viewer.' if overrides
                                else 'Uses current viewer colors, colormap, and scale. Click a geometry type to override.')

    def restore_config(self, state):
        from vorpy.workbench.export_config import config_from_json
        self.config = config_from_json(state)
        self.directory.setText(self.config.directory)
        self.round_to.setValue(self.config.precision)
        self.file_type.setCurrentText(self.config.file_type.upper())
        self.atom_format.setCurrentText(self.config.atom_format.upper())
        self.refresh_color_controls()
        self.refresh_summary()

    def select_preset(self, key):
        self.config.select_preset(key)
        self.refresh_summary()

    def directory_changed(self, value):
        self.config.directory = value
        self.refresh_summary()

    def refresh_summary(self):
        self.size_group.setExclusive(False)
        for key, button in self.size_buttons.items():
            button.setChecked(key == self.config.preset)
        self.size_group.setExclusive(True)
        self.summary_title.setText(f'EXPORT SUMMARY · {PRESET_LABELS[self.config.preset].upper()}')
        self.summary.setText(' • '.join(label for key, label in OPTION_LABELS.items()
                                      if self.config.options.get(key, False)) or 'No export components selected.')
        directory = self.config.output_directory(self._result)
        self.destination.setText('Output: ' + (directory or 'Default output directory (beside source)'))
        self.destination.setToolTip(directory)

    def update_state(self, result, busy):
        self._result, self._busy = result, busy
        ready = result is not None and result.export_group is not None and not result.defaults_stale
        self.export_button.setEnabled(ready and not busy)
        if busy:
            status = 'ⓘ Wait for the current operation to finish.'
        elif result is not None and result.defaults_stale:
            status = 'ⓘ Solve again to export the updated structure.'
        elif not ready:
            status = ('ⓘ Solve this imported result again to enable full export.' if result is not None and result.layers
                      else 'ⓘ Solve the current structure to enable export.')
        else:
            status = f'{result.name} · frame {result.frame_index} · A new folder is created for each export.'
        self.hint.setText(status)
        self.refresh_summary()

    def browse(self):
        directory = QFileDialog.getExistingDirectory(self, 'Choose Export Location',
                                                     self.config.output_directory(self._result))
        if directory:
            self.directory.setText(directory)

    def edit_custom(self):
        dialog = ExportContentsDialog(self.config, self)
        if dialog.exec() == QDialog.Accepted:
            self.config.options = dict(dialog.config.options)
            self.refresh_summary()

    def export(self):
        result = self.result_provider()
        if result is None or result.export_group is None or result.defaults_stale or self._busy:
            return
        if not any(self.config.options.values()):
            QMessageBox.information(self, 'Export', 'Select at least one export component in Customize contents.')
            return
        directory = self.config.output_directory(result)
        if not directory:
            self.browse()
            directory = self.config.output_directory(result)
            if not directory:
                return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.export_button.setEnabled(False)
        error = None
        try:
            destination = export_result(result, directory, self.config.preset, self.config.options,
                                        self.config.precision, self.config.file_type,
                                        atom_format=self.config.atom_format,
                                        color_overrides=self.config.color_overrides)
        except Exception as exc:
            error = str(exc)
        finally:
            QApplication.restoreOverrideCursor()
            self.update_state(self.result_provider(), self._busy)
        if error is not None:
            self.hint.setText('Export failed. Check the output location and try again.')
            QMessageBox.warning(self, 'Export failed', error)
        else:
            self.hint.setText(f'Exported to {destination}')
