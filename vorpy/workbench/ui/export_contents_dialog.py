"""Transactional editor for export contents; Cancel leaves the shared state intact."""
from copy import deepcopy

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QButtonGroup, QCheckBox, QDialog, QDialogButtonBox, QGridLayout,
    QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
)
from vorpy.workbench.export_config import (
    PRESET_LABELS, GENERAL_OPTIONS, ATOMIC_OPTIONS, GEOMETRY_ROWS, GEOMETRY_MODES,
)


def section_label(text):
    label = QLabel(text)
    label.setObjectName('sectionLabel')
    return label


class ExportContentsDialog(QDialog):
    def __init__(self, config, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Export Contents')
        self.initial_config = deepcopy(config)
        self.config = deepcopy(config)
        layout = QVBoxLayout(self)
        layout.setSpacing(10)
        layout.addWidget(section_label('PRESET'))
        presets = QHBoxLayout()
        self.preset_group = QButtonGroup(self)
        self.preset_buttons = {}
        for key, label in PRESET_LABELS.items():
            button = QPushButton(label)
            button.setObjectName('exportPreset')
            button.setCheckable(True)
            self.preset_group.addButton(button)
            self.preset_buttons[key] = button
            button.clicked.connect(lambda checked=False, key=key: self.select_preset(key))
            presets.addWidget(button)
        layout.addLayout(presets)
        self.checks = {}
        for title, options in (('GENERAL', GENERAL_OPTIONS), ('ATOMIC STRUCTURES', ATOMIC_OPTIONS)):
            layout.addWidget(section_label(title))
            row = QGridLayout()
            for index, (key, label) in enumerate(options):
                row.addWidget(self.make_check(key, label), index // 2, index % 2)
            layout.addLayout(row)
        layout.addWidget(section_label('GEOMETRY'))
        matrix = QGridLayout()
        matrix.setHorizontalSpacing(18)
        matrix.setVerticalSpacing(10)
        descriptions = ('All geometry together', 'Group boundary geometry',
                        'Individual geometry files', 'Geometry for each atom cell')
        for column, ((mode, _), description) in enumerate(zip(GEOMETRY_MODES, descriptions), 1):
            label = QLabel(mode)
            label.setToolTip(description)
            matrix.addWidget(label, 0, column, Qt.AlignCenter)
            matrix.setColumnStretch(column, 1)
        for row, (label, key) in enumerate(GEOMETRY_ROWS, 1):
            matrix.addWidget(QLabel(label), row, 0)
            for column, (mode, prefix) in enumerate(GEOMETRY_MODES, 1):
                check = self.make_check(prefix + key, '')
                check.setAccessibleName(f'{mode} {label.lower()}')
                check.setToolTip(f'{mode} {label.lower()}')
                matrix.addWidget(check, row, column, Qt.AlignCenter)
        layout.addLayout(matrix)
        self.count_label = section_label('')
        layout.addWidget(self.count_label)
        self.buttons = QDialogButtonBox(QDialogButtonBox.Reset | QDialogButtonBox.Cancel | QDialogButtonBox.Apply)
        self.buttons.button(QDialogButtonBox.Reset).clicked.connect(self.reset)
        self.buttons.button(QDialogButtonBox.Apply).clicked.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.sync_controls()

    def make_check(self, key, text):
        check = QCheckBox(text)
        self.checks[key] = check
        check.toggled.connect(lambda enabled, key=key: self.set_option(key, enabled))
        return check

    def set_option(self, key, enabled):
        self.config.options[key] = enabled
        self.update_indicator()

    def select_preset(self, key):
        if key != 'Custom':
            self.config.select_preset(key)
            self.sync_controls()
        # Custom retains all current selections until the user edits them.

    def sync_controls(self):
        for key, check in self.checks.items():
            check.blockSignals(True)
            check.setChecked(self.config.options.get(key, False))
            check.blockSignals(False)
        self.update_indicator()

    def update_indicator(self):
        self.preset_buttons[self.config.preset].setChecked(True)
        count = sum(bool(value) for value in self.config.options.values())
        self.count_label.setText(f'{count} export components selected')

    def reset(self):
        self.config = deepcopy(self.initial_config)
        self.sync_controls()
