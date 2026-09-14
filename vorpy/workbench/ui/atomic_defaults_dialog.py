"""Atomic-property editor with a detached draft and explicit inheritance."""
from copy import deepcopy

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QButtonGroup, QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QGridLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QPushButton, QTableWidget, QTableWidgetItem, QTabWidget,
    QVBoxLayout, QWidget,
)
from vorpy.workbench.atomic_defaults import PROPERTIES, ION_RESIDUES, element_defaults, validate_defaults
from vorpy.workbench.ui.atomic_data import (
    _PERIODIC_TABLE, _ELEMENT_NUMBERS, _ELEMENT_NAMES, _PERIODIC_COLORS,
    _element_category, RESIDUE_ATOMS, _RESIDUE_NAMES, ION_NAMES,
)
from vorpy.workbench.ui.panels import scroll_panel

PROPERTY_LABELS = {"radius": "Radius (Å)", "mass": "Mass (Da)", "charge": "Charge (e)"}
PROPERTY_HELP = {
    "radius": "Sphere radius used in the viewer and the next VorPy build.",
    "mass": "Atomic mass in daltons, passed to the next build for mass-dependent properties.",
    "charge": "Signed atomic charge in units of elementary charge. Unset retains the source value. "
              "Charge is stored with the atom; it does not add an electrostatic calculation to VorPy.",
}


class PropertyEditor(QWidget):
    """An unset override is different from an explicit value, including zero charge."""
    def __init__(self, name, inherited, override, changed):
        super().__init__()
        self.setToolTip(PROPERTY_HELP[name])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(3, 2, 3, 2)
        layout.setSpacing(2)
        row = QHBoxLayout()
        row.setSpacing(4)
        self.enabled = QCheckBox("Override")
        self.value = QDoubleSpinBox()
        self.value.setDecimals(4)
        self.value.setRange(-1e6 if name == "charge" else 0.0001, 1e6)
        self.value.setSingleStep(0.1 if name == "charge" else 0.01)
        self.value.setValue(override if override is not None else inherited if inherited is not None else 0 if name == "charge" else 1)
        self.enabled.setChecked(override is not None)
        self.value.setEnabled(override is not None)
        self.enabled.setToolTip("Uncheck to inherit the element or original atom value.")
        self.value.setAccessibleName(PROPERTY_LABELS[name])
        row.addWidget(self.enabled)
        row.addWidget(self.value, 1)
        layout.addLayout(row)
        self.inherited = QLabel(f"Inherited: {inherited:g}" if inherited is not None else "Inherited: source / unset")
        self.inherited.setObjectName("sectionLabel")
        layout.addWidget(self.inherited)

        def toggle(checked):
            self.value.setEnabled(checked)
            if not checked and inherited is not None:
                self.value.setValue(inherited)
            changed(self.value.value() if checked else None)
        self.enabled.toggled.connect(toggle)
        self.value.valueChanged.connect(lambda value: changed(value) if self.enabled.isChecked() else None)


class AtomicDefaultsDialog(QDialog):
    def __init__(self, defaults, atoms=(), parent=None):
        super().__init__(parent)
        self.setWindowTitle("Atomic Defaults")
        self.resize(1280, 760)
        self.setMinimumSize(850, 540)
        self._draft = validate_defaults(defaults)
        self._atoms = list(atoms)
        self._selected_element = "C"
        self._residues = {name: {atom: atom[0] for atom in names} for name, names in RESIDUE_ATOMS.items()}
        # PDB identities disambiguate CA (alpha carbon) from calcium ions.
        for atom in self._atoms:
            self._residues.setdefault(atom.residue_name.upper(), {})[atom.name] = atom.element.upper()
        self._ions = sorted(set(ION_NAMES) | {atom.element.upper() for atom in self._atoms if atom.residue_name.upper() in ION_RESIDUES} | {key.split(':')[1] for key in self._draft if key.startswith('ION:')})
        for key in self._draft:
            if key.startswith('RES:') and len(key.split(':')) == 3:
                _, residue, atom = key.split(':')
                self._residues.setdefault(residue, {}).setdefault(atom, atom[0])
        layout = QVBoxLayout(self)
        description = QLabel("Project atomic defaults · Residue atom → Ion → Element → Original atom value")
        description.setWordWrap(True)
        layout.addWidget(description)
        hint = QLabel("Enable Override to change a property; uncheck it to inherit. Apply updates atoms and marks existing analyses for rebuilding.")
        hint.setWordWrap(True)
        hint.setObjectName("sectionLabel")
        layout.addWidget(hint)
        self.tabs = QTabWidget()
        self.tabs.addTab(self._elements_page(), "Elements")
        self.tabs.addTab(self._residues_page(), "Residues")
        self.ion_table = self._property_table()
        self.tabs.addTab(self.ion_table, "Ions")
        self.tabs.currentChanged.connect(self._refresh_tab)
        layout.addWidget(self.tabs, 1)
        footer = QHBoxLayout()
        self.reset_all = QPushButton("Clear project overrides")
        self.reset_all.clicked.connect(self._reset)
        footer.addWidget(self.reset_all)
        self.changes = QLabel()
        footer.addWidget(self.changes, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Save).setText("Apply to project")
        buttons.button(QDialogButtonBox.Save).setObjectName("primaryAction")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        footer.addWidget(buttons)
        layout.addLayout(footer)
        self._select_element("C")
        self._update_count()

    def values(self):
        return deepcopy(self._draft)

    def _update_count(self):
        self.changes.setText(f"{sum(len(values) for values in self._draft.values())} property overrides")

    def _set(self, scope, name, value):
        if value is None:
            self._draft.get(scope, {}).pop(name, None)
            if not self._draft.get(scope):
                self._draft.pop(scope, None)
        else:
            self._draft.setdefault(scope, {})[name] = value
        self._update_count()

    def _inherited(self, scope, element, name):
        value = element_defaults(element)[name]
        if scope != element:
            value = self._draft.get(element, {}).get(name, value)
        if scope.startswith("RES:"):
            _, residue, atom_name = scope.split(':')
            matches = [atom for atom in self._atoms if atom.residue_name.upper() == residue and atom.name.upper() == atom_name]
            source_values = {(atom.source_properties or {key: getattr(atom, key) for key in PROPERTIES})[name] for atom in matches}
            if matches and name not in self._draft.get(element, {}):
                value = next(iter(source_values)) if len(source_values) == 1 else None
            if residue in ION_RESIDUES:
                value = self._draft.get(f"ION:{element}", {}).get(name, value)
        return value

    def _editor(self, scope, element, name):
        return PropertyEditor(name, self._inherited(scope, element, name),
                              self._draft.get(scope, {}).get(name),
                              lambda value: self._set(scope, name, value))

    def _elements_page(self):
        page = QWidget()
        layout = QHBoxLayout(page)
        grid_widget = QWidget()
        grid = QGridLayout(grid_widget)
        grid.setSpacing(3)
        grid.setAlignment(Qt.AlignTop | Qt.AlignLeft)
        self.element_buttons = {}
        self.element_group = QButtonGroup(self)
        self.element_group.setExclusive(True)
        for symbol, (row, column) in _PERIODIC_TABLE.items():
            button = QPushButton(f"{_ELEMENT_NUMBERS[symbol.upper()]}\n{symbol}")
            button.setFixedSize(46, 54)
            button.setCheckable(True)
            self.element_group.addButton(button)
            button.setToolTip(_ELEMENT_NAMES.get(symbol, symbol))
            button.setAccessibleName(f"Edit {_ELEMENT_NAMES.get(symbol, symbol)} atomic defaults")
            color = _PERIODIC_COLORS[_element_category(symbol)]
            button.setStyleSheet(f"QPushButton {{ background: {color}; color: #14202a; border: 1px solid transparent; padding: 2px; font-weight: 600; }} QPushButton:focus, QPushButton:checked {{ border: 3px solid #806df0; }}")
            button.clicked.connect(lambda checked=False, key=symbol.upper(): self._select_element(key))
            grid.addWidget(button, row, column)
            self.element_buttons[symbol.upper()] = button
        layout.addWidget(scroll_panel(grid_widget), 1)
        detail = QWidget()
        detail.setMinimumWidth(270)
        detail.setMaximumWidth(320)
        self.element_form = QVBoxLayout(detail)
        layout.addWidget(detail)
        return page

    def _select_element(self, element):
        self._selected_element = element
        if element in self.element_buttons:
            self.element_buttons[element].setChecked(True)
        while self.element_form.count():
            item = self.element_form.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        title = QLabel(f"{_ELEMENT_NAMES.get(element.title(), element.title())} · {element.title()}")
        title.setObjectName("viewSettingsTitle")
        self.element_form.addWidget(title)
        description = QLabel("Defaults for this element. Unset properties retain the source atom values. Specific residue and ion overrides take priority.")
        description.setWordWrap(True)
        self.element_form.addWidget(description)
        self.element_editors = {}
        for name in PROPERTIES:
            self.element_form.addWidget(QLabel(PROPERTY_LABELS[name]))
            editor = self._editor(element, element, name)
            self.element_editors[name] = editor
            self.element_form.addWidget(editor)
        self.element_form.addStretch()

    def _property_table(self):
        table = QTableWidget(0, 5)
        table.setHorizontalHeaderLabels(["Atom / ion", "Element", *PROPERTY_LABELS.values()])
        table.setAlternatingRowColors(True)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        for column in range(2, 5):
            table.horizontalHeader().setSectionResizeMode(column, QHeaderView.Interactive)
            table.setColumnWidth(column, 230)
        table.horizontalHeader().setMinimumSectionSize(75)
        return table

    def _residues_page(self):
        page = QWidget()
        layout = QHBoxLayout(page)
        sidebar = QWidget()
        side = QVBoxLayout(sidebar)
        sidebar.setMaximumWidth(230)
        search = QLineEdit()
        search.setPlaceholderText("Find residue…")
        side.addWidget(search)
        self.residue_list = QListWidget()
        for residue in sorted(self._residues):
            item = QListWidgetItem(f"{_RESIDUE_NAMES.get(residue, residue)} ({residue})")
            item.setData(Qt.UserRole, residue)
            self.residue_list.addItem(item)
        search.textChanged.connect(lambda text: [self.residue_list.item(i).setHidden(text.lower() not in self.residue_list.item(i).text().lower()) for i in range(self.residue_list.count())])
        side.addWidget(self.residue_list)
        self.residue_table = self._property_table()
        self.residue_list.currentItemChanged.connect(lambda *_: self._populate_residue())
        layout.addWidget(sidebar)
        layout.addWidget(self.residue_table, 1)
        self.residue_list.setCurrentRow(0)
        return page

    def _populate(self, table, rows):
        table.setRowCount(len(rows))
        for row, (label, element, scope) in enumerate(rows):
            table.setRowHeight(row, 66)
            table.setItem(row, 0, QTableWidgetItem(label))
            table.setItem(row, 1, QTableWidgetItem(element))
            for column, name in enumerate(PROPERTIES, 2):
                table.setCellWidget(row, column, self._editor(scope, element, name))

    def _populate_residue(self):
        item = self.residue_list.currentItem()
        if item:
            residue = item.data(Qt.UserRole)
            self._populate(self.residue_table, [(atom, element, f"RES:{residue}:{atom}".upper()) for atom, element in self._residues[residue].items()])

    def _refresh_tab(self, index):
        if index == 0:
            self._select_element(self._selected_element)
        elif index == 1:
            self._populate_residue()
        else:
            self._populate(self.ion_table, [(ion, ion, f"ION:{ion}") for ion in self._ions])

    def _reset(self):
        self._draft.clear()
        self._refresh_tab(self.tabs.currentIndex())
        self._update_count()
