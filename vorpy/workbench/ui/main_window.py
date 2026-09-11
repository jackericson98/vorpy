"""Viewer-forward scientific workbench application shell."""

from __future__ import annotations

from pathlib import Path
from dataclasses import asdict, replace
from uuid import uuid4

from PySide6.QtCore import QPoint, Qt, QThread
from PySide6.QtGui import QAction, QActionGroup, QCloseEvent, QColor, QBrush, QPainter, QPen, QPolygon
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDoubleSpinBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QGridLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSlider,
    QSpinBox,
    QSplitter,
    QStyle,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from vorpy.workbench.domain import AnalysisResult, Atom
from vorpy.workbench.project import (
    PROJECT_SUFFIX,
    AtomKey,
    GroupDefinition,
    InterfaceDefinition,
    Project,
    StructureSource,
    load_project,
    save_project,
    result_from_json,
    result_to_json,
)
from vorpy.workbench.services.result_directory import load_result_directory
from vorpy.workbench.services.info_parser import Measurement, NetworkSummary, parse_group_info, parse_network_summary
from vorpy.workbench.services.structure_loader import DISPLAY_RADII, load_pdb
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings
from vorpy.workbench.ui.molecular_view import (
    ION_RESIDUES,
    WATER_RESIDUES,
    MolecularView,
)
from vorpy.workbench.workers.solve_worker import SolveWorker

DEFAULT_DATA_DIRECTORY = Path(__file__).resolve().parents[2] / "data"

NETWORK_LAYER_OPTIONS = (
    ("edges", "Edges", "#72bde4"),
    ("vertices", "Vertices", "#efb84f"),
    ("surfaces", "Surfaces", "#4f9fcf"),
    ("shell_edges", "Shell edges", "#42d6c7"),
    ("shell_vertices", "Shell vertices", "#f29f67"),
    ("shell_surfaces", "Shell surfaces", "#806df0"),
)


_PERIODIC_TABLE = {
    "H": (0, 0), "He": (0, 17),
    "Li": (1, 0), "Be": (1, 1), "B": (1, 12), "C": (1, 13), "N": (1, 14), "O": (1, 15), "F": (1, 16), "Ne": (1, 17),
    "Na": (2, 0), "Mg": (2, 1), "Al": (2, 12), "Si": (2, 13), "P": (2, 14), "S": (2, 15), "Cl": (2, 16), "Ar": (2, 17),
    "K": (3, 0), "Ca": (3, 1), "Sc": (3, 2), "Ti": (3, 3), "V": (3, 4), "Cr": (3, 5), "Mn": (3, 6), "Fe": (3, 7), "Co": (3, 8), "Ni": (3, 9), "Cu": (3, 10), "Zn": (3, 11), "Ga": (3, 12), "Ge": (3, 13), "As": (3, 14), "Se": (3, 15), "Br": (3, 16), "Kr": (3, 17),
    "Rb": (4, 0), "Sr": (4, 1), "Y": (4, 2), "Zr": (4, 3), "Nb": (4, 4), "Mo": (4, 5), "Tc": (4, 6), "Ru": (4, 7), "Rh": (4, 8), "Pd": (4, 9), "Ag": (4, 10), "Cd": (4, 11), "In": (4, 12), "Sn": (4, 13), "Sb": (4, 14), "Te": (4, 15), "I": (4, 16), "Xe": (4, 17),
    "Cs": (5, 0), "Ba": (5, 1), "La": (5, 2), "Hf": (5, 3), "Ta": (5, 4), "W": (5, 5), "Re": (5, 6), "Os": (5, 7), "Ir": (5, 8), "Pt": (5, 9), "Au": (5, 10), "Hg": (5, 11), "Tl": (5, 12), "Pb": (5, 13), "Bi": (5, 14), "Po": (5, 15), "At": (5, 16), "Rn": (5, 17),
    "Fr": (6, 0), "Ra": (6, 1), "Ac": (6, 2), "Rf": (6, 3), "Db": (6, 4), "Sg": (6, 5), "Bh": (6, 6), "Hs": (6, 7), "Mt": (6, 8), "Ds": (6, 9), "Rg": (6, 10), "Cn": (6, 11), "Nh": (6, 12), "Fl": (6, 13), "Mc": (6, 14), "Lv": (6, 15), "Ts": (6, 16), "Og": (6, 17),
}


RESIDUE_ATOMS = {
    "ALA": ["N", "CA", "C", "O", "CB"], "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"], "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"], "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"], "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"], "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"], "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"], "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"], "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"], "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"], "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
    "A": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N9"],
    "C": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1"],
    "G": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N9"],
    "T": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1", "C5M"],
    "U": ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N1"],
}
ION_NAMES = ("LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "ZN", "FE", "CL")
_PERIODIC_COLORS = {"alkali": "#f7c6c7", "alkaline": "#f2d3a1", "transition": "#f4e3b2", "post": "#c8e6c9", "metalloid": "#b2dfdb", "nonmetal": "#bbdefb", "halogen": "#d1c4e9", "noble": "#e1bee7", "lanthanide": "#ffe0b2", "actinide": "#ffccbc"}

def _element_category(symbol: str) -> str:
    if symbol in {"H", "C", "N", "O", "P", "S", "Se"}: return "nonmetal"
    if symbol in {"F", "Cl", "Br", "I", "At", "Ts"}: return "halogen"
    if symbol in {"He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og"}: return "noble"
    if symbol in {"Li", "Na", "K", "Rb", "Cs", "Fr"}: return "alkali"
    if symbol in {"Be", "Mg", "Ca", "Sr", "Ba", "Ra"}: return "alkaline"
    if symbol in {"B", "Si", "Ge", "As", "Sb", "Te", "Po"}: return "metalloid"
    if symbol in {"Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi", "Nh", "Fl", "Mc", "Lv"}: return "post"
    if symbol in {"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"}: return "lanthanide"
    if symbol in {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"}: return "actinide"
    return "transition"

for _column, _symbol in enumerate(("Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"), 2):
    _PERIODIC_TABLE[_symbol] = (7, _column)
for _column, _symbol in enumerate(("Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"), 2):
    _PERIODIC_TABLE[_symbol] = (8, _column)
_ELEMENT_ORDER = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
_ELEMENT_NUMBERS = {symbol.upper(): index for index, symbol in enumerate(_ELEMENT_ORDER, 1)}
_ELEMENT_NAMES = {
    "H": "Hydrogen", "He": "Helium", "Li": "Lithium", "Be": "Beryllium", "B": "Boron", "C": "Carbon", "N": "Nitrogen", "O": "Oxygen", "F": "Fluorine", "Ne": "Neon",
    "Na": "Sodium", "Mg": "Magnesium", "Al": "Aluminium", "Si": "Silicon", "P": "Phosphorus", "S": "Sulfur", "Cl": "Chlorine", "Ar": "Argon",
}
_ELEMENT_MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "P": 30.974, "S": 32.06, "Cl": 35.45, "Na": 22.990, "Mg": 24.305, "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Zn": 65.38, "Br": 79.904, "I": 126.904}
_RESIDUE_NAMES = {
    "ALA": "Alanine", "ARG": "Arginine", "ASN": "Asparagine", "ASP": "Aspartic acid", "CYS": "Cysteine", "GLN": "Glutamine", "GLU": "Glutamic acid", "GLY": "Glycine", "HIS": "Histidine", "ILE": "Isoleucine", "LEU": "Leucine", "LYS": "Lysine", "MET": "Methionine", "PHE": "Phenylalanine", "PRO": "Proline", "SER": "Serine", "THR": "Threonine", "TRP": "Tryptophan", "TYR": "Tyrosine", "VAL": "Valine", "A": "Adenine", "C": "Cytosine", "G": "Guanine", "T": "Thymine", "U": "Uracil"
}

class ResidueDiagram(QWidget):
    """Clean schematic of a residue's standard biochemical connectivity."""
    def __init__(self, parent=None):
        super().__init__(parent); self.residue = ""; self.setMinimumSize(280, 250)
    def set_atoms(self, residue: str, atoms: list[str]) -> None:
        self.residue = residue; self.update()
    def paintEvent(self, event) -> None:
        painter = QPainter(self); painter.setRenderHint(QPainter.Antialiasing)
        painter.setPen(QPen(QColor("#c7d2df"), 2)); painter.drawText(10, 22, f"{_RESIDUE_NAMES.get(self.residue, self.residue)} ({self.residue})")
        def atom(x, y, text, color="#4f9fcf", radius=19):
            painter.setBrush(QBrush(QColor(color))); painter.setPen(QPen(QColor("#dce8f2"), 1)); painter.drawEllipse(x-radius, y-radius, radius*2, radius*2); painter.drawText(x-radius, y-7, radius*2, 18, Qt.AlignCenter, text)
        def bond(a, b, double=False):
            painter.setPen(QPen(QColor("#71849a"), 3)); painter.drawLine(*a, *b)
            if double: painter.drawLine(a[0], a[1]+5, b[0], b[1]+5)
        if self.residue in {"A", "C", "G", "T", "U"}:
            # Pentose ring with phosphate on the 5' side and the base attached at C1'.
            ring = {"C4'": (125, 125), "O4'": (165, 95), "C1'": (210, 115), "C2'": (205, 165), "C3'": (155, 180)}
            for a, b in (("C4'", "O4'"), ("O4'", "C1'"), ("C1'", "C2'"), ("C2'", "C3'"), ("C3'", "C4'")): bond(ring[a], ring[b])
            for name, pos in ring.items(): atom(*pos, name, "#56b893" if name == "O4'" else "#4f9fcf", 16)
            bond((125, 125), (70, 125)); atom(48, 125, "P", "#e5b94f", 17); bond((48, 125), (35, 88)); bond((48, 125), (35, 162)); atom(25, 82, "O", "#d66d6d", 13); atom(25, 168, "O", "#d66d6d", 13)
            base_name = "N9" if self.residue in {"A", "G"} else "N1"; bond(ring["C1'"], (255, 115))
            base = [(255, 115), (295, 90), (330, 115), (320, 155), (275, 165), (245, 145)]
            painter.setPen(QPen(QColor("#71849a"), 3)); painter.drawPolygon(QPolygon([QPoint(*point) for point in base])); atom(255, 115, base_name, "#d9778a", 15); painter.setPen(QPen(QColor("#a9b7c5"), 1)); painter.drawText(245, 205, 100, 18, Qt.AlignCenter, f"base {self.residue}")
        else:
            # Shared peptide backbone; the residue-specific side chain branches from Cα.
            n, ca, c, o = (50, 125), (135, 125), (220, 125), (285, 85)
            bond(n, ca); bond(ca, c); bond(c, o, True); atom(*n, "N"); atom(*ca, "Cα", "#62a9d8"); atom(*c, "C"); atom(*o, "O", "#d66d6d")
            side = [name for name in RESIDUE_ATOMS.get(self.residue, []) if name not in {"N", "CA", "C", "O"}]
            side = side[:10]
            points = {}
            for index, name in enumerate(side):
                points[name] = (135 + 45 * (index % 3), 185 + 42 * (index // 3))
                parent = ca if index == 0 else points[side[index - 1]]; bond(parent, points[name])
            for name, pos in points.items(): atom(*pos, name, "#a88bd8" if name.startswith("O") or name.startswith("N") else "#4f9fcf", 15)
            painter.setPen(QPen(QColor("#a9b7c5"), 1)); painter.drawText(85, 245, 180, 18, Qt.AlignCenter, "residue-specific side chain")


class AtomicRadiiDialog(QDialog):
    """Three-tab editor for element, residue-atom, and ion radii."""
    def __init__(self, radii: dict[str, float], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Atomic radii")
        self.setMinimumSize(1150, 560)
        self._radii = {key.upper(): float(value) for key, value in radii.items()}
        self._fields: dict[str, QDoubleSpinBox] = {}
        self._element_values: dict[str, float] = {}
        self._baseline: dict[str, float] = {}
        self._residue_fields: dict[str, QDoubleSpinBox] = {}
        self._ion_fields: dict[str, QDoubleSpinBox] = {}
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Set radii in Å. Residue-atom overrides take precedence over element defaults."))
        tabs = QTabWidget()
        tabs.addTab(self._build_defaults_tab(), "Defaults")
        tabs.addTab(self._build_residue_tab(), "Residues")
        tabs.addTab(self._build_ions_tab(), "Ions")
        layout.addWidget(tabs, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept); buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _spin(self, key: str, value: float) -> QDoubleSpinBox:
        field = QDoubleSpinBox(); field.setRange(0.01, 5.0); field.setDecimals(3); field.setSingleStep(0.01); field.setValue(value)
        self._baseline[key] = value
        return field

    def _build_defaults_tab(self) -> QWidget:
        panel = QWidget(); layout = QVBoxLayout(panel)
        legend = QLabel("Click an element to edit its display radius. Colors identify periodic-table families."); legend.setWordWrap(True); layout.addWidget(legend)
        grid = QGridLayout(); grid.setSpacing(3)
        for symbol, (row, column) in _PERIODIC_TABLE.items():
            card = QFrame(); card.setFixedSize(62, 78); card.setCursor(Qt.PointingHandCursor); card.setStyleSheet(f"QFrame {{ background: {_PERIODIC_COLORS[_element_category(symbol)]}; border: none; border-radius: 3px; }} QLabel {{ color: #17202a; }}")
            card_layout = QVBoxLayout(card); card_layout.setContentsMargins(2, 2, 2, 2); card_layout.setSpacing(0)
            number = _ELEMENT_NUMBERS[symbol.upper()]; mass = _ELEMENT_MASSES.get(symbol, float(number * 2))
            top = QLabel(f"{number}  {mass:g}"); top.setAlignment(Qt.AlignCenter); top.setStyleSheet("font-size: 7px;"); card_layout.addWidget(top)
            name = QLabel(_ELEMENT_NAMES.get(symbol, symbol)); name.setAlignment(Qt.AlignCenter); name.setWordWrap(True); name.setFixedHeight(19); name.setStyleSheet("font-size: 7px;"); card_layout.addWidget(name)
            abbrev = QLabel(symbol); abbrev.setAlignment(Qt.AlignCenter); abbrev.setStyleSheet("font-size: 16px; font-weight: 700;"); card_layout.addWidget(abbrev, 1)
            radius = QLabel(f"{self._radii.get(symbol.upper(), 0.36):.2f} Å"); radius.setAlignment(Qt.AlignCenter); radius.setStyleSheet("font-size: 8px;"); card_layout.addWidget(radius)
            card.mousePressEvent = lambda event, element=symbol, label=radius: self._edit_element(element, label)
            grid.addWidget(card, row, column)
        layout.addLayout(grid); return panel

    def _edit_element(self, symbol: str, label: QLabel) -> None:
        key = symbol.upper(); current = self._element_values.get(key, self._radii.get(key, 0.36))
        value, accepted = QInputDialog.getDouble(self, f"{_ELEMENT_NAMES.get(symbol, symbol)} radius", "Radius (Å)", current, 0.01, 5.0, 3)
        if accepted:
            self._element_values[key] = value; label.setText(f"{value:.2f} Å")

    def _build_residue_tab(self) -> QWidget:
        panel = QWidget(); layout = QHBoxLayout(panel)
        names = QListWidget(); names.setMaximumWidth(180)
        for abbreviation in RESIDUE_ATOMS:
            item = QListWidgetItem(f"{_RESIDUE_NAMES.get(abbreviation, abbreviation)} ({abbreviation})"); item.setData(Qt.UserRole, abbreviation); names.addItem(item)
        diagram = ResidueDiagram(); table = QTableWidget(0, 2); table.setHorizontalHeaderLabels(["Atom name", "Radius (Å)"]); table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(names); layout.addWidget(diagram, 1); layout.addWidget(table, 1)
        self._residue_table, self._residue_diagram = table, diagram
        names.currentItemChanged.connect(lambda item, _previous: self._select_residue_editor(item.data(Qt.UserRole) if item else "")); names.setCurrentRow(0)
        return panel

    def _select_residue_editor(self, residue: str) -> None:
        atoms = RESIDUE_ATOMS.get(residue, []); self._residue_table.setRowCount(len(atoms)); self._residue_diagram.set_atoms(residue, atoms)
        for row, atom_name in enumerate(atoms):
            key = f"RES:{residue}:{atom_name}"; self._residue_table.setItem(row, 0, QTableWidgetItem(atom_name))
            field = self._spin(key, self._radii.get(key, self._radii.get(atom_name[:2].upper(), self._radii.get(atom_name[:1].upper(), 0.36)))); self._residue_table.setCellWidget(row, 1, field); self._residue_fields[key] = field

    def _build_ions_tab(self) -> QWidget:
        panel = QWidget(); layout = QVBoxLayout(panel); table = QTableWidget(len(ION_NAMES), 2); table.setHorizontalHeaderLabels(["Ion", "Radius (Å)"]); table.horizontalHeader().setStretchLastSection(True)
        for row, ion in enumerate(ION_NAMES):
            table.setItem(row, 0, QTableWidgetItem(ion)); key = f"ION:{ion}"; field = self._spin(key, self._radii.get(key, self._radii.get(ion, 0.4))); table.setCellWidget(row, 1, field); self._ion_fields[key] = field
        layout.addWidget(table); return panel

    def values(self) -> dict[str, float]:
        values = dict(self._element_values)
        values.update({key: field.value() for key, field in {**self._residue_fields, **self._ion_fields}.items() if abs(field.value() - self._baseline[key]) > 1e-9})
        return values


class MetricCard(QFrame):
    """Compact numerical readout used in the analysis tray."""

    def __init__(self, title: str, accent: str, parent=None):
        super().__init__(parent)
        self.setObjectName("metricCard")
        self.setStyleSheet(f"QFrame#metricCard {{ border-top-color: {accent}; }}")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 9, 12, 9)
        layout.setSpacing(2)
        label = QLabel(title)
        label.setObjectName("metricTitle")
        self.value = QLabel("—")
        self.value.setObjectName("metricValue")
        layout.addWidget(label)
        layout.addWidget(self.value)


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("VorPy")
        self.resize(1480, 900)
        self.setMinimumSize(1050, 680)
        self.backend = VorPyBackend()
        self.project = Project()
        self.radius_overrides: dict[str, float] = {}
        self.project_file: Path | None = None
        self._project_dirty = False
        self._loading_project = False
        self.source: Path | None = None
        self.current_result: AnalysisResult | None = None
        self._loaded_results: dict[Path, AnalysisResult] = {}
        self._structure_states: dict[Path, tuple[set[int], dict[str, tuple[int, ...]], dict[str, tuple[str, str]]]] = {}
        self._thread: QThread | None = None
        self._worker: SolveWorker | None = None
        self._running_selection: set[int] = set()
        self._groups: dict[str, tuple[int, ...]] = {}
        self._interfaces: dict[str, tuple[str, str]] = {}
        self._selection_entries: list[tuple[str, tuple[int, ...]]] = []

        self.viewer = MolecularView(self)
        self.viewer.selected_atom.connect(self._show_selected_atom)
        self.viewer.selected_residue.connect(self._show_selected_residue)
        self.viewer.selected_chain.connect(self._show_selected_chain)
        self.viewer.selected_molecule.connect(self._show_selected_molecule)
        self.viewer.selection_cleared.connect(self._clear_selection_from_viewer)
        self._build_actions()
        self._build_menu_and_toolbar()
        self._build_workspace()
        self._build_status()
        self.statusBar().showMessage("Ready — load a structure and select a residue")
        self.select_residue_action.setChecked(True)

    def _build_actions(self) -> None:
        style = self.style()
        self.new_project_action = QAction("New project", self)
        self.new_project_action.setShortcut("Ctrl+N")
        self.new_project_action.triggered.connect(self.new_project)
        self.open_project_action = QAction("Open project…", self)
        self.open_project_action.setShortcut("Ctrl+Shift+P")
        self.open_project_action.triggered.connect(self.open_project)
        self.save_project_action = QAction("Save project", self)
        self.save_project_action.setShortcut("Ctrl+S")
        self.save_project_action.triggered.connect(self.save_project)
        self.save_project_as_action = QAction("Save project as…", self)
        self.save_project_as_action.setShortcut("Ctrl+Shift+S")
        self.save_project_as_action.triggered.connect(self.save_project_as)
        self.open_action = QAction(
            style.standardIcon(QStyle.SP_DialogOpenButton), "Open PDB", self
        )
        self.open_action.setShortcut("Ctrl+O")
        self.open_action.triggered.connect(self.open_structure)
        self.open_result_action = QAction(
            style.standardIcon(QStyle.SP_DirOpenIcon), "Open VorPy output", self
        )
        self.open_result_action.setShortcut("Ctrl+Shift+O")
        self.open_result_action.triggered.connect(self.open_result_directory)
        self.solve_action = QAction(
            style.standardIcon(QStyle.SP_MediaPlay), "Run VorPy analysis", self
        )
        self.solve_action.setShortcut("Ctrl+R")
        self.solve_action.triggered.connect(self.solve)
        self.cancel_action = QAction(
            style.standardIcon(QStyle.SP_MediaStop), "Cancel", self
        )
        self.cancel_action.setEnabled(False)
        self.cancel_action.triggered.connect(self.cancel_analysis)
        self.fit_action = QAction(
            style.standardIcon(QStyle.SP_BrowserReload), "Fit view", self
        )
        self.fit_action.setShortcut("F")
        self.fit_action.triggered.connect(self._fit_view)
        self.screenshot_action = QAction(
            style.standardIcon(QStyle.SP_DialogSaveButton), "Screenshot", self
        )
        self.screenshot_action.triggered.connect(self.save_screenshot)

        self.selection_actions = QActionGroup(self)
        self.selection_actions.setExclusive(True)
        self.select_atom_action = QAction("Select atom", self, checkable=True)
        self.select_atom_action.setShortcut("A")
        self.select_residue_action = QAction("Select residue", self, checkable=True)
        self.select_chain_action = QAction("Select chain", self, checkable=True)
        self.select_chain_action.setShortcut("C")
        self.select_molecule_action = QAction("Select molecule", self, checkable=True)
        self.select_molecule_action.setShortcut("M")
        self.selection_actions.addAction(self.select_atom_action)
        self.selection_actions.addAction(self.select_residue_action)
        self.selection_actions.addAction(self.select_chain_action)
        self.selection_actions.addAction(self.select_molecule_action)
        self.select_atom_action.toggled.connect(
            lambda checked: self._set_selection_mode("atom" if checked else None)
        )
        self.select_residue_action.toggled.connect(
            lambda checked: self._set_selection_mode("residue" if checked else None)
        )
        self.select_chain_action.toggled.connect(
            lambda checked: self._set_selection_mode("chain" if checked else None)
        )
        self.select_molecule_action.toggled.connect(
            lambda checked: self._set_selection_mode("molecule" if checked else None)
        )
        self.reset_selection_action = QAction("Reset selection…", self)
        self.reset_selection_action.setShortcut("R")
        self.reset_selection_action.triggered.connect(self.reset_selection)

    def _build_menu_and_toolbar(self) -> None:
        file_menu = self.menuBar().addMenu("File")
        file_menu.addActions([self.new_project_action, self.open_project_action])
        file_menu.addActions([self.save_project_action, self.save_project_as_action])
        file_menu.addSeparator()
        file_menu.addActions([self.open_action, self.open_result_action])
        file_menu.addSeparator()
        file_menu.addAction(self.screenshot_action)
        file_menu.addSeparator()
        file_menu.addAction("Exit", self.close)
        selection_menu = self.menuBar().addMenu("Selection")
        selection_menu.addActions(
            [
                self.select_atom_action,
                self.select_residue_action,
                self.select_chain_action,
                self.select_molecule_action,
            ]
        )
        selection_menu.addSeparator()
        selection_menu.addAction(self.reset_selection_action)
        self.menuBar().addMenu("Analysis").addAction(self.solve_action)
        self.menuBar().addMenu("View").addAction(self.fit_action)
        self.menuBar().addMenu("Help")

    def _build_workspace(self) -> None:
        root = QWidget()
        layout = QVBoxLayout(root)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        workflow = QWidget()
        self.workflow_panel = workflow
        workflow.setObjectName("workflowPanel")
        workflow.setMinimumWidth(350)
        workflow.setMaximumWidth(430)
        workflow_layout = QVBoxLayout(workflow)
        workflow_layout.setContentsMargins(8, 8, 5, 5)
        workflow_layout.setSpacing(8)
        self.workflow_tabs = QTabWidget()
        self.workflow_tabs.addTab(self._build_structure_tab(), "Structure")
        self.workflow_tabs.addTab(self._build_selection_tab(), "Selection")
        self.workflow_tabs.addTab(self._build_groups_tab(), "Groups")
        self.workflow_tabs.addTab(self._build_interfaces_tab(), "Interfaces")
        workflow_layout.addWidget(self.workflow_tabs, 1)
        vertical = QSplitter(Qt.Vertical)
        self.workspace_splitter = vertical
        upper = QSplitter(Qt.Horizontal)
        self.upper_workspace = upper
        upper.addWidget(workflow)
        upper.addWidget(self._build_viewport())
        upper.addWidget(self._build_inspector())
        upper.setStretchFactor(0, 0)
        upper.setStretchFactor(1, 5)
        upper.setStretchFactor(2, 1)
        upper.setSizes([390, 1050, 320])
        vertical.addWidget(upper)
        self.analysis_tray = self._build_analysis_tray()
        vertical.addWidget(self.analysis_tray)
        vertical.setStretchFactor(0, 4)
        vertical.setStretchFactor(1, 2)
        vertical.setSizes([560, 320])
        layout.addWidget(vertical)
        self.setCentralWidget(root)

    def _build_viewport(self) -> QWidget:
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(5, 5, 5, 0)
        layout.setSpacing(5)
        layout.addWidget(self.viewer, 1)
        bar = QFrame()
        bar.setObjectName("viewerBar")
        bar_layout = QHBoxLayout(bar)
        bar_layout.setContentsMargins(12, 6, 12, 6)
        bar_layout.setSpacing(14)
        bar_layout.addStretch()
        self.selection_mode_label = QLabel("Navigation mode")
        self.selection_mode_label.setObjectName("sectionLabel")
        bar_layout.addWidget(self.selection_mode_label)
        layout.addWidget(bar)
        return container

    @staticmethod
    def _visibility_checkbox(label: str, callback, checked: bool = True) -> QCheckBox:
        checkbox = QCheckBox(label)
        checkbox.setChecked(checked)
        checkbox.toggled.connect(callback)
        return checkbox

    def _build_inspector(self) -> QWidget:
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        title = QLabel("View Settings")
        title.setObjectName("viewSettingsTitle")
        layout.addWidget(title)
        self.inspector = QTabWidget()
        self.inspector.setMinimumWidth(290)
        self.inspector.setMaximumWidth(390)
        self.inspector.addTab(self._build_visualization_tab(), "System")
        self.inspector.addTab(self._build_network_tab(), "Network")
        layout.addWidget(self.inspector, 1)
        return container

    def _build_visualization_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        molecule_group = QGroupBox("Molecule")
        molecule_form = QFormLayout(molecule_group)
        self.show_cartoon = self._visibility_checkbox(
            "Cartoon", self.viewer.set_cartoon_visible
        )
        self.show_spheres = self._visibility_checkbox(
            "Spheres",
            lambda checked: self.viewer.set_category_visible("atoms", checked),
            checked=False,
        )
        self.show_sticks = self._visibility_checkbox(
            "Sticks", self.viewer.set_bonds_visible, checked=False
        )
        # Compatibility aliases for code using the former visibility names.
        self.show_atoms = self.show_spheres
        self.show_bonds = self.show_sticks
        representations = QWidget()
        representation_layout = QVBoxLayout(representations)
        representation_layout.setContentsMargins(0, 0, 0, 0)
        for checkbox in (self.show_cartoon, self.show_spheres, self.show_sticks):
            representation_layout.addWidget(checkbox)
        self.molecule_opacity = QSlider(Qt.Horizontal)
        self.molecule_opacity.setRange(5, 100)
        self.molecule_opacity.setValue(50)
        self.molecule_opacity.valueChanged.connect(
            lambda value: self.viewer.set_molecule_opacity(value / 100.0)
        )
        molecule_form.addRow("Representations", representations)
        molecule_form.addRow("Sphere opacity", self.molecule_opacity)
        layout.addWidget(molecule_group)

        water_group = QGroupBox("Water")
        water_form = QFormLayout(water_group)
        self.show_waters = self._visibility_checkbox(
            "Show waters", self.viewer.set_waters_visible, checked=False
        )
        self.water_style = QComboBox()
        self.water_style.addItem("Ball and stick", "ball-and-stick")
        self.water_style.addItem("Sticks only", "sticks")
        self.water_style.addItem("Translucent spheres", "spheres")
        self.water_style.currentIndexChanged.connect(
            lambda: self.viewer.set_water_style(self.water_style.currentData())
        )
        self.water_opacity = QSlider(Qt.Horizontal)
        self.water_opacity.setRange(5, 80)
        self.water_opacity.setValue(25)
        self.water_opacity.valueChanged.connect(
            lambda value: self.viewer.set_water_opacity(value / 100.0)
        )
        water_form.addRow(self.show_waters)
        water_form.addRow("Style", self.water_style)
        water_form.addRow("Sphere opacity", self.water_opacity)
        layout.addWidget(water_group)

        ion_group = QGroupBox("Ions")
        ion_layout = QVBoxLayout(ion_group)
        self.show_ions = self._visibility_checkbox(
            "Show ions", self.viewer.set_ions_visible, checked=False
        )
        ion_layout.addWidget(self.show_ions)
        layout.addWidget(ion_group)
        layout.addStretch()
        return panel

    def _build_structure_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.structure_name = QLabel("No structure loaded")
        self.structure_name.setObjectName("metricValue")
        self.structure_name.setAlignment(Qt.AlignCenter)
        self.structure_summary = QLabel("Browse to load a molecular structure")
        self.structure_summary.setWordWrap(True)
        self.structure_summary.setAlignment(Qt.AlignCenter)
        self.structure_summary.setObjectName("sectionLabel")
        self.structure_summary.setVisible(False)
        layout.addStretch(1)
        layout.addWidget(self.structure_name)
        layout.addStretch(1)
        self.chemical_info = QGroupBox("Chemical information")
        chemical_form = QFormLayout(self.chemical_info)
        self.info_atoms = QLabel("—")
        self.info_bonds = QLabel("—")
        self.info_residues = QLabel("—")
        self.info_chains = QLabel("—")
        self.info_molecules = QLabel("—")
        self.info_waters = QLabel("—")
        self.info_ions = QLabel("—")
        for label, widget in (("Atoms", self.info_atoms), ("Bonds", self.info_bonds),
                              ("Residues", self.info_residues), ("Chains", self.info_chains),
                              ("Molecules", self.info_molecules), ("Waters", self.info_waters),
                              ("Ions", self.info_ions)):
            widget.setObjectName("sectionLabel")
            chemical_form.addRow(label, widget)
        self.chemical_info.setVisible(False)
        layout.addWidget(self.chemical_info)
        self.loaded_structures = QListWidget()
        self.loaded_structures.setAlternatingRowColors(True)
        self.loaded_structures.itemClicked.connect(self._switch_structure_item)
        self.loaded_structures.setToolTip("Loaded structures; click one to switch the rendered scene")
        self.loaded_structures_label = QLabel("Loaded structures")
        self.loaded_structures_label.setVisible(False)
        self.loaded_structures.setVisible(False)
        layout.addWidget(self.loaded_structures_label)
        layout.addWidget(self.loaded_structures, 2)
        structure_buttons = QHBoxLayout()
        self.add_structure_button = QPushButton("Add structure")
        self.delete_structure_button = QPushButton("Delete")
        self.reset_structures_button = QPushButton("Reset all")
        self.add_structure_button.clicked.connect(self.open_structure)
        self.delete_structure_button.clicked.connect(self._delete_current_structure)
        self.reset_structures_button.clicked.connect(self._reset_structures)
        for button in (self.add_structure_button, self.delete_structure_button, self.reset_structures_button):
            structure_buttons.addWidget(button)
        layout.addLayout(structure_buttons)
        return panel

    def _build_selection_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.selection_group = QGroupBox("No selection")
        form = QFormLayout(self.selection_group)
        self.atom_name = QLabel("—")
        self.atom_element = QLabel("—")
        self.atom_residue = QLabel("—")
        self.atom_chain = QLabel("—")
        self.atom_position = QLabel("—")
        self.selection_count = QLabel("—")
        form.addRow("Atom", self.atom_name)
        form.addRow("Element", self.atom_element)
        form.addRow("Residue", self.atom_residue)
        form.addRow("Chain", self.atom_chain)
        form.addRow("Coordinates (Å)", self.atom_position)
        form.addRow("Atoms selected", self.selection_count)
        self.selection_group.setVisible(False)
        layout.addWidget(self.selection_group)
        self.running_selection = QLabel("No atoms selected")
        self.running_selection.setObjectName("sectionLabel")
        self.running_selection.setWordWrap(True)
        self.running_selection.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.running_selection)
        self.structure_search = QLineEdit()
        self.structure_search.setPlaceholderText("Search residue, chain, molecule, or atom…")
        self.structure_search.setClearButtonEnabled(True)
        self.structure_search.textChanged.connect(self._filter_structure_browser)
        layout.addWidget(self.structure_search)
        self.structure_browser_type = QComboBox()
        self.structure_browser_type.addItems(["Atoms", "Residues", "Chains", "Molecules"])
        self.structure_browser_type.currentIndexChanged.connect(self._populate_structure_browser)
        layout.addWidget(self.structure_browser_type)
        self.structure_browser = QListWidget()
        self.structure_browser.setAlternatingRowColors(True)
        self.structure_browser.itemChanged.connect(self._browser_item_changed)
        layout.addWidget(self.structure_browser, 1)
        hint = QLabel("Select atoms, residues, chains, or molecules, then save the active selection as a group.")
        hint.setWordWrap(True)
        hint.setObjectName("sectionLabel")
        layout.addWidget(hint)
        self.make_group_button = QPushButton("Make group")
        self.make_group_button.setEnabled(False)
        self.make_group_button.clicked.connect(self._make_group)
        layout.addWidget(self.make_group_button)
        return panel

    def _build_groups_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        hint = QLabel("Saved groups can be selected to highlight their atoms and used as solve targets.")
        hint.setWordWrap(True)
        hint.setObjectName("sectionLabel")
        layout.addWidget(hint)
        self.groups_list = QListWidget()
        self.groups_list.setAlternatingRowColors(True)
        self.groups_list.itemClicked.connect(self._select_saved_group)
        layout.addWidget(self.groups_list, 1)
        return panel

    def _build_interfaces_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        hint = QLabel("Define an interface between two saved groups for targeted analysis.")
        hint.setWordWrap(True)
        hint.setObjectName("sectionLabel")
        layout.addWidget(hint)
        interface_box = QGroupBox("New interface")
        interface_layout = QFormLayout(interface_box)
        self.interface_group_a = QComboBox()
        self.interface_group_b = QComboBox()
        self.make_interface_button = QPushButton("Make interface")
        self.make_interface_button.clicked.connect(self._make_interface)
        interface_layout.addRow("Group A", self.interface_group_a)
        interface_layout.addRow("Group B", self.interface_group_b)
        interface_layout.addRow(self.make_interface_button)
        layout.addWidget(interface_box)
        self.interfaces_list = QListWidget()
        self.interfaces_list.setAlternatingRowColors(True)
        layout.addWidget(self.interfaces_list, 1)
        return panel

    def _build_solve_section(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        section = QGroupBox("Solve configuration")
        form = QFormLayout(section)
        self.solve_target = QComboBox()
        self.solve_target.addItem("Whole molecule", ("whole", ""))
        self.solve_target.addItem("Active selection", ("selection", ""))
        self.network_type = QComboBox()
        self.network_type.addItems(["Atomic Voronoi", "Power", "Primitive"])
        self.max_vertices = QSpinBox()
        self.max_vertices.setRange(1, 10000)
        self.max_vertices.setValue(40)
        self.box_size = QDoubleSpinBox()
        self.box_size.setRange(1.0, 10.0)
        self.box_size.setSingleStep(0.05)
        self.box_size.setValue(1.25)
        self.surface_resolution = QDoubleSpinBox()
        self.surface_resolution.setRange(0.01, 1.0)
        self.surface_resolution.setSingleStep(0.01)
        self.surface_resolution.setValue(0.2)
        form.addRow("Target", self.solve_target)
        form.addRow("Network", self.network_type)
        max_row = QWidget()
        max_layout = QHBoxLayout(max_row)
        max_layout.setContentsMargins(0, 0, 0, 0)
        max_layout.addWidget(self.max_vertices)
        max_help = QPushButton("?")
        max_help.setFixedWidth(26)
        max_help.setToolTip("Maximum vertices per Voronoi cell. Lower values simplify the network and speed up surface construction.")
        max_layout.addWidget(max_help)
        form.addRow("Max Vert rad", max_row)
        form.addRow("Box size", self.box_size)
        form.addRow("Surface resolution", self.surface_resolution)
        help_button = QPushButton("?  Build setting help")
        help_button.setToolTip("Explain the solve configuration settings")
        help_button.clicked.connect(self._show_build_settings_help)
        form.addRow(help_button)
        layout.addWidget(section)
        self.adjust_radii_button = QPushButton("Adjust atomic radii…")
        self.adjust_radii_button.clicked.connect(self._open_atomic_radii)
        layout.addWidget(self.adjust_radii_button)
        layout.addStretch(1)
        self.solve_network_button = QPushButton("Solve network")
        self.solve_network_button.setObjectName("primaryAction")
        self.solve_network_button.setMinimumHeight(36)
        self.solve_network_button.clicked.connect(self.solve_action.trigger)
        self.solve_action.changed.connect(
            lambda: self.solve_network_button.setEnabled(
                self.solve_action.isEnabled()
            )
        )
        layout.addWidget(self.solve_network_button)
        return panel

    def _show_build_settings_help(self) -> None:
        QMessageBox.information(self, "Build settings",
            "Target chooses the whole structure, active selection, group, or interface.\n\n"
            "Network selects the Voronoi construction. Max Vert rad limits the number of vertices used to represent each cell.\n\n"
            "Box size controls the padding around the group. Surface resolution controls triangle spacing; smaller values create finer surfaces but take longer.")

    def _open_atomic_radii(self) -> None:
        current = dict(DISPLAY_RADII)
        if self.current_result is not None:
            current.update({atom.element.upper(): atom.radius for atom in self.current_result.atoms})
        current.update(self.radius_overrides)
        dialog = AtomicRadiiDialog(current, self)
        if dialog.exec() != QDialog.Accepted:
            return
        self.radius_overrides = dialog.values()
        if self.current_result is not None:
            self.current_result.atoms = [
                replace(atom, radius=self._radius_for_atom(atom))
                for atom in self.current_result.atoms
            ]
            self.viewer.display_result(self.current_result)
            self._set_project_dirty()
        self.statusBar().showMessage("Atomic radii updated")

    def _update_solve_targets(self) -> None:
        current = self.solve_target.currentData()
        self.solve_target.blockSignals(True)
        while self.solve_target.count() > 2:
            self.solve_target.removeItem(2)
        for name in self._groups:
            self.solve_target.addItem(f"Group: {name}", ("group", name))
        for name in self._interfaces:
            self.solve_target.addItem(f"Interface: {name}", ("interface", name))
        if current is not None:
            index = self.solve_target.findData(current)
            self.solve_target.setCurrentIndex(max(0, index))
        self.solve_target.blockSignals(False)

    def _build_network_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.network_layer_checks = {}
        self.network_color_buttons = {}
        self.network_colors = {}

        # Compact three-column layer controls:
        # Show | Object | Color
        layer_table = QWidget()
        layer_grid = QGridLayout(layer_table)
        layer_grid.setContentsMargins(0, 0, 0, 0)
        layer_grid.setHorizontalSpacing(8)
        layer_grid.setVerticalSpacing(2)

        show_header = QLabel("Show")
        object_header = QLabel("Object")
        color_header = QLabel("Color")
        for header in (show_header, object_header, color_header):
            header.setObjectName("sectionLabel")

        layer_grid.addWidget(show_header, 0, 0, Qt.AlignmentFlag.AlignCenter)
        layer_grid.addWidget(object_header, 0, 1, Qt.AlignmentFlag.AlignLeft)
        layer_grid.addWidget(color_header, 0, 2, Qt.AlignmentFlag.AlignCenter)

        for row_index, (key, label, color) in enumerate(
            NETWORK_LAYER_OPTIONS, start=1
        ):
            checkbox = QCheckBox()
            checkbox.setEnabled(False)
            checkbox.setToolTip(f"Show {label.lower()}")
            checkbox.toggled.connect(
                lambda visible, layer_key=key: self._set_network_layer_visible(
                    layer_key, visible
                )
            )
            self.network_layer_checks[key] = checkbox
            self.network_colors[key] = color

            object_label = QLabel(label)
            object_label.setMinimumHeight(22)

            color_button = QPushButton()
            color_button.setFixedSize(22, 22)
            color_button.setEnabled(False)
            color_button.setToolTip(f"Choose solid color for {label.lower()}")
            color_button.clicked.connect(
                lambda _checked=False, layer_key=key: self._choose_network_color(
                    layer_key
                )
            )
            self._set_color_button_swatch(color_button, color)
            self.network_color_buttons[key] = color_button

            layer_grid.addWidget(
                checkbox,
                row_index,
                0,
                Qt.AlignmentFlag.AlignCenter,
            )
            layer_grid.addWidget(object_label, row_index, 1)
            layer_grid.addWidget(
                color_button,
                row_index,
                2,
                Qt.AlignmentFlag.AlignCenter,
            )

        layer_grid.setColumnStretch(0, 0)
        layer_grid.setColumnStretch(1, 1)
        layer_grid.setColumnStretch(2, 0)
        layout.addWidget(layer_table)

        surfaces = QGroupBox("Geometry appearance")
        surface_form = QFormLayout(surfaces)

        self.surface_color_scheme = QComboBox()
        for label, scheme in (
            ("Solid color", "solid"),
            ("Mean curvature", "mean_curvature"),
            ("Gaussian curvature", "gaussian_curvature"),
            ("Integrated mean curvature", "integrated_mean_curvature"),
            ("Integrated Gaussian curvature", "integrated_gaussian_curvature"),
            ("Surface energy", "surface_energy"),
            ("Distance to center", "distance"),
            ("Inside / outside", "inside_outside"),
        ):
            self.surface_color_scheme.addItem(label, scheme)
        self.surface_color_scheme.currentIndexChanged.connect(
            self._set_surface_color_scheme
        )
        surface_form.addRow("Color by", self.surface_color_scheme)

        self.curvature_colormap = QComboBox()
        for label, value in (
            ("Coolwarm", "coolwarm"),
            ("RdBu", "RdBu_r"),
            ("Blue / white / red", "bwr"),
            ("Seismic", "seismic"),
            ("Spectral", "Spectral_r"),
            ("Purple / orange", "PuOr"),
            ("Purple / green", "PRGn"),
            ("Brown / teal", "BrBG"),
            ("Pink / green", "PiYG"),
            ("Viridis", "viridis"),
            ("Cividis", "cividis"),
            ("Plasma", "plasma"),
            ("Inferno", "inferno"),
            ("Magma", "magma"),
            ("Turbo", "turbo"),
            ("Cubehelix", "cubehelix"),
        ):
            self.curvature_colormap.addItem(label, value)
        self.curvature_colormap.currentIndexChanged.connect(
            self._set_curvature_colormap
        )
        surface_form.addRow("Colormap", self.curvature_colormap)

        # Discrete signed-log contrast slider. The current historical
        # Signed log setting corresponds to Log1000 (position 3).
        self.curvature_scale = QSlider(Qt.Horizontal)
        self.curvature_scale.setRange(0, 5)
        self.curvature_scale.setSingleStep(1)
        self.curvature_scale.setPageStep(1)
        self.curvature_scale.setTickPosition(QSlider.TicksBelow)
        self.curvature_scale.setTickInterval(1)
        self.curvature_scale.setValue(3)
        self.curvature_scale.valueChanged.connect(self._set_curvature_scale)

        self.curvature_scale_label = QLabel("Log1000")
        self.curvature_scale_label.setMinimumWidth(68)
        scale_row = QWidget()
        scale_layout = QHBoxLayout(scale_row)
        scale_layout.setContentsMargins(0, 0, 0, 0)
        scale_layout.setSpacing(6)
        scale_layout.addWidget(self.curvature_scale, 1)
        scale_layout.addWidget(self.curvature_scale_label)
        surface_form.addRow("Scale", scale_row)

        self.edge_size = QSlider(Qt.Horizontal)
        self.edge_size.setRange(5, 100)
        self.edge_size.setValue(20)
        self.edge_size.setToolTip("Displayed VorPy edge radius (0.005–0.100 Å)")
        self.edge_size.valueChanged.connect(self._set_edge_size)
        surface_form.addRow("Edge size", self.edge_size)

        self.vertex_size = QSlider(Qt.Horizontal)
        self.vertex_size.setRange(5, 30)
        self.vertex_size.setValue(13)
        self.vertex_size.setToolTip("Displayed VorPy vertex radius (0.05–0.30 Å)")
        self.vertex_size.valueChanged.connect(self._set_vertex_size)
        surface_form.addRow("Vertex size", self.vertex_size)

        self.surface_opacity = QSlider(Qt.Horizontal)
        self.surface_opacity.setRange(0, 100)
        self.surface_opacity.setValue(45)
        self.surface_opacity.valueChanged.connect(self._set_surface_opacity)
        surface_form.addRow("Surface opacity", self.surface_opacity)

        self.surface_color_scheme.setEnabled(False)
        self.curvature_colormap.setEnabled(False)
        self.curvature_scale.setEnabled(False)
        self.edge_size.setEnabled(False)
        self.vertex_size.setEnabled(False)
        self.surface_opacity.setEnabled(False)

        layout.addWidget(surfaces)
        layout.addStretch()
        return panel

    def _build_analysis_tray(self) -> QWidget:
        tray = QSplitter(Qt.Horizontal)
        tray.setObjectName("bottomTray")
        solve = self._build_solve_section()
        solve.setMinimumWidth(350)
        tray.addWidget(solve)
        tabs = QTabWidget()
        self.analysis_tray_tabs = tabs
        overview = QWidget()
        overview_layout = QVBoxLayout(overview)
        cards = QWidget()
        cards_layout = QHBoxLayout(cards)
        cards_layout.setContentsMargins(8, 8, 8, 4)
        self.metric_cards = {}
        for key, title, accent in (
            ("atoms", "Atoms", "#6857d9"), ("residues", "Residues", "#5a9bd5"),
            ("chains", "Chains", "#4cbe9b"), ("vertices", "Vertices", "#d99a32"),
            ("edges", "Edges", "#367bd6"), ("surfaces", "Surfaces", "#9b59db"),
            ("cells", "Complete cells", "#38a873"), ("layers", "Geometry layers", "#d99a32"),
        ):
            card = MetricCard(title, accent)
            self.metric_cards[key] = card
            cards_layout.addWidget(card)
        overview_layout.addWidget(cards)
        self.summary_statement = QLabel("Load a completed network to see its scientific summary.")
        self.summary_statement.setWordWrap(True)
        self.summary_statement.setObjectName("sectionLabel")
        overview_layout.addWidget(self.summary_statement)
        self.summary_status = QLabel("—")
        self.summary_status.setWordWrap(True)
        self.summary_status.setObjectName("sectionLabel")
        overview_layout.addWidget(self.summary_status)
        self.results = QTableWidget(0, 2)
        self.results.setAlternatingRowColors(True)
        self.results.setHorizontalHeaderLabels(["Metric", "Value"])
        self.results.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.results.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        overview_layout.addWidget(self.results, 1)
        tabs.addTab(overview, "Overview")
        self.analysis_section_tables = {}
        for section in ("Composition", "Build Information", "Build Timing",
                         "Voronoi Network", "Group Geometry", "Surface Curvature",
                         "Surface Energy Estimate", "Surface Classification",
                         "Chain Composition", "Residue Composition"):
            columns = (["Chain", "Atoms", "Residues", "Volume", "Boundary SA", "Inter-chain SA", "Solvent-interfacial SA"] if section == "Chain Composition" else ["Chain", "Residue", "Identifier", "Atoms", "Volume", "Boundary SA", "Inter-residue SA", "Solvent-interfacial SA"] if section == "Residue Composition" else ["Value", "Details"])
            table = QTableWidget(0, len(columns))
            table.setAlternatingRowColors(True)
            table.setSortingEnabled(True)
            table.setHorizontalHeaderLabels(columns)
            table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
            for column in range(1, len(columns)): table.horizontalHeader().setSectionResizeMode(column, QHeaderView.Stretch)
            self.analysis_section_tables[section] = table
            tabs.addTab(table, section.replace(" Information", "").replace(" Estimate", ""))
        tray.addWidget(tabs)
        tray.setStretchFactor(0, 0)
        tray.setStretchFactor(1, 1)
        tray.setSizes([430, 1050])
        return tray

    def _fallback_info_sections(self, result: AnalysisResult) -> dict[str, list[tuple[str, str]]]:
        return {
            "Composition": [("Atoms", f"{len(result.atoms):,}"),
                            ("Residues", f"{len(self._residue_groups(result)):,}"),
                            ("Chains", f"{len(self._chain_groups(result)):,}")],
            "Voronoi Network": [("Vertices", str(sum(len(layer.points) for layer in result.layers if layer.kind.lower() in {"vertex", "vertices"}))),
                                 ("Edges", f"{len(result.bonds):,}"),
                                 ("Surfaces", f"{result.surface_count:,}")],
        }

    @staticmethod
    def _format_measurement(value: Measurement | None, decimals: int = 3) -> str:
        if value is None or value.value is None: return "—"
        unit = f" {value.unit}" if value.unit else ""
        return f"{value.value:,.{decimals}f}{unit}"

    def _populate_analysis_sections(self, result: AnalysisResult) -> None:
        summary = result.summary if isinstance(result.summary, NetworkSummary) else None
        atoms = len(result.atoms); residues = len(self._residue_groups(result)); chains = len(self._chain_groups(result))
        network = summary.network if summary else {}
        vertices = int(network.get("vertices", Measurement(None)).value or 0) if network else sum(len(layer.points) for layer in result.layers if layer.kind.lower() in {"vertex", "vertices"})
        edges = int(network.get("edges", Measurement(None)).value or 0) if network else len(result.bonds)
        surfaces = int(network.get("surfaces", Measurement(None)).value or result.surface_count) if network else result.surface_count
        metrics = [("System", summary.system if summary and summary.system else result.name), ("Group", summary.group if summary and summary.group else "Whole system"), ("Atoms", f"{atoms:,}"), ("Residues", f"{residues:,}"), ("Chains", f"{chains:,}"), ("Network", summary.network_type if summary and summary.network_type else "—"), ("Vertices", f"{vertices:,}"), ("Edges", f"{edges:,}"), ("Surfaces", f"{surfaces:,}")]
        if summary:
            metrics += [("Volume", self._format_measurement(summary.geometry.get("volume"))), ("Surface area", self._format_measurement(summary.geometry.get("surface_area"))), ("Representative surface energy", self._format_measurement(summary.energy.get("surf_energy"), 5)), ("Energy / area", self._format_measurement(summary.energy.get("energy_per_area"), 5)), ("Mapped waters", self._format_measurement(summary.classification.get("mapped_surrounding_waters"), 0))]
        if summary:
            group_label = summary.group or "Whole system"
            network_label = summary.network_type or "unknown network"
            volume = self._format_measurement(summary.geometry.get("volume"))
            area = self._format_measurement(summary.geometry.get("surface_area"))
            energy = self._format_measurement(summary.energy.get("surf_energy"), 5)
            waters = self._format_measurement(summary.classification.get("mapped_surrounding_waters"), 0)
            self.summary_statement.setText(f"{group_label} contains {atoms:,} atoms in {residues:,} residues. The {network_label} network contains {vertices:,} vertices, {edges:,} edges, and {surfaces:,} surfaces.\n\nThe group occupies {volume} with {area} of total boundary surface area. Representative surface energy: {energy}. Mapped surrounding waters: {waters}.")
            status = ["Parsed network metadata"]
            if summary.missing_sections: status.append("Missing: " + ", ".join(summary.missing_sections))
            if summary.warnings: status.append("Warnings: " + "; ".join(summary.warnings))
            self.summary_status.setText(" • ".join(status))
        else:
            self.summary_statement.setText(f"{result.name} contains {atoms:,} atoms in {residues:,} residues and {chains:,} chains.")
            self.summary_status.setText("Summary generated from loaded structure metrics; no info.txt metadata was found.")
        self.results.setRowCount(len(metrics))
        self.results.setColumnCount(2); self.results.setHorizontalHeaderLabels(["Metric", "Value"])
        for row, (label, value) in enumerate(metrics): self.results.setItem(row, 0, QTableWidgetItem(label)); self.results.setItem(row, 1, QTableWidgetItem(str(value)))
        sections = result.info_sections or self._fallback_info_sections(result)
        for name, table in self.analysis_section_tables.items():
            if name == "Chain Composition" and summary:
                table.setRowCount(len(summary.chains))
                for row, record in enumerate(summary.chains):
                    values=[record.chain, record.atoms or "—", record.residues or "—", self._format_measurement(record.volume), self._format_measurement(record.boundary_area), self._format_measurement(record.inter_chain_area), self._format_measurement(record.solvent_interfacial_area)]
                    for col, value in enumerate(values): table.setItem(row,col,QTableWidgetItem(str(value)))
            elif name == "Residue Composition" and summary:
                table.setRowCount(len(summary.residues))
                for row, record in enumerate(summary.residues):
                    values=[record.chain, record.name, record.identifier, record.atoms or "—", self._format_measurement(record.volume), self._format_measurement(record.boundary_area), self._format_measurement(record.inter_residue_area), self._format_measurement(record.solvent_interfacial_area)]
                    for col, value in enumerate(values): table.setItem(row,col,QTableWidgetItem(str(value)))
            else:
                rows = sections.get(name, []); table.setRowCount(len(rows))
                for row, (key, value) in enumerate(rows): table.setItem(row, 0, QTableWidgetItem(key)); table.setItem(row, 1, QTableWidgetItem(value))
        values = {"atoms": atoms, "residues": residues, "chains": chains, "vertices": vertices, "edges": edges, "surfaces": surfaces, "cells": result.complete_cells, "layers": len(result.layers)}
        for key, value in values.items():
            if key in self.metric_cards: self.metric_cards[key].value.setText(f"{value:,}")

    def _build_status(self) -> None:
        self.progress_label = QLabel("Ready")
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setFixedWidth(230)
        self.statusBar().addPermanentWidget(self.progress_label)
        self.statusBar().addPermanentWidget(self.progress_bar)

    def new_project(self) -> None:
        if not self._confirm_discard_changes():
            return
        self.project = Project()
        self.project_file = None
        self.source = None
        self.current_result = None
        self._loaded_results.clear()
        self._structure_states.clear()
        self.viewer.clear_result()
        self._running_selection.clear()
        self._groups.clear()
        self._interfaces.clear()
        self._reset_network_controls()
        self.structure_browser.clear()
        self.structure_name.setText("No structure loaded")
        self.structure_summary.setText("Open a PDB structure")
        self._set_project_dirty(False)
        self.statusBar().showMessage("New project ready")

    def open_project(self) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open VorPy Workbench project",
            str(self.project_file.parent if self.project_file else Path.cwd()),
            f"VorPy Workbench projects (*{PROJECT_SUFFIX});;JSON files (*.json)",
        )
        if not filename or not self._confirm_discard_changes():
            return
        try:
            project_file = Path(filename).resolve()
            project = load_project(project_file)
            if project.structure is not None and not project.structure.source_path.exists():
                raise ValueError(
                    f"Structure file is missing: {project.structure.source_path}"
                )
            self.project = project
            self.project_file = project_file
            if project.backend_settings:
                allowed = set(VorPySolveSettings.__dataclass_fields__)
                settings = {key: value for key, value in project.backend_settings.items() if key in allowed}
                self.backend.settings = VorPySolveSettings(**settings)
            self._loading_project = True
            if project.structure is not None:
                self.load_path(project.structure.source_path)
                if project.result_state is not None:
                    restored_result = result_from_json(project.result_state)
                    self._loaded_results[project.structure.source_path] = restored_result
                    self.source = project.structure.source_path
                    self._display_result(restored_result)
                self._restore_project_groups()
                self._interfaces = {
                    interface.name: (interface.group_a, interface.group_b)
                    for interface in project.interfaces
                }
                self._refresh_groups_panel()
                self._refresh_groups_panel()
                self._restore_view_state(project.view_state)
            else:
                self.source = None
                self.current_result = None
                self._loaded_results.clear()
                self._structure_states.clear()
                self.loaded_structures.clear()
                self.viewer.clear_result()
                self._groups.clear()
                self._interfaces.clear()
                self._reset_network_controls()
                self.structure_browser.clear()
            self._set_project_dirty(False)
            self.statusBar().showMessage(f"Opened project {project.name}")
        except Exception as error:  # noqa: BLE001 - project boundary reports failures.
            self._show_error(str(error))
        finally:
            self._loading_project = False

    def save_project(self) -> bool:
        if self.project_file is None:
            return self.save_project_as()
        try:
            self._sync_project_state()
            save_project(self.project, self.project_file)
            self._set_project_dirty(False)
            self.statusBar().showMessage(f"Saved project to {self.project_file}")
            return True
        except Exception as error:  # noqa: BLE001 - project boundary reports failures.
            self._show_error(str(error))
            return False

    def save_project_as(self) -> bool:
        suggested = self.project_file or Path.cwd() / (
            self.project.name.replace(" ", "_") + PROJECT_SUFFIX
        )
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save VorPy Workbench project",
            str(suggested),
            f"VorPy Workbench projects (*{PROJECT_SUFFIX})",
        )
        if not filename:
            return False
        destination = Path(filename)
        if not str(destination).endswith(PROJECT_SUFFIX):
            destination = Path(str(destination) + PROJECT_SUFFIX)
        self.project_file = destination.resolve()
        if self.project.name == "Untitled Project":
            self.project.name = destination.name.removesuffix(PROJECT_SUFFIX)
        return self.save_project()

    def _sync_project_state(self) -> None:
        if self.source is not None:
            existing = self.project.structure
            if existing is None or existing.source_path != self.source.resolve():
                self.project.structure = StructureSource.from_path(
                    self.source, self.current_result.name if self.current_result else None
                )
        result = self.current_result
        existing_groups = {group.name: group for group in self.project.groups}
        groups = []
        if result is not None:
            atoms_by_index = {atom.index: atom for atom in result.atoms}
            for name, indices in self._groups.items():
                atoms = [atoms_by_index[index] for index in indices if index in atoms_by_index]
                previous = existing_groups.get(name)
                group = GroupDefinition.create(name, atoms)
                if previous is not None:
                    group.id = previous.id
                    group.color = previous.color
                    group.description = previous.description
                groups.append(group)
        self.project.groups = groups
        existing_interfaces = {interface.name: interface for interface in self.project.interfaces}
        self.project.interfaces = []
        for name, (group_a, group_b) in self._interfaces.items():
            previous = existing_interfaces.get(name)
            self.project.interfaces.append(
                InterfaceDefinition(
                    id=previous.id if previous else str(uuid4()),
                    name=name,
                    group_a=group_a,
                    group_b=group_b,
                )
            )
        self.project.result_state = result_to_json(result) if result is not None else None
        self.project.view_state = self._view_state_to_json()
        self.project.backend_settings = asdict(self.backend.settings)

    def _view_state_to_json(self) -> dict:
        mode = "residue"
        for action, value in ((
            (self.select_atom_action, "atom"),
            (self.select_residue_action, "residue"),
            (self.select_chain_action, "chain"),
            (self.select_molecule_action, "molecule"),
        )):
            if action.isChecked():
                mode = value
                break
        return {
            "show_cartoon": self.show_cartoon.isChecked(),
            "show_spheres": self.show_spheres.isChecked(),
            "show_sticks": self.show_sticks.isChecked(),
            "show_waters": self.show_waters.isChecked(),
            "water_style": self.water_style.currentData(),
            "show_ions": self.show_ions.isChecked(),
            "molecule_opacity": self.molecule_opacity.value(),
            "water_opacity": self.water_opacity.value(),
            "surface_opacity": self.surface_opacity.value(),
            "surface_color_scheme": self.surface_color_scheme.currentData(),
            "curvature_colormap": self.curvature_colormap.currentData(),
            "curvature_scale": self._curvature_scale_mode(),
            "edge_size": self.edge_size.value(),
            "vertex_size": self.vertex_size.value(),
            "selection_mode": mode,
            "running_selection": sorted(self._running_selection),
            "depth_clip_fraction": getattr(self.viewer, "_depth_clip_fraction", 0.0),
            "radius_overrides": dict(self.radius_overrides),
        }

    def _restore_view_state(self, state: dict) -> None:
        if not state:
            return
        self.radius_overrides = {str(key): float(value) for key, value in state.get("radius_overrides", {}).items()}
        for widget, key in ((
            (self.show_cartoon, "show_cartoon"),
            (self.show_spheres, "show_spheres"),
            (self.show_sticks, "show_sticks"),
            (self.show_waters, "show_waters"),
            (self.show_ions, "show_ions"),
        )):
            if key in state:
                widget.setChecked(bool(state[key]))
        if state.get("water_style") is not None:
            index = self.water_style.findData(state["water_style"])
            if index >= 0:
                self.water_style.setCurrentIndex(index)
        for widget, key in ((
            (self.molecule_opacity, "molecule_opacity"),
            (self.water_opacity, "water_opacity"),
            (self.surface_opacity, "surface_opacity"),
        )):
            if key in state:
                widget.setValue(int(state[key]))
        index = self.surface_color_scheme.findData(
            state.get("surface_color_scheme", "solid")
        )
        if index >= 0:
            self.surface_color_scheme.setCurrentIndex(index)
        if state.get("curvature_colormap") is not None:
            index = self.curvature_colormap.findData(state["curvature_colormap"])
            if index >= 0:
                self.curvature_colormap.setCurrentIndex(index)
        if state.get("curvature_scale") is not None:
            scale_positions = {
                "linear": 0,
                "log10": 1,
                "log100": 2,
                "signed_log": 3,   # compatibility with older projects
                "log1000": 3,
                "log10000": 4,
                "log100000": 5,
            }
            self.curvature_scale.setValue(
                scale_positions.get(str(state["curvature_scale"]), 3)
            )
        if state.get("edge_size") is not None:
            self.edge_size.setValue(int(state["edge_size"]))
        if state.get("vertex_size") is not None:
            self.vertex_size.setValue(int(state["vertex_size"]))
        actions = {
            "atom": self.select_atom_action,
            "residue": self.select_residue_action,
            "chain": self.select_chain_action,
            "molecule": self.select_molecule_action,
        }
        actions.get(state.get("selection_mode", "residue"), self.select_residue_action).setChecked(True)
        self._running_selection = {int(index) for index in state.get("running_selection", [])}
        self._update_running_selection()
        if self.current_result is not None and "radius_overrides" in state:
            self.current_result.atoms = [
                replace(atom, radius=self._radius_for_atom(atom))
                for atom in self.current_result.atoms
            ]
            self.viewer.display_result(self.current_result)
        self.viewer.set_depth_clipping_fraction(float(state.get("depth_clip_fraction", 0.0)))

    def _restore_project_groups(self) -> None:
        result = self.current_result
        if result is None:
            return
        indices_by_key: dict[AtomKey, list[int]] = {}
        for atom in result.atoms:
            indices_by_key.setdefault(AtomKey.from_atom(atom), []).append(atom.index)
        self._groups = {
            group.name: tuple(
                index
                for key in group.atom_keys
                for index in indices_by_key.get(key, ())
            )
            for group in self.project.groups
        }
        self._refresh_groups_panel()

    def _refresh_groups_panel(self) -> None:
        self.groups_list.clear()
        self.interface_group_a.clear()
        self.interface_group_b.clear()
        names = list(self._groups)
        self.interface_group_a.addItems(names)
        self.interface_group_b.addItems(names)
        for name, indices in self._groups.items():
            item = QListWidgetItem(f"{name} ({len(indices):,} atoms)")
            item.setData(Qt.UserRole, name)
            self.groups_list.addItem(item)
        self.interfaces_list.clear()
        for name, (group_a, group_b) in self._interfaces.items():
            self.interfaces_list.addItem(f"{name}: {group_a} ↔ {group_b}")
        self.make_interface_button.setEnabled(len(names) >= 2)
        self._update_solve_targets()

    def _select_saved_group(self, item: QListWidgetItem) -> None:
        name = item.data(Qt.UserRole)
        indices = self._groups.get(name, ())
        self._running_selection = set(indices)
        self._update_running_selection()
        self.selection_group.setTitle(f"Group: {name}")
        self.statusBar().showMessage(f"Selected group {name} ({len(indices):,} atoms)")

    def _make_interface(self) -> None:
        group_a = self.interface_group_a.currentText()
        group_b = self.interface_group_b.currentText()
        if not group_a or not group_b or group_a == group_b:
            QMessageBox.information(self, "Make interface", "Choose two different groups.")
            return
        name, accepted = QInputDialog.getText(
            self, "Make interface", "Interface name", text=f"{group_a}-{group_b}"
        )
        if not accepted or not name.strip():
            return
        self._interfaces[name.strip()] = (group_a, group_b)
        self._sync_project_state()
        self._set_project_dirty()
        self._refresh_groups_panel()
        self.statusBar().showMessage(f"Created interface {name.strip()}")

    def _confirm_discard_changes(self) -> bool:
        if not self._project_dirty:
            return True
        choice = QMessageBox.question(
            self,
            "Unsaved project",
            "Save changes to the current project?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
            QMessageBox.Save,
        )
        if choice == QMessageBox.Cancel:
            return False
        if choice == QMessageBox.Save:
            return self.save_project()
        return True

    def _set_project_dirty(self, dirty: bool = True) -> None:
        self._project_dirty = dirty
        self._update_window_title()

    def _update_window_title(self) -> None:
        subject = self.current_result.name if self.current_result else self.project.name
        marker = " *" if self._project_dirty else ""
        self.setWindowTitle(f"VorPy — {subject}{marker}")

    def closeEvent(self, event: QCloseEvent) -> None:
        if self._confirm_discard_changes():
            event.accept()
        else:
            event.ignore()

    def _update_chemical_info(self, result: AnalysisResult | None) -> None:
        if result is None:
            self.chemical_info.setVisible(False)
            return
        residues = self._residue_groups(result)
        chains = self._chain_groups(result)
        molecules = self._molecule_groups(result)
        waters = sum(1 for atom in result.atoms if atom.residue_name.upper() in WATER_RESIDUES)
        ions = sum(1 for atom in result.atoms if atom.residue_name.upper() in ION_RESIDUES)
        values = ((self.info_atoms, len(result.atoms)), (self.info_bonds, len(result.bonds)),
                  (self.info_residues, len(residues)), (self.info_chains, len(chains)),
                  (self.info_molecules, len(molecules)), (self.info_waters, waters),
                  (self.info_ions, ions))
        for widget, value in values:
            widget.setText(f"{value:,}")
        self.chemical_info.setVisible(True)

    def _delete_current_structure(self) -> None:
        if self.source is None:
            return
        if QMessageBox.question(self, "Delete structure", "Delete the current structure from this session?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            return
        source = self.source
        self._save_current_structure_state()
        self._loaded_results.pop(source, None)
        self._structure_states.pop(source, None)
        self.source = None
        self.current_result = None
        self._running_selection.clear()
        self._groups.clear()
        self._interfaces.clear()
        if self._loaded_results:
            next_source = next(iter(self._loaded_results))
            self.source = next_source
            self._display_result(self._loaded_results[next_source])
            self._restore_structure_state(next_source)
        else:
            self.viewer.clear_result()
            self.structure_name.setText("No structure loaded")
            self.structure_summary.setText("Browse to load a molecular structure")
            self._update_chemical_info(None)
            self._populate_structure_browser()
        self._populate_loaded_structures()
        self._set_project_dirty()

    def _reset_structures(self) -> None:
        if not self._loaded_results:
            return
        if QMessageBox.question(self, "Reset structures", "Remove all loaded structures and selections?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No) != QMessageBox.Yes:
            return
        self._loaded_results.clear()
        self._structure_states.clear()
        self.source = None
        self.current_result = None
        self._running_selection.clear()
        self._groups.clear()
        self._interfaces.clear()
        self.viewer.clear_result()
        self.structure_name.setText("No structure loaded")
        self.structure_summary.setText("Browse to load a molecular structure")
        self._update_chemical_info(None)
        self._populate_loaded_structures()
        self._populate_structure_browser()
        self._set_project_dirty()

    def _save_current_structure_state(self) -> None:
        if self.source is None or self.current_result is None:
            return
        self._structure_states[self.source] = (
            set(self._running_selection),
            dict(self._groups),
            dict(self._interfaces),
        )

    def _restore_structure_state(self, source: Path) -> None:
        selection, groups, interfaces = self._structure_states.get(
            source, (set(), {}, {})
        )
        self._running_selection = set(selection)
        self._groups = dict(groups)
        self._interfaces = dict(interfaces)
        if self.current_result is not None:
            self._populate_structure_browser()
            self._update_running_selection()

    def _populate_loaded_structures(self) -> None:
        self.loaded_structures.blockSignals(True)
        self.loaded_structures.clear()
        for source, result in self._loaded_results.items():
            item = QListWidgetItem(result.name)
            item.setToolTip(str(source))
            item.setData(Qt.UserRole, source)
            self.loaded_structures.addItem(item)
            if source == self.source:
                self.loaded_structures.setCurrentItem(item)
        has_structures = bool(self._loaded_results)
        self.loaded_structures_label.setVisible(has_structures)
        self.loaded_structures.setVisible(has_structures)
        self.loaded_structures.blockSignals(False)

    def _switch_structure_item(self, item: QListWidgetItem) -> None:
        source = Path(item.data(Qt.UserRole))
        if source == self.source:
            return
        self._save_current_structure_state()
        result = self._loaded_results.get(source)
        if result is None:
            return
        self.source = source
        self._display_result(result)
        self._restore_structure_state(source)
        self._populate_loaded_structures()
        self._set_project_dirty()

    def open_structure(self) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open PDB structure",
            str(DEFAULT_DATA_DIRECTORY),
            "PDB structures (*.pdb);;All files (*)",
        )
        if not filename:
            return
        self.load_path(Path(filename))

    def load_path(self, path: Path) -> None:
        """Load a structure, retaining prior structures for later switching."""
        try:
            source = path.expanduser().resolve()
            if source in self._loaded_results:
                self._switch_structure_item(
                    next(
                        self.loaded_structures.item(row)
                        for row in range(self.loaded_structures.count())
                        if self.loaded_structures.item(row).data(Qt.UserRole) == source
                    )
                )
                return
            self._save_current_structure_state()
            result = load_result_directory(source) if source.is_dir() else load_pdb(source)
            if source.is_dir():
                info_path = source / "info.txt"
                if info_path.exists():
                    result.info_sections = parse_group_info(info_path)
                    result.summary = parse_network_summary(info_path)
            self.source = source
            self._loaded_results[source] = result
            self._display_result(result)
            self._restore_structure_state(source)
            self._populate_loaded_structures()
            if not self._loading_project:
                self.project.structure = StructureSource.from_path(
                    source, result.name
                )
                self.project.groups.clear()
                self._set_project_dirty()
        except Exception as error:  # noqa: BLE001 - GUI boundary reports loader failures.
            self._show_error(str(error))

    def open_result_directory(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self, "Open completed VorPy output"
        )
        if not directory:
            return
        self.load_path(Path(directory))

    def _solve_target_indices(self) -> tuple[int, ...] | None:
        target, name = self.solve_target.currentData() or ("whole", "")
        if target == "whole":
            return None
        if target == "selection":
            if not self._running_selection:
                raise ValueError("Choose atoms in Structure or select a group before solving.")
            return tuple(sorted(self._running_selection))
        if target == "group":
            indices = self._groups.get(name, ())
            if not indices:
                raise ValueError(f"Group {name!r} is empty.")
            return tuple(indices)
        if target == "interface":
            group_names = self._interfaces.get(name)
            if group_names is None:
                raise ValueError(f"Interface {name!r} is unavailable.")
            indices = set(self._groups.get(group_names[0], ()))
            indices.update(self._groups.get(group_names[1], ()))
            if not indices:
                raise ValueError(f"Interface {name!r} has no atoms.")
            return tuple(sorted(indices))
        raise ValueError(f"Unknown solve target: {target}")

    def solve(self) -> None:
        if self._thread is not None:
            return
        try:
            selected_indices = self._solve_target_indices()
        except ValueError as error:
            self._show_error(str(error))
            return
        network_types = {0: "aw", 1: "pow", 2: "prm"}
        self.backend = VorPyBackend(
            VorPySolveSettings(
                network_type=network_types[self.network_type.currentIndex()],
                max_vertices=self.max_vertices.value(),
                box_size=self.box_size.value(),
                surface_resolution=self.surface_resolution.value(),
                build_surfaces=True,
                build_vertices=True,
                build_edges=True,
            )
        )
        self.solve_action.setEnabled(False)
        self.cancel_action.setEnabled(True)
        self.progress_bar.setValue(0)
        self._thread = QThread(self)
        self._worker = SolveWorker(self.backend, self.source, selected_indices)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self._show_progress)
        self._worker.completed.connect(self._display_result)
        self._worker.failed.connect(self._show_error)
        self._worker.cancelled.connect(self._show_cancelled)
        self._worker.finished.connect(self._thread.quit)
        self._worker.finished.connect(self._worker.deleteLater)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.finished.connect(self._thread_finished)
        self._thread.start()

    def cancel_analysis(self) -> None:
        if self._worker is not None:
            self._worker.cancel()
            self.progress_label.setText("Cancelling…")
            self.statusBar().showMessage("Cancelling analysis…")

    def _fit_view(self) -> None:
        self.viewer.reset_depth_clipping()
        self.viewer.plotter.reset_camera()

    def save_screenshot(self) -> None:
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save viewer screenshot", "vorpy_view.png", "PNG images (*.png)"
        )
        if filename:
            try:
                scale, accepted = QInputDialog.getInt(
                    self,
                    "Screenshot resolution",
                    "Resolution multiplier",
                    value=2,
                    minValue=1,
                    maxValue=8,
                )
                if not accepted:
                    return
                self.viewer.save_screenshot(filename, scale)
                self.statusBar().showMessage(
                    f"Screenshot saved to {filename} at {scale}× resolution"
                )
            except Exception as error:  # noqa: BLE001 - VTK writers expose varied exceptions.
                self._show_error(str(error))

    def _radius_for_atom(self, atom: Atom) -> float:
        residue = atom.residue_name.upper()
        specific = self.radius_overrides.get(f"RES:{residue}:{atom.name}")
        if specific is not None:
            return specific
        if residue in ION_RESIDUES:
            ion_radius = self.radius_overrides.get(f"ION:{atom.element.upper()}")
            if ion_radius is not None:
                return ion_radius
        return self.radius_overrides.get(atom.element.upper(), atom.radius)

    def _display_result(self, result: AnalysisResult) -> None:
        previous = self.current_result
        same_structure = (
            previous is not None
            and previous.source == result.source
            and len(previous.atoms) == len(result.atoms)
        )
        if self.radius_overrides:
            result.atoms = [replace(atom, radius=self._radius_for_atom(atom)) for atom in result.atoms]
        self.current_result = result
        main_atom_count = sum(
            not (atom.residue_name.upper() in WATER_RESIDUES or atom.residue_name.upper() in ION_RESIDUES)
            for atom in result.atoms
        )
        self.show_cartoon.setChecked(True)
        self.show_sticks.setChecked(main_atom_count <= 1500)
        if self.source is not None:
            self._loaded_results[self.source] = result
        self.viewer.display_result(result)
        self._update_window_title()
        self.structure_name.setText(result.name)
        residues = self._residue_groups(result)
        chains = self._chain_groups(result)
        self.structure_summary.setText(
            f"{len(result.atoms):,} atoms  •  {len(result.bonds):,} bonds  •  "
            f"{len(residues):,} residues  •  {len(chains):,} chains"
        )
        self._update_chemical_info(result)
        if not same_structure:
            self._running_selection.clear()
            self._groups.clear()
        else:
            valid_indices = {atom.index for atom in result.atoms}
            self._running_selection.intersection_update(valid_indices)
            self._groups = {
                name: tuple(index for index in indices if index in valid_indices)
                for name, indices in self._groups.items()
            }
        self._refresh_groups_panel()
        self._populate_network_controls(result)
        self._populate_structure_browser()
        self._update_running_selection()
        self._populate_analysis_sections(result)
        self.progress_label.setText("Complete")
        self.progress_bar.setValue(100)
        self.statusBar().showMessage(f"{result.name} ready")

    @staticmethod
    def _network_layer_key(layer) -> str | None:
        kind = layer.kind.lower().strip()
        kind = {"vertex": "vertices", "surface": "surfaces", "edge": "edges"}.get(
            kind, kind
        )
        if kind not in {"edges", "vertices", "surfaces"}:
            return None
        return f"shell_{kind}" if "shell" in layer.name.lower() else kind

    def _network_layers(self, key: str):
        if self.current_result is None:
            return []
        return [
            layer
            for layer in self.current_result.layers
            if self._network_layer_key(layer) == key
        ]

    @staticmethod
    def _set_color_button_swatch(button: QPushButton, color: str) -> None:
        button.setStyleSheet(
            f"QPushButton {{ border-left: 18px solid {color}; padding-left: 7px; }}"
        )

    def _reset_network_controls(self) -> None:
        for key, checkbox in self.network_layer_checks.items():
            checkbox.blockSignals(True)
            checkbox.setChecked(False)
            checkbox.setEnabled(False)
            checkbox.blockSignals(False)
            self.network_color_buttons[key].setEnabled(False)
        self.surface_color_scheme.setEnabled(False)
        self.curvature_colormap.setEnabled(False)
        self.curvature_scale.setEnabled(False)
        self.edge_size.setEnabled(False)
        self.vertex_size.setEnabled(False)
        self.surface_opacity.setEnabled(False)

    def _populate_network_controls(self, result: AnalysisResult) -> None:
        for key, checkbox in self.network_layer_checks.items():
            layers = [
                layer for layer in result.layers if self._network_layer_key(layer) == key
            ]
            checkbox.blockSignals(True)
            checkbox.setEnabled(bool(layers))
            checkbox.setChecked(any(layer.visible for layer in layers))
            checkbox.blockSignals(False)
            self.network_color_buttons[key].setEnabled(bool(layers))
            if layers:
                self.network_colors[key] = layers[0].color
                self._set_color_button_swatch(
                    self.network_color_buttons[key], layers[0].color
                )
        surface_layers = [layer for layer in result.layers if self._network_layer_key(layer) in {"surfaces", "shell_surfaces"}]
        scalar_layers = [layer for layer in result.layers if layer.cell_scalars]
        self.surface_opacity.blockSignals(True)
        if surface_layers:
            self.surface_opacity.setValue(round(surface_layers[0].opacity * 100))
        has_scalar_colors = bool(scalar_layers)
        self.surface_color_scheme.setEnabled(has_scalar_colors)
        self.curvature_colormap.setEnabled(has_scalar_colors)
        self.curvature_scale.setEnabled(has_scalar_colors)
        if has_scalar_colors:
            self._set_surface_color_scheme(self.surface_color_scheme.currentIndex())
            self._set_curvature_colormap(self.curvature_colormap.currentIndex())
            self._set_curvature_scale(self.curvature_scale.value())
        else:
            self.surface_color_scheme.blockSignals(True)
            self.surface_color_scheme.setCurrentIndex(0)
            self.surface_color_scheme.blockSignals(False)
        self.surface_opacity.setEnabled(bool(surface_layers))
        self.surface_opacity.blockSignals(False)
        self.edge_size.setEnabled(
            bool(self._network_layers("edges") or self._network_layers("shell_edges"))
        )
        self.vertex_size.setEnabled(
            bool(
                self._network_layers("vertices")
                or self._network_layers("shell_vertices")
            )
        )

    @staticmethod
    def _residue_groups(result: AnalysisResult) -> list[tuple[str, tuple[int, ...]]]:
        groups: dict[tuple[str, str, str], list[int]] = {}
        for atom in result.atoms:
            key = (atom.chain, atom.residue_sequence, atom.residue_name)
            groups.setdefault(key, []).append(atom.index)
        return [
            (
                (
                    f"{residue or 'Unknown'} {sequence or '?'}"
                    f"{f' · chain {chain}' if chain else ''} ({len(indices)} atoms)"
                ),
                tuple(indices),
            )
            for (chain, sequence, residue), indices in groups.items()
        ]

    @staticmethod
    def _chain_groups(result: AnalysisResult) -> list[tuple[str, tuple[int, ...]]]:
        groups: dict[tuple[str, str, str], list[int]] = {}
        for atom in result.atoms:
            if atom.chain:
                key = ("chain", atom.chain, "")
            elif atom.residue_name.strip().upper() in WATER_RESIDUES:
                key = ("water", atom.residue_sequence, atom.residue_name)
            else:
                key = ("unassigned", "", "")
            groups.setdefault(key, []).append(atom.index)
        entries = []
        for (kind, identity, residue), indices in groups.items():
            if kind == "chain":
                label = f"Chain {identity}"
            elif kind == "unassigned":
                label = "No chain · non-water atoms"
            else:
                label = f"No chain · {residue or 'Unknown'} {identity or '?'}"
            atom_label = "atom" if len(indices) == 1 else "atoms"
            entries.append((f"{label} ({len(indices)} {atom_label})", tuple(indices)))
        return entries

    @staticmethod
    def _molecule_groups(result: AnalysisResult) -> list[tuple[str, tuple[int, ...]]]:
        neighbors = {atom.index: set() for atom in result.atoms}
        for bond in result.bonds:
            neighbors[bond.atom_a].add(bond.atom_b)
            neighbors[bond.atom_b].add(bond.atom_a)
        components: list[tuple[str, tuple[int, ...]]] = []
        unseen = set(neighbors)
        while unseen:
            first = min(unseen)
            pending = [first]
            component: set[int] = set()
            while pending:
                atom_index = pending.pop()
                if atom_index in component:
                    continue
                component.add(atom_index)
                pending.extend(neighbors[atom_index] - component)
            unseen -= component
            indices = tuple(sorted(component))
            components.append(
                (f"Molecule {len(components) + 1} ({len(indices)} atoms)", indices)
            )
        return components

    def _populate_structure_browser(self) -> None:
        self.structure_browser.blockSignals(True)
        self.structure_browser.clear()
        result = self.current_result
        if result is None:
            self._selection_entries = []
        else:
            category = self.structure_browser_type.currentText()
            if category == "Atoms":
                self._selection_entries = [
                    (
                        (
                            f"#{atom.serial} · {atom.name} · "
                            f"{atom.residue_name} {atom.residue_sequence}"
                            f"{f' · chain {atom.chain}' if atom.chain else ''}"
                        ),
                        (atom.index,),
                    )
                    for atom in result.atoms
                ]
            elif category == "Residues":
                self._selection_entries = self._residue_groups(result)
            elif category == "Chains":
                self._selection_entries = self._chain_groups(result)
            else:
                self._selection_entries = self._molecule_groups(result)

            for label, atom_indices in self._selection_entries:
                item = QListWidgetItem(label)
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                selected_count = len(self._running_selection.intersection(atom_indices))
                if selected_count == len(atom_indices):
                    state = Qt.Checked
                elif selected_count:
                    state = Qt.PartiallyChecked
                else:
                    state = Qt.Unchecked
                item.setCheckState(state)
                item.setData(Qt.UserRole, atom_indices)
                self.structure_browser.addItem(item)
        self.structure_browser.blockSignals(False)
        self._filter_structure_browser(self.structure_search.text())

    def _filter_structure_browser(self, query: str) -> None:
        query = query.strip().casefold()
        for row in range(self.structure_browser.count()):
            item = self.structure_browser.item(row)
            item.setHidden(bool(query) and query not in item.text().casefold())

    def _browser_item_changed(self, item: QListWidgetItem) -> None:
        atom_indices = set(item.data(Qt.UserRole) or ())
        if item.checkState() == Qt.Checked:
            self._running_selection.update(atom_indices)
        else:
            self._running_selection.difference_update(atom_indices)
        self._update_running_selection()
        self._populate_structure_browser()

    def _update_running_selection(self) -> None:
        result = self.current_result
        if result is None or not self._running_selection:
            self.viewer.set_group_selection([])
            self.running_selection.setText("No atoms selected")
            self.make_group_button.setEnabled(False)
            return
        atoms = [result.atoms[index] for index in sorted(self._running_selection)]
        self.viewer.set_group_selection(atoms)
        residues = {
            (atom.residue_name, atom.residue_sequence, atom.chain) for atom in atoms
        }
        chains = {atom.chain or "(blank)" for atom in atoms}
        residue_preview = ", ".join(
            f"{name} {sequence}{f' ({chain})' if chain else ''}"
            for name, sequence, chain in sorted(residues)[:4]
        )
        if len(residues) > 4:
            residue_preview += f", +{len(residues) - 4} more"
        self.running_selection.setText(
            f"{len(atoms):,} atoms · {len(residues):,} residues · "
            f"{len(chains):,} chains\n{residue_preview}"
        )
        self.make_group_button.setEnabled(True)

    def reset_selection(self) -> None:
        choice = QMessageBox.question(
            self,
            "Reset selection",
            "Reset Selection?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if choice != QMessageBox.Yes:
            return
        self._clear_selection("Selection reset")

    def _clear_selection_from_viewer(self) -> None:
        self._clear_selection("Selection cleared")

    def _clear_selection(self, status: str) -> None:
        self._running_selection.clear()
        self._update_running_selection()
        self._populate_structure_browser()
        self.selection_group.setTitle("No selection")
        self.atom_name.setText("—")
        self.atom_element.setText("—")
        self.atom_residue.setText("—")
        self.atom_chain.setText("—")
        self.atom_position.setText("—")
        self.selection_count.setText("—")
        self.statusBar().showMessage(status)

    def _make_group(self) -> None:
        if not self._running_selection or self.current_result is None:
            return
        default_name = f"Group {len(self._groups) + 1}"
        name, accepted = QInputDialog.getText(
            self, "Make group", "Group name", text=default_name
        )
        name = name.strip()
        if not accepted or not name:
            return
        self._groups[name] = tuple(sorted(self._running_selection))
        self._sync_project_state()
        self._set_project_dirty()
        self._refresh_groups_panel()
        self.statusBar().showMessage(
            f"Created {name} with {len(self._running_selection):,} atoms"
        )

    def _set_network_layer_visible(self, key: str, visible: bool) -> None:
        for layer in self._network_layers(key):
            layer.visible = visible
            self.viewer.set_layer_visible(layer.name, visible)

    def _set_surface_color_scheme(self, _index: int) -> None:
        scheme = self.surface_color_scheme.currentData()
        eligible = [layer for layer in (self.current_result.layers if self.current_result else []) if scheme == "solid" or scheme in layer.cell_scalars]
        if scheme != "solid" and not any(layer.visible for layer in eligible):
            shell = next((layer for layer in eligible if "shell" in layer.name.lower()), None)
            preferred = shell or (eligible[0] if eligible else None)
            if preferred is not None:
                key = self._network_layer_key(preferred)
                if key in self.network_layer_checks:
                    self.network_layer_checks[key].setChecked(True)
        for layer in (self.current_result.layers if self.current_result else []):
            if scheme == "solid" or scheme in layer.cell_scalars:
                layer.color_scheme = scheme
                self.viewer.set_layer_color_scheme(layer.name, scheme)

    def _set_curvature_colormap(self, _index: int) -> None:
        color_map = self.curvature_colormap.currentData() or "coolwarm"
        for layer in (self.current_result.layers if self.current_result else []):
            layer.color_map = color_map
            setter = getattr(self.viewer, "set_layer_colormap", None)
            if setter is not None:
                setter(layer.name, color_map)

    @staticmethod
    def _curvature_scale_options() -> tuple[tuple[str, str], ...]:
        return (
            ("Linear", "linear"),
            ("Log10", "log10"),
            ("Log100", "log100"),
            ("Log1000", "log1000"),
            ("Log10000", "log10000"),
            ("Log100000", "log100000"),
        )

    def _curvature_scale_mode(self) -> str:
        options = self._curvature_scale_options()
        index = max(0, min(int(self.curvature_scale.value()), len(options) - 1))
        return options[index][1]

    def _set_curvature_scale(self, _value: int) -> None:
        options = self._curvature_scale_options()
        index = max(0, min(int(self.curvature_scale.value()), len(options) - 1))
        label, scale = options[index]
        self.curvature_scale_label.setText(label)
        self.curvature_scale.setToolTip(
            "Linear" if scale == "linear" else f"Signed logarithmic compression ({label})"
        )

        for layer in (self.current_result.layers if self.current_result else []):
            layer.scale_mode = scale
            setter = getattr(self.viewer, "set_layer_scale_mode", None)
            if setter is not None:
                setter(layer.name, scale)

    def _set_edge_size(self, value: int) -> None:
        setter = getattr(self.viewer, "set_edge_radius", None)
        if setter is not None:
            setter(float(value) / 1000.0)

    def _set_vertex_size(self, value: int) -> None:
        setter = getattr(self.viewer, "set_vertex_radius", None)
        if setter is not None:
            setter(float(value) / 100.0)

    def _set_surface_opacity(self, value: int) -> None:
        opacity = value / 100.0
        for key in ("surfaces", "shell_surfaces"):
            for layer in self._network_layers(key):
                layer.opacity = opacity
                self.viewer.set_layer_opacity(layer.name, opacity)


    def _choose_network_color(self, key: str) -> None:
        color = QColorDialog.getColor(
            QColor(self.network_colors[key]), self, "Choose network color"
        )
        if color.isValid():
            # Manual per-layer colors are the solid-color representation.
            self.surface_color_scheme.setCurrentIndex(
                self.surface_color_scheme.findData("solid")
            )
            color_name = color.name()
            self.network_colors[key] = color_name
            self._set_color_button_swatch(self.network_color_buttons[key], color_name)
            for layer in self._network_layers(key):
                layer.color = color_name
                self.viewer.set_layer_color(layer.name, color_name)

    def _apply_picked_atoms(self, atoms: list[Atom], additive: bool) -> str:
        indices = {atom.index for atom in atoms}
        if additive:
            if indices.issubset(self._running_selection):
                self._running_selection.difference_update(indices)
                action = "Removed"
            else:
                self._running_selection.update(indices)
                action = "Added"
        else:
            self._running_selection = indices
            action = "Selected"
        self._update_running_selection()
        self._populate_structure_browser()
        return action

    def _show_selected_atom(self, atom: Atom, additive: bool = False) -> None:
        action = self._apply_picked_atoms([atom], additive)
        self.selection_group.setTitle("Selected atoms" if additive else "Selected atom")
        self.atom_name.setText(f"{atom.name} (#{atom.serial})")
        self.atom_element.setText(atom.element)
        self.atom_residue.setText(
            f"{atom.residue_name} {atom.residue_sequence}".strip() or "—"
        )
        self.atom_chain.setText(atom.chain or "—")
        self.atom_position.setText(", ".join(f"{value:.3f}" for value in atom.position))
        self.selection_count.setText(str(len(self._running_selection)))
        self.statusBar().showMessage(
            f"{action} {atom.name}, "
            f"{atom.residue_name} {atom.residue_sequence} "
            f"({len(self._running_selection)} atoms total)"
        )

    def _show_selected_residue(
        self, atoms: list[Atom], additive: bool = False
    ) -> None:
        if not atoms:
            return
        # Mouse residue picks and Structure-browser residue picks share the
        # same temporary selection, highlight, and group-creation path.
        action = self._apply_picked_atoms(atoms, additive)
        atom = atoms[0]
        self.selection_group.setTitle(
            "Selected residues" if additive else "Selected residue"
        )
        self.atom_name.setText("Multiple")
        self.atom_element.setText("—")
        self.atom_residue.setText(
            f"{atom.residue_name} {atom.residue_sequence}".strip()
        )
        self.atom_chain.setText(atom.chain or "—")
        self.atom_position.setText("—")
        self.selection_count.setText(str(len(self._running_selection)))
        self.statusBar().showMessage(
            f"{action} {atom.residue_name} "
            f"{atom.residue_sequence}, chain {atom.chain or '—'} "
            f"({len(self._running_selection)} atoms total)"
        )

    def _show_selected_chain(self, atoms: list[Atom], additive: bool = False) -> None:
        if not atoms:
            return
        atom = atoms[0]
        chain = atom.chain or "—"
        description = (
            f"Chain {chain}"
            if atom.chain
            else f"No chain · {atom.residue_name or 'Unknown'} "
            f"{atom.residue_sequence or '?'}"
        )
        self._show_selected_collection(
            atoms, "chain", description, chain, additive
        )

    def _show_selected_molecule(
        self, atoms: list[Atom], additive: bool = False
    ) -> None:
        if not atoms:
            return
        residues = {
            (atom.chain, atom.residue_sequence, atom.residue_name) for atom in atoms
        }
        self._show_selected_collection(
            atoms,
            "molecule",
            f"Molecule ({len(residues)} residues)",
            atoms[0].chain or "—",
            additive,
        )

    def _show_selected_collection(
        self,
        atoms: list[Atom],
        kind: str,
        description: str,
        chain: str,
        additive: bool = False,
    ) -> None:
        action = self._apply_picked_atoms(atoms, additive)
        self.selection_group.setTitle(
            f"Selected {kind}s" if additive else f"Selected {kind}"
        )
        self.atom_name.setText("Multiple")
        self.atom_element.setText("—")
        self.atom_residue.setText(description)
        self.atom_chain.setText(chain)
        self.atom_position.setText("—")
        self.selection_count.setText(str(len(self._running_selection)))
        self.statusBar().showMessage(
            f"{action} {description} "
            f"({len(self._running_selection)} atoms total)"
        )

    def _set_selection_mode(self, mode: str | None) -> None:
        active = self.selection_actions.checkedAction()
        if mode is None and active is not None:
            return
        self.viewer.set_selection_mode(mode)
        self.selection_mode_label.setText(
            f"{mode.title()} selection active" if mode else "Navigation mode"
        )
        if mode is not None:
            self.inspector.setCurrentIndex(1)
            self.statusBar().showMessage(
                f"{mode.title()} selection active — Shift-click toggles selections"
            )

    def _show_progress(self, label: str, value: int) -> None:
        self.progress_label.setText(label)
        self.progress_bar.setValue(value)
        self.statusBar().showMessage(label)

    def _show_error(self, message: str) -> None:
        QMessageBox.critical(self, "Operation failed", message)
        self.statusBar().showMessage(message)

    def _show_cancelled(self) -> None:
        self.progress_label.setText("Cancelled")
        self.statusBar().showMessage("Analysis cancelled")

    def _thread_finished(self) -> None:
        self._thread = None
        self._worker = None
        self.solve_action.setEnabled(True)
        self.cancel_action.setEnabled(False)
