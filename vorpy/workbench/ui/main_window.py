"""Viewer-forward scientific workbench application shell."""

from __future__ import annotations

from pathlib import Path
from uuid import uuid4

from PySide6.QtCore import Qt, QThread
from PySide6.QtGui import QAction, QActionGroup, QCloseEvent, QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
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
    QTreeWidget,
    QTreeWidgetItem,
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
)
from vorpy.workbench.services.result_directory import load_result_directory
from vorpy.workbench.services.structure_loader import load_pdb
from vorpy.workbench.services.vorpy_backend import VorPyBackend, VorPySolveSettings
from vorpy.workbench.ui.molecular_view import (
    ION_RESIDUES,
    WATER_RESIDUES,
    MolecularView,
)
from vorpy.workbench.workers.solve_worker import SolveWorker

DEFAULT_DATA_DIRECTORY = Path(__file__).resolve().parents[2] / "data"


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
    LAYER_ROLE = Qt.UserRole + 1

    def __init__(self):
        super().__init__()
        self.setWindowTitle("VorPy")
        self.resize(1480, 900)
        self.setMinimumSize(1050, 680)
        self.backend = VorPyBackend()
        self.project = Project()
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
        self._build_actions()
        self._build_menu_and_toolbar()
        self._build_workspace()
        self._build_status()
        self.statusBar().showMessage(
            "Ready — load a structure, then choose atom or residue selection"
        )

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
        self.fit_action.triggered.connect(self.viewer.plotter.reset_camera)
        self.screenshot_action = QAction(
            style.standardIcon(QStyle.SP_DialogSaveButton), "Screenshot", self
        )
        self.screenshot_action.triggered.connect(self.save_screenshot)

        self.selection_actions = QActionGroup(self)
        self.selection_actions.setExclusive(True)
        self.select_atom_action = QAction("Select atom", self, checkable=True)
        self.select_atom_action.setShortcut("A")
        self.select_residue_action = QAction("Select residue", self, checkable=True)
        self.select_residue_action.setShortcut("R")
        self.selection_actions.addAction(self.select_atom_action)
        self.selection_actions.addAction(self.select_residue_action)
        self.select_atom_action.toggled.connect(
            lambda checked: self._set_selection_mode("atom" if checked else None)
        )
        self.select_residue_action.toggled.connect(
            lambda checked: self._set_selection_mode("residue" if checked else None)
        )

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
        selection_menu.addActions([self.select_atom_action, self.select_residue_action])
        self.menuBar().addMenu("Analysis").addAction(self.solve_action)
        self.menuBar().addMenu("View").addAction(self.fit_action)
        self.menuBar().addMenu("Help")

    def _build_workspace(self) -> None:
        root = QWidget()
        layout = QHBoxLayout(root)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        workflow = QWidget()
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
        layout.addWidget(workflow)

        vertical = QSplitter(Qt.Vertical)
        upper = QSplitter(Qt.Horizontal)
        upper.addWidget(self._build_viewport())
        upper.addWidget(self._build_inspector())
        upper.setStretchFactor(0, 5)
        upper.setStretchFactor(1, 1)
        upper.setSizes([1050, 320])
        vertical.addWidget(upper)
        vertical.addWidget(self._build_analysis_tray())
        vertical.setStretchFactor(0, 5)
        vertical.setStretchFactor(1, 1)
        vertical.setSizes([610, 270])
        layout.addWidget(vertical, 1)
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
        self.inspector = QTabWidget()
        self.inspector.setMinimumWidth(290)
        self.inspector.setMaximumWidth(390)
        self.inspector.addTab(self._build_visualization_tab(), "Visual")
        self.inspector.addTab(self._build_layers_tab(), "Layers")
        return self.inspector

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
        self.mesh_format = QComboBox()
        self.mesh_format.addItems(["OFF", "PLY", "VTP"])
        form.addRow("Target", self.solve_target)
        form.addRow("Network", self.network_type)
        form.addRow("Mesh format", self.mesh_format)
        outputs = QHBoxLayout()
        self.build_vertices = QCheckBox("Vertices")
        self.build_edges = QCheckBox("Edges")
        self.build_surfaces = QCheckBox("Surfaces")
        for checkbox in (self.build_vertices, self.build_edges, self.build_surfaces):
            checkbox.setChecked(True)
            outputs.addWidget(checkbox)
        form.addRow("Build", outputs)
        self.build_settings_button = QPushButton("Detailed build settings…")
        self.build_settings_button.clicked.connect(self._open_build_settings)
        form.addRow(self.build_settings_button)
        layout.addWidget(section)
        return panel

    def _open_build_settings(self) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("VorPy build settings")
        dialog.setModal(True)
        form = QFormLayout(dialog)
        network = QComboBox()
        network.addItems(["Atomic Voronoi", "Power", "Primitive"])
        network.setCurrentIndex(self.network_type.currentIndex())
        max_vertices = QSpinBox()
        max_vertices.setRange(1, 10000)
        max_vertices.setValue(self.max_vertices.value())
        box_size = QDoubleSpinBox()
        box_size.setRange(1.0, 10.0)
        box_size.setSingleStep(0.05)
        box_size.setValue(self.box_size.value())
        surface_resolution = QDoubleSpinBox()
        surface_resolution.setRange(0.01, 1.0)
        surface_resolution.setSingleStep(0.01)
        surface_resolution.setValue(self.surface_resolution.value())
        form.addRow("Network", network)
        form.addRow("Maximum vertices", max_vertices)
        form.addRow("Box size", box_size)
        form.addRow("Surface resolution", surface_resolution)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        form.addRow(buttons)
        if dialog.exec() == QDialog.Accepted:
            self.network_type.setCurrentIndex(network.currentIndex())
            self.max_vertices.setValue(max_vertices.value())
            self.box_size.setValue(box_size.value())
            self.surface_resolution.setValue(surface_resolution.value())

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

    def _build_layers_tab(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        self.layer_tree = QTreeWidget()
        self.layer_tree.setHeaderLabel("VorPy geometry")
        self.layer_tree.itemChanged.connect(self._tree_item_changed)
        layout.addWidget(self.layer_tree, 1)
        styling = QGroupBox("Selected layer")
        styling_form = QFormLayout(styling)
        self.layer_opacity = QSlider(Qt.Horizontal)
        self.layer_opacity.setRange(0, 100)
        self.layer_opacity.setValue(70)
        self.layer_opacity.valueChanged.connect(self._change_selected_layer_opacity)
        self.layer_color = QPushButton("Choose color…")
        self.layer_color.clicked.connect(self._change_selected_layer_color)
        styling_form.addRow("Opacity", self.layer_opacity)
        styling_form.addRow("Color", self.layer_color)
        layout.addWidget(styling)
        return panel

    def _build_analysis_tray(self) -> QWidget:
        tray = QSplitter(Qt.Horizontal)
        tray.setObjectName("bottomTray")
        solve = self._build_solve_section()
        solve.setMinimumWidth(300)
        tray.addWidget(solve)
        tabs = QTabWidget()
        self.analysis_tray_tabs = tabs
        analysis = QWidget()
        cards_layout = QHBoxLayout(analysis)
        cards_layout.setContentsMargins(8, 8, 8, 8)
        self.metric_cards = {}
        for key, title, accent in ((
            "atoms", "Atoms", "#6857d9"),
            ("bonds", "Bonds", "#367bd6"),
            ("cells", "Complete cells", "#38a873"),
            ("surfaces", "Surfaces", "#9b59db"),
            ("layers", "Geometry layers", "#d99a32"),
        ):
            card = MetricCard(title, accent)
            self.metric_cards[key] = card
            cards_layout.addWidget(card)
        self.results = QTableWidget(0, 2)
        self.results.setAlternatingRowColors(True)
        self.results.setHorizontalHeaderLabels(["Metric", "Value"])
        self.results.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.results.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeToContents
        )
        tabs.addTab(analysis, "Analysis")
        tabs.addTab(self.results, "Results")
        tray.addWidget(tabs)
        tray.setStretchFactor(0, 3)
        tray.setStretchFactor(1, 7)
        tray.setSizes([360, 840])
        return tray

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
        self.layer_tree.clear()
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
            self._loading_project = True
            if project.structure is not None:
                self.load_path(project.structure.source_path)
                self._restore_project_groups()
                self._interfaces = {
                    interface.name: (interface.group_a, interface.group_b)
                    for interface in project.interfaces
                }
                self._refresh_groups_panel()
            else:
                self.source = None
                self.current_result = None
                self._loaded_results.clear()
                self._structure_states.clear()
                self.loaded_structures.clear()
                self.viewer.clear_result()
                self._groups.clear()
                self._interfaces.clear()
                self.layer_tree.clear()
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
                build_surfaces=self.build_surfaces.isChecked(),
                build_vertices=self.build_vertices.isChecked(),
                build_edges=self.build_edges.isChecked(),
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

    def _display_result(self, result: AnalysisResult) -> None:
        previous = self.current_result
        same_structure = (
            previous is not None
            and previous.source == result.source
            and len(previous.atoms) == len(result.atoms)
        )
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
        self._populate_layer_tree(result)
        self._populate_structure_browser()
        self._update_running_selection()
        metrics = [
            ("Structure", result.name),
            ("Atoms", f"{len(result.atoms):,}"),
            ("Bonds", f"{len(result.bonds):,}"),
            ("Geometry layers", str(len(result.layers))),
            ("Complete cells", f"{result.complete_cells:,}"),
            ("Surfaces", f"{result.surface_count:,}"),
        ]
        self.results.setRowCount(len(metrics))
        for row, (label, value) in enumerate(metrics):
            self.results.setItem(row, 0, QTableWidgetItem(label))
            self.results.setItem(row, 1, QTableWidgetItem(value))
        values = {
            "atoms": len(result.atoms),
            "bonds": len(result.bonds),
            "cells": result.complete_cells,
            "surfaces": result.surface_count,
            "layers": len(result.layers),
        }
        for key, value in values.items():
            self.metric_cards[key].value.setText(f"{value:,}")
        self.progress_label.setText("Complete")
        self.progress_bar.setValue(100)
        self.statusBar().showMessage(f"{result.name} ready")

    def _populate_layer_tree(self, result: AnalysisResult) -> None:
        self.layer_tree.blockSignals(True)
        self.layer_tree.clear()
        for layer in result.layers:
            item = QTreeWidgetItem([layer.name])
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(0, Qt.Checked if layer.visible else Qt.Unchecked)
            item.setData(0, self.LAYER_ROLE, layer.name)
            self.layer_tree.addTopLevelItem(item)
        self.layer_tree.blockSignals(False)

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
        groups: dict[str, list[int]] = {}
        for atom in result.atoms:
            groups.setdefault(atom.chain or "(blank)", []).append(atom.index)
        return [
            (f"Chain {chain} ({len(indices)} atoms)", tuple(indices))
            for chain, indices in groups.items()
        ]

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

    def _tree_item_changed(self, item: QTreeWidgetItem, column: int) -> None:
        layer_name = item.data(column, self.LAYER_ROLE)
        if layer_name:
            self.viewer.set_layer_visible(
                layer_name, item.checkState(column) == Qt.Checked
            )

    def _change_selected_layer_opacity(self, value: int) -> None:
        item = self.layer_tree.currentItem()
        if item is not None:
            layer_name = item.data(0, self.LAYER_ROLE)
            if layer_name:
                self.viewer.set_layer_opacity(layer_name, value / 100.0)

    def _change_selected_layer_color(self) -> None:
        item = self.layer_tree.currentItem()
        if item is None:
            return
        layer_name = item.data(0, self.LAYER_ROLE)
        if not layer_name:
            return
        color = QColorDialog.getColor(QColor("#6857d9"), self, "Choose layer color")
        if color.isValid():
            self.viewer.set_layer_color(layer_name, color.name())

    def _show_selected_atom(self, atom: Atom) -> None:
        self._running_selection = {atom.index}
        self._update_running_selection()
        self._populate_structure_browser()
        self.selection_group.setTitle("Selected atom")
        self.atom_name.setText(f"{atom.name} (#{atom.serial})")
        self.atom_element.setText(atom.element)
        self.atom_residue.setText(
            f"{atom.residue_name} {atom.residue_sequence}".strip() or "—"
        )
        self.atom_chain.setText(atom.chain or "—")
        self.atom_position.setText(", ".join(f"{value:.3f}" for value in atom.position))
        self.selection_count.setText("1")
        self.statusBar().showMessage(
            f"Selected {atom.name}, {atom.residue_name} {atom.residue_sequence}"
        )

    def _show_selected_residue(self, atoms: list[Atom]) -> None:
        if not atoms:
            return
        # Mouse residue picks and Structure-browser residue picks share the
        # same temporary selection, highlight, and group-creation path.
        self._running_selection = {atom.index for atom in atoms}
        self._update_running_selection()
        self._populate_structure_browser()
        atom = atoms[0]
        self.selection_group.setTitle("Selected residue")
        self.atom_name.setText("Multiple")
        self.atom_element.setText("—")
        self.atom_residue.setText(
            f"{atom.residue_name} {atom.residue_sequence}".strip()
        )
        self.atom_chain.setText(atom.chain or "—")
        self.atom_position.setText("—")
        self.selection_count.setText(str(len(atoms)))
        self.statusBar().showMessage(
            f"Selected {atom.residue_name} {atom.residue_sequence}, chain "
            f"{atom.chain or '—'} ({len(atoms)} atoms)"
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
                f"{mode.title()} selection active — each click replaces the selection"
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
