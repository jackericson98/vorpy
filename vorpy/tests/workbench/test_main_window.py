import os
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QApplication, QToolBar, QWidget
from vtkmodules.vtkCommonDataModel import vtkPlane

from vorpy.workbench.domain import AnalysisResult, Atom, Bond, GeometryLayer
from vorpy.workbench.ui import main_window
from vorpy.workbench.ui.molecular_view import MolecularView


class PlotterStub:
    def reset_camera(self):
        pass

    def screenshot(self, _filename):
        pass


class ViewerStub(QWidget):
    loading_progress = Signal(str, int)
    selected_atom = Signal(object, bool)
    selected_residue = Signal(object, bool)
    selected_chain = Signal(object, bool)
    selected_molecule = Signal(object, bool)
    selection_cleared = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.plotter = PlotterStub()
        self.calls = []
        self.result = None

    def display_result(self, result, **kwargs):
        self.result = result
        self.calls.append(("display", kwargs.get("preserve_camera", False)))

    def clear_result(self):
        self.result = None
        self.calls.append(("clear-result",))

    def set_selection_mode(self, mode):
        self.calls.append(("selection", mode))

    def set_cartoon_visible(self, visible):
        self.calls.append(("cartoon", visible))

    def set_category_visible(self, category, visible):
        self.calls.append((category, visible))

    def set_bonds_visible(self, visible):
        self.calls.append(("bonds", visible))

    def set_group_selection(self, atoms):
        self.calls.append(("group-selection", tuple(atom.index for atom in atoms)))

    def set_molecule_opacity(self, opacity):
        self.calls.append(("molecule-opacity", opacity))

    def set_waters_visible(self, visible):
        self.calls.append(("waters", visible))

    def set_water_style(self, style):
        self.calls.append(("water-style", style))

    def set_water_opacity(self, opacity):
        self.calls.append(("water-opacity", opacity))

    def set_ions_visible(self, visible):
        self.calls.append(("ions", visible))

    def set_ion_opacity(self, opacity):
        self.calls.append(("ion-opacity", opacity))

    def set_layer_visible(self, name, visible):
        self.calls.append((f"layer:{name}", visible))

    def set_layer_opacity(self, name, opacity):
        self.calls.append((f"opacity:{name}", opacity))

    def set_layer_color(self, name, color):
        self.calls.append((f"color:{name}", color))

    def set_layer_color_scheme(self, name, scheme):
        self.calls.append((f"scheme:{name}", scheme))

    def reset_depth_clipping(self):
        self.calls.append(("reset-depth",))

    def set_depth_clipping_fraction(self, fraction):
        self.calls.append(("depth-clipping", fraction))

    def save_screenshot(self, filename, scale=1):
        self.calls.append(("screenshot", filename, scale))


def make_window(monkeypatch):
    QApplication.instance() or QApplication([])
    monkeypatch.setattr(main_window, "MolecularView", ViewerStub)
    return main_window.MainWindow()


def sample_result():
    atoms = [
        Atom(0, 1, "N", "N", (0.0, 1.0, 2.0), "GLY", "7", "A"),
        Atom(1, 2, "CA", "C", (1.0, 1.0, 2.0), "GLY", "7", "A"),
    ]
    return AnalysisResult(
        source=None,
        name="tiny",
        atoms=atoms,
        bonds=[Bond(0, 1)],
        layers=[GeometryLayer("edges", "edges", visible=True)],
        complete_cells=2,
        surface_count=3,
    )


def result_metric(window, label):
    """Return an overview metric by label rather than fragile row position."""
    for row in range(window.results.rowCount()):
        label_item = window.results.item(row, 0)
        if label_item is not None and label_item.text() == label:
            value_item = window.results.item(row, 1)
            return value_item.text() if value_item is not None else None
    raise AssertionError(f"Overview metric {label!r} was not found")


def test_action_state_and_visibility_controls(monkeypatch):
    window = make_window(monkeypatch)

    assert window.select_residue_action.isChecked()
    assert ("selection", "residue") in window.viewer.calls
    assert not window.select_residue_action.shortcut().toString()
    assert window.reset_selection_action.shortcut().toString() == "R"

    assert [window.workflow_tabs.tabText(i) for i in range(window.workflow_tabs.count())] == [
        "Structure", "Selection", "Groups", "Interfaces"
    ]
    assert window.inspector.tabText(0) == "System"
    assert window.inspector.tabText(1) == "Network"

    assert not window.solve_action.isEnabled()
    assert window.inspector.currentIndex() == 0
    assert not window.cancel_action.isEnabled()
    assert window.show_cartoon.isChecked()
    assert not window.show_spheres.isChecked()
    assert not window.show_sticks.isChecked()
    assert not window.show_waters.isChecked()
    assert window.show_ions.isChecked()
    assert window.water_style.currentData() == "ball-and-stick"
    assert window.water_opacity.value() == 25
    assert window.molecule_opacity.value() == 50
    assert window.ion_opacity.value() == 75

    window.show_spheres.setChecked(True)
    window.show_sticks.setChecked(True)
    assert ("atoms", True) in window.viewer.calls
    assert ("bonds", True) in window.viewer.calls

    window.molecule_opacity.setValue(65)
    assert ("molecule-opacity", 0.65) in window.viewer.calls

    window.show_waters.setChecked(True)
    window.show_waters.setChecked(False)
    window.show_ions.setChecked(True)
    window.show_ions.setChecked(False)
    assert ("waters", False) in window.viewer.calls
    assert ("ions", False) in window.viewer.calls

    window.water_style.setCurrentIndex(2)
    window.water_opacity.setValue(35)
    assert ("water-style", "spheres") in window.viewer.calls
    assert ("water-opacity", 0.35) in window.viewer.calls
    window.ion_opacity.setValue(40)
    assert ("ion-opacity", 0.4) in window.viewer.calls

    window.select_atom_action.setChecked(True)
    window.select_residue_action.setChecked(True)
    assert not window.select_atom_action.isChecked()
    assert window.select_residue_action.isChecked()
    assert window.viewer.calls[-1] == ("selection", "residue")
    assert window.selection_mode_label.text() == "Residue selection active"

    window.select_chain_action.setChecked(True)
    assert window.viewer.calls[-1] == ("selection", "chain")
    assert window.selection_mode_label.text() == "Chain selection active"

    window.select_molecule_action.setChecked(True)
    assert window.viewer.calls[-1] == ("selection", "molecule")
    assert window.selection_mode_label.text() == "Molecule selection active"


def test_bottom_tray_and_small_molecule_representation_defaults(monkeypatch):
    window = make_window(monkeypatch)
    assert window.workspace_splitter.orientation() == Qt.Vertical
    assert window.workspace_splitter.widget(0) is window.upper_workspace
    assert window.workspace_splitter.widget(1) is window.lower_workspace
    assert window.upper_workspace.widget(0) is window.viewer_panel
    assert window.upper_workspace.widget(1) is window.view_settings_panel
    assert window.workflow_panel.parentWidget() is window.lower_workspace
    assert window.analysis_tray.parentWidget() is window.lower_workspace
    assert window.max_vertices.isVisibleTo(window)
    assert window.box_size.isVisibleTo(window)
    assert window.surface_resolution.isVisibleTo(window)
    assert window.solve_network_button.text() == "Solve network"
    assert window.solve_network_button.objectName() == "primaryAction"
    assert not hasattr(window, "build_settings_button")
    assert not hasattr(window, "build_vertices")
    assert not hasattr(window, "build_edges")
    assert not hasattr(window, "build_surfaces")
    window.solve_action.setEnabled(False)
    assert not window.solve_network_button.isEnabled()
    window.solve_action.setEnabled(True)
    assert window.solve_network_button.isEnabled()
    assert [window.analysis_tray_tabs.tabText(i) for i in range(window.analysis_tray_tabs.count())] == [
        "Overview", "Composition", "Build", "Timing", "Network / Topology",
        "Group Geometry", "Surface Curvature", "Surface Energy",
        "Surface Classification", "Chain Composition", "Residue Composition", "Interfaces",
    ]
    assert window.solve_target.count() == 2
    assert not window.findChildren(QToolBar)

    small = sample_result()
    window._display_result(small)
    assert window.show_cartoon.isChecked()
    assert window.show_sticks.isChecked()

    large = sample_result()
    large.atoms.extend(
        Atom(i, i + 1, "C", "C", (float(i), 0.0, 0.0), "LIG", str(i), "A")
        for i in range(2, 1502)
    )
    window._display_result(large)
    assert not window.show_sticks.isChecked()


def test_result_and_selection_populate_inspector(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()

    window._display_result(result)
    assert window.viewer.result is result
    assert window.structure_name.text() == "tiny"
    assert window.structure_summary.text() == (
        "2 atoms  •  1 bonds  •  1 residues  •  1 chains"
    )
    assert window.metric_cards["atoms"].value.text() == "2"
    assert window.metric_cards["cells"].value.text() == "2"
    assert result_metric(window, "Surfaces") == "3"

    assert window.network_layer_checks["edges"].isEnabled()
    assert window.network_layer_checks["edges"].isChecked()
    assert not window.network_layer_checks["vertices"].isEnabled()
    window.network_layer_checks["edges"].setChecked(False)
    assert ("layer:edges", False) in window.viewer.calls

    window._show_selected_atom(result.atoms[1])
    assert window.selection_group.title() == "Selected atom"
    assert window.atom_name.text() == "CA (#2)"
    assert window.atom_residue.text() == "GLY 7"
    assert window.atom_position.text() == "1.000, 1.000, 2.000"

    window._show_selected_residue(result.atoms)
    assert window.selection_group.title() == "Selected residue"
    assert window.selection_count.text() == "2"
    assert window._running_selection == {0, 1}
    assert ("group-selection", (0, 1)) in window.viewer.calls


def test_network_controls_manage_categories_color_and_surface_opacity(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    result.layers = [
        GeometryLayer("network edges", "edges", color="#111111"),
        GeometryLayer("network vertices", "vertices", color="#222222"),
        GeometryLayer("network surfaces", "surfaces", color="#333333", opacity=0.6, cell_scalars={"mean_curvature": object()}),
        GeometryLayer("shell edges", "edges", color="#444444"),
        GeometryLayer("shell vertices", "vertices", color="#555555"),
        GeometryLayer("shell surfaces", "surfaces", color="#666666", opacity=0.6, cell_scalars={"mean_curvature": object()}),
    ]

    window._display_result(result)

    assert set(window.network_layer_checks) == {
        "edges", "vertices", "surfaces",
        "shell_edges", "shell_vertices", "shell_surfaces",
    }
    assert all(control.isEnabled() for control in window.network_layer_checks.values())
    assert window.surface_opacity.isEnabled()
    assert window.surface_opacity.value() == 60
    assert window.surface_color_scheme.isEnabled()
    window.surface_color_scheme.setCurrentIndex(
        window.surface_color_scheme.findData("mean_curvature")
    )
    assert ("scheme:network surfaces", "mean_curvature") in window.viewer.calls
    assert ("scheme:shell surfaces", "mean_curvature") in window.viewer.calls

    window.network_layer_checks["shell_edges"].setChecked(False)
    assert ("layer:shell edges", False) in window.viewer.calls

    window.surface_opacity.setValue(35)
    assert ("opacity:network surfaces", 0.35) in window.viewer.calls
    assert ("opacity:shell surfaces", 0.35) in window.viewer.calls

    monkeypatch.setattr(
        main_window.QColorDialog,
        "getColor",
        lambda *args, **kwargs: QColor("#abcdef"),
    )
    window._choose_network_color("vertices")
    assert ("color:network vertices", "#abcdef") in window.viewer.calls
    assert result.layers[1].color == "#abcdef"


def test_shift_selection_toggles_and_normal_selection_replaces(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    window._display_result(result)

    window._show_selected_atom(result.atoms[0])
    window._show_selected_atom(result.atoms[1], additive=True)
    assert window._running_selection == {0, 1}
    assert window.selection_count.text() == "2"
    assert window.selection_group.title() == "Selected atoms"

    window._show_selected_atom(result.atoms[1], additive=True)
    assert window._running_selection == {0}
    assert window.selection_count.text() == "1"
    assert window.statusBar().currentMessage().startswith("Removed CA")

    window._show_selected_atom(result.atoms[1], additive=True)
    assert window._running_selection == {0, 1}

    window._show_selected_atom(result.atoms[0])
    assert window._running_selection == {0}
    assert window.selection_count.text() == "1"


def test_reset_selection_requires_confirmation(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    window._display_result(result)
    window._show_selected_residue(result.atoms)

    monkeypatch.setattr(
        main_window.QMessageBox,
        "question",
        lambda *_args, **_kwargs: main_window.QMessageBox.No,
    )
    window.reset_selection()
    assert window._running_selection == {0, 1}

    monkeypatch.setattr(
        main_window.QMessageBox,
        "question",
        lambda *_args, **_kwargs: main_window.QMessageBox.Yes,
    )
    window.reset_selection()
    assert window._running_selection == set()
    assert window.running_selection.text() == "No atoms selected"
    assert window.selection_group.title() == "No selection"
    assert window.statusBar().currentMessage() == "Selection reset"
    assert ("group-selection", ()) in window.viewer.calls


def test_clicking_viewer_background_clears_selection(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    window._display_result(result)
    window._show_selected_residue(result.atoms)

    window.viewer.selection_cleared.emit()

    assert window._running_selection == set()
    assert window.running_selection.text() == "No atoms selected"
    assert window.selection_group.title() == "No selection"
    assert window.statusBar().currentMessage() == "Selection cleared"
    assert ("group-selection", ()) in window.viewer.calls


def test_empty_viewer_pick_emits_selection_cleared():
    events = []
    viewer = SimpleNamespace(
        _pending_pick=None,
        _selection_additive=False,
        _result=sample_result(),
        _selection_mode="residue",
        _positions=[(0.0, 0.0, 0.0)],
        _selection_dragged=False,
        _clear_pick_highlight=lambda: events.append("highlight-cleared"),
        selection_cleared=SimpleNamespace(
            emit=lambda: events.append("selection-cleared")
        ),
    )

    MolecularView._apply_pending_pick(viewer)

    assert events == ["highlight-cleared", "selection-cleared"]


def test_depth_clipping_advances_and_restores_all_rendered_actors():
    class Mapper:
        def __init__(self):
            self.planes = []

        def RemoveAllClippingPlanes(self):
            self.planes.clear()

        def AddClippingPlane(self, plane):
            self.planes.append(plane)

    mapper = Mapper()
    actor = SimpleNamespace(
        GetMapper=lambda: mapper,
        GetBounds=lambda: (-2.0, 2.0, -1.0, 1.0, -3.0, 3.0),
    )
    label = SimpleNamespace(
        setText=lambda text: None,
        show=lambda: None,
        hide=lambda: None,
    )
    viewer = SimpleNamespace(
        _depth_clip_fraction=0.0,
        _depth_clip_plane=vtkPlane(),
        depth_clip_label=label,
        plotter=SimpleNamespace(
            camera=SimpleNamespace(
                GetDirectionOfProjection=lambda: (0.0, 0.0, -1.0)
            ),
            render=lambda: None,
        ),
        _rendered_actors=lambda: [actor],
    )
    viewer._apply_depth_clip_to_actor = lambda item: (
        MolecularView._apply_depth_clip_to_actor(viewer, item)
    )
    viewer._depth_projection_range = lambda direction: (
        MolecularView._depth_projection_range(viewer, direction)
    )
    viewer._apply_depth_clipping = lambda render=True: (
        MolecularView._apply_depth_clipping(viewer, render)
    )

    MolecularView.adjust_depth_clipping(viewer, 5)
    assert viewer._depth_clip_fraction == 0.2
    assert mapper.planes == [viewer._depth_clip_plane]

    MolecularView.adjust_depth_clipping(viewer, -5)
    assert viewer._depth_clip_fraction == 0.0
    assert mapper.planes == []


def test_ion_classification_uses_residue_identity():
    alpha_carbon = Atom(0, 1, "CA", "C", (0.0, 0.0, 0.0), "GLY", "7", "A")
    calcium = Atom(1, 2, "CA", "CA", (1.0, 0.0, 0.0), "CA", "8", "A")

    assert not MolecularView._is_ion(alpha_carbon)
    assert MolecularView._is_ion(calcium)


def test_cartoon_builds_nucleic_trace_with_terminal_sugar_fallback():
    phosphate = Atom(0, 1, "P", "P", (0.0, 0.0, 0.0), "DA", "1", "A")
    sugar_same_residue = Atom(1, 2, "C4'", "C", (1.0, 0.0, 0.0), "DA", "1", "A")
    terminal_sugar = Atom(2, 3, "C4'", "C", (6.0, 0.0, 0.0), "DT", "2", "A")
    segments = []

    viewer = SimpleNamespace(
        _is_water=MolecularView._is_water,
        _is_ion=MolecularView._is_ion,
        _residue_sort_key=MolecularView._residue_sort_key,
        _add_cartoon_segment=lambda atoms, color: segments.append((atoms, color)),
    )
    viewer._add_cartoon_chain = lambda atoms, color, maximum_gap: (
        MolecularView._add_cartoon_chain(viewer, atoms, color, maximum_gap)
    )

    MolecularView._add_cartoon(
        viewer, [phosphate, sugar_same_residue, terminal_sugar]
    )

    assert len(segments) == 1
    assert [atom.serial for atom in segments[0][0]] == [1, 3]


def test_screenshot_uses_selected_resolution(monkeypatch):
    window = make_window(monkeypatch)
    monkeypatch.setattr(
        main_window.QFileDialog,
        "getSaveFileName",
        lambda *args: ("/tmp/view.png", "PNG images (*.png)"),
    )
    monkeypatch.setattr(
        main_window.QInputDialog,
        "getInt",
        lambda *args, **kwargs: (4, True),
    )

    window.save_screenshot()

    assert ("screenshot", "/tmp/view.png", 4) in window.viewer.calls


def test_screenshot_temporarily_hides_orientation_axes():
    class AxesStub:
        visible = True

        def GetVisibility(self):
            return self.visible

        def SetVisibility(self, visible):
            self.visible = visible

    class PlotterStub:
        def __init__(self, axes):
            self.axes = axes
            self.captured = None
            self.rendered = False

        def screenshot(self, filename, scale):
            self.captured = (filename, scale, self.axes.visible)

        def render(self):
            self.rendered = True

    axes = AxesStub()
    plotter = PlotterStub(axes)
    viewer = SimpleNamespace(_axes_actor=axes, plotter=plotter)

    MolecularView.save_screenshot(viewer, "high-resolution.png", 6)

    assert plotter.captured == ("high-resolution.png", 6, False)
    assert axes.visible
    assert plotter.rendered


def test_structure_browser_builds_running_selection_and_group(monkeypatch):
    window = make_window(monkeypatch)
    window._display_result(sample_result())

    assert window.structure_browser_type.currentText() == "Atoms"
    assert window.structure_browser.count() == 2
    assert not window.make_group_button.isEnabled()

    window.structure_browser.item(0).setCheckState(Qt.Checked)
    assert window._running_selection == {0}
    assert ("group-selection", (0,)) in window.viewer.calls
    assert window.running_selection.text().startswith("1 atoms · 1 residues · 1 chains")
    assert window.make_group_button.isEnabled()

    window.structure_browser_type.setCurrentText("Residues")
    assert window.structure_browser.count() == 1
    assert window.structure_browser.item(0).checkState() == Qt.PartiallyChecked

    window.structure_search.setText("gly 7")
    assert not window.structure_browser.item(0).isHidden()
    window.structure_search.setText("missing")
    assert window.structure_browser.item(0).isHidden()

    monkeypatch.setattr(
        main_window.QInputDialog,
        "getText",
        lambda *args, **kwargs: ("Active site", True),
    )
    window._make_group()

    assert window._groups["Active site"] == (0,)
    assert window.groups_list.count() == 1
    assert window.groups_list.item(0).text() == "Active site (1 atoms)"
    assert window.inspector.tabText(0) == "System"
    assert window.inspector.tabText(1) == "Network"
    assert window.groups_list.item(0).text() == "Active site (1 atoms)"

    window._select_saved_group(window.groups_list.item(0))
    assert window._running_selection == {0}
    assert ("group-selection", (0,)) in window.viewer.calls


def test_interface_section_connects_two_saved_groups(monkeypatch):
    window = make_window(monkeypatch)
    window._display_result(sample_result())
    window._groups = {"A": (0,), "B": (1,)}
    window._refresh_groups_panel()
    window.interface_group_a.setCurrentText("A")
    window.interface_group_b.setCurrentText("B")
    monkeypatch.setattr(
        main_window.QInputDialog, "getText", lambda *args, **kwargs: ("Contact", True)
    )

    window._make_interface()

    assert window._interfaces == {"Contact": ("A", "B")}
    assert window.interfaces_list.item(0).text() == "Contact: A ↔ B"


def test_multiple_structures_switch_with_independent_selection_and_statistics(monkeypatch, tmp_path):
    first = tmp_path / "first.pdb"
    second = tmp_path / "second.pdb"
    first.write_text("END\n", encoding="utf-8")
    second.write_text("END\n", encoding="utf-8")
    first_result = sample_result()
    first_result.source = first
    second_result = sample_result()
    second_result.source = second
    second_result.name = "second"
    second_result.complete_cells = 99
    monkeypatch.setattr(
        main_window, "load_pdb", lambda path, **kwargs: {first: first_result, second: second_result}[path]
    )
    window = make_window(monkeypatch)

    window.load_path(first)
    window._running_selection = {1}
    window._groups = {"First group": (1,)}
    window.load_path(second)
    assert window.source == second.resolve()
    assert window._running_selection == set()
    assert window.metric_cards["cells"].value.text() == "99"

    window._switch_structure_item(window.loaded_structures.item(0))
    assert window.source == first.resolve()
    assert window._running_selection == {1}
    assert window._groups == {"First group": (1,)}
    assert window.metric_cards["cells"].value.text() == "2"


def test_molecule_browser_uses_bond_connected_components(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    result.atoms.append(Atom(2, 3, "O", "O", (4.0, 0.0, 0.0), "HOH", "8", ""))
    window._display_result(result)

    window.structure_browser_type.setCurrentText("Molecules")

    assert window.structure_browser.count() == 2
    assert window.structure_browser.item(0).data(Qt.UserRole) == (0, 1)
    assert window.structure_browser.item(1).data(Qt.UserRole) == (2,)


def test_viewer_molecule_selection_uses_bond_connected_component():
    result = sample_result()
    isolated = Atom(2, 3, "O", "O", (4.0, 0.0, 0.0), "HOH", "8", "")
    result.atoms.append(isolated)
    viewer = SimpleNamespace(_result=result)

    bonded = MolecularView._molecule_atoms(viewer, result.atoms[0])
    separate = MolecularView._molecule_atoms(viewer, isolated)

    assert [atom.index for atom in bonded] == [0, 1]
    assert [atom.index for atom in separate] == [2]


def test_blank_chain_selection_keeps_water_residues_separate():
    atoms = [
        Atom(0, 1, "CA", "C", (0.0, 0.0, 0.0), "GLY", "1", ""),
        Atom(1, 2, "CB", "C", (0.5, 0.0, 0.0), "ALA", "2", ""),
        Atom(2, 3, "OW", "O", (4.0, 0.0, 0.0), "SOL", "47", ""),
        Atom(3, 4, "HW1", "H", (4.5, 0.0, 0.0), "SOL", "47", ""),
        Atom(4, 5, "OW", "O", (8.0, 0.0, 0.0), "SOL", "48", ""),
    ]
    result = AnalysisResult(source=None, name="waters", atoms=atoms)
    viewer = SimpleNamespace(
        _result=result,
        _residue_key=MolecularView._residue_key,
        _is_water=MolecularView._is_water,
    )

    selected_protein = MolecularView._chain_atoms(viewer, atoms[0])
    selected_water = MolecularView._chain_atoms(viewer, atoms[2])
    groups = main_window.MainWindow._chain_groups(result)

    assert [atom.index for atom in selected_protein] == [0, 1]
    assert [atom.index for atom in selected_water] == [2, 3]
    assert groups == [
        ("No chain · non-water atoms (2 atoms)", (0, 1)),
        ("No chain · SOL 47 (2 atoms)", (2, 3)),
        ("No chain · SOL 48 (1 atom)", (4,)),
    ]


def test_groups_survive_analysis_result_for_same_structure(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    window._display_result(result)
    window._running_selection = {0, 1}
    window._groups["Protein"] = (0, 1)

    replacement = sample_result()
    window._display_result(replacement)

    assert window._running_selection == {0, 1}
    assert window._groups == {"Protein": (0, 1)}


def test_group_highlight_replaces_actor_without_resetting_camera(monkeypatch):
    class CloudStub:
        def __setitem__(self, _key, _value):
            pass

        def glyph(self, **_kwargs):
            return "highlight-mesh"

    class GroupPlotterStub:
        def __init__(self):
            self.removed = []
            self.mesh_kwargs = None
            self.rendered = False

        def remove_actor(self, name, render):
            self.removed.append((name, render))

        def add_mesh(self, _mesh, **kwargs):
            self.mesh_kwargs = kwargs

        def render(self):
            self.rendered = True

    monkeypatch.setattr(
        "vorpy.workbench.ui.molecular_view.pv.PolyData", lambda _points: CloudStub()
    )
    monkeypatch.setattr(
        "vorpy.workbench.ui.molecular_view.pv.Sphere", lambda **_kwargs: "sphere"
    )
    plotter = GroupPlotterStub()
    viewer = SimpleNamespace(plotter=plotter)

    MolecularView.set_group_selection(viewer, sample_result().atoms)

    assert plotter.removed == [("group-selection", False)]
    assert plotter.mesh_kwargs["name"] == "group-selection"
    assert plotter.mesh_kwargs["reset_camera"] is False
    assert plotter.mesh_kwargs["pickable"] is False
    assert plotter.rendered


def test_open_pdb_starts_in_packaged_data_directory(monkeypatch):
    window = make_window(monkeypatch)
    captured = {}

    def choose_file(_parent, _title, directory, _filters):
        captured["directory"] = directory
        return "", ""

    monkeypatch.setattr(main_window.QFileDialog, "getOpenFileName", choose_file)

    window.open_structure()

    assert captured["directory"] == str(main_window.DEFAULT_DATA_DIRECTORY)
    assert main_window.DEFAULT_DATA_DIRECTORY.name == "data"
    assert main_window.DEFAULT_DATA_DIRECTORY.parent.name == "vorpy"


def test_project_save_and_open_restores_groups(monkeypatch, tmp_path):
    pdb = tmp_path / "tiny.pdb"
    pdb.write_text(
        "ATOM      1  N   GLY A   7       0.000   0.000   0.000  1.00  0.00           N  \n"
        "ATOM      2  CA  GLY A   7       1.450   0.000   0.000  1.00  0.00           C  \n"
        "CONECT    1    2\nEND\n",
        encoding="utf-8",
    )
    project_file = tmp_path / "saved.vpyworkbench.json"
    window = make_window(monkeypatch)
    window.load_path(pdb)
    window._groups["Backbone"] = (1,)
    window._set_project_dirty()
    window.project_file = project_file

    assert window.save_project()
    assert not window._project_dirty

    restored = make_window(monkeypatch)
    monkeypatch.setattr(
        main_window.QFileDialog,
        "getOpenFileName",
        lambda *args: (str(project_file), ""),
    )
    restored.open_project()

    assert restored.source == pdb.resolve()
    assert restored._groups == {"Backbone": (1,)}
    assert not restored._project_dirty
    assert restored.project_file == project_file.resolve()


def test_solve_target_resolves_group_and_interface_membership(monkeypatch):
    window = make_window(monkeypatch)
    window._display_result(sample_result())
    window._groups = {"A": (0,), "B": (1,)}
    window._interfaces = {"Contact": ("A", "B")}
    window._refresh_groups_panel()

    window.solve_target.setCurrentIndex(window.solve_target.findText("Group: B"))
    assert window._solve_target_indices() == (1,)
    window.solve_target.setCurrentIndex(window.solve_target.findText("Interface: Contact"))
    assert window._solve_target_indices() == (0, 1)
    window.solve_target.setCurrentIndex(window.solve_target.findText("Whole molecule"))
    assert window._solve_target_indices() is None


def test_workbench_empty_loaded_busy_and_reset_states(monkeypatch):
    window = make_window(monkeypatch)
    assert not window.solve_network_button.isEnabled()
    assert not window.system_display.isEnabled()
    assert window.analysis_empty.isVisibleTo(window)
    assert not window.make_interface_button.isEnabled()
    assert window.lower_workspace.widget(1).isAncestorOf(window.solve_network_button)
    assert "0 selected" in window.state_summary.text()

    result = sample_result()
    window._display_result(result)
    assert window.solve_network_button.isEnabled()
    assert not window.analysis_empty.isVisibleTo(window)
    assert window.network_layer_checks["edges"].isEnabled()
    window._show_selected_atom(result.atoms[0])
    assert "1 selected" in window.state_summary.text()
    window.inspector.setCurrentIndex(0)
    window.select_atom_action.setChecked(True)
    assert window.inspector.currentIndex() == 0

    window._thread = object()
    window._refresh_ui_state()
    assert not window.solve_network_button.isEnabled()
    assert not window.open_action.isEnabled()
    assert window.progress_bar.isVisibleTo(window)
    window._thread_finished()
    assert window.solve_network_button.isEnabled()

    monkeypatch.setattr(window, "_confirm_discard_changes", lambda: True)
    window.new_project()
    assert not window.solve_network_button.isEnabled()
    assert window.results.rowCount() == 0
    assert window.analysis_empty.isVisibleTo(window)
    assert not window.network_layer_checks["edges"].isEnabled()
    assert "0 selected" in window.state_summary.text()


def test_layout_state_round_trip_and_results_collapse(monkeypatch):
    window = make_window(monkeypatch)
    window.show()
    QApplication.processEvents()
    window.upper_workspace.setSizes([1000, 400])
    window.lower_workspace.setSizes([450, 350, 600])
    window.analysis_tray.set_expanded(False)
    state = window._view_state_to_json()
    assert not state["results_expanded"]
    window.analysis_tray.set_expanded(True)
    window._restore_view_state(state)
    assert not window.analysis_tray.content.isVisibleTo(window)
    assert window.upper_workspace.sizes() == state["panel_sizes"]
    assert window.lower_workspace.sizes() == state["bottom_panel_sizes"]
    window._focus_results()
    assert window.analysis_tray.content.isVisibleTo(window)
    window.close()


def test_dragged_pick_preserves_selection():
    import numpy as np
    emitted = []
    viewer = SimpleNamespace(
        _pending_pick=np.array([0., 1., 2.]), _selection_additive=False,
        _result=sample_result(), _selection_mode="atom",
        _positions=np.array([[0., 1., 2.], [1., 1., 2.]]),
        _selection_dragged=True,
        selected_atom=SimpleNamespace(emit=lambda *args: emitted.append(args)),
        selection_cleared=SimpleNamespace(emit=lambda: emitted.append("cleared")),
    )
    MolecularView._apply_pending_pick(viewer)
    assert emitted == []


def test_atomic_defaults_update_loaded_atoms_and_mark_results_outdated(monkeypatch):
    window = make_window(monkeypatch)
    window._display_result(sample_result())
    window.atomic_defaults = {'C': {'radius': 2.0, 'mass': 13., 'charge': -0.25}}
    window._apply_project_atomic_defaults()
    carbon = window.current_result.atoms[1]
    assert (carbon.radius, carbon.mass, carbon.charge) == (2., 13., -0.25)
    assert window.current_result.defaults_stale
    assert 'outdated' in window.state_summary.text()
    assert window.solve_hint.isVisibleTo(window)
    state = window._view_state_to_json()
    assert state['atomic_defaults']['C']['charge'] == -0.25
    window.atomic_defaults = {}
    window._apply_project_atomic_defaults()
    assert window.current_result.atoms[1].radius == 0.35
    window._restore_view_state(state)
    assert window.current_result.atoms[1].mass == 13.
    monkeypatch.setattr(window, '_confirm_discard_changes', lambda: True)
    window.new_project()
    assert window.atomic_defaults == {}


def test_filtered_browser_switch_preserves_selection(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    from dataclasses import replace
    result.atoms[0] = replace(result.atoms[0], residue_name="ILE")
    window._display_result(result)
    window.structure_browser.item(0).setCheckState(Qt.Checked)
    window.structure_search.setText("ILE")
    for category in ("Residues", "Chains", "Molecules", "Atoms"):
        window.structure_browser_type.setCurrentText(category)
        assert window._running_selection == {0}
        for row in range(window.structure_browser.count()):
            item = window.structure_browser.item(row)
            assert item.isHidden() == ("ile" not in item.text().casefold())
    window.close()


def test_atom_rendering_instances_spheres():
    from collections import defaultdict
    from vtkmodules.vtkRenderingCore import vtkGlyph3DMapper
    class Plotter:
        def add_actor(self, actor, **kwargs):
            pass
    viewer = SimpleNamespace(
        plotter=Plotter(), _actors=defaultdict(list), _water_opacity=0.3,
        _molecule_opacity=1.0, _waters_visible=False, _water_style="spheres",
        _ions_visible=True, _spheres_visible=True,
        _is_water=MolecularView._is_water, _is_ion=MolecularView._is_ion,
    )
    MolecularView._add_atoms(viewer, sample_result().atoms)
    actors = viewer._actors["atoms"]
    assert len(actors) == 2
    for actor in actors:
        mapper = actor.GetMapper()
        assert isinstance(mapper, vtkGlyph3DMapper)
        assert mapper.GetInput().GetNumberOfPoints() == 1
        assert mapper.GetScaleFactor() == 1.0


def test_loading_keeps_event_loop_alive_and_recovers_from_error(monkeypatch, tmp_path):
    from threading import Event, get_ident
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication
    window = make_window(monkeypatch)
    source = tmp_path / 'slow.pdb'
    source.write_text('END\n')
    gui_thread = get_ident()
    heartbeat = Event()
    observed = []

    def read(path, progress, preview=None):
        assert get_ident() != gui_thread
        progress('Reading test structure', 25)
        assert heartbeat.wait(3), 'GUI event loop stopped during reading'
        return sample_result()

    def tick():
        if window._loading_structure:
            observed.append(window.progress_bar.isVisibleTo(window))
            heartbeat.set()

    timer = QTimer(window)
    timer.timeout.connect(tick)
    timer.start(5)
    monkeypatch.setattr(main_window, 'load_pdb', read)
    window.load_path(source)
    timer.stop()
    assert heartbeat.is_set() and all(observed)
    assert window._load_dialog is None
    assert not window._loading_structure
    assert window.current_result.name == 'tiny'
    previous = window.current_result
    errors = []
    monkeypatch.setattr(window, '_show_error', errors.append)
    def fail(path, progress, preview=None):
        raise ValueError('Invalid structure')
    monkeypatch.setattr(main_window, 'load_pdb', fail)
    window.load_path(tmp_path / 'invalid.pdb')
    assert errors == ['Invalid structure']
    assert window.current_result is previous
    assert window.open_action.isEnabled()
    assert not window._loading_structure
    window._project_dirty = False
    window.close()


def test_molecule_groups_cover_disconnected_atoms_and_cycles():
    result = sample_result()
    result.atoms = [Atom(i, i + 1, 'O', 'O', (float(i), 0., 0.)) for i in range(10000)]
    result.bonds = [Bond(0, 1), Bond(1, 2), Bond(2, 0), Bond(5, 7)]
    groups = main_window.MainWindow._molecule_groups(result)
    assert groups[:4] == [('Molecule 1 (3 atoms)', (0, 1, 2)),
                          ('Molecule 2 (1 atoms)', (3,)),
                          ('Molecule 3 (1 atoms)', (4,)),
                          ('Molecule 4 (2 atoms)', (5, 7))]
    assert sorted(i for label, indices in groups for i in indices) == list(range(10000))


def test_hidden_bonds_are_built_once_when_enabled():
    from collections import defaultdict
    from types import MethodType
    import numpy as np
    import pyvista as pv
    class Plotter:
        def __init__(self):
            self.added = []
        def add_mesh(self, mesh, **kwargs):
            self.added.append(mesh)
            return pv.Actor()
        def render(self):
            pass
    result = sample_result()
    view = SimpleNamespace(plotter=Plotter(), _actors=defaultdict(list),
                           _pending_bond_groups={}, _ions_visible=False,
                           _waters_visible=False, _water_style='ball-and-stick',
                           _apply_depth_clip_to_actor=lambda actor: None)
    view._add_bond_group = MethodType(MolecularView._add_bond_group, view)
    view._ensure_bond_group = MethodType(MolecularView._ensure_bond_group, view)
    view._add_bond_group(np.asarray([atom.position for atom in result.atoms]),
                         result.atoms, result.bonds, 'bonds', False)
    assert not view.plotter.added
    assert 'bonds' in view._pending_bond_groups
    MolecularView.set_bonds_visible(view, True)
    assert len(view.plotter.added) == 2
    assert 'bonds' not in view._pending_bond_groups
    MolecularView.set_bonds_visible(view, False)
    MolecularView.set_bonds_visible(view, True)
    assert len(view.plotter.added) == 2
    assert all(actor.GetVisibility() for actor in view._actors['bonds'])


def test_preview_is_interactive_and_remaps_selection_on_completion(monkeypatch, tmp_path):
    from dataclasses import replace
    from threading import Event
    from PySide6.QtCore import QTimer
    from vorpy.workbench.domain import AnalysisResult
    source = tmp_path / 'solvated.pdb'
    source.write_text('END\n')
    primary = sample_result()
    primary.source = source
    primary.layers = []
    primary.complete_cells = 0
    full = AnalysisResult(source, 'full', atoms=[
        Atom(0, 10, 'O', 'O', (10., 0., 0.), 'HOH'),
        replace(primary.atoms[0], index=1),
        Atom(2, 11, 'O', 'O', (20., 0., 0.), 'HOH'),
        replace(primary.atoms[1], index=3),
    ], bonds=[Bond(1, 3)])
    interacted = Event()
    window = make_window(monkeypatch)
    observations = []
    def read(path, progress, preview):
        preview(primary)
        assert interacted.wait(3)
        return full
    def interact():
        if window._molecule_preview and window._load_dialog is None and not interacted.is_set():
            observations.append((window.solve_action.isEnabled(), window.save_project_action.isEnabled(),
                                 window.structure_browser.isEnabled(), source not in window._loaded_results))
            window.structure_browser.item(1).setCheckState(Qt.Checked)
            window._groups = {'Chosen': (0, 1)}
            window.show_sticks.setChecked(False)
            interacted.set()
    timer = QTimer(window)
    timer.timeout.connect(interact)
    timer.start(5)
    monkeypatch.setattr(main_window, 'load_pdb', read)
    window.load_path(source)
    timer.stop()
    assert observations == [(False, False, True, True)]
    assert window._running_selection == {3}
    assert window._groups == {'Chosen': (1, 3)}
    assert ('display', True) in window.viewer.calls
    assert not window.show_sticks.isChecked()
    assert window.current_result is full
    assert window._loaded_results[source] is full
    assert window.solve_action.isEnabled() and window.save_project_action.isEnabled()
    assert not window._molecule_preview
    window._project_dirty = False
    window.close()


def test_failed_background_load_restores_previous_structure(monkeypatch, tmp_path):
    from threading import Event
    from PySide6.QtCore import QTimer
    window = make_window(monkeypatch)
    old = sample_result()
    old.source = tmp_path / 'old.pdb'
    window.source = old.source
    window._display_result(old)
    window._running_selection = {1}
    window._groups = {'Old group': (1,)}
    source = tmp_path / 'new.pdb'
    source.write_text('END\n')
    primary = sample_result()
    primary.source = source
    shown = Event()
    def read(path, progress, preview):
        preview(primary)
        assert shown.wait(3)
        raise ValueError('Solvent read failed')
    timer = QTimer(window)
    timer.timeout.connect(lambda: shown.set() if window._molecule_preview and window._load_dialog is None else None)
    timer.start(5)
    errors = []
    monkeypatch.setattr(main_window, 'load_pdb', read)
    monkeypatch.setattr(window, '_show_error', errors.append)
    window.load_path(source)
    timer.stop()
    assert errors == ['Solvent read failed']
    assert window.source == old.source and window.current_result is old
    assert window._running_selection == {1}
    assert window._groups == {'Old group': (1,)}
    assert source not in window._loaded_results
    assert not window._molecule_preview and not window._loading_structure
    window._project_dirty = False
    window.close()


def test_virtual_browser_filtered_checkbox_maps_to_source_atom(monkeypatch):
    window = make_window(monkeypatch)
    window._display_result(sample_result())
    window.structure_search.setText('CA')
    proxy = window.structure_browser.model()
    assert proxy.rowCount() == 1
    assert proxy.setData(proxy.index(0, 0), Qt.Checked, Qt.CheckStateRole)
    assert window._running_selection == {1}
    window.structure_search.clear()
    assert proxy.rowCount() == 2
    assert window.structure_browser.item(0).checkState() == Qt.Unchecked
    assert window.structure_browser.item(1).checkState() == Qt.Checked
    window._running_selection = {0}
    window._update_running_selection()
    assert window.structure_browser.item(0).checkState() == Qt.Checked
    assert window.structure_browser.item(1).checkState() == Qt.Unchecked
    window.close()


def test_preview_display_failure_releases_worker_and_restores_empty_state(monkeypatch, tmp_path):
    window = make_window(monkeypatch)
    source = tmp_path / 'preview.pdb'
    source.write_text('END\n')
    primary = sample_result()
    primary.source = source
    def read(path, progress, preview):
        preview(primary)
        return primary
    original_display = window._display_result
    def display(result, **kwargs):
        if kwargs.get('preview'):
            raise ValueError('Preview could not be displayed')
        original_display(result, **kwargs)
    errors = []
    monkeypatch.setattr(window, '_display_result', display)
    monkeypatch.setattr(window, '_show_error', errors.append)
    monkeypatch.setattr(main_window, 'load_pdb', read)
    window.load_path(source)
    assert errors == ['Preview could not be displayed']
    assert window.current_result is None and window.source is None
    assert window._load_worker is None
    assert not window._loading_structure and not window._molecule_preview
    assert source not in window._loaded_results
    window.close()


def test_browser_selection_populates_information_for_every_category(monkeypatch):
    window = make_window(monkeypatch)
    window._display_result(sample_result())
    for category in ('Atoms', 'Residues', 'Chains', 'Molecules'):
        window._clear_selection('Reset for test')
        window.structure_browser_type.setCurrentText(category)
        window.structure_browser.item(0).setCheckState(Qt.Checked)
        assert window.selection_group.title() != 'No selection'
        assert window.atom_residue.text() == 'GLY 7'
        assert window.atom_chain.text() == 'A'
        assert window.atom_name.text() != '—'
        assert window.atom_element.text() != '—'
        assert window.atom_position.text() != '—'
        assert 'Å' in window.atom_properties.text()
        window.structure_browser.item(0).setCheckState(Qt.Unchecked)
        assert window.selection_group.title() == 'No selection'
        assert window.atom_residue.text() == '—'
        assert window.atom_position.text() == '—'
    window.close()


def test_information_tracks_remaining_atoms_after_additive_removal(monkeypatch):
    from dataclasses import replace
    window = make_window(monkeypatch)
    result = sample_result()
    result.atoms[0] = replace(result.atoms[0], radius=1.5, mass=14., charge=0.)
    result.atoms[1] = replace(result.atoms[1], residue_name='ILE', residue_sequence='8',
                              chain='B', radius=1.8, mass=12., charge=-0.5)
    window._display_result(result)
    window._show_selected_atom(result.atoms[0])
    window._show_selected_atom(result.atoms[1], additive=True)
    assert window.atom_name.text() == 'CA, N'
    assert window.atom_element.text() == 'C, N'
    assert window.atom_residue.text() == 'GLY 7, ILE 8'
    assert window.atom_chain.text() == 'A, B'
    assert window.atom_position.text() == 'Center: 0.500, 1.000, 2.000'
    assert window.atom_properties.text() == '1.5 to 1.8 Å · 12 to 14 Da · -0.5 to +0 e'
    window._show_selected_atom(result.atoms[1], additive=True)
    assert window.atom_name.text() == 'N (#1)'
    assert window.atom_residue.text() == 'GLY 7'
    assert window.atom_chain.text() == 'A'
    assert window.atom_properties.text() == '1.5 Å · 14 Da · +0 e'
    window._show_selected_atom(result.atoms[0], additive=True)
    assert window.selection_group.title() == 'No selection'
    assert window.atom_name.text() == '—'
    assert window.atom_properties.text() == '—'
    window.close()


def test_saved_group_selection_refreshes_all_information(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    window._display_result(result)
    window._show_selected_atom(result.atoms[0])
    window._groups = {'Carbon': (1,)}
    window._refresh_groups_panel()
    window._select_saved_group(window.groups_list.item(0))
    assert window.selection_group.title() == 'Group: Carbon'
    assert window.atom_name.text() == 'CA (#2)'
    assert window.atom_element.text() == 'C'
    assert window.atom_residue.text() == 'GLY 7'
    assert window.atom_position.text() == '1.000, 1.000, 2.000'
    window.close()
