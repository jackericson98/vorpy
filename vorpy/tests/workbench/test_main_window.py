import os
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QApplication, QToolBar, QWidget

from vorpy.workbench.domain import AnalysisResult, Atom, Bond, GeometryLayer
from vorpy.workbench.ui import main_window
from vorpy.workbench.ui.molecular_view import MolecularView


class PlotterStub:
    def reset_camera(self):
        pass

    def screenshot(self, _filename):
        pass


class ViewerStub(QWidget):
    selected_atom = Signal(object)
    selected_residue = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.plotter = PlotterStub()
        self.calls = []
        self.result = None

    def display_result(self, result):
        self.result = result

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

    def set_layer_visible(self, name, visible):
        self.calls.append((f"layer:{name}", visible))

    def set_layer_opacity(self, name, opacity):
        self.calls.append((f"opacity:{name}", opacity))

    def set_layer_color(self, name, color):
        self.calls.append((f"color:{name}", color))

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


def test_action_state_and_visibility_controls(monkeypatch):
    window = make_window(monkeypatch)

    assert [window.workflow_tabs.tabText(i) for i in range(window.workflow_tabs.count())] == [
        "Structure", "Selection", "Groups", "Interfaces"
    ]
    assert window.inspector.tabText(0) == "Visual"
    assert window.inspector.tabText(1) == "Layers"

    assert window.solve_action.isEnabled()
    assert not window.cancel_action.isEnabled()
    assert window.show_cartoon.isChecked()
    assert not window.show_spheres.isChecked()
    assert not window.show_sticks.isChecked()
    assert not window.show_waters.isChecked()
    assert not window.show_ions.isChecked()
    assert window.water_style.currentData() == "ball-and-stick"
    assert window.water_opacity.value() == 25
    assert window.molecule_opacity.value() == 50

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

    window.select_atom_action.setChecked(True)
    window.select_residue_action.setChecked(True)
    assert not window.select_atom_action.isChecked()
    assert window.select_residue_action.isChecked()
    assert window.viewer.calls[-1] == ("selection", "residue")
    assert window.selection_mode_label.text() == "Residue selection active"


def test_bottom_tray_and_small_molecule_representation_defaults(monkeypatch):
    window = make_window(monkeypatch)
    assert [window.analysis_tray_tabs.tabText(i) for i in range(window.analysis_tray_tabs.count())] == [
        "Analysis", "Results"
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
    assert window.results.item(5, 1).text() == "3"

    layer_item = window.layer_tree.topLevelItem(0)
    layer_item.setCheckState(0, Qt.Unchecked)
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
    assert window.inspector.tabText(0) == "Visual"
    assert window.inspector.tabText(1) == "Layers"
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
        main_window, "load_pdb", lambda path: {first: first_result, second: second_result}[path]
    )
    window = make_window(monkeypatch)

    window.load_path(first)
    window._running_selection = {1}
    window._groups = {"First group": (1,)}
    window.load_path(second)
    assert window.source == second.resolve()
    assert window._running_selection == set()
    assert window.results.item(4, 1).text() == "99"

    window._switch_structure_item(window.loaded_structures.item(0))
    assert window.source == first.resolve()
    assert window._running_selection == {1}
    assert window._groups == {"First group": (1,)}
    assert window.results.item(4, 1).text() == "2"


def test_molecule_browser_uses_bond_connected_components(monkeypatch):
    window = make_window(monkeypatch)
    result = sample_result()
    result.atoms.append(Atom(2, 3, "O", "O", (4.0, 0.0, 0.0), "HOH", "8", ""))
    window._display_result(result)

    window.structure_browser_type.setCurrentText("Molecules")

    assert window.structure_browser.count() == 2
    assert window.structure_browser.item(0).data(Qt.UserRole) == (0, 1)
    assert window.structure_browser.item(1).data(Qt.UserRole) == (2,)


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
