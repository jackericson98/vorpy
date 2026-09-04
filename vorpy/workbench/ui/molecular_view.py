"""Scalable molecular-context and VorPy-result rendering."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pyvista as pv
from PySide6.QtCore import QEvent, Qt, QTimer, Signal
from PySide6.QtWidgets import QVBoxLayout, QWidget
from pyvistaqt import QtInteractor

from vorpy.workbench.domain import AnalysisResult, Atom, GeometryLayer

ELEMENT_COLORS = {
    "C": "#aeb7c2",
    "N": "#3979d4",
    "O": "#d94a4a",
    "H": "#f3f3f3",
    "S": "#d4ad32",
    "P": "#e38a36",
    "F": "#62c96f",
    "CL": "#62c96f",
}
VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "MG": 1.73,
    "CA": 2.31,
    "ZN": 1.39,
    "FE": 1.56,
}
WATER_RESIDUES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE"}
ION_RESIDUES = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "ZN",
    "CD",
    "FE",
    "FE2",
    "FE3",
    "MN",
    "CU",
    "CU1",
    "CU2",
    "CO",
    "NI",
    "AL",
    "F",
    "CL",
    "BR",
    "I",
    "IOD",
    "SO4",
    "PO4",
    "NH4",
}
CARTOON_COLORS = ("#6f7ee8", "#39a88e", "#d27a43", "#9c68cf", "#cf5f7b")


class MolecularView(QWidget):
    selected_atom = Signal(object)
    selected_residue = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor)
        self.plotter.interactor.installEventFilter(self)
        self.plotter.set_background("#111820")
        self._result: AnalysisResult | None = None
        self._positions = np.empty((0, 3))
        self._actors: dict[str, list[object]] = defaultdict(list)
        self._axes_actor = None
        self._layer_definitions: dict[str, GeometryLayer] = {}
        self._selection_mode: str | None = None
        self._cartoon_visible = True
        self._spheres_visible = False
        self._molecule_opacity = 0.5
        self._waters_visible = False
        self._water_style = "ball-and-stick"
        self._water_opacity = 0.25
        self._ions_visible = False
        self._bonds_visible = False
        self._press_position: tuple[float, float] | None = None
        self._selection_dragged = False
        self._pending_pick = None

    def clear_result(self) -> None:
        """Clear molecular state and picking when closing the active project."""
        try:
            self.plotter.disable_picking()
        except (AttributeError, RuntimeError):
            pass
        self.plotter.clear()
        self._actors.clear()
        self._layer_definitions.clear()
        self._result = None
        self._positions = np.empty((0, 3))
        self.plotter.render()

    def display_result(self, result: AnalysisResult) -> None:
        # Picking survives plotter.clear(), so explicitly remove the previous
        # observer before installing a picker for the newly loaded result.
        try:
            self.plotter.disable_picking()
        except (AttributeError, RuntimeError):
            pass
        self.plotter.clear()
        self._actors.clear()
        self._layer_definitions = {layer.name: layer for layer in result.layers}
        self._result = result
        self._positions = np.asarray(
            [atom.position for atom in result.atoms], dtype=float
        )
        self._add_cartoon(result.atoms)
        self._add_atoms(result.atoms)
        self._add_bonds(result)
        for layer in result.layers:
            if layer.visible:
                self._add_layer(layer)
        self._axes_actor = self.plotter.add_axes()
        self.plotter.reset_camera()
        self.plotter.render()
        if self._selection_mode is not None:
            self.set_selection_mode(self._selection_mode)

    def _add_atoms(self, atoms: list[Atom]) -> None:
        grouped: dict[tuple[str, str], list[Atom]] = defaultdict(list)
        for atom in atoms:
            if self._is_water(atom):
                category = "waters"
            elif self._is_ion(atom):
                category = "ions"
            else:
                category = "atoms"
            grouped[(atom.element.upper(), category)].append(atom)
        unit_sphere = pv.Sphere(radius=1.0, theta_resolution=14, phi_resolution=14)
        for (element, category), element_atoms in grouped.items():
            cloud = pv.PolyData(np.asarray([atom.position for atom in element_atoms]))
            if category == "waters":
                cloud["radius"] = np.asarray(
                    [atom.radius * 1.25 for atom in element_atoms]
                )
            else:
                cloud["radius"] = np.asarray(
                    [VDW_RADII.get(atom.element.upper(), 1.7) for atom in element_atoms]
                )
            glyphs = cloud.glyph(scale="radius", orient=False, geom=unit_sphere)
            opacity = (
                self._water_opacity if category == "waters" else self._molecule_opacity
            )
            actor = self.plotter.add_mesh(
                glyphs,
                color=ELEMENT_COLORS.get(element, "#bf8bd8"),
                opacity=opacity,
                smooth_shading=True,
                name=f"{category}-{element}",
            )
            visible = (
                self._waters_visible and self._water_style == "ball-and-stick"
                if category == "waters"
                else (
                    self._ions_visible if category == "ions" else self._spheres_visible
                )
            )
            actor.SetVisibility(visible)
            self._actors[category].append(actor)
            if category == "waters":
                vdw_cloud = pv.PolyData(
                    np.asarray([atom.position for atom in element_atoms])
                )
                vdw_cloud["radius"] = np.asarray(
                    [VDW_RADII.get(atom.element.upper(), 1.7) for atom in element_atoms]
                )
                vdw_actor = self.plotter.add_mesh(
                    vdw_cloud.glyph(scale="radius", orient=False, geom=unit_sphere),
                    color=ELEMENT_COLORS.get(element, "#bf8bd8"),
                    opacity=self._water_opacity,
                    smooth_shading=True,
                    name=f"water-spheres-{element}",
                )
                vdw_actor.SetVisibility(
                    self._waters_visible and self._water_style == "spheres"
                )
                self._actors["water_spheres"].append(vdw_actor)

    def _add_cartoon(self, atoms: list[Atom]) -> None:
        """Render protein alpha-carbon and nucleic phosphate/sugar traces."""
        protein_chains: dict[str, list[Atom]] = defaultdict(list)
        nucleotide_atoms: dict[tuple[str, str, str], dict[str, Atom]] = defaultdict(dict)
        for atom in atoms:
            if self._is_water(atom) or self._is_ion(atom):
                continue
            atom_name = atom.name.strip().upper().replace("*", "'")
            if atom_name == "CA":
                protein_chains[atom.chain].append(atom)
            elif atom_name in {"P", "C4'"}:
                residue = (atom.chain, atom.residue_sequence, atom.residue_name)
                nucleotide_atoms[residue][atom_name] = atom

        nucleic_chains: dict[str, list[Atom]] = defaultdict(list)
        for residue_atoms in nucleotide_atoms.values():
            # P atoms make the conventional nucleic backbone trace. C4' keeps
            # terminal nucleotides (which commonly omit P) in the cartoon.
            anchor = residue_atoms.get("P") or residue_atoms.get("C4'")
            if anchor is not None:
                nucleic_chains[anchor.chain].append(anchor)

        color_index = 0
        for chain_atoms in protein_chains.values():
            self._add_cartoon_chain(chain_atoms, color_index, maximum_gap=5.0)
            color_index += 1
        for chain_atoms in nucleic_chains.values():
            self._add_cartoon_chain(chain_atoms, color_index, maximum_gap=10.0)
            color_index += 1

    def _add_cartoon_chain(
        self,
        chain_atoms: list[Atom],
        color_index: int,
        maximum_gap: float,
    ) -> None:
        """Split one ordered backbone into continuous drawable segments."""
        chain_atoms.sort(key=self._residue_sort_key)
        segment: list[Atom] = []
        for atom in chain_atoms:
            if segment:
                distance = np.linalg.norm(
                    np.asarray(atom.position) - np.asarray(segment[-1].position)
                )
                if distance > maximum_gap:
                    self._add_cartoon_segment(segment, color_index)
                    segment = []
            segment.append(atom)
        self._add_cartoon_segment(segment, color_index)

    def _add_cartoon_segment(self, atoms: list[Atom], color_index: int) -> None:
        if len(atoms) < 2:
            return
        points = np.asarray([atom.position for atom in atoms], dtype=float)
        if len(points) == 2:
            centerline = pv.Line(points[0], points[1], resolution=8)
        else:
            centerline = pv.Spline(points, n_points=max(12, len(points) * 8))
        actor = self.plotter.add_mesh(
            centerline.tube(radius=0.3, n_sides=16),
            color=CARTOON_COLORS[color_index % len(CARTOON_COLORS)],
            smooth_shading=True,
            name=f"cartoon-{color_index}-{len(self._actors['cartoon'])}",
        )
        actor.SetVisibility(self._cartoon_visible)
        self._actors["cartoon"].append(actor)

    def _add_bonds(self, result: AnalysisResult) -> None:
        if not result.bonds:
            return
        points = np.asarray([atom.position for atom in result.atoms], dtype=float)
        regular_bonds = []
        water_bonds = []
        ion_bonds = []
        for bond in result.bonds:
            if self._is_water(result.atoms[bond.atom_a]) or self._is_water(
                result.atoms[bond.atom_b]
            ):
                target = water_bonds
            elif self._is_ion(result.atoms[bond.atom_a]) or self._is_ion(
                result.atoms[bond.atom_b]
            ):
                target = ion_bonds
            else:
                target = regular_bonds
            target.append(bond)
        self._add_bond_group(
            points, result.atoms, regular_bonds, "bonds", self._bonds_visible
        )
        self._add_bond_group(
            points,
            result.atoms,
            water_bonds,
            "water_bonds",
            self._waters_visible and self._water_style in {"sticks", "ball-and-stick"},
        )
        self._add_bond_group(
            points,
            result.atoms,
            ion_bonds,
            "ion_bonds",
            self._bonds_visible and self._ions_visible,
        )

    def _add_bond_group(
        self,
        points: np.ndarray,
        atoms: list[Atom],
        bonds: list,
        category: str,
        visible: bool,
    ) -> None:
        if not bonds:
            return
        half_bonds: dict[str, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(list)
        for bond in bonds:
            point_a = points[bond.atom_a]
            point_b = points[bond.atom_b]
            midpoint = (point_a + point_b) / 2.0
            element_a = atoms[bond.atom_a].element.upper()
            element_b = atoms[bond.atom_b].element.upper()
            half_bonds[element_a].append((point_a, midpoint))
            half_bonds[element_b].append((midpoint, point_b))

        for element, segments in half_bonds.items():
            line_points = np.asarray(
                [point for segment in segments for point in segment], dtype=float
            )
            cells = np.asarray(
                [[2, index, index + 1] for index in range(0, len(line_points), 2)],
                dtype=np.int64,
            ).ravel()
            lines = pv.PolyData(line_points)
            lines.lines = cells
            actor = self.plotter.add_mesh(
                lines.tube(radius=0.1, n_sides=10),
                color=ELEMENT_COLORS.get(element, "#bf8bd8"),
                smooth_shading=True,
                name=f"{category}-{element}",
            )
            actor.SetVisibility(visible)
            self._actors[category].append(actor)

    def _add_layer(self, layer: GeometryLayer) -> None:
        try:
            if layer.source_path is not None:
                mesh = pv.read(layer.source_path)
            elif layer.lines is not None:
                cells = np.asarray(
                    [[2, a, b] for a, b in layer.lines], dtype=np.int64
                ).ravel()
                mesh = pv.PolyData(np.asarray(layer.points, dtype=float))
                mesh.lines = cells
                mesh = mesh.tube(radius=0.045, n_sides=6)
            else:
                cloud = pv.PolyData(np.asarray(layer.points, dtype=float))
                cloud["radius"] = np.full(len(layer.points), 0.13)
                mesh = cloud.glyph(
                    scale="radius", orient=False, geom=pv.Icosahedron(radius=1.0)
                )
            actor = self.plotter.add_mesh(
                mesh,
                color=layer.color,
                opacity=layer.opacity,
                name=f"layer-{layer.name}",
                show_edges=False,
            )
            actor.SetVisibility(layer.visible)
            self._actors[layer.name].append(actor)
        except Exception as error:  # noqa: BLE001 - mesh readers expose varied exceptions.
            layer.visible = False
            layer.name = f"{layer.name} [load failed: {error}]"

    def _picked_point(self, point) -> None:
        # PyVista reports a point at the start of a left-button gesture. Defer
        # changing selection until release so a camera drag cannot replace it.
        if point is not None:
            self._pending_pick = np.asarray(point, dtype=float)

    def _apply_pending_pick(self) -> None:
        point = self._pending_pick
        self._pending_pick = None
        if (
            self._result is None
            or self._selection_mode is None
            or len(self._positions) == 0
            or point is None
            or self._selection_dragged
        ):
            return
        selectable = [
            index
            for index, atom in enumerate(self._result.atoms)
            if (self._waters_visible or not self._is_water(atom))
            and (self._ions_visible or not self._is_ion(atom))
        ]
        if not selectable:
            return
        distances = np.linalg.norm(
            self._positions[selectable] - np.asarray(point), axis=1
        )
        atom = self._result.atoms[selectable[int(np.argmin(distances))]]
        if self._selection_mode == "residue":
            key = self._residue_key(atom)
            residue_atoms = [
                item for item in self._result.atoms if self._residue_key(item) == key
            ]
            self._highlight_atoms(residue_atoms, "selected-residue")
            self.selected_residue.emit(residue_atoms)
        else:
            self._highlight_atoms([atom], "selected-atom")
            self.selected_atom.emit(atom)

    def eventFilter(self, watched, event) -> bool:
        if watched is self.plotter.interactor and self._selection_mode is not None:
            if (
                event.type() == QEvent.Type.MouseButtonPress
                and event.button() == Qt.MouseButton.LeftButton
            ):
                position = event.position()
                self._press_position = (position.x(), position.y())
                self._selection_dragged = False
                self._pending_pick = None
            elif (
                event.type() == QEvent.Type.MouseMove
                and self._press_position is not None
            ):
                position = event.position()
                dx = position.x() - self._press_position[0]
                dy = position.y() - self._press_position[1]
                if dx * dx + dy * dy > 25:
                    self._selection_dragged = True
            elif (
                event.type() == QEvent.Type.MouseButtonRelease
                and event.button() == Qt.MouseButton.LeftButton
                and self._press_position is not None
            ):
                self._press_position = None
                QTimer.singleShot(0, self._apply_pending_pick)
        return super().eventFilter(watched, event)

    def _highlight_atoms(self, atoms: list[Atom], name: str) -> None:
        self.plotter.remove_actor("selected-atom", render=False)
        self.plotter.remove_actor("selected-residue", render=False)
        cloud = pv.PolyData(np.asarray([atom.position for atom in atoms], dtype=float))
        cloud["radius"] = np.asarray([atom.radius * 1.22 for atom in atoms])
        highlight = cloud.glyph(
            scale="radius",
            orient=False,
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
        )
        self.plotter.add_mesh(
            highlight,
            name=name,
            style="wireframe",
            color="#ffd45a",
            line_width=3,
            pickable=False,
            reset_camera=False,
        )

    def set_group_selection(self, atoms: list[Atom]) -> None:
        """Highlight atoms accumulated by the Structure group builder."""
        self.plotter.remove_actor("group-selection", render=False)
        if not atoms:
            self.plotter.render()
            return
        cloud = pv.PolyData(np.asarray([atom.position for atom in atoms], dtype=float))
        cloud["radius"] = np.asarray([atom.radius * 1.35 for atom in atoms])
        highlight = cloud.glyph(
            scale="radius",
            orient=False,
            geom=pv.Sphere(radius=1.0, theta_resolution=16, phi_resolution=16),
        )
        self.plotter.add_mesh(
            highlight,
            name="group-selection",
            style="wireframe",
            color="#42d6c7",
            line_width=3,
            pickable=False,
            reset_camera=False,
        )
        self.plotter.render()

    def set_selection_mode(self, mode: str | None) -> None:
        try:
            self.plotter.disable_picking()
        except (AttributeError, RuntimeError):
            pass
        self._selection_mode = mode
        if mode is not None and self._result is not None:
            self.plotter.enable_point_picking(
                callback=self._picked_point,
                show_message=False,
                show_point=False,
                pickable_window=False,
                left_clicking=True,
            )

    @staticmethod
    def _is_water(atom: Atom) -> bool:
        return atom.residue_name.strip().upper() in WATER_RESIDUES

    @staticmethod
    def _is_ion(atom: Atom) -> bool:
        name = atom.residue_name.strip().upper().rstrip("+-")
        return name in ION_RESIDUES

    @staticmethod
    def _residue_key(atom: Atom) -> tuple[str, str, str]:
        return atom.chain, atom.residue_sequence, atom.residue_name

    @staticmethod
    def _residue_sort_key(atom: Atom) -> tuple[int, str]:
        sequence = atom.residue_sequence.strip()
        digits = "".join(
            character
            for character in sequence
            if character.isdigit() or character == "-"
        )
        try:
            number = int(digits)
        except ValueError:
            number = 0
        return number, sequence

    def set_category_visible(self, category: str, visible: bool) -> None:
        if category == "atoms":
            self._spheres_visible = visible
        for actor in self._actors.get(category, ()):
            actor.SetVisibility(visible)
        self.plotter.render()

    def set_molecule_opacity(self, opacity: float) -> None:
        self._molecule_opacity = opacity
        for actor in self._actors.get("atoms", ()):
            actor.GetProperty().SetOpacity(opacity)
        self.plotter.render()

    def set_cartoon_visible(self, visible: bool) -> None:
        self._cartoon_visible = visible
        for actor in self._actors.get("cartoon", ()):
            actor.SetVisibility(visible)
        self.plotter.render()

    def set_bonds_visible(self, visible: bool) -> None:
        self._bonds_visible = visible
        for actor in self._actors.get("bonds", ()):
            actor.SetVisibility(visible)
        for actor in self._actors.get("water_bonds", ()):
            actor.SetVisibility(
                self._waters_visible
                and self._water_style in {"sticks", "ball-and-stick"}
            )
        for actor in self._actors.get("ion_bonds", ()):
            actor.SetVisibility(visible and self._ions_visible)
        self.plotter.render()

    def set_waters_visible(self, visible: bool) -> None:
        self._waters_visible = visible
        for actor in self._actors.get("waters", ()):
            actor.SetVisibility(visible and self._water_style == "ball-and-stick")
        for actor in self._actors.get("water_spheres", ()):
            actor.SetVisibility(visible and self._water_style == "spheres")
        for actor in self._actors.get("water_bonds", ()):
            actor.SetVisibility(
                visible and self._water_style in {"sticks", "ball-and-stick"}
            )
        self.plotter.render()

    def set_water_style(self, style: str) -> None:
        self._water_style = style
        self.set_waters_visible(self._waters_visible)

    def set_water_opacity(self, opacity: float) -> None:
        self._water_opacity = opacity
        for category in ("waters", "water_spheres"):
            for actor in self._actors.get(category, ()):
                actor.GetProperty().SetOpacity(opacity)
        self.plotter.render()

    def set_ions_visible(self, visible: bool) -> None:
        self._ions_visible = visible
        for actor in self._actors.get("ions", ()):
            actor.SetVisibility(visible)
        for actor in self._actors.get("ion_bonds", ()):
            actor.SetVisibility(visible and self._bonds_visible)
        self.plotter.render()

    def set_layer_visible(self, name: str, visible: bool) -> None:
        if visible and not self._actors.get(name):
            layer = self._layer_definitions.get(name)
            if layer is not None:
                self._add_layer(layer)
        for actor in self._actors.get(name, ()):
            actor.SetVisibility(visible)
        self.plotter.render()

    def set_layer_opacity(self, name: str, opacity: float) -> None:
        for actor in self._actors.get(name, ()):
            actor.GetProperty().SetOpacity(opacity)
        self.plotter.render()

    def set_layer_color(self, name: str, color: str) -> None:
        for actor in self._actors.get(name, ()):
            actor.GetProperty().SetColor(pv.Color(color).float_rgb)
        layer = self._layer_definitions.get(name)
        if layer is not None:
            layer.color = color
        self.plotter.render()

    def save_screenshot(self, filename: str, scale: int = 1) -> None:
        axes_visible = (
            bool(self._axes_actor.GetVisibility())
            if self._axes_actor is not None
            else False
        )
        if self._axes_actor is not None:
            self._axes_actor.SetVisibility(False)
        try:
            self.plotter.screenshot(filename, scale=scale)
        finally:
            if self._axes_actor is not None:
                self._axes_actor.SetVisibility(axes_visible)
            self.plotter.render()

    def closeEvent(self, event) -> None:
        self.plotter.close()
        super().closeEvent(event)
