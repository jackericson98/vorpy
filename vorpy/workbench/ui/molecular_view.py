"""Scalable molecular-context and VorPy-result rendering."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pyvista as pv
from PySide6.QtCore import QEvent, Qt, QTimer, Signal
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget
from pyvistaqt import QtInteractor
from vtkmodules.vtkCommonDataModel import vtkPlane

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
    selected_atom = Signal(object, bool)
    selected_residue = Signal(object, bool)
    selected_chain = Signal(object, bool)
    selected_molecule = Signal(object, bool)
    selection_cleared = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.plotter = QtInteractor(self)
        layout.addWidget(self.plotter.interactor, 1)
        self.plotter.interactor.installEventFilter(self)
        self.depth_clip_label = QLabel("Depth cut: off")
        self.depth_clip_label.setObjectName("sectionLabel")
        self.depth_clip_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.depth_clip_label.hide()
        layout.addWidget(self.depth_clip_label)
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
        self._selection_additive = False
        self._pending_pick = None
        self._depth_clip_fraction = 0.0
        self._depth_clip_plane = vtkPlane()
        # Display radii for VorPy geometry only. Atom spheres remain 1:1 with
        # the computational radii stored on Atom.radius.
        self._edge_radius = 0.020
        self._vertex_radius = 0.130

    def clear_result(self) -> None:
        """Clear molecular state and picking when closing the active project."""
        try:
            self.plotter.disable_picking()
        except (AttributeError, RuntimeError):
            pass
        self.reset_depth_clipping(render=False)
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
        self.reset_depth_clipping(render=False)
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
        self._refresh_scalar_bar()
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
            # Render the exact radius used by the VorPy computation.
            # The unit sphere is scaled 1:1, so a 1.70 Å atom is displayed
            # with a 1.70 Å radius in the molecular coordinate system.
            cloud["radius"] = np.asarray(
                [atom.radius for atom in element_atoms],
                dtype=float,
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
                    [atom.radius for atom in element_atoms],
                    dtype=float,
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

    def _layer_scalar_limit(self, scheme: str, interpretation: str) -> float:
        values = []
        for layer in self._layer_definitions.values():
            if layer.interpretation != interpretation or scheme not in layer.cell_scalars:
                continue
            array = np.asarray(layer.cell_scalars[scheme], dtype=float).ravel()
            finite = array[np.isfinite(array)]
            if len(finite):
                values.extend(np.abs(finite))
        if not values:
            return 1.0
        limit = float(np.percentile(np.asarray(values, dtype=float), 98.0))
        return limit if np.isfinite(limit) and limit > 0.0 else 1.0

    def _display_scalars(self, layer: GeometryLayer, raw: np.ndarray) -> tuple[np.ndarray, tuple[float, float]]:
        values = np.asarray(raw, dtype=float).reshape(-1)
        finite = values[np.isfinite(values)]
        replacement = float(np.median(finite)) if len(finite) else 0.0
        values = np.nan_to_num(values, nan=replacement, posinf=replacement, neginf=replacement)
        if layer.color_scheme == "inside_outside":
            return np.clip(values, 0.0, 1.0), (0.0, 1.0)

        if layer.color_scheme == "distance":
            low, high = np.percentile(values, (2.0, 98.0)) if len(values) else (0.0, 1.0)
            if np.isclose(low, high):
                padding = max(abs(float(low)) * 0.01, 1e-9)
                low, high = low - padding, high + padding
            return values, (float(low), float(high))

        limit = self._layer_scalar_limit(layer.color_scheme, layer.interpretation)
        ratio = np.clip(values / limit, -1.0, 1.0)
        if layer.interpretation == "magnitude":
            ratio = np.abs(ratio)
        scale_gains = {
            "log10": 10.0,
            "log100": 100.0,
            "signed_log": 1000.0,  # compatibility with older projects
            "log1000": 1000.0,
            "log10000": 10000.0,
            "log100000": 100000.0,
        }
        gain = scale_gains.get(layer.scale_mode)
        if gain is not None:
            ratio = (
                np.sign(ratio)
                * np.log1p(gain * np.abs(ratio))
                / np.log1p(gain)
            )
        if layer.interpretation == "magnitude":
            ratio = 0.5 + 0.5 * np.abs(ratio)
            return ratio, (0.0, 1.0)
        return ratio, (-1.0, 1.0)

    @staticmethod
    def _physical_line_ranges(lines: np.ndarray) -> list[tuple[int, int]]:
        """Return contiguous segment ranges for disconnected polyline edges."""
        lines = np.asarray(lines, dtype=np.int64).reshape((-1, 2))
        if len(lines) == 0:
            return []
        starts = [0]
        for index in range(1, len(lines)):
            if int(lines[index, 0]) != int(lines[index - 1, 1]):
                starts.append(index)
        starts.append(len(lines))
        return list(zip(starts[:-1], starts[1:]))

    def _attach_line_scalars(
        self,
        line_mesh: pv.PolyData,
        lines: np.ndarray,
        values: np.ndarray,
        layer_name: str,
        scheme: str,
    ) -> str:
        """Attach edge scalars robustly and return PyVista preference.

        Accepted backend representations:
          * one value per rendered point,
          * one value per rendered segment,
          * one value per physical polyline edge.
        """
        values = np.asarray(values, dtype=float).reshape(-1)

        if len(values) == line_mesh.n_points:
            line_mesh.point_data["layer_values"] = values
            return "point"

        if len(values) == line_mesh.n_cells:
            line_mesh.cell_data["layer_values"] = values
            converted = line_mesh.cell_data_to_point_data(pass_cell_data=False)
            line_mesh.point_data.clear()
            line_mesh.copy_from(converted)
            return "point"

        ranges = self._physical_line_ranges(lines)
        if len(values) == len(ranges):
            expanded = np.empty(line_mesh.n_cells, dtype=float)
            for value, (start, stop) in zip(values, ranges):
                expanded[start:stop] = value
            line_mesh.cell_data["layer_values"] = expanded
            converted = line_mesh.cell_data_to_point_data(pass_cell_data=False)
            line_mesh.point_data.clear()
            line_mesh.copy_from(converted)
            return "point"

        raise ValueError(
            f"{scheme} has {len(values)} values for layer {layer_name!r}; "
            f"expected {line_mesh.n_points} points, {line_mesh.n_cells} "
            f"segments, or {len(ranges)} physical polyline edges"
        )

    @staticmethod
    def _scheme_title(layer: GeometryLayer) -> str:
        titles = {
            "mean_curvature": "Mean Curvature",
            "gaussian_curvature": "Gaussian Curvature",
            "integrated_mean_curvature": "Integrated Mean Curvature",
            "integrated_gaussian_curvature": "Integrated Gaussian Curvature",
            "surface_energy": "Surface Energy",
            "distance": "Distance to Center",
            "inside_outside": "Inside / Outside",
        }
        title = titles.get(
            layer.color_scheme,
            layer.color_scheme.replace("_", " ").title(),
        )
        if layer.color_scheme in {
            "mean_curvature",
            "gaussian_curvature",
            "integrated_mean_curvature",
            "integrated_gaussian_curvature",
        }:
            scale_labels = {
                "linear": "linear",
                "signed_log": "Log1000",
                "log10": "Log10",
                "log100": "Log100",
                "log1000": "Log1000",
                "log10000": "Log10000",
                "log100000": "Log100000",
            }
            scale = scale_labels.get(layer.scale_mode, layer.scale_mode)
            return f"{title} · {scale} normalized"
        return title

    def _refresh_scalar_bar(self) -> None:
        """Show exactly one scalar bar for the active scalar visualization.

        PyVista raises StopIteration when remove_scalar_bar() is called with no
        registered scalar bars, so cleanup must explicitly tolerate the empty
        state. This method is called during initial result loading before any
        scalar bar necessarily exists.
        """
        try:
            scalar_bars = getattr(self.plotter, "scalar_bars", None)
            actors = getattr(scalar_bars, "_scalar_bar_actors", None)
            if actors:
                self.plotter.remove_scalar_bar()
        except (AttributeError, KeyError, RuntimeError, StopIteration):
            pass

        for layer in self._layer_definitions.values():
            if not layer.visible or layer.color_scheme == "solid":
                continue
            actors = self._actors.get(layer.name, ())
            actor = next((item for item in actors if item is not None), None)
            if actor is None:
                continue
            try:
                mapper = actor.GetMapper()
                self.plotter.add_scalar_bar(
                    title=self._scheme_title(layer),
                    mapper=mapper,
                    n_labels=2 if layer.color_scheme == "inside_outside" else 5,
                    vertical=False,
                    position_x=0.27,
                    position_y=0.025,
                    width=0.46,
                    height=0.055,
                    color="white",
                    title_font_size=10,
                    label_font_size=8,
                    fmt="%.2g",
                )
                return
            except (AttributeError, RuntimeError, TypeError):
                continue

    def _add_layer(self, layer: GeometryLayer) -> None:
        try:
            preference = "cell"
            scalar_values = layer.cell_scalars.get(layer.color_scheme)
            if layer.source_path is not None:
                mesh = pv.read(layer.source_path)
            elif layer.faces is not None:
                points = np.asarray(layer.points, dtype=float).reshape((-1, 3))
                faces = np.asarray(layer.faces, dtype=np.int64).reshape((-1, 3))
                face_cells = np.column_stack((np.full(len(faces), 3, dtype=np.int64), faces)).ravel()
                mesh = pv.PolyData(points, faces=face_cells)
                if scalar_values is not None:
                    values, _ = self._display_scalars(layer, scalar_values)

                    if len(values) != mesh.n_cells:
                        raise ValueError(f"{layer.color_scheme} has {len(values)} values for {mesh.n_cells} surface triangles")
                    mesh.cell_data["layer_values"] = values
            elif layer.lines is not None:
                cells = np.asarray(
                    [[2, a, b] for a, b in layer.lines],
                    dtype=np.int64,
                ).ravel()
                # Build a line-only PolyData. Calling pv.PolyData(points)
                # creates one vertex cell per point; adding lines afterward
                # leaves those vertex cells in the dataset, so n_cells becomes
                # n_points + n_lines. That made correctly sized edge scalar
                # arrays appear too short. Construct the dataset explicitly
                # without vertex cells instead.
                line_mesh = pv.PolyData()
                line_mesh.points = np.asarray(layer.points, dtype=float)
                line_mesh.lines = cells

                if scalar_values is not None:
                    values, _ = self._display_scalars(layer, scalar_values)
                    if len(values) not in {
                        line_mesh.n_points,
                        line_mesh.n_lines,
                        len(self._physical_line_ranges(layer.lines)),
                    }:
                        print(
                            f"[EDGE SCALAR DEBUG] {layer.name}: "
                            f"points={line_mesh.n_points}, "
                            f"lines={line_mesh.n_lines}, "
                            f"cells={line_mesh.n_cells}, "
                            f"physical_edges={len(self._physical_line_ranges(layer.lines))}, "
                            f"values={len(values)}, "
                            f"scheme={layer.color_scheme}"
                        )
                    preference = self._attach_line_scalars(
                        line_mesh,
                        layer.lines,
                        values,
                        layer.name,
                        layer.color_scheme,
                    )

                mesh = line_mesh.tube(
                    radius=self._edge_radius,
                    n_sides=10,
                    capping=True,
                )

                if scalar_values is not None and "layer_values" not in mesh.array_names:
                    raise ValueError(
                        f"Scalar data were lost while tubing layer {layer.name!r}"
                    )
            else:
                cloud = pv.PolyData(np.asarray(layer.points, dtype=float))
                cloud["radius"] = np.full(
                    len(layer.points),
                    self._vertex_radius,
                    dtype=float,
                )
                if scalar_values is not None:
                    values, _ = self._display_scalars(layer, scalar_values)
                    if len(values) != cloud.n_points:
                        raise ValueError(f"{layer.color_scheme} has {len(values)} values for {cloud.n_points} vertices")
                    cloud.point_data["layer_values"] = values
                    preference = "point"
                mesh = cloud.glyph(scale="radius", orient=False, geom=pv.Icosahedron(radius=1.0))

            mesh_options = {
                "opacity": layer.opacity,
                "name": f"layer-{layer.name}",
                "show_edges": False,
            }
            if scalar_values is not None and "layer_values" in mesh.array_names:
                _, clim = self._display_scalars(layer, scalar_values)
                mesh_options.update(
                    scalars="layer_values",
                    preference=preference,
                    cmap=layer.color_map,
                    clim=clim,
                    interpolate_before_map=False,
                    lighting=layer.kind.lower() != "surfaces",
                    smooth_shading=layer.kind.lower() in {"edges", "vertices"},
                    show_scalar_bar=False,
                )
            else:
                mesh_options["color"] = layer.color
            actor = self.plotter.add_mesh(
                mesh,
                reset_camera=False,
                **mesh_options,
            )
            actor.SetVisibility(layer.visible)
            self._actors[layer.name].append(actor)
            self._apply_depth_clip_to_actor(actor)
        except Exception as error:  # noqa: BLE001
            print(f"[VorPy Workbench] Failed to render layer {layer.name!r}: {error}")

    def _picked_point(self, point) -> None:
        # PyVista reports a point at the start of a left-button gesture. Defer
        # changing selection until release so a camera drag cannot replace it.
        if point is not None:
            self._pending_pick = np.asarray(point, dtype=float)

    def _apply_pending_pick(self) -> None:
        point = self._pending_pick
        self._pending_pick = None
        additive = self._selection_additive
        self._selection_additive = False
        if (
            self._result is None
            or self._selection_mode is None
            or len(self._positions) == 0
            or self._selection_dragged
        ):
            return
        if point is None:
            self._clear_pick_highlight()
            self.selection_cleared.emit()
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
            selected_atoms = [
                item for item in self._result.atoms if self._residue_key(item) == key
            ]
            signal = self.selected_residue
        elif self._selection_mode == "chain":
            selected_atoms = self._chain_atoms(atom)
            signal = self.selected_chain
        elif self._selection_mode == "molecule":
            selected_atoms = self._molecule_atoms(atom)
            signal = self.selected_molecule
        else:
            self._highlight_atoms([atom], "selected-atom")
            self.selected_atom.emit(atom, additive)
            return
        self._highlight_atoms(selected_atoms, f"selected-{self._selection_mode}")
        signal.emit(selected_atoms, additive)

    def eventFilter(self, watched, event) -> bool:
        if (
            watched is self.plotter.interactor
            and event.type() == QEvent.Type.Wheel
            and self._result is not None
        ):
            steps = event.angleDelta().y() / 120.0
            if steps:
                self.adjust_depth_clipping(steps)
                event.accept()
                return True
        if watched is self.plotter.interactor and self._selection_mode is not None:
            if (
                event.type() == QEvent.Type.MouseButtonPress
                and event.button() == Qt.MouseButton.LeftButton
            ):
                position = event.position()
                self._press_position = (position.x(), position.y())
                self._selection_dragged = False
                self._selection_additive = bool(
                    event.modifiers() & Qt.KeyboardModifier.ShiftModifier
                )
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

    def _rendered_actors(self) -> list[object]:
        return [
            actor
            for actors in self._actors.values()
            for actor in actors
            if actor is not None
        ]

    def _apply_depth_clip_to_actor(self, actor) -> None:
        mapper = actor.GetMapper()
        mapper.RemoveAllClippingPlanes()
        if self._depth_clip_fraction > 0.0:
            mapper.AddClippingPlane(self._depth_clip_plane)

    def _depth_projection_range(
        self, direction: np.ndarray
    ) -> tuple[float, float] | None:
        projections = []
        for actor in self._rendered_actors():
            bounds = actor.GetBounds()
            if bounds is None or len(bounds) != 6:
                continue
            for x in (bounds[0], bounds[1]):
                for y in (bounds[2], bounds[3]):
                    for z in (bounds[4], bounds[5]):
                        projections.append(float(np.dot((x, y, z), direction)))
        if not projections:
            return None
        return min(projections), max(projections)

    def adjust_depth_clipping(self, wheel_steps: float) -> None:
        self._depth_clip_fraction = float(
            np.clip(self._depth_clip_fraction + 0.04 * wheel_steps, 0.0, 0.98)
        )
        self._apply_depth_clipping()

    def _apply_depth_clipping(self, render: bool = True) -> None:
        if self._depth_clip_fraction > 0.0:
            direction = np.asarray(
                self.plotter.camera.GetDirectionOfProjection(), dtype=float
            )
            magnitude = np.linalg.norm(direction)
            if magnitude:
                direction /= magnitude
                projection_range = self._depth_projection_range(direction)
                if projection_range is not None:
                    near, far = projection_range
                    cut = near + self._depth_clip_fraction * (far - near)
                    self._depth_clip_plane.SetNormal(direction)
                    self._depth_clip_plane.SetOrigin(direction * cut)

        for actor in self._rendered_actors():
            self._apply_depth_clip_to_actor(actor)

        if self._depth_clip_fraction > 0.0:
            self.depth_clip_label.setText(
                f"Depth cut: {self._depth_clip_fraction:.0%} · scroll out to restore"
            )
            self.depth_clip_label.show()
        else:
            self.depth_clip_label.hide()
        if render:
            self.plotter.render()

    def set_depth_clipping_fraction(self, fraction: float) -> None:
        self._depth_clip_fraction = float(np.clip(fraction, 0.0, 0.98))
        self._apply_depth_clipping()

    def reset_depth_clipping(self, render: bool = True) -> None:
        self._depth_clip_fraction = 0.0
        for actor in self._rendered_actors():
            actor.GetMapper().RemoveAllClippingPlanes()
        self.depth_clip_label.hide()
        if render:
            self.plotter.render()

    def _highlight_atoms(self, atoms: list[Atom], name: str) -> None:
        self._clear_pick_highlight()
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

    def _clear_pick_highlight(self) -> None:
        for selection_name in ("atom", "residue", "chain", "molecule"):
            self.plotter.remove_actor(f"selected-{selection_name}", render=False)

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

    def _molecule_atoms(self, atom: Atom) -> list[Atom]:
        """Return the bond-connected component containing ``atom``."""
        if self._result is None:
            return []
        neighbors = {item.index: set() for item in self._result.atoms}
        for bond in self._result.bonds:
            neighbors[bond.atom_a].add(bond.atom_b)
            neighbors[bond.atom_b].add(bond.atom_a)
        pending = [atom.index]
        selected: set[int] = set()
        while pending:
            atom_index = pending.pop()
            if atom_index in selected:
                continue
            selected.add(atom_index)
            pending.extend(neighbors[atom_index] - selected)
        return [item for item in self._result.atoms if item.index in selected]

    def _chain_atoms(self, atom: Atom) -> list[Atom]:
        """Select a chain while keeping blank-chain waters independent."""
        if self._result is None:
            return []
        if atom.chain:
            return [item for item in self._result.atoms if item.chain == atom.chain]
        if self._is_water(atom):
            residue = self._residue_key(atom)
            return [
                item
                for item in self._result.atoms
                if self._residue_key(item) == residue
            ]
        return [
            item
            for item in self._result.atoms
            if not item.chain and not self._is_water(item)
        ]

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
        layer = self._layer_definitions.get(name)
        if layer is not None:
            layer.visible = visible
        if visible and not self._actors.get(name):
            if layer is not None:
                self._add_layer(layer)
        for actor in self._actors.get(name, ()):
            actor.SetVisibility(visible)
        self._refresh_scalar_bar()
        self.plotter.render()

    def set_layer_opacity(self, name: str, opacity: float) -> None:
        layer = self._layer_definitions.get(name)
        if layer is not None:
            layer.opacity = opacity
        for actor in self._actors.get(name, ()):
            actor.GetProperty().SetOpacity(opacity)
        self.plotter.render()

    def _redraw_layer(self, layer: GeometryLayer) -> None:
        """Redraw one geometry layer without changing the current camera view."""
        camera_position = self.plotter.camera_position
        parallel_scale = self.plotter.camera.parallel_scale
        view_angle = self.plotter.camera.view_angle

        self.plotter.remove_actor(f"layer-{layer.name}", render=False)
        self._actors[layer.name].clear()

        if layer.visible:
            self._add_layer(layer)

        # Re-adding an actor can cause PyVista/VTK to recompute the camera.
        # Restore the exact user view so size/color changes feel live rather
        # than like a new scene load.
        self.plotter.camera_position = camera_position
        self.plotter.camera.parallel_scale = parallel_scale
        self.plotter.camera.view_angle = view_angle

        self._refresh_scalar_bar()
        self._apply_depth_clipping(render=False)
        self.plotter.render()

    def set_layer_color_scheme(self, name: str, scheme: str) -> None:
        layer = self._layer_definitions.get(name)
        if layer is None:
            return
        layer.color_scheme = scheme if scheme == "solid" or scheme in layer.cell_scalars else "solid"
        self._redraw_layer(layer)

    def set_layer_colormap(self, name: str, color_map: str) -> None:
        layer = self._layer_definitions.get(name)
        if layer is None:
            return
        layer.color_map = str(color_map)
        if layer.color_scheme != "solid":
            self._redraw_layer(layer)

    def set_layer_scale_mode(self, name: str, scale_mode: str) -> None:
        layer = self._layer_definitions.get(name)
        if layer is None:
            return
        allowed = {
            "linear", "log10", "log100", "log1000", "log10000", "log100000"
        }
        # Older saved projects used the generic name ``signed_log``.
        if scale_mode == "signed_log":
            scale_mode = "log1000"
        layer.scale_mode = scale_mode if scale_mode in allowed else "log1000"
        if layer.color_scheme != "solid":
            self._redraw_layer(layer)

    def set_edge_radius(self, radius: float) -> None:
        """Set the displayed VorPy edge-tube radius and redraw edge layers."""
        radius = float(radius)
        if not np.isfinite(radius) or radius <= 0.0:
            return
        if np.isclose(radius, self._edge_radius):
            return
        self._edge_radius = radius
        for layer in list(self._layer_definitions.values()):
            if layer.kind.lower() == "edges" and layer.visible:
                self._redraw_layer(layer)

    def set_vertex_radius(self, radius: float) -> None:
        """Set the displayed VorPy vertex-glyph radius and redraw vertex layers."""
        radius = float(radius)
        if not np.isfinite(radius) or radius <= 0.0:
            return
        if np.isclose(radius, self._vertex_radius):
            return
        self._vertex_radius = radius
        for layer in list(self._layer_definitions.values()):
            if layer.kind.lower() == "vertices" and layer.visible:
                self._redraw_layer(layer)

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
