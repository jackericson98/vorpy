"""Production adapter from the Analysis Studio to VorPy's solver objects."""

from __future__ import annotations

import time
import tempfile
import shutil
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from vorpy.src.group import Group
from vorpy.src.interface import Interface
from vorpy.src.network.views import NetworkViewError
from vorpy.src.results import (
    GeometryLineage,
    InterfaceResult,
    ResultsDependency,
    ResultIdentity,
    Provenance,
    SharedNetworkIdentity,
    adapt_interface_analysis,
    make_geometry_lineage,
)
from vorpy.src.system import System
from vorpy.src.boundary import (
    BoundaryConfig,
    BoundaryMode,
    SolventMode,
    DEFAULT_PROBE_RADIUS,
    DEFAULT_BOUNDARY_GENERATOR_RADIUS,
    DEFAULT_SHELL_OFFSET,
    DEFAULT_SHELL_SPACING,
    normalize_retained_system_indices,
    normalize_solvent_mode,
    prepare_network_geometry,
)
from vorpy.workbench.domain import AnalysisResult, Atom, GeometryLayer
from vorpy.workbench.atomic_defaults import apply_system_defaults
from vorpy.workbench.services.backend import CancellationCheck, ProgressCallback
from vorpy.workbench.services.structure_loader import load_pdb
from vorpy.workbench.services.trajectory import read_frame_lines, index_pdb_frames
from vorpy.src.output.curvature_colors import (
    component_value,
    mean_vertex_display_value,
)


@dataclass(frozen=True)
class VorPySolveSettings:
    network_type: str = "aw"
    max_vertices: int = 5
    box_size: float = 1.25
    surface_resolution: float = 0.2
    build_surfaces: bool = True
    build_vertices: bool = True
    build_edges: bool = True
    atomic_defaults: dict[str, dict[str, float]] = field(default_factory=dict)
    boundary_mode: str | None = None
    probe_radius: float = DEFAULT_PROBE_RADIUS
    shell_spacing: float = DEFAULT_SHELL_SPACING
    shell_offset: float = DEFAULT_SHELL_OFFSET
    # Optional NEXUS handoff; values are System-global atom IDs. A mapping
    # preserves Workbench group labels; the tuple form stays compatible with
    # positional callers.
    typed_group_selections: (
        Mapping[str, tuple[int, ...]] | tuple[tuple[int, ...], ...]
    ) = ()
    # Producer identity is separate from the archive object's persisted ID.
    source_network_id: str | None = None
    retained_system_indices: tuple[int, ...] = ()
    # Scientific generator-universe contract shared with the CLI. Boundary
    # mode remains a separate numerical boundary setting.
    solvent_mode: str = SolventMode.MOLECULE_ONLY.value

    def __post_init__(self):
        object.__setattr__(
            self,
            "solvent_mode",
            normalize_solvent_mode(self.solvent_mode).value,
        )
        object.__setattr__(
            self,
            "retained_system_indices",
            normalize_retained_system_indices(self.retained_system_indices),
        )


def _interface_indices(value, label):
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{label} must be a non-empty sequence of System indices")
    try:
        indices = tuple(value)
    except TypeError as error:
        raise ValueError(
            f"{label} must be a non-empty sequence of System indices"
        ) from error
    if not indices:
        raise ValueError(f"{label} must not be empty")
    if any(
        isinstance(index, bool) or not isinstance(index, int) or index < 0
        for index in indices
    ):
        raise ValueError(f"{label} must contain non-negative integer System indices")
    if tuple(sorted(set(indices))) != indices:
        raise ValueError(f"{label} must be sorted and contain no duplicate indices")
    return indices


@dataclass(frozen=True)
class InterfaceSolveSelection:
    """Validated, ordered Group A/Group B input for a dedicated interface.

    The selections are System-global indices.  Their order is semantic (A is
    the first physical side), while each side is required to be canonical so
    the same target cannot acquire multiple cache or identity spellings.
    """

    interface_id: str
    group_a_id: str
    group_a_system_indices: tuple[int, ...]
    group_b_id: str
    group_b_system_indices: tuple[int, ...]

    def __post_init__(self):
        for name in ("interface_id", "group_a_id", "group_b_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty stable identifier")
        if self.group_a_id == self.group_b_id:
            raise ValueError("Group A and Group B identifiers must be distinct")
        group_a = _interface_indices(
            self.group_a_system_indices, "group_a_system_indices"
        )
        group_b = _interface_indices(
            self.group_b_system_indices, "group_b_system_indices"
        )
        overlap = sorted(set(group_a).intersection(group_b))
        if overlap:
            raise ValueError(
                "Group A and Group B selections overlap: "
                f"{overlap[:5]}"
            )
        object.__setattr__(self, "group_a_system_indices", group_a)
        object.__setattr__(self, "group_b_system_indices", group_b)


@dataclass(frozen=True)
class InterfaceResultHandoff:
    """Typed producer evidence returned with one live interface Result."""

    result: InterfaceResult
    identity: SharedNetworkIdentity
    dependency: ResultsDependency
    geometry_lineage: GeometryLineage


class _ProgressBridge:
    def __init__(self, progress: ProgressCallback, is_cancelled: CancellationCheck):
        self._progress = progress
        self._is_cancelled = is_cancelled

    def update_progress(
        self, process: str, value: float, network: str | None = None
    ) -> None:
        if self._is_cancelled():
            raise RuntimeError("Analysis cancelled")
        label = f"{network}: {process}" if network else process
        self._progress(label, round(value))


class VorPyBackend:
    """Run VorPy and convert its output to the viewer's stable data model."""

    supports_interface_results = True

    def __init__(self, settings: VorPySolveSettings | None = None):
        self.settings = settings or VorPySolveSettings()

    def make_interface_network(
        self, system: System, selection: InterfaceSolveSelection
    ) -> Interface:
        """Construct a genuine ordered A/B interface network from a System.

        This is intentionally a lower-level producer seam.  It does not infer
        a Workbench Result or adapt ATLAS provenance; callers must supply that
        separately once the Results dependency/context contract is selected.
        """
        if not isinstance(selection, InterfaceSolveSelection):
            raise TypeError("selection must be InterfaceSolveSelection")
        source_network_id = self.settings.source_network_id
        if source_network_id is not None and (
            not isinstance(source_network_id, str) or not source_network_id.strip()
        ):
            raise ValueError(
                "source_network_id must be null or a non-empty stable identifier"
            )

        system_mode = normalize_solvent_mode(
            getattr(system, "solvent_mode", None)
        ).value
        if system_mode != self.settings.solvent_mode:
            raise ValueError(
                "Interface system solvent_mode does not match backend settings"
            )
        system_retained = normalize_retained_system_indices(
            getattr(system, "retained_system_indices", ())
        )
        if system_retained != self.settings.retained_system_indices:
            raise ValueError(
                "Interface system retained_system_indices do not match backend settings"
            )
        existing_group_ids = {
            getattr(group, "group_id", None) for group in getattr(system, "groups", ())
        }
        if {selection.group_a_id, selection.group_b_id}.intersection(existing_group_ids):
            raise ValueError("Interface group identifiers already exist in the System")
        if any(
            getattr(interface, "interface_id", None) == selection.interface_id
            for interface in getattr(system, "ifaces", ())
        ):
            raise ValueError("Interface identifier already exists in the System")

        geometry = prepare_network_geometry(system, self.settings.solvent_mode)
        available = {
            int(index) for index in geometry.system_indices if index is not None
        }
        for label, indices in (
            ("Group A", selection.group_a_system_indices),
            ("Group B", selection.group_b_system_indices),
        ):
            missing = sorted(set(indices).difference(available))
            if missing:
                raise ValueError(
                    f"{label} contains generators excluded by solvent_mode "
                    f"{self.settings.solvent_mode}: {missing[:5]}"
                )

        group_kwargs = dict(
            sys=system,
            net_type=self.settings.network_type,
            max_vert=self.settings.max_vertices,
            box_size=self.settings.box_size,
            surf_res=self.settings.surface_resolution,
            print_metrics=False,
            make_net=False,
            build_net=False,
            mode="interface",
            solvent_mode=self.settings.solvent_mode,
        )
        group_a = Group(
            name=selection.group_a_id,
            atoms=list(selection.group_a_system_indices),
            group_id=selection.group_a_id,
            **group_kwargs,
        )
        group_b = Group(
            name=selection.group_b_id,
            atoms=list(selection.group_b_system_indices),
            group_id=selection.group_b_id,
            **group_kwargs,
        )

        # In-memory Workbench systems have no output directory.  Keep this
        # constructor bounded without changing the normal CLI output policy.
        files = getattr(system, "files", {}) or {}
        compact_before = getattr(system, "_compact_interface_workflow", None)
        if files.get("dir") is None:
            system._compact_interface_workflow = True
        try:
            interface = Interface(
                sys=system,
                group1=group_a,
                group2=group_b,
                interface_id=selection.interface_id,
            )
            interface.make_net()
            interface.net.source_network_id = source_network_id
        finally:
            if compact_before is None:
                try:
                    del system._compact_interface_workflow
                except AttributeError:
                    pass
            else:
                system._compact_interface_workflow = compact_before

        system.groups.extend((group_a, group_b))
        system.ifaces.append(interface)
        return interface

    def make_interface_result_handoff(
        self,
        interface: Interface,
        selection: InterfaceSolveSelection,
        analysis,
        *,
        system_id: str,
        result_network_id: str,
        radius_configuration_id: str,
        selection_id: str,
        group_a_selection_id: str,
        group_b_selection_id: str,
        orientation: str,
        source_identity_status: str,
        source_revision: str,
        environment: str,
        completeness_proven: bool,
    ) -> InterfaceResultHandoff:
        """Adapt one built ordered interface into the typed Results seam.

        Every identity that is not present on the solved network is explicit
        input.  This method rejects provisional producer identity rather than
        presenting calculated interface quantities as authoritative.
        """
        if not isinstance(interface, Interface):
            raise TypeError("interface must be a VorPy Interface")
        if not isinstance(selection, InterfaceSolveSelection):
            raise TypeError("selection must be InterfaceSolveSelection")
        if analysis is None:
            raise ValueError("A built interface analysis is required for typed Results")
        if source_identity_status != "AUTHORITATIVE":
            raise ValueError(
                "Typed live interface Results require authoritative source identity"
            )
        if orientation not in {"A -> B", "B -> A"}:
            raise ValueError("Interface orientation must be A -> B or B -> A")
        for label, value in (
            ("system_id", system_id),
            ("result_network_id", result_network_id),
            ("radius_configuration_id", radius_configuration_id),
            ("selection_id", selection_id),
            ("group_a_selection_id", group_a_selection_id),
            ("group_b_selection_id", group_b_selection_id),
            ("source_revision", source_revision),
            ("environment", environment),
        ):
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{label} must be a non-empty stable identifier")
        if not isinstance(completeness_proven, bool):
            raise TypeError("completeness_proven must be boolean")

        network = getattr(interface, "net", None)
        if network is None or getattr(network, "network_mode", None) != "interface":
            raise ValueError("Typed interface Results require a dedicated interface network")
        if interface.interface_id != selection.interface_id:
            raise ValueError("Interface identity does not match InterfaceSolveSelection")
        if tuple(sorted(interface.group1.ball_ndxs)) != selection.group_a_system_indices:
            raise ValueError("Interface Group A membership does not match selection")
        if tuple(sorted(interface.group2.ball_ndxs)) != selection.group_b_system_indices:
            raise ValueError("Interface Group B membership does not match selection")

        iface_groups = getattr(network, "iface_grps", None)
        generator_system_indices = getattr(network, "generator_system_indices", None)
        if (
            not isinstance(iface_groups, tuple)
            or len(iface_groups) != 2
            or generator_system_indices is None
        ):
            raise ValueError(
                "dedicated interface network is missing distinct generator-side evidence"
            )
        local_by_system = {
            int(system_index): local_index
            for local_index, system_index in enumerate(generator_system_indices)
            if system_index is not None
        }
        expected_a = {
            local_by_system[index]
            for index in selection.group_a_system_indices
            if index in local_by_system
        }
        expected_b = {
            local_by_system[index]
            for index in selection.group_b_system_indices
            if index in local_by_system
        }
        if (
            set(iface_groups[0]) != expected_a
            or set(iface_groups[1]) != expected_b
            or set(getattr(network, "group", ()) or ()) != expected_a | expected_b
        ):
            raise ValueError("Interface network membership is not the ordered A/B selection")

        source_network_id = getattr(network, "source_network_id", None)
        if not isinstance(source_network_id, str) or not source_network_id.strip():
            raise ValueError("Typed interface Results require explicit source_network_id")
        if result_network_id == source_network_id:
            raise ValueError("Result network identity must remain distinct from producer identity")

        context = network.get_solve_context_key()
        partition = {"aw": "aw", "pow": "power"}.get(context.scheme)
        if partition is None:
            raise ValueError(
                f"Typed interface Results do not support geometry scheme {context.scheme!r}"
            )
        solvent_mode = getattr(network.solvent_mode, "value", network.solvent_mode)
        if not isinstance(solvent_mode, str) or not solvent_mode:
            raise ValueError("Interface network is missing solvent_mode provenance")

        shared_identity = SharedNetworkIdentity(
            source_network_id=source_network_id,
            source_identity_status=source_identity_status,
            solve_universe_id=context.solve_universe_id,
            frame_id=context.frame_id,
            analysis_context_kind="GLOBAL_NETWORK",
            geometry_scheme=context.scheme,
            partition_scheme=partition,
            scope_mode=context.scope_mode,
            solvent_mode=solvent_mode,
            solvent_generator_digest=getattr(
                network, "solvent_generator_digest", None
            ),
            solvent_cells_complete=getattr(
                network, "solvent_cells_complete", None
            ),
            retained_system_indices=tuple(
                getattr(network, "retained_system_indices", ()) or ()
            ),
            retained_generator_ids=tuple(
                getattr(network, "retained_generator_ids", ()) or ()
            ),
            retained_generator_digest=getattr(
                network, "retained_generator_digest", None
            ),
            selection_id=selection_id,
            group_A_selection_id=group_a_selection_id,
            group_B_selection_id=group_b_selection_id,
            orientation=orientation,
        )
        dependency = ResultsDependency(
            identity=shared_identity,
            source_id=selection.interface_id,
            source_revision=source_revision,
        )

        from vorpy.src.network.aw_geometry_reuse import build_geometry_lineage_metadata

        lineage_metadata = build_geometry_lineage_metadata(
            network,
            selection_scope={
                "selection_kind": "interface",
                "interface_id": selection.interface_id,
                "group_A_id": selection.group_a_id,
                "group_B_id": selection.group_b_id,
                "group_A_system_indices": selection.group_a_system_indices,
                "group_B_system_indices": selection.group_b_system_indices,
                "selection_id": selection_id,
            },
            orientation={
                "convention": (
                    "physical_interface_A_to_B"
                    if orientation == "A -> B"
                    else "physical_interface_B_to_A"
                )
            },
            completeness={"completeness_proven": completeness_proven},
        )
        lineage_metadata["source_interface_id"] = selection.interface_id
        geometry_lineage = make_geometry_lineage(
            lineage_metadata,
            solver_evidence={
                "schema_version": 1,
                "enabled": True,
                "solve_key": context,
                "completeness": {
                    "completeness_proven": completeness_proven,
                },
            },
        )
        if geometry_lineage.origin.value != "DIRECT_SOLVE":
            raise ValueError("Explicit interface solve evidence did not produce DIRECT_SOLVE lineage")

        result_identity = ResultIdentity(
            result_kind="interface",
            system=system_id,
            frame=context.frame_id,
            network_id=result_network_id,
            network_scope="interface",
            environment=environment,
            partition=partition,
            representation="voronoi",
            radius_configuration_id=radius_configuration_id,
            interface_id=selection.interface_id,
            group_A=selection.group_a_id,
            group_B=selection.group_b_id,
            side="shared",
            orientation=(
                "group1 -> group2"
                if orientation == "A -> B"
                else "group2 -> group1"
            ),
            source_network_id=source_network_id,
            source_identity_status=source_identity_status,
            source_interface_id=selection.interface_id,
            solve_universe_id=context.solve_universe_id,
            analysis_context_kind="GLOBAL_NETWORK",
            selection_id=selection_id,
            solvent_mode=solvent_mode,
        )
        provenance = Provenance(
            "FORGE ordered interface producer",
            source_revision=source_revision,
            details={
                "shared_network_identity": shared_identity.payload(),
                "dependency_digest": dependency.dependency_digest,
                "source_solvent_provenance": network.solvent_provenance,
                "source_retained_generator_provenance": (
                    network.retained_generator_provenance
                ),
            },
        )
        result = adapt_interface_analysis(
            analysis,
            result_identity,
            provenance=provenance,
            geometry_lineage=geometry_lineage,
        )
        return InterfaceResultHandoff(
            result=result,
            identity=shared_identity,
            dependency=dependency,
            geometry_lineage=geometry_lineage,
        )

    def solve(
        self,
        source: Path | None,
        progress: ProgressCallback,
        is_cancelled: CancellationCheck,
        selected_indices: tuple[int, ...] | None = None,
        interface_selection: InterfaceSolveSelection | None = None,
        frame_index: int = 1,
        frame_ranges: tuple[tuple[int, int], ...] | None = None,
    ) -> AnalysisResult:
        if source is None:
            raise ValueError("Load a structure before running analysis")
        if source.suffix.lower() != ".pdb":
            raise ValueError("The integrated solver currently accepts PDB structures")

        started = time.perf_counter()
        progress("Reading structure", 0)
        ranges = frame_ranges or index_pdb_frames(source)
        if not 1 <= frame_index <= len(ranges):
            raise ValueError(f"Frame {frame_index} is outside 1–{len(ranges)}")
        frame_count = len(ranges)
        directory = tempfile.mkdtemp(prefix="vorpy_solve_frame_")
        solve_source = source
        if frame_count > 1:
            solve_source = Path(directory) / f"{source.stem}_frame_{frame_index:04d}.pdb"
            solve_source.write_text("".join(read_frame_lines(source, ranges[frame_index - 1])))
        system_options = {
            "file": str(solve_source),
            "gui": _ProgressBridge(progress, is_cancelled),
            "print_actions": False,
            "solvent_mode": self.settings.solvent_mode,
            "retained_system_indices": self.settings.retained_system_indices,
            "boundary_config": BoundaryConfig(
                self.settings.boundary_mode,
                self.settings.probe_radius,
                self.settings.shell_spacing,
                self.settings.shell_offset,
            ),
        }
        if frame_count > 1:
            system_options["output_directory"] = directory
        system = System(**system_options)
        # The temporary solve input contains one model, but progress should
        # identify its position in the original trajectory.
        system.frame_index = frame_index
        system.frame_count = frame_count
        system.loaded_frame_count = 1
        loaded_structure = load_pdb(solve_source)
        if loaded_structure.atoms:
            from vorpy.src.boundary import has_explicit_solvent
            solvent_present = has_explicit_solvent(loaded_structure.atoms)
        else:
            solvent_present = False
        resolved_boundary = BoundaryConfig(
            self.settings.boundary_mode,
            self.settings.probe_radius,
            self.settings.shell_spacing,
            self.settings.shell_offset,
        )
        from vorpy.src.boundary import resolve_boundary
        resolved_boundary = resolve_boundary(resolved_boundary, solvent_present)
        system.boundary_config = resolved_boundary
        if is_cancelled():
            raise RuntimeError("Analysis cancelled")

        apply_system_defaults(system, self.settings.atomic_defaults)
        interface = None
        if interface_selection is not None:
            interface = self.make_interface_network(system, interface_selection)
            network_owner = interface
        else:
            group = Group(
                system,
                name=source.stem,
                atoms=list(selected_indices or ()),
                net_type=self.settings.network_type,
                max_vert=self.settings.max_vertices,
                box_size=self.settings.box_size,
                surf_res=self.settings.surface_resolution,
                print_metrics=False,
            )
            if group not in system.groups:
                system.groups.append(group)
            network_owner = group

        system.start_run()
        if interface is not None:
            interface.build(analyze_waters=False)
        else:
            group.build()
        network = network_owner.net
        shell_surfs = shell_edges = shell_verts = ()
        if interface is None and all(
            getattr(group.net, name, None) is not None
            for name in ("surfs", "edges", "verts")
        ):
            group.get_layers(max_layers=1, group_resids=False)
            if group.layer_surfs:
                shell_surfs = group.layer_surfs[0]
            if group.layer_edges:
                shell_edges = group.layer_edges[0]
            if group.layer_verts:
                shell_verts = group.layer_verts[0]
        system.finish_run()
        if is_cancelled():
            raise RuntimeError("Analysis cancelled")

        result = loaded_structure
        boundary_info = [
            ("Mode", {BoundaryMode.NONE: "No boundary", BoundaryMode.SHELL: "Virtual solvent shell",
                      BoundaryMode.EXPLICIT: "Explicit solvent"}[resolved_boundary.mode]),
            ("Explicit solvent detected", "Yes" if solvent_present else "No"),
        ]
        if resolved_boundary.mode is BoundaryMode.SHELL:
            boundary_info.extend([
                ("Probe radius", f"{resolved_boundary.probe_radius:.3f} Å"),
                ("Shell offset", f"{resolved_boundary.shell_offset:.3f} Å"),
                ("Shell spacing", f"{resolved_boundary.shell_spacing:.3f} Å"),
                ("Boundary generator radius", f"{DEFAULT_BOUNDARY_GENERATOR_RADIUS:.3f} Å"),
                ("Boundary generators", str(len(system.boundary_generators))),
            ])
        result.info_sections["Boundary"] = boundary_info
        result.export_group = network_owner
        system.gui = None
        result.source = source
        result.frame_index = frame_index
        result.frame_count = frame_count
        result.frame_ranges = tuple(ranges)
        result.atoms = _merge_atoms(result.atoms, _atoms_from_system(system))
        result.layers = _layers_from_network(
            network,
            shell_surfs=shell_surfs,
            shell_edges=shell_edges,
            shell_verts=shell_verts,
        )
        requested_layers = {"vertices", "edges", "surfaces"}
        if not self.settings.build_vertices:
            requested_layers.discard("vertices")
        if not self.settings.build_edges:
            requested_layers.discard("edges")
        if not self.settings.build_surfaces:
            requested_layers.discard("surfaces")
        result.layers = [
            layer for layer in result.layers
            if layer.kind.lower() in requested_layers
        ]
        complete = network.balls.get("complete")
        result.complete_cells = (
            int(sum(bool(value) for value in complete.iloc[:len(system.balls)]))
            if complete is not None else 0
        )
        if network.surfs is None:
            result.surface_count = 0
        else:
            boundary_indices = set(getattr(network, "boundary_indices", ()))
            if boundary_indices and "balls" in network.surfs:
                result.surface_count = sum(
                    not bool({int(value) for value in row} & boundary_indices)
                    or bool({int(value) for value in row} - boundary_indices)
                    for row in network.surfs["balls"]
                )
            else:
                result.surface_count = len(network.surfs)
        if interface is not None:
            result.typed_results = self._interface_typed_results(
                interface, interface_selection, source, frame_index
            )
            result.typed_selection_system_ids = {
                interface_selection.group_a_id: interface_selection.group_a_system_indices,
                interface_selection.group_b_id: interface_selection.group_b_system_indices,
            }
        elif self.settings.typed_group_selections:
            from types import SimpleNamespace

            from vorpy.src.group.shared_view import attach_group_view_results

            selections = self.settings.typed_group_selections
            if isinstance(selections, Mapping):
                selection_items = tuple(
                    (str(group_id), tuple(indices))
                    for group_id, indices in selections.items()
                )
            else:
                selection_items = tuple(
                    (f"{source.stem}:group-{index}", tuple(selection))
                    for index, selection in enumerate(selections)
                )
            requested_groups = tuple(
                SimpleNamespace(ball_ndxs=selection, group_id=group_id)
                for group_id, selection in selection_items
            )
            try:
                attach_group_view_results(
                    result,
                    network,
                    requested_groups,
                    system=source.stem,
                    radius_configuration_id="vorpy-backend-radii-v1",
                    source_network_id=self.settings.source_network_id,
                )
            except NetworkViewError as error:
                # Missing producer identity is a valid unresolved Workbench
                # state; every other producer failure must remain visible.
                if (
                    self.settings.source_network_id is not None
                    or getattr(group.net, "source_network_id", None) is not None
                    or "producer-attested source_network_id is required" not in str(error)
                ):
                    raise
                result.typed_results = ()
            result.typed_selection_system_ids = dict(selection_items)
        result.elapsed_seconds = time.perf_counter() - started
        progress("Preparing viewer", 100)
        shutil.rmtree(directory, ignore_errors=True)
        return result

    def _interface_typed_results(
        self,
        interface: Interface,
        selection: InterfaceSolveSelection,
        source: Path,
        frame_index: int,
    ) -> tuple[object, ...]:
        """Adapt one explicit FORGE interface solve into ATLAS Results.

        Missing producer evidence remains an empty typed handoff; the viewer
        then presents the solved geometry as unresolved rather than certifying
        an interface from labels or the union of A/B atoms.
        """
        from vorpy.src.network.views import NetworkViewError
        from vorpy.src.results.model import stable_id

        try:
            network = interface.net
            context = network.get_solve_context_key()
            source_network_id = getattr(network, "source_network_id", None)
            if not isinstance(source_network_id, str) or not source_network_id.strip():
                return ()
            view = network.make_view(
                iface_grps=(
                    selection.group_a_system_indices,
                    selection.group_b_system_indices,
                ),
                source_network_id=source_network_id,
            )
            group_a_selection_id = stable_id(
                "view",
                {"group_id": selection.group_a_id, "system_ids": selection.group_a_system_indices},
            )
            group_b_selection_id = stable_id(
                "view",
                {"group_id": selection.group_b_id, "system_ids": selection.group_b_system_indices},
            )
            result_network_id = getattr(network, "archive_id", None) or stable_id(
                "archive", {"source": str(source), "interface_id": selection.interface_id}
            )
            environment = (
                "solvent_competing"
                if getattr(network, "solvent_mode", None) != SolventMode.MOLECULE_ONLY.value
                else "dry"
            )
            handoff = self.make_interface_result_handoff(
                interface,
                selection,
                interface.geometry_analysis,
                system_id=source.stem,
                result_network_id=result_network_id,
                radius_configuration_id="vorpy-backend-radii-v1",
                selection_id=view.selection_id,
                group_a_selection_id=group_a_selection_id,
                group_b_selection_id=group_b_selection_id,
                orientation="A -> B",
                source_identity_status="AUTHORITATIVE",
                source_revision=f"interface:{selection.interface_id}",
                environment=environment,
                completeness_proven=False,
            )
            return (handoff.result,)
        except (NetworkViewError, TypeError, ValueError, KeyError):
            return ()


def _numeric_charge(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _atoms_from_system(system: System) -> list[Atom]:
    atoms: list[Atom] = []
    for index, row in system.balls.reset_index(drop=True).iterrows():
        atoms.append(
            Atom(
                index=index,
                serial=int(row.get("num", index)) + 1,
                name=str(row.get("name", "")),
                element=str(row.get("element", "")),
                position=tuple(float(value) for value in row["loc"]),
                residue_name=str(row.get("res_name", row.get("residue", ""))),
                residue_sequence=str(row.get("res_seq", "")),
                chain=str(row.get("chain_name", row.get("chain", ""))),
                radius=float(row["rad"]),
                mass=float(row["mass"]) if row.get("mass") is not None else None,
                charge=_numeric_charge(row.get("charge")),
            )
        )
    return atoms


def _merge_atoms(display_atoms: list[Atom], solved_atoms: list[Atom]) -> list[Atom]:
    """Keep reliable PDB identities while adopting solved coordinates."""
    if len(display_atoms) != len(solved_atoms):
        return display_atoms
    return [
        Atom(
            index=display.index,
            serial=display.serial,
            name=display.name,
            element=display.element,
            position=solved.position,
            residue_name=display.residue_name,
            residue_sequence=display.residue_sequence,
            chain=display.chain,
            # The solved System is the source of truth for geometry:
            # preserve PDB identity fields, but keep the exact radius that
            # VorPy used to build the network.
            radius=solved.radius,
            mass=solved.mass,
            charge=solved.charge if solved.charge is not None else display.charge,
            source_properties={"radius": display.radius, "mass": display.mass, "charge": display.charge},
        )
        for display, solved in zip(display_atoms, solved_atoms, strict=True)
    ]


def _surface_geometry(
    surfaces, balls, target_cells=None, interpretation="magnitude"
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    points: list[np.ndarray] = []
    faces: list[np.ndarray] = []
    scalar_values = {
        "gaussian_curvature": [],
        "mean_curvature": [],
        "integrated_mean_curvature": [],
        "integrated_gaussian_curvature": [],
        "surface_energy": [],
        "distance": [],
        "inside_outside": [],
    }
    offset = 0
    for _, surface in surfaces.iterrows():
        surface_points = np.asarray(surface["points"], dtype=float).reshape((-1, 3))
        surface_faces = np.asarray(surface["tris"], dtype=np.int64).reshape((-1, 3))
        face_count = len(surface_faces)
        points.extend(surface_points)
        faces.extend(surface_faces + offset)
        offset += len(surface_points)

        for scheme, column in (
            ("gaussian_curvature", "gauss_tri_curvs"),
            ("mean_curvature", "mean_tri_curvs"),
        ):
            values = np.asarray(surface.get(column, []), dtype=float).ravel()
            if len(values) != face_count:
                fallback = "gauss_curv" if scheme == "gaussian_curvature" else "mean_curv"
                values = np.full(face_count, float(surface.get(fallback, 0.0) or 0.0))
            scalar_values[scheme].extend(values)

        mean_values = np.asarray(scalar_values["mean_curvature"][-face_count:], dtype=float)
        scalar_values["surface_energy"].extend(2.0 * np.square(mean_values))

        for display_name, core_name in (
            ("integrated_mean_curvature", "int_mean_curv"),
            ("integrated_gaussian_curvature", "int_gauss_curv"),
        ):
            value = component_value(
                surface, "surface", core_name, target_cells, mode=interpretation
            )
            scalar_values[display_name].extend(
                np.full(face_count, 0.0 if value is None else float(value))
            )

        center = np.asarray(surface.get("loc", surface_points.mean(axis=0)), dtype=float)
        point_distances = np.linalg.norm(surface_points - center, axis=1)
        scalar_values["distance"].extend(np.max(point_distances[surface_faces], axis=1))

        inside = np.zeros(len(surface_points), dtype=bool)
        defining = surface.get("balls", [])
        if balls is not None and len(defining):
            matches = balls
            if "num" in balls:
                matches = balls.loc[balls["num"] == int(defining[0])]
            elif int(defining[0]) in balls.index:
                matches = balls.loc[[int(defining[0])]]
            if len(matches):
                atom = matches.iloc[0]
                inside = np.linalg.norm(
                    surface_points - np.asarray(atom["loc"], dtype=float), axis=1
                ) < float(atom["rad"])
        scalar_values["inside_outside"].extend(
            np.all(inside[surface_faces], axis=1).astype(float)
        )

    return (
        np.asarray(points, dtype=float).reshape((-1, 3)),
        np.asarray(faces, dtype=np.int64).reshape((-1, 3)),
        {key: np.asarray(values, dtype=float) for key, values in scalar_values.items()},
    )


def _edge_geometry(
    edges, target_cells=None, interpretation="magnitude"
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    points: list[np.ndarray] = []
    lines: list[tuple[int, int]] = []
    scalars = {
        "integrated_mean_curvature": [],
        "integrated_gaussian_curvature": [],
    }
    for _, edge in edges.iterrows():
        edge_points = np.asarray(edge["points"], dtype=float)
        start = len(points)
        points.extend(edge_points)
        segment_count = max(len(edge_points) - 1, 0)
        lines.extend((start + index, start + index + 1) for index in range(segment_count))
        for display_name, core_name in (
            ("integrated_mean_curvature", "int_mean_curv"),
            ("integrated_gaussian_curvature", "int_gauss_curv"),
        ):
            value = component_value(
                edge, "edge", core_name, target_cells, mode=interpretation
            )
            scalars[display_name].extend(
                np.full(segment_count, 0.0 if value is None else float(value))
            )
    line_array = np.asarray(lines, dtype=np.int64).reshape((-1, 2))
    scalar_arrays = {
        key: np.asarray(values, dtype=float)
        for key, values in scalars.items()
    }
    for name, values in scalar_arrays.items():
        if len(values) != len(line_array):
            raise ValueError(
                f"Edge scalar {name!r} has {len(values)} values for "
                f"{len(line_array)} rendered line segments"
            )

    return (
        np.asarray(points, dtype=float).reshape((-1, 3)),
        line_array,
        scalar_arrays,
    )


def _vertex_geometry(
    network, vertices, target_cells=None, interpretation="magnitude"
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    points = np.asarray(list(vertices["loc"]), dtype=float).reshape((-1, 3))
    mean_values, gauss_values = [], []
    for _, vertex in vertices.iterrows():
        mean_value = mean_vertex_display_value(
            network, vertex, target_cells, mode=interpretation
        )
        gauss_value = component_value(
            vertex, "vertex", "int_gauss_curv", target_cells, mode=interpretation
        )
        mean_values.append(0.0 if mean_value is None else float(mean_value))
        gauss_values.append(0.0 if gauss_value is None else float(gauss_value))
    return points, {
        "integrated_mean_curvature": np.asarray(mean_values, dtype=float),
        "integrated_gaussian_curvature": np.asarray(gauss_values, dtype=float),
    }


def _layers_from_network(
    network,
    shell_surfs=(),
    shell_edges=(),
    shell_verts=(),
) -> list[GeometryLayer]:
    layers: list[GeometryLayer] = []
    raw_group = getattr(network, "group", None)
    target_cells = () if raw_group is None else tuple(int(value) for value in raw_group)
    boundary_indices = set(getattr(network, "boundary_indices", ()))

    def split_rows(frame):
        """Separate molecular, boundary-facing and shell-only topology rows."""
        if frame is None or frame.empty or not boundary_indices or "balls" not in frame:
            return frame, frame.iloc[0:0] if frame is not None else frame
        molecular, boundary = [], []
        for index, row in frame.iterrows():
            try:
                adjacent = {int(value) for value in row["balls"]}
            except (TypeError, ValueError):
                molecular.append(index)
                continue
            artificial = adjacent & boundary_indices
            if not artificial:
                molecular.append(index)
            elif adjacent - boundary_indices:
                boundary.append(index)
            # Rows surrounded only by generated sites are internal shell geometry.
        return frame.loc[molecular], frame.loc[boundary]

    physical_edges, boundary_edges = split_rows(network.edges)
    physical_vertices, boundary_vertices = split_rows(network.verts)
    physical_surfaces, boundary_surfaces = split_rows(network.surfs)

    if network.edges is not None and "points" in network.edges:
        points, lines, scalars = _edge_geometry(
            physical_edges, target_cells, interpretation="magnitude"
        )
        if not physical_edges.empty:
            layers.append(GeometryLayer(
                "Voronoi edges", "edges", points, lines, color="#55a9d9",
                cell_scalars=scalars, interpretation="magnitude",
            ))
        shell_edge_indices = network.edges.iloc[list(shell_edges)].index
        shell_edge_rows = physical_edges.loc[physical_edges.index.intersection(shell_edge_indices)]
        if not shell_edge_rows.empty:
            points, lines, scalars = _edge_geometry(
                shell_edge_rows, target_cells, interpretation="boundary"
            )
            layers.append(GeometryLayer(
                "Voronoi shell edges", "edges", points, lines, color="#42d6c7",
                cell_scalars=scalars, interpretation="boundary",
            ))
        if not boundary_edges.empty:
            points, lines, scalars = _edge_geometry(
                boundary_edges, target_cells, interpretation="boundary"
            )
            layers.append(GeometryLayer(
                "Virtual boundary edges", "edges", points, lines, color="#9467bd",
                cell_scalars=scalars, interpretation="boundary", visible=False,
            ))

    if network.verts is not None and "loc" in network.verts:
        points, scalars = _vertex_geometry(
            network, physical_vertices, target_cells, interpretation="magnitude"
        )
        if not physical_vertices.empty:
            layers.append(GeometryLayer(
                "Voronoi vertices", "vertices", points, color="#efb84f",
                cell_scalars=scalars, interpretation="magnitude",
            ))
        shell_vertex_indices = network.verts.iloc[list(shell_verts)].index
        shell_vertex_rows = physical_vertices.loc[physical_vertices.index.intersection(shell_vertex_indices)]
        if not shell_vertex_rows.empty:
            points, scalars = _vertex_geometry(
                network, shell_vertex_rows, target_cells, interpretation="boundary"
            )
            layers.append(GeometryLayer(
                "Voronoi shell vertices", "vertices", points, color="#f29f67",
                cell_scalars=scalars, interpretation="boundary",
            ))
        if not boundary_vertices.empty:
            points, scalars = _vertex_geometry(
                network, boundary_vertices, target_cells, interpretation="boundary"
            )
            layers.append(GeometryLayer(
                "Virtual boundary vertices", "vertices", points, color="#a98bd4",
                cell_scalars=scalars, interpretation="boundary", visible=False,
            ))

    if network.surfs is not None and "points" in network.surfs and "tris" in network.surfs:
        points, faces, scalars = _surface_geometry(
            physical_surfaces, getattr(network, "balls", None), target_cells,
            interpretation="magnitude",
        )
        if not physical_surfaces.empty:
            layers.append(GeometryLayer(
                "Voronoi surfaces", "surfaces", points=points, faces=faces,
                color="#4f9fcf", cell_scalars=scalars, opacity=0.45,
                visible=False, interpretation="magnitude",
            ))
        shell_surface_indices = network.surfs.iloc[list(shell_surfs)].index
        shell_surface_rows = physical_surfaces.loc[physical_surfaces.index.intersection(shell_surface_indices)]
        if not shell_surface_rows.empty:
            points, faces, scalars = _surface_geometry(
                shell_surface_rows, getattr(network, "balls", None), target_cells,
                interpretation="boundary",
            )
            layers.append(GeometryLayer(
                "Voronoi shell surfaces", "surfaces", points=points, faces=faces,
                color="#806df0", cell_scalars=scalars, opacity=0.45,
                interpretation="boundary",
            ))
        if not boundary_surfaces.empty:
            points, faces, scalars = _surface_geometry(
                boundary_surfaces, getattr(network, "balls", None), target_cells,
                interpretation="boundary",
            )
            layers.append(GeometryLayer(
                "Virtual boundary surfaces", "surfaces", points=points, faces=faces,
                color="#9467bd", cell_scalars=scalars, opacity=0.45,
                visible=False, interpretation="boundary",
            ))
    return layers
