import os
from time import perf_counter
import numpy as np
from copy import deepcopy
from vorpy.src.group import Group
from vorpy.src.network import Network
from vorpy.src.boundary import prepare_network_geometry
from vorpy.src.interface.export import interface_exports
from vorpy.src.interface.water import analyze_interface_waters
from vorpy.src.interface.water import build_buried_water_groups
from vorpy.src.output.layout import interface_folder_name, unique_output_directory


class Interface:
    """
    Represents the shared interface between two groups.

    Each side owns a partial Group whose network contains only the
    vertices, edges, and surfaces relevant to the interface.
    """

    def __init__(self, sys, group1, group2=None, name=None, interface_id=None):
        self.sys = sys

        # Original group definitions. These are references, not copies.
        self.group1 = group1
        self.group1_name = None
        self._surrounding_interface_requested = group2 is None
        if group2 is None:
            group2 = self._make_surrounding_group(group1)
        self.group2 = group2
        self.interface_id = (
                interface_id
                or self.make_interface_id(self.group1, self.group2)
        )
        self.settings = self._get_settings()
        self.dir = None
        self.name = name or self._make_name()

        # The interface owns its own network.
        self.net = None
        # Full-system mode replaces the dedicated network with a read-only
        # source view while retaining explicit provenance on the record.
        self.source_network = None
        self.source_view = None
        self.geometry_source = None
        self.geometry_analysis = None

        # Cached index collections defining the two interface sides.
        self.group1_indices = None
        self.group2_indices = None

        # Partial groups created specifically for this interface
        self.partial_group1 = None
        self.partial_group2 = None

        self._register_with_groups()

        self.water_groups = []
        self.buried_water_groups = []
        self.water_geometries = []
        self.water_topology = None
        self._buried_water_exports_complete = False

    def _make_name(self):
        second_name = "surrounding" if self.group2 is None else self.group2.name
        self.name = f"{self.group1.name}_{second_name}_interface"
        return self.name

    def set_dir(self):
        root = self.sys.files['dir']
        occupied = [getattr(interface, 'dir', None)
                    for interface in getattr(self.sys, 'ifaces', ()) or ()
                    if getattr(interface, 'dir', None)]
        self.dir = str(unique_output_directory(
            root, interface_folder_name(self.group1, self.group2), occupied,
        ))
        os.makedirs(self.dir, exist_ok=False)

    def _get_settings(self):
        """
        Return independent settings for the interface network.

        Group 1 provides the initial construction settings, but the dictionary
        is copied so interface-specific changes do not affect the source group.
        """
        return deepcopy(self.group1.settings)

    def _make_surrounding_group(self, source_group):
        """
        Create a concrete Group containing the spatially restricted
        surrounding balls for source_group.
        """

        surrounding_indices = self.get_surrounding_indices(
            source_group.ball_ndxs
        )

        surrounding_group = Group(
            sys=self.sys,
            name=f"{source_group.name}_surrounding",
            settings=deepcopy(source_group.settings),
            make_net=False,
            build_net=False,
            mode="interface",
        )

        # The constructor may populate a default selection, so explicitly
        # replace it with the spatially restricted surrounding selection.
        surrounding_group.ball_ndxs = sorted(
            set(int(index) for index in surrounding_indices)
        )

        surrounding_group.group_id = (
            f"{source_group.group_id}__surrounding"
        )

        return surrounding_group

    def make_net(self, verts=None):
        """
        Create the interface-specific Network without building its topology.

        The Network uses the complete system geometry but restricts vertex
        discovery to vertices involving balls from both interface sides.
        """
        if self.dir is None and not getattr(self.sys, '_compact_interface_workflow', False):
            self.set_dir()

        self.group1_indices = set(self.group1.ball_ndxs)

        if self.group2 is None:
            self.group2_indices = self.get_surrounding_indices(self.group1_indices)
            self.group2_name = "surrounding"
        else:
            self.group2_indices = set(self.group2.ball_ndxs)
            self.group2_name = self.group2.name

        geometry = prepare_network_geometry(
            self.sys,
            getattr(self.group1, "solvent_mode", None),
        )
        system_to_local = {
            int(system_index): local_index
            for local_index, system_index in enumerate(geometry.system_indices)
            if system_index is not None
        }
        missing_group1 = sorted(set(self.group1_indices) - set(system_to_local))
        if missing_group1:
            raise ValueError(
                "Interface group 1 contains generators excluded by solvent_mode: "
                f"{missing_group1[:5]}"
            )
        # A default molecule-only interface may be requested against the
        # surrounding System selection; excluded physical solvent is simply
        # absent from that surrounding side. Explicit solvent modes retain it.
        group1_local = {
            system_to_local[index] for index in self.group1_indices
        }
        group2_local = {
            system_to_local[index]
            for index in self.group2_indices
            if index in system_to_local
        }
        interface_local_indices = sorted(group1_local | group2_local)
        self.net = Network(
            # Retain the selected solve universe so surrounding generators can
            # participate in geometric validity checks.
            locs=geometry.locations,
            rads=geometry.radii,
            masses=geometry.masses,

            # Balls whose interface topology is being represented.
            group=interface_local_indices,

            # Original side membership used during interface vertex filtering.
            iface_grps=(group1_local, group2_local),

            group_name=self.name,
            settings=self.settings,
            sort_balls=False,
            verts=verts,
            system=self.sys,
            generator_system_indices=geometry.system_indices,
            generator_stable_ids=geometry.stable_ids,
            solvent_local_indices=geometry.solvent_local_indices,
            solvent_generator_ids=geometry.solvent_generator_ids,
            solvent_mode=geometry.solvent_mode,
            retained_system_indices=geometry.retained_system_indices,
            retained_generator_ids=geometry.retained_generator_ids,
        )
        self.net.boundary_indices = geometry.boundary_indices
        self.net.boundary_config = geometry.boundary_config
        self.net.boundary_mode = getattr(self.sys, "boundary_mode", None)
        self.net.boundary_generators = geometry.boundary_generators
        if geometry.boundary_indices:
            self.net.balls["is_boundary_generator"] = False
            self.net.balls.loc[list(geometry.boundary_indices), "is_boundary_generator"] = True
        self.net.refresh_solvent_provenance()

        self._update_group_metadata(
            network_created=True,
            built=False,
        )

        return self.net

    def build(self, *, analyze_waters=True, calculate_curvature=None):
        """
        Create and build this interface's dedicated network.
        """
        if self.net is None:
            self.make_net()

        if calculate_curvature is None:
            calculate_curvature = getattr(self.sys, '_interface_curvature_requested', True)
        if calculate_curvature:
            # Preserve the historical no-argument call for legacy Network
            # implementations and test doubles.
            self.net.build()
        else:
            self.net.build(calculate_curvature=False)
        self._record_interface_geometry_evidence()
        self.geometry_source = 'DIRECT_SOLVE'
        self.geometry_analysis = None
        self._geometry_analysis_key = None
        self._geometry_source_key = None

        # Make completed topology available to later pairwise interface builds.
        self.sys.cache_interface_geometry(self)

        # Cache the solved geometry independently of later buried-water solves.
        analysis_started = perf_counter()
        update_progress = getattr(self.sys, 'update_progress', None)
        if update_progress is not None:
            update_progress(
                process='Interface analysis | Preparing cached physical interface',
                progress=0.0,
                network=self.name,
            )
        self.analyze_geometry(
            calculate_curvature=calculate_curvature
        )
        self._interface_stage_timings = getattr(self, '_interface_stage_timings', {})
        self._interface_stage_timings['physical_interface'] = perf_counter() - analysis_started

        # Water classification and buried-water solves are optional.  Keep the
        # direct API's historical default, while the normal interface CLI
        # opts out unless explicitly asked for them.
        if not analyze_waters:
            self.water_geometries = []
            self.water_topology = None
            self.buried_water_groups = []
            self.water_groups = self.buried_water_groups
            self._interface_stage_timings['interface_water'] = 0.0
            self._interface_stage_timings['buried_water'] = 0.0
            self._update_group_metadata(network_created=True, built=True)
            return

        # Stage 1: classify all strict interface waters from the completed
        # interface topology.  This pass is topology-only and does not build
        # separate water networks.
        water_started = perf_counter()
        if update_progress is not None:
            update_progress(
                process='Interface analysis | Classifying interface waters',
                progress=0.0,
                network=self.name,
            )
        self.water_geometries = analyze_interface_waters(
            iface=self,
        )
        self._interface_stage_timings['interface_water'] = perf_counter() - water_started

        topology = self.water_topology or {}
        waters = topology.get("waters", {})
        class_counts = {"buried": 0, "semi_buried": 0, "peripheral": 0}
        for water in waters.values():
            burial_class = water.get("burial_class")
            if burial_class in class_counts:
                class_counts[burial_class] += 1

        buried_sets = topology.get("buried_cycle_sets", [])
        non_buried_count = (
            class_counts["peripheral"] + class_counts["semi_buried"]
        )
        # Keep water classification in the interface info export; normal CLI output uses the
        # single system-wide progress stream rather than standalone debug prints.
        if hasattr(self.sys, "update_progress"):
            self.sys.update_progress(
                process=(
                    f"Interface waters classified: {len(waters)} total, "
                    f"{non_buried_count} non-buried, "
                    f"{class_counts['buried']} buried, "
                    f"{len(buried_sets)} buried groups"
                ),
                progress=100.0,
                network=self.name,
            )

        # Stage 2: solve every CLOSED buried-water set as a normal complete
        # Group against the parent System.  These networks live under
        # <interface>/waters/<buried_group>/ and provide V/A/C/Q/X plus cell
        # geometry using the standard VorPy machinery.
        buried_started = perf_counter()
        if update_progress is not None:
            update_progress(
                process='Interface analysis | Analyzing buried water groups',
                progress=0.0,
                network=self.name,
            )
        self.buried_water_groups = build_buried_water_groups(self)
        self._interface_stage_timings['buried_water'] = perf_counter() - buried_started
        if update_progress is not None:
            update_progress(
                process=(
                    'Interface analysis | Buried water groups complete '
                    f'({self._interface_stage_timings["buried_water"]:.2f} s)'
                ),
                progress=100.0,
                network=self.name,
            )

        # Backward-compatible alias for code that already expects water_groups.
        self.water_groups = self.buried_water_groups

        self._update_group_metadata(
            network_created=True,
            built=True,
        )

    def build_selective_geometry(self, *, calculate_curvature=False):
        """Build only the selected physical A/B interface geometry.

        The Network keeps every generator in the prepared solve universe as a
        geometric competitor.  Its ``iface_grps`` filter limits constructed
        topology to vertices, edges, and surfaces spanning the two requested
        groups.  This entrypoint deliberately stops before interface analysis
        and water/buried-water work; callers can request those stages later.
        """
        fallback = self.selective_geometry_fallback()
        if fallback is not None:
            return {
                'schema_version': 1,
                'producer': 'FULL_SOURCE_FALLBACK',
                'fallback': fallback,
            }
        if self.net is None:
            self.make_net()
        if getattr(self.net, "iface_grps", None) is None:
            raise ValueError("Selective interface geometry requires two groups.")

        started = perf_counter()
        self.net.build(calculate_curvature=bool(calculate_curvature))
        evidence = self._record_interface_geometry_evidence()
        self.geometry_source = 'DIRECT_SOLVE'
        self.geometry_analysis = None
        self._geometry_analysis_key = None
        self._geometry_source_key = None
        self._interface_stage_timings = getattr(self, '_interface_stage_timings', {})
        self._interface_stage_timings['selective_geometry'] = perf_counter() - started

        # Keep the same cache and group lifecycle as the historical direct
        # builder, without invoking unrelated interface analysis.
        self.sys.cache_interface_geometry(self)
        self._update_group_metadata(network_created=True, built=True)
        return evidence

    def selective_geometry_fallback(self):
        """Return a structured full-source fallback reason, if required.

        Explicit A/B interfaces are supported by the direct producer.  A
        surrounding interface is not: its legacy spatial cutoff is not a
        proof that every competing generator was retained.
        """
        if getattr(self, '_surrounding_interface_requested', False):
            return {
                'required': True,
                'mode': 'FULL_SOURCE_MATERIALIZED_VIEW',
                'reason_code': 'SURROUNDING_COMPETITOR_CLOSURE_UNPROVEN',
                'reason': (
                    'The surrounding-group cutoff is not a proof of the full '
                    'competitor universe.'
                ),
            }

        network = getattr(self, 'net', None)
        settings = (
            getattr(self, 'settings', None)
            or getattr(network, 'settings', None)
            or {}
        )
        scheme = str(settings.get('net_type', 'aw')).strip().lower()
        if scheme not in {'prm', 'pow', 'aw'}:
            return {
                'required': True,
                'mode': 'FULL_SOURCE_MATERIALIZED_VIEW',
                'reason_code': 'UNSUPPORTED_INTERFACE_SCHEME',
                'reason': f'No selective producer is registered for {scheme!r}.',
            }

        if network is not None:
            groups = getattr(network, 'iface_grps', None)
            if groups is None or len(groups) != 2:
                return {
                    'required': True,
                    'mode': 'FULL_SOURCE_MATERIALIZED_VIEW',
                    'reason_code': 'EXPLICIT_INTERFACE_GROUPS_UNAVAILABLE',
                    'reason': 'The direct network lacks two explicit interface sides.',
                }
        return None

    def _record_interface_geometry_evidence(self):
        """Record auditable local support evidence for a direct solve.

        A direct interface network intentionally does not contain complete
        competing cells.  Therefore local finite incidence is reported
        separately from cell completeness; no incomplete cell is promoted.
        The payload is stored in archive provenance so a checkpoint can be
        inspected without rebuilding geometry.
        """
        net = self.net
        iface_groups = getattr(net, 'iface_grps', None)
        if iface_groups is None or len(iface_groups) != 2:
            raise ValueError('Interface geometry evidence requires two groups.')

        group_a = set(int(value) for value in iface_groups[0])
        group_b = set(int(value) for value in iface_groups[1])
        balls = net.balls
        verts = net.verts
        edges = net.edges
        surfs = net.surfs
        ball_count = len(balls)
        vertex_count = len(verts)
        edge_count = len(edges)

        system_indices = tuple(getattr(net, 'generator_system_indices', ()) or ())
        stable_ids = tuple(getattr(net, 'generator_stable_ids', ()) or ())

        def generator_record(index):
            index = int(index)
            system_id = (
                system_indices[index]
                if 0 <= index < len(system_indices) else None
            )
            stable_id = (
                stable_ids[index]
                if 0 <= index < len(stable_ids) else None
            )
            if stable_id is None and 0 <= index < ball_count and 'stable_id' in balls.columns:
                value = balls.iloc[index].get('stable_id')
                stable_id = None if value is None else str(value)
            return {
                'local_index': index,
                'system_id': None if system_id is None else int(system_id),
                'stable_id': stable_id,
            }

        def references(row, field, limit):
            value = row.get(field, ())
            if value is None:
                return (), False
            try:
                values = tuple(int(item) for item in value)
            except (TypeError, ValueError, OverflowError):
                return (), False
            return values, all(0 <= item < limit for item in values)

        def finite_geometry(value):
            if value is None:
                return False
            try:
                values = np.asarray(value, dtype=float)
            except (TypeError, ValueError):
                return False
            return bool(values.size and np.all(np.isfinite(values)))

        def finite_area(value):
            try:
                return bool(np.isfinite(float(value)))
            except (TypeError, ValueError, OverflowError):
                return False

        complete_column = 'complete' in balls.columns

        def cell_status(dependencies):
            dependencies = tuple(sorted(set(int(value) for value in dependencies)))
            invalid = [value for value in dependencies if value < 0 or value >= ball_count]
            incomplete = []
            if complete_column:
                incomplete = [
                    value for value in dependencies
                    if value not in invalid and not bool(balls.iloc[value].get('complete'))
                ]
            else:
                incomplete = list(dependencies)
            return dependencies, invalid, incomplete

        def ordered_pairs(ball_indices):
            pairs = []
            for left in sorted(set(ball_indices).intersection(group_a)):
                for right in sorted(set(ball_indices).intersection(group_b)):
                    left_id = generator_record(left)['system_id']
                    right_id = generator_record(right)['system_id']
                    pairs.append((left_id, right_id))
            return tuple(pairs)

        selected_surface_positions = []
        surface_records = []
        support_generators = set()
        unsupported_surface_reasons = []
        for position, row in surfs.iterrows():
            position = int(position)
            ball_refs, balls_valid = references(row, 'balls', ball_count)
            crosses = bool(set(ball_refs).intersection(group_a)) and bool(
                set(ball_refs).intersection(group_b)
            )
            if not crosses:
                continue

            edge_refs, edges_valid = references(row, 'edges', edge_count)
            vertex_refs, vertices_valid = references(row, 'verts', vertex_count)
            reasons = []
            dependency_generators = set(ball_refs)
            if not balls_valid:
                reasons.append('invalid_surface_generator_incidence')
            if not edges_valid:
                reasons.append('invalid_surface_edge_incidence')
            if not vertices_valid:
                reasons.append('invalid_surface_vertex_incidence')
            surface_vertices = set(vertex_refs)
            for edge_index in edge_refs:
                edge_row = edges.iloc[edge_index] if edges_valid else None
                edge_balls, edge_balls_valid = references(
                    edge_row, 'balls', ball_count
                ) if edge_row is not None else ((), False)
                edge_vertices, edge_vertices_valid = references(
                    edge_row, 'verts', vertex_count
                ) if edge_row is not None else ((), False)
                dependency_generators.update(edge_balls)
                if not edge_balls_valid or not edge_vertices_valid:
                    reasons.append(f'invalid_edge_incidence:{edge_index}')
                if not set(edge_vertices).issubset(surface_vertices):
                    reasons.append(f'open_surface_edge:{edge_index}')
            for vertex_index in vertex_refs:
                vertex_row = verts.iloc[vertex_index] if vertices_valid else None
                vertex_balls, vertex_balls_valid = references(
                    vertex_row, 'balls', ball_count
                ) if vertex_row is not None else ((), False)
                dependency_generators.update(vertex_balls)
                if not vertex_balls_valid:
                    reasons.append(f'invalid_vertex_incidence:{vertex_index}')
                if vertex_row is not None and not finite_geometry(vertex_row.get('loc')):
                    reasons.append(f'nonfinite_vertex_geometry:{vertex_index}')

            dependencies, invalid_cells, incomplete_cells = cell_status(
                dependency_generators
            )
            support_generators.update(dependencies)
            if invalid_cells:
                reasons.append('invalid_generator_dependency')
            if incomplete_cells:
                reasons.append('generator_cell_completeness_not_proven')
            bounded = finite_geometry(row.get('points')) and finite_area(row.get('sa'))
            if not bounded:
                reasons.append('nonfinite_surface_geometry')
            supported = not any(
                reason.startswith(('invalid_', 'open_', 'nonfinite_'))
                for reason in reasons
            )
            complete = supported and not incomplete_cells
            if supported:
                selected_surface_positions.append(position)
            else:
                unsupported_surface_reasons.extend(reasons)
            surface_records.append({
                'source_id': position,
                'topology_index': position,
                'generator_local_indices': tuple(ball_refs),
                'generator_system_ids': tuple(
                    generator_record(value)['system_id'] for value in ball_refs
                ),
                'generator_stable_ids': tuple(
                    generator_record(value)['stable_id'] for value in ball_refs
                ),
                'generator_pairs_A_to_B': ordered_pairs(ball_refs),
                'incidence': {'edges': tuple(edge_refs), 'vertices': tuple(vertex_refs)},
                'supported': bool(supported),
                'complete': bool(complete),
                'bounded': bool(bounded),
                'support_generator_local_indices': dependencies,
                'incomplete_generator_local_indices': tuple(incomplete_cells),
                'reasons': tuple(sorted(set(reasons))),
                'area': float(row.get('sa')) if bounded else None,
            })

        selected_edges = set()
        selected_vertices = set()
        for position in selected_surface_positions:
            row = surfs.iloc[position]
            edge_refs, _ = references(row, 'edges', edge_count)
            vertex_refs, _ = references(row, 'verts', vertex_count)
            selected_edges.update(edge_refs)
            selected_vertices.update(vertex_refs)

        edge_records = []
        for position in sorted(selected_edges):
            row = edges.iloc[position]
            ball_refs, balls_valid = references(row, 'balls', ball_count)
            vertex_refs, vertices_valid = references(row, 'verts', vertex_count)
            dependencies, invalid_cells, incomplete_cells = cell_status(ball_refs)
            reasons = []
            if not balls_valid:
                reasons.append('invalid_edge_generator_incidence')
            if not vertices_valid:
                reasons.append('invalid_edge_vertex_incidence')
            if not finite_geometry(row.get('points')):
                reasons.append('nonfinite_edge_geometry')
            support_generators.update(dependencies)
            reasons.extend(
                'generator_cell_completeness_not_proven' for _ in incomplete_cells
            )
            supported = not any(reason.startswith(('invalid_', 'nonfinite_')) for reason in reasons)
            edge_records.append({
                'source_id': int(position),
                'topology_index': int(position),
                'generator_local_indices': tuple(ball_refs),
                'generator_system_ids': tuple(
                    generator_record(value)['system_id'] for value in ball_refs
                ),
                'generator_stable_ids': tuple(
                    generator_record(value)['stable_id'] for value in ball_refs
                ),
                'incidence': {'vertices': tuple(vertex_refs), 'surfaces': tuple(
                    references(row, 'surfs', len(surfs))[0]
                )},
                'supported': bool(supported),
                'complete': bool(supported and not incomplete_cells),
                'bounded': finite_geometry(row.get('points')),
                'support_generator_local_indices': dependencies,
                'incomplete_generator_local_indices': tuple(incomplete_cells),
                'reasons': tuple(sorted(set(reasons))),
            })

        vertex_records = []
        for position in sorted(selected_vertices):
            row = verts.iloc[position]
            ball_refs, balls_valid = references(row, 'balls', ball_count)
            dependencies, invalid_cells, incomplete_cells = cell_status(ball_refs)
            reasons = []
            if not balls_valid:
                reasons.append('invalid_vertex_generator_incidence')
            if not finite_geometry(row.get('loc')):
                reasons.append('nonfinite_vertex_geometry')
            support_generators.update(dependencies)
            supported = not any(reason.startswith(('invalid_', 'nonfinite_')) for reason in reasons)
            vertex_records.append({
                'source_id': int(position),
                'topology_index': int(position),
                'generator_local_indices': tuple(ball_refs),
                'generator_system_ids': tuple(
                    generator_record(value)['system_id'] for value in ball_refs
                ),
                'generator_stable_ids': tuple(
                    generator_record(value)['stable_id'] for value in ball_refs
                ),
                'incidence': {
                    'edges': references(row, 'edges', edge_count)[0],
                    'surfaces': references(row, 'surfs', len(surfs))[0],
                },
                'supported': bool(supported),
                'complete': bool(supported and not incomplete_cells),
                'bounded': finite_geometry(row.get('loc')),
                'support_generator_local_indices': dependencies,
                'incomplete_generator_local_indices': tuple(incomplete_cells),
                'reasons': tuple(sorted(set(reasons))),
            })

        all_features = surface_records + edge_records + vertex_records
        feature_complete = all(bool(record['complete']) for record in all_features)
        feature_supported = all(bool(record['supported']) for record in all_features)
        completeness_reasons = list(unsupported_surface_reasons)
        if not complete_column:
            completeness_reasons.append('generator_cell_completeness_not_recorded')
        elif not feature_complete:
            completeness_reasons.append('dedicated_interface_cell_closure_not_proven')
        if not selected_surface_positions:
            completeness_reasons.append('no_supported_interface_surface')

        def count(name, fallback):
            value = getattr(net, name, None)
            return fallback if not isinstance(value, dict) else value

        context_id = None
        try:
            context_id = net.get_solve_context_key().digest
        except (AttributeError, TypeError, ValueError):
            completeness_reasons.append('solve_context_unavailable')

        evidence = {
            'schema_version': 1,
            'producer': 'DIRECT_SOLVE',
            'source_network_id': getattr(net, 'source_network_id', None),
            'solve_context_id': context_id,
            'selection': {
                'group_a_local_indices': tuple(sorted(group_a)),
                'group_b_local_indices': tuple(sorted(group_b)),
                'group_a_system_ids': tuple(
                    generator_record(value)['system_id'] for value in sorted(group_a)
                ),
                'group_b_system_ids': tuple(
                    generator_record(value)['system_id'] for value in sorted(group_b)
                ),
                'orientation': 'A -> B',
            },
            'competitor_universe': {
                'generator_count': ball_count,
                'physical_generator_count': sum(
                    value is not None for value in system_indices
                ),
                'boundary_indices': tuple(getattr(net, 'boundary_indices', ()) or ()),
                'solvent_mode': getattr(net, 'solvent_mode', None),
                'generator_system_indices': system_indices,
                'generator_stable_ids': stable_ids,
            },
            'support_generator_local_indices': tuple(sorted(support_generators)),
            'support_generators': tuple(
                generator_record(value) for value in sorted(support_generators)
            ),
            'counts': {
                'constructed_vertices': len(verts),
                'constructed_edges': len(edges),
                'constructed_surfaces': len(surfs),
                'selected_surfaces': len(surface_records),
                'supported_surfaces': len(selected_surface_positions),
                'selected_edges': len(edge_records),
                'selected_vertices': len(vertex_records),
                'surface_counts': count('surface_counts', {}),
                'build_net_counts': count('build_net_counts', {}),
            },
            'completeness': {
                'proven': bool(feature_supported and feature_complete and not completeness_reasons),
                'scope': 'direct_interface_local_geometry',
                'feature_support_proven': bool(feature_supported),
                'feature_completeness_proven': bool(feature_complete),
                'reasons': tuple(sorted(set(completeness_reasons))),
            },
            'surfaces': tuple(surface_records),
            'edges': tuple(edge_records),
            'vertices': tuple(vertex_records),
        }
        net.interface_geometry_evidence = evidence
        provenance = dict(getattr(net, 'archive_provenance', {}) or {})
        provenance['interface_geometry_evidence'] = evidence
        provenance['geometry_source'] = 'DIRECT_SOLVE'
        net.archive_provenance = provenance
        return evidence

    def analyze_geometry(self, alpha=None, *, refresh=False, calculate_curvature=True):
        """Return the cached alpha-selected physical interface analysis."""
        from vorpy.src.interface.geometry_analysis import analyze_interface_geometry
        if self.net is None:
            raise ValueError('Build the interface network before analyzing geometry')
        return analyze_interface_geometry(
            self, alpha=alpha, refresh=refresh,
            calculate_curvature=calculate_curvature,
        )

    def _register_with_groups(self):
        self.group1.register_interface(
            interface_id=self.interface_id,
            interface=self,
            other_group=self.group2,
            side="group1",
            interface_name=self.name,
        )

        if self.group2 is not None:
            self.group2.register_interface(
                interface_id=self.interface_id,
                interface=self,
                other_group=self.group1,
                side="group2",
                interface_name=self.name,
            )

    def make_interface_id(self, group1, group2=None):
        group1_id = group1.group_id

        if group2 is None:
            return f"{group1_id}__surrounding"

        group2_id = group2.group_id

        # Sorting makes A-B and B-A resolve to the same interface.
        first, second = sorted((group1_id, group2_id))
        return f"{first}__{second}"

    def get_surrounding_indices(self, group1_indices, surr_dist=5):
        """
        Return system balls close enough to group 1 to participate in a
        vertex whose radius does not exceed the network maximum.

        A candidate surrounding ball j is retained when at least one group-1
        ball i satisfies:

            distance(i, j) <= radius_i + radius_j + 2 * max_vert
        """

        group1_indices = set(int(i) for i in group1_indices)

        all_locations = np.asarray(
            self.sys.balls["loc"].tolist(),
            dtype=float,
        )
        all_radii = np.asarray(
            self.sys.balls["rad"].to_numpy(),
            dtype=float,
        )

        surrounding_indices = set()

        for group1_index in group1_indices:
            group1_location = all_locations[group1_index]
            group1_radius = all_radii[group1_index]

            distances = np.linalg.norm(
                all_locations - group1_location,
                axis=1,
            )

            cutoffs = (
                    group1_radius
                    + surr_dist
            )

            candidate_indices = np.flatnonzero(
                distances <= cutoffs
            )

            surrounding_indices.update(
                int(index)
                for index in candidate_indices
                if int(index) not in group1_indices
            )

        return surrounding_indices

    def _update_group_metadata(self, network_created=None, built=None):
        groups = [self.group1]

        if self.group2 is not None:
            groups.append(self.group2)

        for group in groups:
            metadata = group.interface_metadata[self.interface_id]

            if network_created is not None:
                metadata["network_created"] = network_created

            if built is not None:
                metadata["built"] = built

    def export(self, all_=False, atoms=False, surfs=False, edges=False, verts=False, logs=False, info=False,
               group_info=False, round_to=3, dual=False, buried_water=None):
        interface_exports(iface=self, all_=all_, atoms=atoms, surfs=surfs, edges=edges, verts=verts, logs=logs,
                          info=info, group_info=group_info, round_to=round_to, dual=dual,
                          buried_water=buried_water)
