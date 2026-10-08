import os
from time import perf_counter
import numpy as np
from copy import deepcopy
from vorpy.src.group import Group
from vorpy.src.network import Network
from vorpy.src.boundary import network_geometry
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

        interface_indices = sorted(self.group1_indices | self.group2_indices)

        locs, rads, masses, boundary_indices, boundary_config = network_geometry(self.sys)
        self.net = Network(
            # Retain the complete system geometry so surrounding balls can
            # participate in geometric validity checks.
            locs=locs,
            rads=rads,
            masses=masses,

            # Balls whose interface topology is being represented.
            group=interface_indices,

            # Original side membership used during interface vertex filtering.
            iface_grps=(self.group1_indices, self.group2_indices),

            group_name=self.name,
            settings=self.settings,
            sort_balls=False,
            verts=verts,
            system=self.sys,
        )
        self.net.boundary_indices = boundary_indices
        self.net.boundary_config = boundary_config
        self.net.boundary_mode = getattr(self.sys, "boundary_mode", None)
        if boundary_indices:
            self.net.balls["is_boundary_generator"] = False
            self.net.balls.loc[list(boundary_indices), "is_boundary_generator"] = True

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
        # Keep water classification in info.txt; normal CLI output uses the
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
