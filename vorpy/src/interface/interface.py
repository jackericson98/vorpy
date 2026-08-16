import os
import numpy as np
from copy import deepcopy
from vorpy.src.group import Group
from vorpy.src.network import Network
from vorpy.src.interface.export import interface_exports
from vorpy.src.interface.water import analyze_interface_waters
from vorpy.src.interface.water import build_interface_water_groups


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

        # Cached index collections defining the two interface sides.
        self.group1_indices = None
        self.group2_indices = None

        # Partial groups created specifically for this interface
        self.partial_group1 = None
        self.partial_group2 = None

        self._register_with_groups()

        self.water_groups = []
        self.water_geometries = []

    def _make_name(self):
        second_name = "surrounding" if self.group2 is None else self.group2.name
        self.name = f"{self.group1.name}_{second_name}_interface"
        return self.name

    def set_dir(self):

        i = 1
        my_dir = self.sys.files['dir'] + "/" + self.name
        first = True
        while os.path.exists(my_dir):
            if first:
                my_dir += "__"
                first = False
            my_dir = my_dir[:-(1 + len(str(i)))] + '_' + str(i)
            i += 1
        self.dir = my_dir
        os.mkdir(self.dir)

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
        if self.dir is None:
            self.set_dir()

        self.group1_indices = set(self.group1.ball_ndxs)

        if self.group2 is None:
            self.group2_indices = self.get_surrounding_indices(self.group1_indices)
            self.group2_name = "surrounding"
        else:
            self.group2_indices = set(self.group2.ball_ndxs)
            self.group2_name = self.group2.name

        interface_indices = sorted(self.group1_indices | self.group2_indices)

        self.net = Network(
            # Retain the complete system geometry so surrounding balls can
            # participate in geometric validity checks.
            locs=self.sys.balls["loc"],
            rads=self.sys.balls["rad"],
            masses=self.sys.balls["mass"],

            # Balls whose interface topology is being represented.
            group=interface_indices,

            # Original side membership used during interface vertex filtering.
            iface_grps=(self.group1_indices, self.group2_indices),

            group_name=self.name,
            settings=self.settings,
            sort_balls=True,
            verts=verts,
        )

        self._update_group_metadata(
            network_created=True,
            built=False,
        )

        return self.net

    def build(self):
        """
        Create and build this interface's dedicated network.
        """
        if self.net is None:
            self.make_net()

        self.net.build()

        # Make completed topology available to later pairwise interface builds.
        self.sys.cache_interface_geometry(self)

        # Identify all waters touching the completed interface.
        self.water_geometries = analyze_interface_waters(
            iface=self,
        )

        print("\n[INTERFACE WATER DETECTION]")
        print(f"touching waters: {len(self.water_geometries)}")

        for water_geometry in self.water_geometries[:10]:
            residue = water_geometry["residue"]
            interface_geometry = water_geometry["interface_geometry"]

            print(
                f"  {residue.name} {residue.seq}: "
                f"group1_surfs="
                f"{interface_geometry['group1']['surface_count']}, "
                f"group2_surfs="
                f"{interface_geometry['group2']['surface_count']}, "
                f"bridging="
                f"{interface_geometry['bridging_water']}"
            )

        # ----------------------------------------------------------
        # TEMPORARY TEST: fully build only the first touching water.
        # ----------------------------------------------------------
        test_water_geometries = self.water_geometries[:1]

        self.water_groups = build_interface_water_groups(
            iface=self,
            water_geometries=test_water_geometries,
        )

        if self.water_groups:
            water_group = self.water_groups[0]

            print("\n[WATER GROUP TEST]")
            print(f"name: {water_group.name}")
            print(f"balls: {water_group.ball_ndxs}")
            print(f"net: {water_group.net}")

            print(
                "verts:",
                len(water_group.net.verts)
                if water_group.net is not None
                   and water_group.net.verts is not None
                else 0,
            )

            print(
                "edges:",
                len(water_group.net.edges)
                if water_group.net is not None
                   and water_group.net.edges is not None
                else 0,
            )

            print(
                "surfs:",
                len(water_group.net.surfs)
                if water_group.net is not None
                   and water_group.net.surfs is not None
                else 0,
            )

            print(
                "water group added to system:",
                water_group in (self.sys.groups or []),
            )

        self._update_group_metadata(
            network_created=True,
            built=True,
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
               group_info=False, round_to=3):
        interface_exports(iface=self, all_=all_, atoms=atoms, surfs=surfs, edges=edges, verts=verts, logs=logs,
                          info=info, group_info=group_info, round_to=round_to)