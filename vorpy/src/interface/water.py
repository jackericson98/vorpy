"""
Interface-water topology analysis.

This module deliberately defines an *interface water* from the retained
Interface topology rather than from a radial solvent cutoff.

For a two-sided interface, every retained direct interface surface is defined
by one Group-1 atom and one Group-2 atom.  A retained interface edge is defined
by three atoms.  Therefore an edge whose defining atoms contain Group 1,
Group 2, and a solvent-water atom identifies a water that directly interrupts
or terminates the Group-1 <-> Group-2 surface patch.

Burial is then treated as a topology problem.  Water-touching interface edges
are connected through their retained vertices.  Waters whose edges are joined
through shared vertices form a linked water component.  The edge/vertex graph
for each component is inspected for closed cycles and open ends.

The first implementation is intentionally diagnostic-first.  It reports the
raw graph quantities and a provisional topology class rather than embedding a
chemistry-specific empirical cutoff.
"""

from collections import defaultdict, deque
import os

import numpy as np

from vorpy.src.group import Group


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def _safe_list(value):
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return []


def _is_missing(value):
    if value is None:
        return True
    try:
        return bool(np.isnan(value))
    except (TypeError, ValueError):
        return False


def _residue_label(residue):
    name = getattr(residue, "name", "SOL")
    seq = getattr(residue, "seq", None)
    return f"{name} {seq}" if seq is not None else str(name)


def _residue_key(residue):
    """Stable key for one solvent residue based on parent-system atom IDs."""
    atoms = []
    for value in _safe_list(getattr(residue, "atoms", [])):
        try:
            atoms.append(int(value))
        except (TypeError, ValueError):
            pass
    if atoms:
        return tuple(sorted(atoms))
    return ("object", id(residue))


def _residue_debug_label(residue):
    """Human-readable residue label with atom IDs for unambiguous diagnostics."""
    label = _residue_label(residue)
    atoms = []
    for value in _safe_list(getattr(residue, "atoms", [])):
        try:
            atoms.append(int(value))
        except (TypeError, ValueError):
            pass
    if not atoms:
        return label
    atoms = sorted(atoms)
    if len(atoms) == 1:
        atom_text = str(atoms[0])
    elif atoms == list(range(atoms[0], atoms[-1] + 1)):
        atom_text = f"{atoms[0]}-{atoms[-1]}"
    else:
        atom_text = ",".join(str(v) for v in atoms)
    return f"{label} [atoms {atom_text}]"


# ---------------------------------------------------------------------------
# Interface topology/system mapping
# ---------------------------------------------------------------------------


def _build_interface_index_context(iface):
    """
    Resolve network TOPOLOGY identifiers back to parent SYSTEM atom indices.

    Interface edge/vertex ``balls`` use network topology IDs.  Group membership
    and solvent-residue membership use parent-system IDs.  ``system_num`` is
    authoritative when present.  The fallback topology==system convention is
    retained only for older/full-system networks.
    """
    net = iface.net
    topology_to_system = {}
    system_to_topology = {}

    if net is not None and net.balls is not None:
        columns = list(net.balls.columns)
        num_idx = columns.index("num") if "num" in columns else None
        system_idx = columns.index("system_num") if "system_num" in columns else None

        for row_pos, row in enumerate(net.balls.itertuples(index=False, name=None)):
            try:
                topology_id = int(row[num_idx]) if num_idx is not None else int(row_pos)
            except (TypeError, ValueError):
                topology_id = int(row_pos)

            system_id = topology_id
            if system_idx is not None:
                value = row[system_idx]
                if not _is_missing(value):
                    try:
                        system_id = int(value)
                    except (TypeError, ValueError):
                        pass

            topology_to_system[topology_id] = system_id
            system_to_topology.setdefault(system_id, topology_id)

    group1_system = set(int(v) for v in getattr(iface, "group1_indices", set()))
    group2_system = set(int(v) for v in getattr(iface, "group2_indices", set()))

    # Build water atom -> residue lookup once.
    water_by_system = {}
    solvent_residues = []
    sol = getattr(iface.sys, "sol", None)
    if sol is not None:
        solvent_residues = getattr(sol, "residues", None) or []

    for residue in solvent_residues:
        for system_id in _safe_list(getattr(residue, "atoms", [])):
            try:
                water_by_system[int(system_id)] = residue
            except (TypeError, ValueError):
                pass

    system_ball_count = len(iface.sys.balls) if getattr(iface.sys, "balls", None) is not None else 0

    def topology_to_parent(topology_id):
        """Map one topology ID; diagnostic fallback supports older networks."""
        topology_id = int(topology_id)
        mapped = topology_to_system.get(topology_id)
        if mapped is not None:
            return mapped, "mapped"

        # Older VorPy interface networks often retained parent-system numbering
        # for outside defining atoms even when those atoms were not represented
        # as rows in net.balls.  Keep this as an explicit, countable fallback.
        if 0 <= topology_id < system_ball_count:
            return topology_id, "identity_fallback"

        return None, "unmapped"

    return {
        "topology_to_system": topology_to_system,
        "system_to_topology": system_to_topology,
        "topology_to_parent": topology_to_parent,
        "group1_system": group1_system,
        "group2_system": group2_system,
        "water_by_system": water_by_system,
        "solvent_residues": solvent_residues,
    }


# ---------------------------------------------------------------------------
# Interfacial-water discovery from interface edges
# ---------------------------------------------------------------------------


def _classify_interface_edge(iface, edge_index, edge, context):
    """
    Classify one retained edge by its defining atoms.

    The strict interfacial-water signature is:
        >=1 defining atom in Group 1
        >=1 defining atom in Group 2
        >=1 defining atom belonging to a solvent-water residue

    For a regular Voronoi edge there are three defining balls, so the expected
    common case is exactly one atom from each category.
    """
    record = {
        "edge_index": int(edge_index),
        "topology_balls": [],
        "system_balls": [],
        "mapping_modes": [],
        "group1_atoms": [],
        "group2_atoms": [],
        "water_residues": [],
        "water_system_atoms": [],
        "other_system_atoms": [],
        "vertex_indices": [],
        "is_interface_edge": False,
        "is_interfacial_water_edge": False,
        "mapping_failed": False,
    }

    try:
        topology_balls = [int(v) for v in edge["balls"]]
    except (KeyError, TypeError, ValueError):
        record["mapping_failed"] = True
        return record

    record["topology_balls"] = topology_balls

    residue_seen = set()
    for topology_id in topology_balls:
        system_id, mode = context["topology_to_parent"](topology_id)
        record["mapping_modes"].append(mode)
        record["system_balls"].append(system_id)

        if system_id is None:
            record["mapping_failed"] = True
            continue

        if system_id in context["group1_system"]:
            record["group1_atoms"].append(system_id)
        if system_id in context["group2_system"]:
            record["group2_atoms"].append(system_id)

        residue = context["water_by_system"].get(system_id)
        if residue is not None:
            key = _residue_key(residue)
            if key not in residue_seen:
                record["water_residues"].append(residue)
                residue_seen.add(key)
            record["water_system_atoms"].append(system_id)
        elif (
            system_id not in context["group1_system"]
            and system_id not in context["group2_system"]
        ):
            record["other_system_atoms"].append(system_id)

    record["vertex_indices"] = sorted(
        int(v) for v in _safe_list(edge.get("verts", []))
    )
    record["is_interface_edge"] = bool(record["group1_atoms"] and record["group2_atoms"])
    record["is_interfacial_water_edge"] = bool(
        record["is_interface_edge"] and record["water_residues"]
    )
    return record


def discover_interfacial_water_edges(iface):
    """Return strict Group1-Group2-water edge records and discovery QC."""
    context = _build_interface_index_context(iface)
    net = iface.net

    records = []
    qc = {
        "total_edges": 0,
        "group1_group2_edges": 0,
        "interfacial_water_edges": 0,
        "unique_interfacial_waters": 0,
        "identity_fallback_uses": 0,
        "unmapped_ball_references": 0,
        "edges_with_mapping_failure": 0,
        "edges_with_multiple_waters": 0,
    }

    if net is None or net.edges is None:
        return records, qc, context

    qc["total_edges"] = len(net.edges)
    water_keys = set()

    for edge_index, edge in net.edges.iterrows():
        record = _classify_interface_edge(iface, edge_index, edge, context)

        qc["identity_fallback_uses"] += sum(
            1 for mode in record["mapping_modes"] if mode == "identity_fallback"
        )
        qc["unmapped_ball_references"] += sum(
            1 for mode in record["mapping_modes"] if mode == "unmapped"
        )
        if record["mapping_failed"]:
            qc["edges_with_mapping_failure"] += 1

        if record["is_interface_edge"]:
            qc["group1_group2_edges"] += 1

        if not record["is_interfacial_water_edge"]:
            continue

        records.append(record)
        qc["interfacial_water_edges"] += 1
        if len(record["water_residues"]) > 1:
            qc["edges_with_multiple_waters"] += 1
        for residue in record["water_residues"]:
            water_keys.add(_residue_key(residue))

    qc["unique_interfacial_waters"] = len(water_keys)
    return records, qc, context


# ---------------------------------------------------------------------------
# Water linkage and edge<->vertex graph analysis
# ---------------------------------------------------------------------------


def _build_water_records(edge_records):
    waters = {}

    for edge_record in edge_records:
        for residue in edge_record["water_residues"]:
            key = _residue_key(residue)
            entry = waters.setdefault(
                key,
                {
                    "key": key,
                    "residue": residue,
                    "label": _residue_label(residue),
                    "debug_label": _residue_debug_label(residue),
                    "edge_indices": set(),
                    "vertex_indices": set(),
                    "linked_water_keys": set(),
                    "component_id": None,
                    "topology_class": None,
                },
            )
            entry["edge_indices"].add(edge_record["edge_index"])
            entry["vertex_indices"].update(edge_record["vertex_indices"])

    return waters


def _link_waters_through_vertices(waters):
    """
    Link interfacial waters only through shared retained interface vertices.

    This prevents traversal from escaping into the bulk solvent graph: the
    universe is already restricted to waters identified from retained
    Group1-Group2-water edges.
    """
    vertex_to_water_keys = defaultdict(set)

    for key, water in waters.items():
        for vertex_index in water["vertex_indices"]:
            vertex_to_water_keys[int(vertex_index)].add(key)

    for keys in vertex_to_water_keys.values():
        if len(keys) < 2:
            continue
        keys = set(keys)
        for key in keys:
            waters[key]["linked_water_keys"].update(keys - {key})

    return vertex_to_water_keys


def _connected_water_components(waters):
    components = []
    unseen = set(waters)

    while unseen:
        start = next(iter(unseen))
        queue = deque([start])
        unseen.remove(start)
        component = set()

        while queue:
            key = queue.popleft()
            component.add(key)
            for nbr in waters[key]["linked_water_keys"]:
                if nbr in unseen:
                    unseen.remove(nbr)
                    queue.append(nbr)

        components.append(component)

    return components


def _edge_vertex_graph_for_component(iface, component_keys, waters):
    """
    Build the selected interface-edge graph for one linked water component.

    Each selected VorPy edge contributes a connection between all of its
    retained endpoint vertices.  The common case is two endpoint vertices.
    Degenerate/multi-vertex edges are retained diagnostically.
    """
    net = iface.net
    edge_indices = set()
    for key in component_keys:
        edge_indices.update(waters[key]["edge_indices"])

    vertex_to_edges = defaultdict(set)
    edge_to_vertices = {}

    for edge_index in sorted(edge_indices):
        try:
            edge = net.edges.loc[edge_index]
        except (KeyError, IndexError):
            continue

        verts = sorted(set(int(v) for v in _safe_list(edge.get("verts", []))))
        edge_to_vertices[int(edge_index)] = verts
        for vertex_index in verts:
            vertex_to_edges[vertex_index].add(int(edge_index))

    # Edge adjacency follows the user's edge -> vertex -> edge traversal.
    edge_adjacency = {edge_index: set() for edge_index in edge_to_vertices}
    for incident_edges in vertex_to_edges.values():
        for edge_index in incident_edges:
            edge_adjacency[edge_index].update(incident_edges - {edge_index})

    return {
        "edge_indices": sorted(edge_indices),
        "edge_to_vertices": edge_to_vertices,
        "vertex_to_edges": vertex_to_edges,
        "edge_adjacency": edge_adjacency,
    }


def _build_bipartite_adjacency(graph):
    """Build the literal edge-node <-> vertex-node graph requested for loop QC."""
    adjacency = defaultdict(set)
    for edge_index, vertices in graph["edge_to_vertices"].items():
        edge_node = ("edge", int(edge_index))
        adjacency.setdefault(edge_node, set())
        for vertex_index in vertices:
            vertex_node = ("vertex", int(vertex_index))
            adjacency[edge_node].add(vertex_node)
            adjacency[vertex_node].add(edge_node)
    return dict(adjacency)


def _generic_graph_components(adjacency):
    unseen = set(adjacency)
    components = []
    while unseen:
        start = next(iter(unseen))
        stack = [start]
        unseen.remove(start)
        found = set()
        while stack:
            node = stack.pop()
            found.add(node)
            for nbr in adjacency.get(node, ()):
                if nbr in unseen:
                    unseen.remove(nbr)
                    stack.append(nbr)
        components.append(found)
    return components


def _find_one_bipartite_cycle(adjacency):
    """Return one exact edge->vertex->edge cycle in the bipartite graph."""
    visited = set()

    def dfs(node, parent, path, position):
        visited.add(node)
        position[node] = len(path)
        path.append(node)

        for nbr in adjacency.get(node, ()):
            if nbr == parent:
                continue
            if nbr in position:
                start = position[nbr]
                return path[start:] + [nbr]
            if nbr not in visited:
                cycle = dfs(nbr, node, path, position)
                if cycle is not None:
                    return cycle

        path.pop()
        position.pop(node, None)
        return None

    for start in adjacency:
        if start in visited:
            continue
        cycle = dfs(start, None, [], {})
        if cycle is not None:
            return cycle
    return None


def _format_bipartite_cycle(cycle):
    if not cycle:
        return None
    return " -> ".join(
        (f"E{value}" if kind == "edge" else f"V{value}")
        for kind, value in cycle
    )




def _find_bridges(adjacency):
    """Return undirected graph links that are bridges using Tarjan DFS."""
    timer = [0]
    discovery = {}
    low = {}
    bridges = set()

    def canon(a, b):
        return frozenset((a, b))

    def dfs(node, parent=None):
        discovery[node] = timer[0]
        low[node] = timer[0]
        timer[0] += 1

        for nbr in adjacency.get(node, ()):
            if nbr == parent:
                continue
            if nbr not in discovery:
                dfs(nbr, node)
                low[node] = min(low[node], low[nbr])
                if low[nbr] > discovery[node]:
                    bridges.add(canon(node, nbr))
            else:
                low[node] = min(low[node], discovery[nbr])

    for node in adjacency:
        if node not in discovery:
            dfs(node)

    return bridges


def _cycle_core_edges(graph, bipartite):
    """Return VorPy edge IDs that participate in at least one exact graph cycle."""
    bridges = _find_bridges(bipartite) if bipartite else set()
    cycle_edges = set()

    for edge_index, vertices in graph["edge_to_vertices"].items():
        edge_node = ("edge", int(edge_index))
        incidence_links = [
            frozenset((edge_node, ("vertex", int(vertex_index))))
            for vertex_index in vertices
        ]
        # A regular VorPy edge lies on a graph cycle iff both endpoint
        # incidences are non-bridges.  Degenerate edges remain diagnostic only.
        if len(incidence_links) == 2 and all(link not in bridges for link in incidence_links):
            cycle_edges.add(int(edge_index))

    return cycle_edges, bridges


def _cycle_edge_components(graph, cycle_edge_indices):
    """Connected components formed only by cycle-participating VorPy edges."""
    cycle_edge_indices = set(int(v) for v in cycle_edge_indices)
    if not cycle_edge_indices:
        return []

    vertex_to_cycle_edges = defaultdict(set)
    for edge_index in cycle_edge_indices:
        for vertex_index in graph["edge_to_vertices"].get(edge_index, []):
            vertex_to_cycle_edges[int(vertex_index)].add(edge_index)

    adjacency = {edge_index: set() for edge_index in cycle_edge_indices}
    for incident in vertex_to_cycle_edges.values():
        for edge_index in incident:
            adjacency[edge_index].update(incident - {edge_index})

    unseen = set(adjacency)
    components = []
    while unseen:
        start = next(iter(unseen))
        stack = [start]
        unseen.remove(start)
        found = set()
        while stack:
            edge_index = stack.pop()
            found.add(edge_index)
            for nbr in adjacency[edge_index]:
                if nbr in unseen:
                    unseen.remove(nbr)
                    stack.append(nbr)
        components.append(found)
    return components

def _analyze_component_graph(iface, component_keys, waters):
    graph = _edge_vertex_graph_for_component(iface, component_keys, waters)
    vertex_to_edges = graph["vertex_to_edges"]

    vertex_degrees = {
        int(vertex): len(edges)
        for vertex, edges in vertex_to_edges.items()
    }
    open_vertices = sorted(v for v, degree in vertex_degrees.items() if degree == 1)
    branch_vertices = sorted(v for v, degree in vertex_degrees.items() if degree > 2)
    isolated_vertices = sorted(v for v, degree in vertex_degrees.items() if degree == 0)

    nonbinary_edges = sorted(
        edge_index
        for edge_index, verts in graph["edge_to_vertices"].items()
        if len(verts) != 2
    )

    # The burial test uses the literal bipartite topology:
    #   edge -> vertex -> edge -> vertex -> ...
    # This avoids false cycles caused by projecting several edges meeting at
    # one branch vertex into an edge-adjacency graph.
    bipartite = _build_bipartite_adjacency(graph)
    bipartite_components = _generic_graph_components(bipartite) if bipartite else []
    incidence_count = sum(len(nbrs) for node, nbrs in bipartite.items() if node[0] == "edge")
    node_count = len(bipartite)
    graph_component_count = len(bipartite_components)
    cycle_rank = incidence_count - node_count + graph_component_count if bipartite else 0

    explicit_cycle_nodes = _find_one_bipartite_cycle(bipartite)
    explicit_cycle = _format_bipartite_cycle(explicit_cycle_nodes)
    has_cycle = cycle_rank > 0
    cycle_edge_indices, bridge_links = _cycle_core_edges(graph, bipartite)
    cycle_components = _cycle_edge_components(graph, cycle_edge_indices)

    # Component class is retained only as a coarse diagnostic.  Final per-water
    # burial uses local membership in the cycle core below.
    if has_cycle and len(open_vertices) == 0:
        topology_class = "closed_cycle_component"
    elif has_cycle:
        topology_class = "mixed_cycle_open_component"
    else:
        topology_class = "open_component"

    graph.update(
        {
            "vertex_degrees": vertex_degrees,
            "open_vertices": open_vertices,
            "branch_vertices": branch_vertices,
            "isolated_vertices": isolated_vertices,
            "nonbinary_edges": nonbinary_edges,
            "edge_graph_component_count": graph_component_count,
            "cycle_rank": cycle_rank,
            "has_cycle": bool(has_cycle),
            "cycle_edge_indices": sorted(cycle_edge_indices),
            "cycle_edge_count": len(cycle_edge_indices),
            "bridge_incidence_count": len(bridge_links),
            "cycle_edge_components": [sorted(v) for v in cycle_components],
            "explicit_cycle": explicit_cycle,
            "explicit_cycle_nodes": explicit_cycle_nodes,
            "bipartite_adjacency": bipartite,
            "topology_class": topology_class,
        }
    )
    return graph


def _classify_cycle_set_group(iface, cycle_set, waters):
    """Classify one connected cycle-core set at the GROUP level."""
    member_keys = set(cycle_set.get("water_keys", set()))
    group_graph = _analyze_component_graph(iface, member_keys, waters)
    member_edges = set(int(v) for v in group_graph.get("edge_indices", []))
    cycle_edges = set(int(v) for v in cycle_set.get("edge_indices", []))
    extra_edges = sorted(member_edges - cycle_edges)
    open_vertices = list(group_graph.get("open_vertices", []))

    is_buried_group = bool(cycle_edges) and not extra_edges and not open_vertices
    return {
        "burial_class": "buried" if is_buried_group else "semi_buried",
        "group_graph": group_graph,
        "extra_noncycle_edge_indices": extra_edges,
        "open_vertex_indices": open_vertices,
    }


def analyze_interface_water_topology(iface):
    """
    Discover and classify linked sets of interface waters from edge topology.

    Returns a dictionary suitable for info-file diagnostics and later machine
    analysis.  No full per-water Group builds are needed for this pass.
    """
    edge_records, discovery_qc, context = discover_interfacial_water_edges(iface)
    waters = _build_water_records(edge_records)
    vertex_to_water_keys = _link_waters_through_vertices(waters)
    component_keys = _connected_water_components(waters)

    components = []
    for component_id, keys in enumerate(component_keys, start=1):
        graph = _analyze_component_graph(iface, keys, waters)
        labels = sorted(waters[key]["label"] for key in keys)

        component = {
            "component_id": component_id,
            "water_keys": set(keys),
            "water_labels": labels,
            "water_count": len(keys),
            **graph,
        }
        components.append(component)

        for key in keys:
            waters[key]["component_id"] = component_id
            waters[key]["topology_class"] = graph["topology_class"]

    # Build exact cycle sets from each coarse linked component.  Peripheral
    # chains can remain in the same coarse component without inheriting burial.
    cycle_sets = []
    cycle_set_id = 1
    for component in components:
        for cycle_edges in component.get("cycle_edge_components", []):
            cycle_edges = set(int(v) for v in cycle_edges)
            member_keys = set()
            for key in component["water_keys"]:
                if set(waters[key]["edge_indices"]) & cycle_edges:
                    member_keys.add(key)
            if not member_keys:
                continue
            cycle_sets.append({
                "cycle_set_id": cycle_set_id,
                "component_id": component["component_id"],
                "edge_indices": sorted(cycle_edges),
                "water_keys": member_keys,
                "water_labels": sorted(waters[key]["label"] for key in member_keys),
                "water_debug_labels": sorted(waters[key]["debug_label"] for key in member_keys),
                "water_count": len(member_keys),
            })
            cycle_set_id += 1

    # Classify each cycle set as a GROUP.  A true buried group must be closed
    # as a whole: every retained interface edge belonging to its member waters
    # must lie in this cycle-core set, and the combined member-water graph must
    # have no open vertices.  This preserves multi-water closed holes while
    # excluding semi-buried waters that contain a local loop plus open branches.
    buried_cycle_sets = []
    for cycle_set in cycle_sets:
        group_classification = _classify_cycle_set_group(iface, cycle_set, waters)
        cycle_set.update(group_classification)
        if cycle_set["burial_class"] == "buried":
            buried_cycle_sets.append(cycle_set)

    # Per-water burial is LOCAL: only the water's own edges that lie in a
    # cycle-core count.  This prevents an open peripheral water from inheriting
    # the classification of a distant loop in the same connected component.
    cycle_edge_union = set()
    edge_to_cycle_sets = defaultdict(set)
    for cycle_set in cycle_sets:
        for edge_index in cycle_set["edge_indices"]:
            cycle_edge_union.add(int(edge_index))
            edge_to_cycle_sets[int(edge_index)].add(cycle_set["cycle_set_id"])

    for key, water in waters.items():
        individual_graph = _analyze_component_graph(iface, {key}, waters)
        water["individual_graph"] = individual_graph
        water["linked_water_labels"] = sorted(
            waters[nbr]["debug_label"] for nbr in water["linked_water_keys"]
        )

        edge_indices = set(int(v) for v in water["edge_indices"])
        cycle_edges = edge_indices & cycle_edge_union
        water["cycle_edge_indices"] = sorted(cycle_edges)
        water["cycle_edge_count"] = len(cycle_edges)
        water["cycle_edge_fraction"] = (
            len(cycle_edges) / len(edge_indices) if edge_indices else 0.0
        )
        water["cycle_set_ids"] = sorted({
            cycle_id
            for edge_index in cycle_edges
            for cycle_id in edge_to_cycle_sets.get(edge_index, ())
        })

        open_count = len(individual_graph.get("open_vertices", []))
        if not cycle_edges:
            burial_class = "peripheral"
        elif len(cycle_edges) == len(edge_indices) and open_count == 0:
            burial_class = "buried"
        else:
            burial_class = "semi_buried"
        water["burial_class"] = burial_class

    return {
        "edge_records": edge_records,
        "discovery_qc": discovery_qc,
        "waters": waters,
        "components": components,
        "cycle_sets": cycle_sets,
        "buried_cycle_sets": buried_cycle_sets,
        "vertex_to_water_keys": vertex_to_water_keys,
        "context": context,
    }


# ---------------------------------------------------------------------------
# Compatibility wrappers used by Interface.build()
# ---------------------------------------------------------------------------


def get_interface_water_residues(iface):
    """
    Return waters discovered from strict Group1-Group2-water interface edges.

    This replaces the old surface-endpoint test, which necessarily returned
    zero for a network containing only direct Group1-Group2 surfaces.
    """
    topology = analyze_interface_water_topology(iface)
    residues = [entry["residue"] for entry in topology["waters"].values()]
    return sorted(residues, key=lambda r: (str(getattr(r, "name", "")), getattr(r, "seq", -1)))


def get_water_interface_geometry(iface, residue, topology_analysis=None):
    """Return retained edge/vertex topology associated with one interface water."""
    if topology_analysis is None:
        topology_analysis = analyze_interface_water_topology(iface)

    key = _residue_key(residue)
    water = topology_analysis["waters"].get(key)
    if water is None:
        return None

    edge_indices = sorted(water["edge_indices"])
    vertex_indices = sorted(water["vertex_indices"])
    net = iface.net

    edge_table = net.edges.loc[edge_indices].copy() if edge_indices else net.edges.iloc[0:0].copy()
    vertex_table = net.verts.loc[vertex_indices].copy() if vertex_indices else net.verts.iloc[0:0].copy()

    return {
        "residue": residue,
        "water_balls": sorted(int(v) for v in _safe_list(getattr(residue, "atoms", []))),
        "surface_indices": [],
        "edge_indices": edge_indices,
        "vertex_indices": vertex_indices,
        "surfaces": net.surfs.iloc[0:0].copy() if net.surfs is not None else None,
        "edges": edge_table,
        "vertices": vertex_table,
        "surface_area": 0.0,
        "contact_area": 0.0,
        "retained_geometry": None,
        "component_id": water["component_id"],
        "topology_class": water["topology_class"],
        "linked_water_labels": list(water.get("linked_water_labels", [])),
        "individual_graph": water.get("individual_graph"),
        "topology_analysis": topology_analysis,
    }


def analyze_interface_waters(iface):
    """
    Analyze all interface waters once and retain the global topology result.
    """
    topology = analyze_interface_water_topology(iface)
    iface.water_topology = topology

    geometries = []
    for water in topology["waters"].values():
        geometry = get_water_interface_geometry(
            iface=iface,
            residue=water["residue"],
            topology_analysis=topology,
        )
        if geometry is not None:
            geometries.append(geometry)

    geometries.sort(
        key=lambda g: (
            str(getattr(g["residue"], "name", "")),
            getattr(g["residue"], "seq", -1),
        )
    )
    return geometries


# ---------------------------------------------------------------------------
# Optional complete-water Group calculations retained for later geometry work
# ---------------------------------------------------------------------------


def get_water_group_volume(water_group):
    """Return the complete volume assigned to the water Group."""
    value = getattr(water_group, "vol", None)
    if value is not None:
        try:
            return float(value)
        except (TypeError, ValueError):
            pass

    system_balls = water_group.sys.balls
    water_indices = list(water_group.ball_ndxs)
    if "vol" in system_balls.columns:
        try:
            return float(system_balls.loc[water_indices, "vol"].dropna().sum())
        except (KeyError, TypeError, IndexError):
            return float(system_balls.iloc[water_indices]["vol"].dropna().sum())
    return None


def summarize_water_group(water_group):
    """Summarize the complete analyzed boundary geometry of one water Group."""
    net = water_group.net
    if net is None:
        return None

    area = float(getattr(water_group, "sa", 0.0) or 0.0)
    oriented_c = float(getattr(water_group, "oriented_int_mean_curv", 0.0) or 0.0)
    q = float(getattr(water_group, "int_mean_curv_sq", 0.0) or 0.0)
    x = float(getattr(water_group, "int_gauss_curv", 0.0) or 0.0)

    return {
        "group": water_group,
        "group_name": water_group.name,
        "ball_indices": list(water_group.ball_ndxs),
        "vertex_count": len(net.verts) if net.verts is not None else 0,
        "edge_count": len(net.edges) if net.edges is not None else 0,
        "surface_count": len(net.surfs) if net.surfs is not None else 0,
        "surface_area": area,
        "volume": get_water_group_volume(water_group),
        "oriented_int_mean_curv": oriented_c,
        "int_mean_curv_sq": q,
        "int_gauss_curv": x,
        "surf_energy": 2.0 * q,
        "convex_sa": float(getattr(water_group, "convex_sa", 0.0) or 0.0),
        "concave_sa": float(getattr(water_group, "concave_sa", 0.0) or 0.0),
        "flat_sa": float(getattr(water_group, "flat_sa", 0.0) or 0.0),
        "failed_curved_sa": float(getattr(water_group, "curved_unoriented_sa", 0.0) or 0.0),
        "vertices": net.verts,
        "edges": net.edges,
        "surfaces": net.surfs,
    }


def _buried_group_name(group_index, residues):
    """Return a short deterministic name for one buried-water network."""
    return f"buried_{group_index}"


def _buried_group_color(group_index):
    """Return a distinct matplotlib color map for a buried-water shell."""
    color_maps = (
        "Reds",
        "Blues",
        "Greens",
        "Purples",
        "Oranges",
        "Greys",
        "YlGn",
        "PuRd",
    )
    return color_maps[(group_index - 1) % len(color_maps)]


def _cycle_set_residues(cycle_set, topology):
    """Resolve one cycle set's stable water keys back to residue objects."""
    waters = topology.get("waters", {})
    residues = []
    for key in cycle_set.get("water_keys", set()):
        water = waters.get(key)
        if water is not None and water.get("residue") is not None:
            residues.append(water["residue"])
    return sorted(
        residues,
        key=lambda r: (
            str(getattr(r, "name", "")),
            getattr(r, "seq", -1),
            _residue_key(r),
        ),
    )


def build_buried_water_groups(iface):
    """
    Build and analyze a normal VorPy Group for every CLOSED buried water set.

    Interface topology determines membership only.  Each selected residue set is
    then solved through the standard Group/Network path against the complete
    parent system, so V/A/C/Q/X and exported cell geometry have the same meaning
    as for any other complete Group.
    """
    topology = getattr(iface, "water_topology", None)
    if not topology:
        return []

    buried_sets = list(topology.get("buried_cycle_sets", []))
    waters_root = os.path.join(iface.dir, "waters")
    os.makedirs(waters_root, exist_ok=True)

    water_groups = []
    total = len(buried_sets)

    for group_index, cycle_set in enumerate(buried_sets, start=1):
        residues = _cycle_set_residues(cycle_set, topology)
        if not residues:
            continue

        group_name = _buried_group_name(group_index, residues)
        group_dir = os.path.join(waters_root, group_name)
        os.makedirs(group_dir, exist_ok=True)

        labels = [_residue_label(residue) for residue in residues]

        water_group = Group(
            sys=iface.sys,
            name=group_name,
            residues=residues,
            settings=iface.group1.settings.copy(),
            make_net=True,
            build_net=False,
            output_directory=group_dir,
            group_id=f"{iface.interface_id}__water_{group_index:03d}",
        )
        # Buried-water networks are subordinate interface work.  Keep their
        # detailed Network progress on the system-wide progress line rather than
        # presenting them as an unrelated top-level Group solve.
        if getattr(water_group, "net", None) is not None:
            water_group.net.progress_network_name = iface.name
            water_group.net.progress_process_prefix = (
                f"Interface water {group_index}/{total}: {group_name}"
            )
            water_group.net.completion_kind = "interface water network"


        water_group.build()
        water_group.get_info()

        # Give every buried-water shell a visibly different color map while
        # retaining the normal curvature-based coloring within that shell.
        shell_color_map = _buried_group_color(group_index)
        water_group.settings["surf_col"] = shell_color_map
        if getattr(water_group, "net", None) is not None:
            water_group.net.settings["surf_col"] = shell_color_map

        summary = summarize_water_group(water_group)
        cycle_set["solved_group_index"] = group_index
        cycle_set["solved_group_name"] = group_name
        cycle_set["solved_group_directory"] = group_dir
        cycle_set["full_group_geometry"] = summary

        # Lightweight metadata on the Group itself makes export/reporting simple
        # without coupling Group core code to Interface topology.
        water_group.parent_interface = iface
        water_group.buried_water_metadata = {
            "group_index": group_index,
            "cycle_set_id": cycle_set.get("cycle_set_id"),
            "residues": residues,
            "water_labels": labels,
            "cycle_edge_indices": list(cycle_set.get("edge_indices", [])),
            "burial_class": "buried",
            "shell_color_map": shell_color_map,
        }
        water_groups.append(water_group)

    return water_groups


def build_interface_water_groups(iface, water_geometries=None):
    """Backward-compatible wrapper: production now builds buried sets only."""
    return build_buried_water_groups(iface)
