import time


def _is_missing(value):
    """Return True for None/NaN-like scalar values without importing NumPy."""
    if value is None:
        return True
    try:
        return value != value
    except Exception:
        return False


def _build_index_maps(grp):
    """
    Build explicit mappings between the two atom-index spaces used by VorPy.

    Group selections and PDB exports use SYSTEM indices (grp.sys.balls rows).
    Network topology (surface['balls'], edge['balls'], vertex['balls']) uses
    TOPOLOGY identifiers stored in net.balls['num'].

    net.balls['system_num'], when present, is the authoritative mapping back
    to the parent System.  Older/full-system networks often have identical
    values, but this function does not assume that.
    """
    net = grp.net

    topology_to_row = {}
    topology_to_system = {}
    system_to_topology = {}

    for row_pos, (_, atom) in enumerate(net.balls.iterrows()):
        try:
            topology_id = int(atom.get("num", row_pos))
        except (TypeError, ValueError):
            topology_id = int(row_pos)

        system_id = None
        if "system_num" in net.balls.columns:
            value = atom.get("system_num", None)
            if not _is_missing(value):
                try:
                    system_id = int(value)
                except (TypeError, ValueError):
                    system_id = None

        if system_id is None:
            # Backward-compatible full-system convention.
            system_id = topology_id

        topology_to_row[topology_id] = row_pos
        topology_to_system[topology_id] = system_id

        # Prefer the first topology record if malformed data duplicates a
        # system index.  Diagnostics below will make this visible.
        system_to_topology.setdefault(system_id, topology_id)

    return topology_to_row, topology_to_system, system_to_topology


def _safe_list(value):
    if value is None:
        return []
    try:
        return list(value)
    except TypeError:
        return []


def _collect_shell_vertices(net, shell_surfs, shell_edges):
    """
    Collect shell vertices redundantly from three topology relationships.

    The old implementation trusted only surface['verts'].  The shell vertex
    export could therefore be incomplete when a surface's cached vertex list
    was incomplete even though its edge or reverse vertex adjacency was sound.

    The returned shell is the union of:
      1. surface['verts']
      2. endpoints of surface['edges']
      3. vertices whose reverse 'surfs' adjacency references a shell surface

    """
    surf_vert_set = set()
    edge_vert_set = set()
    reverse_vert_set = set()

    for surf_ndx in shell_surfs:
        try:
            surf = net.surfs.iloc[int(surf_ndx)]
        except (IndexError, TypeError, ValueError):
            continue

        for vert_ndx in _safe_list(surf.get("verts", [])):
            try:
                surf_vert_set.add(int(vert_ndx))
            except (TypeError, ValueError):
                pass

    for edge_ndx in shell_edges:
        try:
            edge = net.edges.iloc[int(edge_ndx)]
        except (IndexError, TypeError, ValueError):
            continue

        for vert_ndx in _safe_list(edge.get("verts", [])):
            try:
                edge_vert_set.add(int(vert_ndx))
            except (TypeError, ValueError):
                pass

    shell_surf_set = set(int(_) for _ in shell_surfs)
    if "surfs" in net.verts.columns:
        for vert_ndx, vert in net.verts.iterrows():
            try:
                attached = set(int(_) for _ in _safe_list(vert.get("surfs", [])))
            except (TypeError, ValueError):
                continue
            if attached & shell_surf_set:
                reverse_vert_set.add(int(vert_ndx))

    shell_vert_set = surf_vert_set | edge_vert_set | reverse_vert_set
    return (
        sorted(shell_vert_set),
        surf_vert_set,
        edge_vert_set,
        reverse_vert_set,
    )


def get_layers(grp, max_layers=50, group_resids=True, build_surfs=True):
    """
    Identify successive topology layers surrounding a Group.

    Index contract
    --------------
    grp.ball_ndxs / grp.layer_atoms
        SYSTEM atom indices.  These are safe to pass to write_pdb(..., sys=grp.sys).

    grp.layer_net_atoms / surface['balls'] / edge['balls'] / vertex['balls']
        NETWORK TOPOLOGY atom identifiers.

    grp.layer_surfs / grp.layer_edges / grp.layer_verts
        Row indices into net.surfs / net.edges / net.verts.

    This separation prevents a network topology index from being accidentally
    exported as a parent-system atom index.
    """
    del build_surfs  # retained in the public API for backward compatibility

    net = grp.net
    t0 = time.perf_counter()

    topology_to_row, topology_to_system, system_to_topology = _build_index_maps(grp)

    # ------------------------------------------------------------------
    # Initial group/body layer: Group membership is SYSTEM-indexed.
    # ------------------------------------------------------------------
    body_system_atoms = [int(_) for _ in grp.ball_ndxs]
    body_topology_atoms = []

    for system_id in body_system_atoms:
        topology_id = system_to_topology.get(system_id)
        if topology_id is None:
            continue
        body_topology_atoms.append(topology_id)

    grp.layer_body_topology = list(body_topology_atoms)

    grp.layer_atoms = [list(body_system_atoms), []]
    grp.layer_net_atoms = [list(body_topology_atoms), []]
    grp.layer_surfs = [[]]
    grp.layer_verts = [[]]
    grp.layer_edges = [[]]
    grp.layer_info = [[0, 0]]

    layer_system_sets = [set(body_system_atoms), set()]
    layer_topology_sets = [set(body_topology_atoms), set()]
    previous_surf_set = set()

    counter = 0

    while counter < max_layers:
        current_topology_atoms = grp.layer_net_atoms[-2]
        current_topology_set = layer_topology_sets[-2]

        next_system_atoms = grp.layer_atoms[-1]
        next_system_set = layer_system_sets[-1]

        next_topology_atoms = grp.layer_net_atoms[-1]
        next_topology_set = layer_topology_sets[-1]

        current_surfs = grp.layer_surfs[-1]
        current_surf_set = set()

        current_edges = grp.layer_edges[-1]
        current_edge_set = set()

        # --------------------------------------------------------------
        # Find surfaces crossing out of the current atom layer.
        # Surface ball IDs are TOPOLOGY IDs.
        # --------------------------------------------------------------
        for topology_id in current_topology_atoms:
            row_pos = topology_to_row.get(int(topology_id))
            if row_pos is None:
                continue

            atom = net.balls.iloc[row_pos]

            for surf_ndx in _safe_list(atom.get("surfs", [])):
                try:
                    surf_ndx = int(surf_ndx)
                except (TypeError, ValueError):
                    continue

                if surf_ndx in current_surf_set or surf_ndx in previous_surf_set:
                    continue

                try:
                    surf = net.surfs.iloc[surf_ndx]
                    ball0, ball1 = [int(_) for _ in surf["balls"]]
                except (IndexError, KeyError, TypeError, ValueError):
                    continue

                ball0_inside = ball0 in current_topology_set
                ball1_inside = ball1 in current_topology_set

                # Internal surface: not part of this layer boundary.
                if ball0_inside and ball1_inside:
                    continue

                # Ignore a malformed/stale adjacency that does not actually
                # touch the current topology layer.
                if not ball0_inside and not ball1_inside:
                    continue

                current_surfs.append(surf_ndx)
                current_surf_set.add(surf_ndx)

                for edge_ndx in _safe_list(surf.get("edges", [])):
                    try:
                        edge_ndx = int(edge_ndx)
                    except (TypeError, ValueError):
                        continue
                    if edge_ndx not in current_edge_set:
                        current_edge_set.add(edge_ndx)
                        current_edges.append(edge_ndx)

                # Atom across the shell surface.
                outside_topology = ball1 if ball0_inside else ball0
                outside_system = topology_to_system.get(outside_topology)

                if outside_system is None:
                    continue

                if outside_topology not in next_topology_set:
                    next_topology_set.add(outside_topology)
                    next_topology_atoms.append(outside_topology)

                if outside_system not in next_system_set:
                    next_system_set.add(outside_system)
                    next_system_atoms.append(outside_system)

        # --------------------------------------------------------------
        # Keep surrounding residues intact, but keep the two index spaces
        # synchronized. Residue atom membership belongs to the SYSTEM.
        # --------------------------------------------------------------
        if group_resids:
            atoms_to_expand = next_system_atoms[:]

            for system_id in atoms_to_expand:
                if system_id < 0 or system_id >= len(grp.sys.balls):
                    continue

                sys_atom = grp.sys.balls.iloc[system_id]
                residue = sys_atom.get("res", None)

                if residue is None:
                    continue

                for residue_system_id in _safe_list(getattr(residue, "atoms", [])):
                    try:
                        residue_system_id = int(residue_system_id)
                    except (TypeError, ValueError):
                        continue

                    residue_topology_id = system_to_topology.get(residue_system_id)
                    if residue_topology_id is None:
                        continue

                    if residue_system_id not in next_system_set:
                        next_system_set.add(residue_system_id)
                        next_system_atoms.append(residue_system_id)

                    if residue_topology_id not in next_topology_set:
                        next_topology_set.add(residue_topology_id)
                        next_topology_atoms.append(residue_topology_id)

        # --------------------------------------------------------------
        # Reconstruct shell vertices redundantly.
        # --------------------------------------------------------------
        shell_verts, _, _, _ = _collect_shell_vertices(
            net,
            current_surfs,
            current_edges,
        )
        grp.layer_verts[-1].extend(shell_verts)

        # --------------------------------------------------------------
        # Layer volume: network geometry using mapped topology rows.
        # --------------------------------------------------------------
        layer_vol = 0.0
        for topology_id in next_topology_atoms:
            row_pos = topology_to_row.get(int(topology_id))
            if row_pos is None:
                continue
            try:
                value = float(net.balls.iloc[row_pos]["vol"])
            except (KeyError, TypeError, ValueError):
                continue
            if value == value:
                layer_vol += value

        # --------------------------------------------------------------
        # Layer shell surface area.
        # --------------------------------------------------------------
        layer_sa = 0.0
        for surf_ndx in current_surfs:
            try:
                value = float(net.surfs.iloc[int(surf_ndx)]["sa"])
            except (IndexError, KeyError, TypeError, ValueError):
                continue
            if value == value:
                layer_sa += value

        grp.layer_info[-1][0] = layer_vol
        grp.layer_info[-1][1] = layer_sa

        # Stop if no shell surfaces were found.
        if len(current_surfs) == 0:
            grp.layer_surfs.pop()
            grp.layer_edges.pop()
            grp.layer_verts.pop()
            grp.layer_info.pop()
            break

        previous_surf_set = current_surf_set

        grp.layer_surfs.append([])
        grp.layer_atoms.append([])
        grp.layer_net_atoms.append([])
        grp.layer_edges.append([])
        grp.layer_verts.append([])
        grp.layer_info.append([0, 0])

        layer_system_sets.append(set())
        layer_topology_sets.append(set())

        counter += 1

    grp.layer_time = time.perf_counter() - t0
