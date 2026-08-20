import time


def get_layers(grp, max_layers=50, group_resids=True, build_surfs=True):
    """
    Identify successive layers of atoms, surfaces, edges, and vertices
    surrounding a Group.

    Results are stored in:
        grp.layer_atoms
        grp.layer_surfs
        grp.layer_edges
        grp.layer_verts
        grp.layer_info
    """

    net = grp.net
    t0 = time.perf_counter()

    # Initial group/body layer
    body_atoms = list(grp.ball_ndxs)
    body_set = set(body_atoms)

    grp.layer_atoms = [body_atoms, []]
    grp.layer_surfs = [[]]
    grp.layer_verts = [[]]
    grp.layer_edges = [[]]
    grp.layer_info = [[0, 0]]

    # Sets used internally for fast membership testing
    layer_atom_sets = [body_set, set()]
    previous_surf_set = set()

    counter = 0

    while counter < max_layers:
        current_atoms = grp.layer_atoms[-2]
        current_atom_set = layer_atom_sets[-2]

        next_atoms = grp.layer_atoms[-1]
        next_atom_set = layer_atom_sets[-1]

        current_surfs = grp.layer_surfs[-1]
        current_surf_set = set()

        current_verts = grp.layer_verts[-1]
        current_vert_set = set()

        current_edges = grp.layer_edges[-1]
        current_edge_set = set()

        # --------------------------------------------------------------
        # Find surfaces crossing out of the current atom layer
        # --------------------------------------------------------------

        for atom_ndx in current_atoms:
            atom = net.balls.iloc[atom_ndx]

            for surf_ndx in atom['surfs']:

                # Already processed in this or previous layer
                if surf_ndx in current_surf_set or surf_ndx in previous_surf_set:
                    continue

                surf = net.surfs.iloc[surf_ndx]
                ball0, ball1 = surf['balls']

                ball0_inside = ball0 in current_atom_set
                ball1_inside = ball1 in current_atom_set

                # Surface is completely internal to the current layer
                if ball0_inside and ball1_inside:
                    continue

                current_surfs.append(surf_ndx)

                current_surf_set.add(surf_ndx)

                # Vertices
                for vert_ndx in surf['verts']:
                    if vert_ndx not in current_vert_set:
                        current_vert_set.add(vert_ndx)
                        current_verts.append(vert_ndx)

                # Edges
                for edge_ndx in surf['edges']:
                    if edge_ndx not in current_edge_set:
                        current_edge_set.add(edge_ndx)
                        current_edges.append(edge_ndx)

                # Atom across the surface
                if ball0_inside and not ball1_inside:
                    if ball1 not in next_atom_set:
                        next_atom_set.add(ball1)
                        next_atoms.append(ball1)

                elif ball1_inside and not ball0_inside:
                    if ball0 not in next_atom_set:
                        next_atom_set.add(ball0)
                        next_atoms.append(ball0)

        # --------------------------------------------------------------
        # Keep residues intact
        # --------------------------------------------------------------

        if group_resids:
            # Iterate over a copy because this loop can add atoms.
            atoms_to_expand = next_atoms[:]

            for atom_ndx in atoms_to_expand:
                atom = net.balls.iloc[atom_ndx]

                if 'res' not in atom or atom['res'] is None:
                    continue

                for resid_atom in atom['res'].atoms:
                    if resid_atom not in next_atom_set:
                        next_atom_set.add(resid_atom)
                        next_atoms.append(resid_atom)

        # --------------------------------------------------------------
        # Layer volume
        # --------------------------------------------------------------

        layer_vol = 0.0

        for atom_ndx in next_atoms:
            layer_vol += net.balls.iloc[atom_ndx]['vol']

        # --------------------------------------------------------------
        # Layer surface area
        # --------------------------------------------------------------

        layer_sa = 0.0

        for surf_ndx in current_surfs:
            layer_sa += net.surfs['sa'][surf_ndx]

        grp.layer_info[-1][0] = layer_vol
        grp.layer_info[-1][1] = layer_sa

        # --------------------------------------------------------------
        # Stop if no surfaces were found
        # --------------------------------------------------------------

        if len(current_surfs) == 0:
            grp.layer_surfs.pop()
            break

        # --------------------------------------------------------------
        # Prepare next layer
        # --------------------------------------------------------------

        previous_surf_set = current_surf_set

        grp.layer_surfs.append([])
        grp.layer_atoms.append([])
        grp.layer_edges.append([])
        grp.layer_verts.append([])
        grp.layer_info.append([0, 0])

        layer_atom_sets.append(set())

        counter += 1
    print(
        f"LAYER PROFILE | layers={counter} "
        f"atoms={sum(len(x) for x in grp.layer_atoms):,} "
        f"surfs={sum(len(x) for x in grp.layer_surfs):,} "
        f"edges={sum(len(x) for x in grp.layer_edges):,} "
        f"verts={sum(len(x) for x in grp.layer_verts):,} "
        f"time={time.perf_counter() - t0:.3f}s"
    )
