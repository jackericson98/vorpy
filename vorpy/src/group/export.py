import os
from os import path
from vorpy.src.output import write_pdb
from vorpy.src.output import write_atom_cells
from vorpy.src.output import write_logs
from vorpy.src.output import write_surfs
from vorpy.src.output import write_edges
from vorpy.src.output import write_off_verts


def _group_topology_indices(grp):
    """
    Convert grp.ball_ndxs (SYSTEM indices) to network TOPOLOGY identifiers.
    """
    net_balls = getattr(grp.net, "balls", None)
    if net_balls is None or len(net_balls) == 0:
        return list(grp.ball_ndxs)

    system_to_topology = {}

    for row_pos, (_, atom) in enumerate(net_balls.iterrows()):
        try:
            topology_id = int(atom.get("num", row_pos))
        except (TypeError, ValueError):
            topology_id = int(row_pos)

        system_id = None
        if "system_num" in net_balls.columns:
            try:
                value = atom.get("system_num", None)
                if value is not None and not (isinstance(value, float) and value != value):
                    system_id = int(value)
            except (TypeError, ValueError):
                system_id = None

        if system_id is None:
            system_id = topology_id

        system_to_topology.setdefault(system_id, topology_id)

    return [
        system_to_topology[system_id]
        for system_id in (int(_) for _ in grp.ball_ndxs)
        if system_id in system_to_topology
    ]


def export_info(grp, directory=None):
    """
    Export a detailed human-readable summary of a VorPy group.
    """

    import numpy as np

    # --------------------------------------------------------------
    # Directory
    # --------------------------------------------------------------

    if directory is not None and os.path.exists(directory):
        os.chdir(directory)

    os.chdir(grp.dir)

    # Make sure group-level geometry is current
    grp.get_info()

    # --------------------------------------------------------------
    # Helpers
    # --------------------------------------------------------------

    def fmt(value, decimals=5):
        try:
            value = float(value)
            if np.isfinite(value):
                return f"{value:.{decimals}f}"
        except (TypeError, ValueError):
            pass
        return "N/A"

    def get_attr(obj, name, default=None):
        value = getattr(obj, name, default)
        return default if value is None else value

    def safe_sum(values):
        vals = []
        for value in values:
            try:
                value = float(value)
                if np.isfinite(value):
                    vals.append(value)
            except (TypeError, ValueError):
                pass
        return sum(vals)

    # --------------------------------------------------------------
    # Group atoms
    # --------------------------------------------------------------

    ball_ndxs = list(grp.ball_ndxs)

    if len(ball_ndxs) > 0:
        atoms = grp.sys.balls.iloc[ball_ndxs]
    else:
        atoms = grp.sys.balls.iloc[[]]

    # --------------------------------------------------------------
    # Surface-energy estimate
    # --------------------------------------------------------------

    # Prefer the value already calculated and stored on the group.
    # Fall back to the current representative model only for objects
    # created before surf_energy was added.
    int_mean_sq = get_attr(grp, "int_mean_curv_sq", 0.0)
    representative_energy = get_attr(grp, "surf_energy", None)

    if representative_energy is None:
        representative_energy = 2.0 * float(int_mean_sq or 0.0)

    # --------------------------------------------------------------
    # Write file
    # --------------------------------------------------------------

    with open("info.txt", "w", encoding="utf-8") as info:

        # ==========================================================
        # HEADER
        # ==========================================================

        info.write("=" * 72 + "\n")
        info.write(f"VORPY GROUP INFORMATION: {grp.name}\n")
        info.write("=" * 72 + "\n\n")

        info.write(f"System: {grp.sys.name}\n")
        info.write(f"Group:  {grp.name}\n\n")

        # ==========================================================
        # COMPOSITION
        # ==========================================================

        info.write("COMPOSITION\n")
        info.write("-" * 72 + "\n")

        info.write(f"Atoms:     {len(ball_ndxs)}\n")
        info.write(f"Residues:  {len(grp.rsds)}\n")
        info.write(f"Chains:    {len(grp.chns)}\n")

        info.write(
            f"Mass:      {fmt(get_attr(grp, 'mass'))} Da\n"
        )

        info.write("\n")

        # ==========================================================
        # BUILD INFORMATION
        # ==========================================================

        info.write("BUILD INFORMATION\n")
        info.write("-" * 72 + "\n")

        net = grp.net
        settings = getattr(net, "settings", {}) or {}
        metrics = getattr(net, "metrics", {}) or {}

        # Settings used to construct this group's independent network.
        setting_labels = [
            ("net_type", "Network Type"),
            ("surf_res", "Surface Resolution"),
            ("box_size", "Box Size"),
            ("max_vert", "Maximum Allowable Vertex"),
        ]

        for key, label in setting_labels:
            try:
                value = settings.get(key, None)
            except AttributeError:
                try:
                    value = settings[key]
                except (KeyError, TypeError):
                    value = None

            if value is not None:
                info.write(f"{label}: {value}\n")

        info.write("\nBUILD TIMING\n")
        info.write("-" * 72 + "\n")

        # Canonical timing fields currently stored by Network.
        timing_labels = [
            ("vert", "Vertex Time"),
            ("con", "Connection Time"),
            ("surf", "Surface Building Time"),
            ("anal", "Analysis Time"),
            ("tot", "Total Time"),
        ]

        printed_metric_keys = set()

        for key, label in timing_labels:
            try:
                value = metrics.get(key, None)
            except AttributeError:
                try:
                    value = metrics[key]
                except (KeyError, TypeError):
                    value = None

            if value is not None:
                try:
                    seconds = float(value)
                    if np.isfinite(seconds):
                        info.write(f"{label}: {seconds:.5f} s\n")
                        printed_metric_keys.add(key)
                except (TypeError, ValueError):
                    pass

        # Include any additional numeric timing metrics produced by newer
        # builds without hard-coding them into this exporter.
        try:
            metric_items = metrics.items()
        except AttributeError:
            metric_items = []

        extras = []
        ignored_metric_keys = {
            "START", "start", "max_vert", "max_vertex",
            "num_balls", "num_verts", "num_edges", "num_surfs",
        }

        for key, value in metric_items:
            if key in printed_metric_keys or key in ignored_metric_keys:
                continue

            key_lower = str(key).lower()
            if "time" not in key_lower and key_lower not in {
                "setup", "doublet", "edge", "filter", "adjacency",
                "surface", "packaging", "validation"
            }:
                continue

            try:
                seconds = float(value)
            except (TypeError, ValueError):
                continue

            if not np.isfinite(seconds):
                continue

            extras.append((str(key), seconds))

        if extras:
            info.write("\nAdditional Timing Metrics:\n")
            for key, seconds in extras:
                label = key.replace("_", " ").strip().title()
                info.write(f"  {label}: {seconds:.5f} s\n")

        info.write("\n")

        # ==========================================================
        # NETWORK
        # ==========================================================

        info.write("VORONOI NETWORK\n")
        info.write("-" * 72 + "\n")

        info.write(
            f"Vertices:  {len(grp.net.verts):,}\n"
        )

        info.write(
            f"Edges:     {len(grp.net.edges):,}\n"
        )

        info.write(
            f"Surfaces:  {len(grp.net.surfs):,}\n"
        )

        info.write("\n")

        # ==========================================================
        # GEOMETRY
        # ==========================================================

        info.write("GROUP GEOMETRY\n")
        info.write("-" * 72 + "\n")

        info.write(
            f"Volume:                    "
            f"{fmt(grp.vol)} Å³\n"
        )

        info.write(
            f"van der Waals Volume:      "
            f"{fmt(get_attr(grp, 'vdw_vol'))} Å³\n"
        )

        info.write(
            f"Surface Area:              "
            f"{fmt(grp.sa)} Å²\n"
        )

        info.write(
            f"Density:                   "
            f"{fmt(grp.density)}\n"
        )

        # Useful compactness descriptor
        if grp.vol and grp.vol > 0:
            sa_vol = grp.sa / grp.vol
            info.write(
                f"Surface Area / Volume:     "
                f"{fmt(sa_vol)} Å⁻¹\n"
            )

        info.write("\n")

        # ==========================================================
        # CURVATURE
        # ==========================================================

        info.write("SURFACE CURVATURE\n")
        info.write("-" * 72 + "\n")

        info.write(
            f"Integrated Mean Curvature:          "
            f"{fmt(get_attr(grp, 'int_mean_curv'))} Å\n"
        )

        info.write(
            f"Integrated Mean Curvature Squared:  "
            f"{fmt(get_attr(grp, 'int_mean_curv_sq'))}\n"
        )

        info.write(
            f"Integrated Gaussian Curvature:      "
            f"{fmt(get_attr(grp, 'int_gauss_curv'))}\n"
        )

        info.write("\n")

        info.write(
            f"Area-Weighted Mean Curvature:        "
            f"{fmt(get_attr(grp, 'avg_mean_curv'))} Å⁻¹\n"
        )

        info.write(
            f"Area-Weighted Gaussian Curvature:    "
            f"{fmt(get_attr(grp, 'avg_gauss_curv'))} Å⁻²\n"
        )

        info.write("\n")

        # ==========================================================
        # SURFACE ENERGY ESTIMATE
        # ==========================================================

        info.write("SURFACE ENERGY ESTIMATE\n")
        info.write("-" * 72 + "\n")

        info.write(
            f"Representative Surface Energy:       "
            f"{fmt(representative_energy)} kBT\n"
        )

        if grp.sa and grp.sa > 0:
            info.write(
                f"Representative Energy / Area:        "
                f"{fmt(representative_energy / grp.sa)} kBT Å⁻²\n"
            )

        info.write("\n")

        # ==========================================================
        # CHAIN / RESIDUE GEOMETRY
        # ==========================================================

        # PDB identity lives on sys.balls, while solved Voronoi geometry
        # and topology live on grp.net.balls / grp.net.surfs.
        #
        # Surface-area rule:
        #   - A surface between two atoms in the same residue is internal
        #     to that residue and contributes 0 to residue surface area.
        #   - A surface between different residues contributes its full area
        #     to each participating residue.
        #   - The identical rule is applied at the chain level.
        #
        # This produces the boundary area of each residue/chain rather than
        # the sum of its atoms' individual surface areas.

        from collections import defaultdict

        net_balls = getattr(grp.net, "balls", None)
        net_surfs = getattr(grp.net, "surfs", None)

        group_system_indices = set(int(_) for _ in ball_ndxs)

        atom_meta = {}
        residue_atoms = defaultdict(list)
        chain_atoms = defaultdict(list)

        if net_balls is not None and len(net_balls) > 0:

            for net_index, net_atom in net_balls.iterrows():

                # Topology tables normally refer to net.balls["num"].
                try:
                    topology_index = int(net_atom.get("num", net_index))
                except (TypeError, ValueError):
                    topology_index = int(net_index)

                # Map the network atom back to sys.balls using the same
                # convention as VorPy's log writer:
                #
                #   1. system_num, when present
                #   2. num, otherwise
                #
                # Do NOT fall back through grp.ball_ndxs here. That list only
                # describes atoms belonging to the group and therefore cannot
                # map surrounding atoms such as solvent.

                sys_index = None

                if "system_num" in net_balls.columns:
                    try:
                        value = net_atom["system_num"]
                        if value is not None and not (
                                isinstance(value, float) and np.isnan(value)
                        ):
                            sys_index = int(value)
                    except (TypeError, ValueError):
                        sys_index = None

                if sys_index is None:
                    try:
                        sys_index = int(net_atom.get("num", net_index))
                    except (TypeError, ValueError):
                        continue

                try:
                    sys_atom = grp.sys.balls.iloc[sys_index]
                except (IndexError, KeyError, TypeError):
                    continue

                chain_name = sys_atom.get(
                    "chain_name",
                    sys_atom.get("chain", "")
                )

                res_name = sys_atom.get(
                    "res_name",
                    sys_atom.get("residue", "")
                )

                res_seq = sys_atom.get(
                    "res_seq",
                    sys_atom.get("residue_sequence", "")
                )

                # Chain is included in the residue key so identical residue
                # numbers on different chains are never merged.
                residue_key = (
                    str(chain_name),
                    res_seq,
                    str(res_name),
                )

                chain_key = str(chain_name)

                in_group = sys_index in group_system_indices

                atom_meta[topology_index] = {
                    "system_index": sys_index,
                    "chain_key": chain_key,
                    "residue_key": residue_key,
                    "res_name": str(res_name).strip().upper(),
                    "in_group": in_group,
                    "vol": net_atom.get("vol", np.nan),
                }

                if in_group:
                    residue_atoms[residue_key].append(topology_index)
                    chain_atoms[chain_key].append(topology_index)

        # ----------------------------------------------------------
        # Volumes: additive over solved atom cells
        # ----------------------------------------------------------

        residue_volumes = defaultdict(float)
        chain_volumes = defaultdict(float)

        for topology_index, meta in atom_meta.items():

            if not meta["in_group"]:
                continue

            try:
                atom_vol = float(meta["vol"])
            except (TypeError, ValueError):
                continue

            if not np.isfinite(atom_vol):
                continue

            residue_volumes[meta["residue_key"]] += atom_vol
            chain_volumes[meta["chain_key"]] += atom_vol

        # ----------------------------------------------------------
        # Surface-area decomposition
        # ----------------------------------------------------------

        # Three distinct surface-area measures are reported:
        #
        # 1. Total Boundary SA
        #    Any surface that is not internal to the same residue/chain.
        #
        # 2. Inter-Residue / Inter-Chain SA
        #    Surfaces between different residues/chains when BOTH atoms
        #    are members of this solved group.
        #
        # 3. Solvent-Interfacial SA
        #    Surfaces where the opposing atom is a recognized water residue.
        #
        # These quantities are intentionally overlapping:
        # solvent-interfacial SA is part of total boundary SA.

        water_names = {
            "SOL", "HOH", "WAT", "H2O",
            "TIP3", "TIP3P", "TIP4", "TIP4P",
            "SPC", "SPCE",
        }

        residue_total_sa = defaultdict(float)
        residue_inter_sa = defaultdict(float)
        residue_sol_sa = defaultdict(float)

        chain_total_sa = defaultdict(float)
        chain_inter_sa = defaultdict(float)
        chain_sol_sa = defaultdict(float)

        if net_surfs is not None and len(net_surfs) > 0:

            for _, surf in net_surfs.iterrows():

                try:
                    ball1, ball2 = surf["balls"]
                    ball1 = int(ball1)
                    ball2 = int(ball2)
                    surf_area = float(surf["sa"])
                except (KeyError, TypeError, ValueError):
                    continue

                if not np.isfinite(surf_area):
                    continue

                meta1 = atom_meta.get(ball1)
                meta2 = atom_meta.get(ball2)

                if meta1 is None or meta2 is None:
                    continue

                ball1_is_water = meta1["res_name"] in water_names
                ball2_is_water = meta2["res_name"] in water_names

                # ==================================================
                # RESIDUE LEVEL
                # ==================================================

                same_group_residue = (
                        meta1["in_group"]
                        and meta2["in_group"]
                        and meta1["residue_key"] == meta2["residue_key"]
                )

                # Total residue boundary area:
                # exclude only surfaces internal to the same residue.
                if not same_group_residue:
                    if meta1["in_group"]:
                        residue_total_sa[meta1["residue_key"]] += surf_area
                    if meta2["in_group"]:
                        residue_total_sa[meta2["residue_key"]] += surf_area

                # Inter-residue area:
                # both atoms must belong to the group and to different residues.
                if (
                        meta1["in_group"]
                        and meta2["in_group"]
                        and meta1["residue_key"] != meta2["residue_key"]
                ):
                    residue_inter_sa[meta1["residue_key"]] += surf_area
                    residue_inter_sa[meta2["residue_key"]] += surf_area

                # Solvent-interfacial residue area.
                if meta1["in_group"] and ball2_is_water:
                    residue_sol_sa[meta1["residue_key"]] += surf_area

                if meta2["in_group"] and ball1_is_water:
                    residue_sol_sa[meta2["residue_key"]] += surf_area

                # ==================================================
                # CHAIN LEVEL
                # ==================================================

                same_group_chain = (
                        meta1["in_group"]
                        and meta2["in_group"]
                        and meta1["chain_key"] == meta2["chain_key"]
                )

                # Total chain boundary area:
                # exclude only surfaces internal to the same chain.
                if not same_group_chain:
                    if meta1["in_group"]:
                        chain_total_sa[meta1["chain_key"]] += surf_area
                    if meta2["in_group"]:
                        chain_total_sa[meta2["chain_key"]] += surf_area

                # Inter-chain area:
                # both atoms must belong to the group and to different chains.
                if (
                        meta1["in_group"]
                        and meta2["in_group"]
                        and meta1["chain_key"] != meta2["chain_key"]
                ):
                    chain_inter_sa[meta1["chain_key"]] += surf_area
                    chain_inter_sa[meta2["chain_key"]] += surf_area

                # Solvent-interfacial chain area.
                if meta1["in_group"] and ball2_is_water:
                    chain_sol_sa[meta1["chain_key"]] += surf_area

                if meta2["in_group"] and ball1_is_water:
                    chain_sol_sa[meta2["chain_key"]] += surf_area

        # ----------------------------------------------------------
        # Surface classification summary
        # ----------------------------------------------------------

        mapped_group_atoms = sum(
            1 for meta in atom_meta.values()
            if meta["in_group"]
        )

        mapped_water_atoms = sum(
            1 for meta in atom_meta.values()
            if meta["res_name"] in water_names
        )

        info.write("SURFACE CLASSIFICATION\n")
        info.write("-" * 72 + "\n")
        info.write(
            f"Mapped Group Atoms:       {mapped_group_atoms:,}\n"
        )
        info.write(
            f"Mapped Surrounding Waters: {mapped_water_atoms:,}\n"
        )
        info.write("\n")

        # ==========================================================
        # CHAINS
        # ==========================================================

        info.write("CHAIN COMPOSITION\n")
        info.write("-" * 72 + "\n")

        if chain_atoms:

            for chain_key, topology_indices in chain_atoms.items():
                # Count unique residues represented by the group atoms in
                # this chain.
                chain_residues = {
                    atom_meta[ndx]["residue_key"]
                    for ndx in topology_indices
                    if ndx in atom_meta
                }

                info.write(
                    f"Chain {chain_key}: "
                    f"{len(topology_indices)} atoms, "
                    f"{len(chain_residues)} residues  "
                    f"Volume: {chain_volumes[chain_key]:.3f} Å³  "
                    f"Total Boundary SA: {chain_total_sa[chain_key]:.3f} Å²  "
                    f"Inter-Chain SA: {chain_inter_sa[chain_key]:.3f} Å²  "
                    f"Solvent-Interfacial SA: {chain_sol_sa[chain_key]:.3f} Å²\n"
                )

        else:
            info.write(
                "Detailed solved chain geometry unavailable.\n"
            )

        info.write("\n")

        # ==========================================================
        # RESIDUES
        # ==========================================================

        info.write("RESIDUE COMPOSITION\n")
        info.write("-" * 72 + "\n")

        if residue_atoms:

            for residue_key, topology_indices in residue_atoms.items():
                chain_name, res_seq, res_name = residue_key

                info.write(
                    f"{str(chain_name):>3}  "
                    f"{str(res_name):>4} "
                    f"{str(res_seq):>5}  "
                    f"{len(topology_indices):>4} atoms  "
                    f"Volume: {residue_volumes[residue_key]:.3f} Å³  "
                    f"Total Boundary SA: {residue_total_sa[residue_key]:.3f} Å²  "
                    f"Inter-Residue SA: {residue_inter_sa[residue_key]:.3f} Å²  "
                    f"Solvent-Interfacial SA: {residue_sol_sa[residue_key]:.3f} Å²\n"
                )

        else:
            info.write(
                "Detailed solved residue geometry unavailable.\n"
            )

        info.write("\n")
        info.write("=" * 72 + "\n")


def group_exports(grp, all_=False, atoms=False, atom_surfs=False, atom_edges=False, atom_verts=False, surfs=False,
                  sep_surfs=False, shell_surfs=False, edges=False, sep_edges=False, shell_edges=False,
                  verts=False, sep_verts=False, shell_verts=False, layers=-1, info=False, surr_atoms=False, logs=False,
                  ext_atoms=False, concave_colors=False, round_to=3, file_type=None):
    """
    Exports various components of a Group object to files based on specified parameters.
    This function provides flexible export options for different aspects of a molecular group,
    including atoms, surfaces, edges, vertices, and surrounding elements.

    Parameters
    ----------
    grp : Group
        The Group object containing the data to be exported
    all_ : bool, optional
        If True, exports all possible components of the group. Default is False
    atoms : bool, optional
        If True, exports a PDB file containing only the atoms of the group. Default is False
    atom_surfs : bool, optional
        If True, exports the surfaces associated with each atom. Default is False
    atom_edges : bool, optional
        If True, exports the edges associated with each atom. Default is False
    atom_verts : bool, optional
        If True, exports the vertices associated with each atom. Default is False
    surfs : bool, optional
        If True, exports all surfaces in the group as a single object. Default is False
    sep_surfs : bool, optional
        If True, exports each surface as a separate file, named by their constituent atoms. Default is False
    shell_surfs : bool, optional
        If True, exports all surfaces for the group's shell. Default is False
    edges : bool, optional
        If True, exports all edges in the group as a single object. Default is False
    sep_edges : bool, optional
        If True, exports each edge as a separate file. Default is False
    shell_edges : bool, optional
        If True, exports all edges for the group's shell. Default is False
    verts : bool, optional
        If True, exports all vertices as a single OFF file. Default is False
    sep_verts : bool, optional
        If True, exports each vertex as a separate file. Default is False
    shell_verts : bool, optional
        If True, exports all vertices for the group's shell. Default is False
    layers : int, optional
        Number of layers to export around the group. If -1, exports all layers. Default is -1
    info : bool, optional
        If True, exports group information to info.txt. Default is False
    surr_atoms : bool, optional
        If True, exports atoms directly surrounding the group with intact residues. Default is False
    logs : bool, optional
        If True, exports log files. Default is False
    ext_atoms : bool, optional
        If True, exports the outermost atoms in the group's shell. Default is False
    concave_colors : bool, optional
        If True, exports the concave colors for the surfaces. Default is False
    round_to : int, optional
    file_type : {'off', 'ply', 'vtp'}, optional
        Geometry format used for surface, edge, and vertex meshes. Default is 'off'.

    Returns
    -------
    None

    Notes
    -----
    - All exports are written to the group's directory (grp.dir)
    - If the directory doesn't exist, it will be created
    - Surface colors and schemes are inherited from network settings if not specified
    - For atom-related exports, a subdirectory 'atoms' is created if needed
    - Layer exports require the group to have calculated layers first

    Examples
    --------
    # Export all components of a group
    >>> group_exports(my_group, all_=True)
    # This will create a comprehensive export of the group, including:
    # - A PDB file of all atoms
    # - All surfaces, edges, and vertices
    # - Individual atom components (surfaces, edges, vertices)
    # - Layer information
    # - Group information
    # - Log files

    # Export only surfaces and edges
    >>> group_exports(my_group, surfs=True, edges=True)
    # This will create:
    # - A single file containing all surfaces
    # - A single file containing all edges
    # Useful for visualization of the group's structure without atom details

    # Export atom-related components
    >>> group_exports(my_group, atoms=True, atom_surfs=True, atom_edges=True)
    # This will create:
    # - A PDB file of the group's atoms
    # - Individual surface files for each atom in the 'atoms' subdirectory
    # - Individual edge files for each atom in the 'atoms' subdirectory
    # Useful for detailed analysis of individual atoms

    # Export surrounding atoms and information
    >>> group_exports(my_group, surr_atoms=True, info=True)
    # This will create:
    # - A file containing atoms surrounding the group
    # - An info.txt file with group statistics
    # Useful for analyzing the group's environment and properties
    """
    # Set the surface colors and scheme
    if file_type is None:
        file_type = getattr(grp, 'settings', {}).get('file_type', 'off')
    file_type = str(file_type).strip().lower().lstrip('.')
    if file_type == 'vtk':
        file_type = 'vtp'
    if file_type not in {'off', 'ply', 'vtp'}:
        raise ValueError("file_type must be 'off', 'ply', or 'vtp'")
    if grp.settings['surf_col'] is None:
        grp.settings['surf_col'] = grp.net.settings['surf_col']
    # Set the surface scheme
    if grp.settings['surf_scheme'] is None:
        grp.settings['surf_scheme'] = grp.net.settings['surf_scheme']
    # Get the surfaces if they haven't been got
    if grp.net.surfs is None or len(grp.net.surfs) == 0:
        return
    # Create the output directory inside the system's directory
    if grp.dir is None:
        i = 1
        my_dir = grp.sys.files['dir'] + "/" + grp.name
        first = True
        while os.path.exists(my_dir):
            if first:
                my_dir += "__"
                first = False
            my_dir = my_dir[:-(1 + len(str(i)))] + '_' + str(i)
            i += 1
        grp.dir = my_dir
        os.mkdir(grp.dir)
    # Go back to the group directory
    os.chdir(grp.dir)
    # Build layer information once if any requested export requires it
    needs_layers = shell_surfs or shell_edges or shell_verts or surr_atoms or ext_atoms or layers > 0 or all_

    if needs_layers and grp.layer_surfs is None:
        grp.get_layers(max_layers=max(1, layers) if layers > 0 else 1)
    # Export the log file first
    if logs or all_:
        write_logs(grp, round_to=round_to)
    # If the user wants to export the atoms for the group
    if atoms or all_:
        if grp.sys.files['base_file'][-3:] == 'txt':
            pass
        else:
            write_pdb(atoms=grp.ball_ndxs, file_name="group_atoms", sys=grp.sys)
    # If the atoms surfaces are selected go for it
    if atom_verts or atom_edges or atom_surfs or all_:
        if not path.exists(grp.dir + '/atoms'):
            os.mkdir(grp.dir + '/atoms')
        write_atom_cells(grp.net, atoms=_group_topology_indices(grp), directory=grp.dir + '/atoms', surfs=atom_surfs or all_,
                         edges=atom_edges or all_, verts=atom_verts or all_, concave_colors=concave_colors,
                         file_type=file_type)
        os.chdir(grp.dir)

    # If the user wants to export the shell for the group
    if shell_surfs or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # noinspection PyUnresolvedReferences
        if grp.layer_surfs is not None and len(grp.layer_surfs) > 0:
            write_surfs(net=grp.net, surfs=grp.layer_surfs[0], file_name="shell_surfs", directory=grp.dir,
                        concave_colors=concave_colors, ref_surfs=_group_topology_indices(grp), universal_max=False,
                        file_type=file_type)
    # If the user wants all of the surfaces in one file
    if surfs or all_:
        write_surfs(grp.net, [i for i in range(len(grp.net.surfs))], 'surfs', file_type=file_type)
    # Separate surfaces
    if sep_surfs or all_:
        # Make the surfaces directory
        if not os.path.exists(grp.dir + '/surfs'):
            os.mkdir(grp.dir + '/surfs')
        # Create the surfaces' files
        for j, my_surf in grp.net.surfs.iterrows():
            write_surfs(grp.net, [j], file_name='b{}_b{}'.format(*my_surf['balls']),
                        directory=grp.dir + '/surfs', file_type=file_type)
    # Shell edges
    if shell_edges or all_:
        if grp.layer_edges is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        write_edges(grp.net, grp.layer_edges[0], file_name="shell_edges", directory=grp.dir,
                    color=grp.settings['edge_col'], file_type=file_type)
    # All one big edge file
    if edges or all_:
        write_edges(grp.net, edges=[i for i in range(len(grp.net.edges))], file_name="edges", directory=grp.dir,
                    color=grp.settings['edge_col'], file_type=file_type)
    # If the separate edges are called
    if sep_edges or all_:
        # Make the edges directory
        if not os.path.exists(grp.dir + '/edges'):
            os.mkdir(grp.dir + '/edges')
        for j, my_edge in grp.net.edges.iterrows():
            write_edges(grp.net, [j], 'b{}_b{}_b{}'.format(*my_edge['balls']),
                        directory=grp.dir + '/edges', file_type=file_type)
    # Run the separate vertices
    if sep_verts:
        # Make the vertices directory
        if not path.exists(grp.dir + '/verts'):
            os.mkdir(grp.dir + "/verts")
        for j, vert in grp.net.verts.iterrows():
            write_off_verts(grp.net, [j], 'b{}_b{}_b{}_b{}'.format(*vert['balls']),
                            directory=grp.dir + "/verts", file_type=file_type)
    # Export all the vertices in one file
    if verts or all_:
        write_off_verts(grp.net, [i for i in range(len(grp.net.verts))], directory=grp.dir, file_name='verts',
                        color=grp.settings['vert_col'], file_type=file_type)
    # Export the shell vertices
    if shell_verts or all_:
        if grp.layer_verts is None:
            grp.get_layers(max_layers=1, build_surfs=False)
        write_off_verts(grp.net, grp.layer_verts[0], file_name="shell_verts", directory=grp.dir,
                        color=grp.settings['vert_col'], file_type=file_type)
    # If the user wants layers
    if layers > 0 or all_:
        # First check to see if the number of layers is greater than 1
        if grp.layer_atoms is None or len(grp.layer_atoms) <= 1:
            grp.get_layers(max_layers=layers)
        # Create the layers directory
        i = 1
        my_dir = os.getcwd() + "/layers"
        while os.path.exists(my_dir):
            if my_dir[-1] == 's':
                my_dir += '__'
            my_dir = my_dir[:-2] + str(i)
            i += 1
        os.mkdir(my_dir)
        os.chdir(my_dir)
        # Create the layer and atoms files
        for i in range(len(grp.layer_surfs)):
            write_pdb(grp.layer_atoms[i + 1], file_name=str(i) + "_atoms", sys=grp.sys)
            write_surfs(grp.net, grp.layer_surfs[i], file_name=str(i) + "_surfs", file_type=file_type)
        # If the user wants info and layers create a layers info file
        if info or all_:
            # Create the information file
            info = open(grp.name + "_layer_info.txt", 'w')
            info.write(grp.name + " body: \n")
            # Go through the layers in the group's layers
            for i in range(len(grp.layer_surfs)):
                info.write("Number of atoms: " + str(len(grp.layer_atoms[i])) + "\n")
                info.write("Volume: " + str(grp.layer_info[i][0]) + "\n")
                info.write("Surface Area: " + str(grp.layer_info[i][1]) + "\n")
            info.close()
        # Change back to the group directory
        os.chdir(grp.dir)
    # If the user wants a full information file on the group
    if info or all_:
        export_info(grp)
    # Surrounding atoms
    if surr_atoms or all_:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        try:
            write_pdb(atoms=grp.layer_atoms[1], file_name="surr_atoms", directory=grp.dir, sys=grp.sys)
        except IndexError:
            pass
    if (ext_atoms or all_) and len(grp.atms) > 15:
        if grp.layer_surfs is None:
            # Get the first layer
            grp.get_layers(max_layers=1)
        # write the surrounding atoms
        write_pdb(sys=grp.sys, atoms=grp.layer_atoms[0], file_name="ext_atoms", directory=grp.dir)
    # Check to see if there is verts file in the system directory
    for file in os.listdir(grp.sys.files['dir']):
        if file.endswith("_verts.txt"):
            source = os.path.join(grp.sys.files["dir"], file)
            destination = os.path.join(grp.dir, file)

            if os.path.exists(destination):
                print(f"Skipping existing vertex file: {destination}")
                continue

            os.rename(source, destination)
    os.chdir("..")
    # Change back to the system directory
    os.chdir(grp.sys.files['dir'])
