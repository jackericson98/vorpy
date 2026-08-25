import os


def export_sys_info(sys):
    """
    Export a top-level manifest for a VorPy System.

    The system info file contains:
      - input/output locations
      - molecular composition derived from the input structure
      - group metadata and build settings
      - interface metadata

    Solved geometric quantities, curvature, surface energy, and build timing
    remain in the individual group/interface output folders.
    """

    balls = getattr(sys, "balls", None)
    files = getattr(sys, "files", {}) or {}

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def first_existing(columns, *names):
        for name in names:
            if name in columns:
                return name
        return None

    def clean(value):
        if value is None:
            return ""
        value = str(value).strip()
        return value if value else "-"

    def object_len(obj, *names):
        for name in names:
            value = getattr(obj, name, None)
            if value is not None:
                try:
                    return len(value)
                except TypeError:
                    pass
        return 0

    def get_mapping_value(mapping, key, default=None):
        if mapping is None:
            return default

        try:
            return mapping.get(key, default)
        except AttributeError:
            try:
                return mapping[key]
            except (KeyError, TypeError):
                return default

    def get_interfaces(system):
        """
        Support the common System interface-container names without making
        the info exporter depend on one exact internal spelling.
        """
        for attr in ("interfaces", "ifaces"):
            value = getattr(system, attr, None)
            if value is not None:
                return value
        return []

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------

    base_file = files.get("base_file", None)
    output_dir = files.get("dir", None)

    input_location = os.path.abspath(base_file) if base_file else "N/A"
    output_location = os.path.abspath(output_dir) if output_dir else os.getcwd()

    output_path = (
        os.path.join(output_dir, sys.name + "_info.txt")
        if output_dir
        else sys.name + "_info.txt"
    )

    # ------------------------------------------------------------------
    # Write file
    # ------------------------------------------------------------------

    with open(output_path, "w", encoding="utf-8") as info:

        info.write("=" * 72 + "\n")
        info.write(f"VORPY SYSTEM INFORMATION: {sys.name}\n")
        info.write("=" * 72 + "\n\n")

        # ==============================================================
        # FILE LOCATIONS
        # ==============================================================

        info.write("FILE LOCATIONS\n")
        info.write("-" * 72 + "\n")
        info.write(f"Input File:       {input_location}\n")
        info.write(f"Output Directory: {output_location}\n")
        info.write("\n")

        # ==============================================================
        # MOLECULAR COMPOSITION
        # ==============================================================

        info.write("MOLECULAR COMPOSITION\n")
        info.write("-" * 72 + "\n")

        if balls is None:
            info.write("Atom information unavailable.\n\n")

        else:
            cols = set(balls.columns)

            res_name_col = first_existing(
                cols,
                "res_name",
                "Residue",
                "residue",
            )

            res_seq_col = first_existing(
                cols,
                "res_seq",
                "Residue Sequence",
                "residue_sequence",
            )

            chain_col = first_existing(
                cols,
                "chain_name",
                "chain",
                "Chain",
            )

            element_col = first_existing(
                cols,
                "element",
                "Element",
            )

            # ----------------------------------------------------------
            # Chemistry-aware classification
            # ----------------------------------------------------------

            water_names = {
                "SOL", "HOH", "WAT", "H2O",
                "TIP3", "TIP3P", "TIP4", "TIP4P",
                "SPC", "SPCE"
            }

            ion_names = {
                "NA", "NA+", "SOD",
                "K", "K+", "POT",
                "CL", "CL-", "CLA",
                "MG", "MG2", "MG2+",
                "CA", "CA2", "CA2+",
                "ZN", "ZN2", "ZN2+",
                "FE", "FE2", "FE3",
                "MN", "MN2",
                "CU", "CU1", "CU2",
                "CO", "NI", "CD",
                "CS", "LI", "RB",
            }

            if res_name_col:
                res_names = (
                    balls[res_name_col]
                    .fillna("")
                    .astype(str)
                    .str.strip()
                    .str.upper()
                )

                water_mask = res_names.isin(water_names)
                ion_mask = res_names.isin(ion_names)
            else:
                water_mask = balls.index.to_series().map(lambda _: False)
                ion_mask = balls.index.to_series().map(lambda _: False)

            molecular_mask = ~(water_mask | ion_mask)

            total_atoms = len(balls)
            molecular_atoms = int(molecular_mask.sum())
            water_atoms = int(water_mask.sum())
            ion_atoms = int(ion_mask.sum())

            # Count residues using chain + sequence + residue name when possible.
            def count_residues(mask):
                if not (res_name_col and res_seq_col):
                    return None

                keys = [res_seq_col, res_name_col]
                if chain_col:
                    keys.insert(0, chain_col)

                return len(
                    balls.loc[mask, keys].drop_duplicates()
                )

            total_residues = count_residues(
                balls.index.to_series().map(lambda _: True)
            )

            nonwater_residues = count_residues(~water_mask)
            molecular_residues = count_residues(molecular_mask)
            water_count = count_residues(water_mask)
            ion_residues = count_residues(ion_mask)

            # Primary chemistry-focused summary
            info.write(f"Molecular Atoms:      {molecular_atoms:,}\n")

            if nonwater_residues is not None:
                info.write(
                    f"Non-Water Residues:    {nonwater_residues:,}\n"
                )

            if water_count is not None:
                info.write(
                    f"Waters:                {water_count:,}\n"
                )

            if water_atoms:
                info.write(
                    f"Water Atoms:            {water_atoms:,}\n"
                )

            if ion_atoms:
                info.write(
                    f"Ions:                   {ion_atoms:,}\n"
                )

            if chain_col:
                chain_count = len(
                    balls.loc[~water_mask, chain_col].drop_duplicates()
                )
                info.write(
                    f"Non-Water Chains:       {chain_count:,}\n"
                )

            info.write("\n")

            # Keep totals available, but subordinate them to the chemically
            # useful composition above.
            info.write("SYSTEM TOTALS\n")
            info.write("-" * 72 + "\n")
            info.write(f"Total Atoms:            {total_atoms:,}\n")

            if total_residues is not None:
                info.write(
                    f"Total Residues:         {total_residues:,}\n"
                )

            info.write("\n")

            # ==========================================================
            # SOLVENT / ION COMPOSITION
            # ==========================================================

            info.write("SOLVENT / ION COMPOSITION\n")
            info.write("-" * 72 + "\n")

            if water_count:
                info.write(
                    f"Water Molecules:        {water_count:,}\n"
                )
                info.write(
                    f"Water Atoms:            {water_atoms:,}\n"
                )

            if ion_atoms:
                ion_counts = (
                    res_names[ion_mask]
                    .value_counts(sort=False)
                )

                for ion_name, count in ion_counts.items():
                    info.write(
                        f"{ion_name:>8}: {int(count):,}\n"
                    )

            if not water_count and not ion_atoms:
                info.write("No recognized solvent or ions.\n")

            info.write("\n")

            # ==========================================================
            # ELEMENT COMPOSITION
            # ==========================================================

            if element_col:

                info.write("ELEMENT COMPOSITION\n")
                info.write("-" * 72 + "\n")

                counts = (
                    balls[element_col]
                    .fillna("")
                    .astype(str)
                    .str.strip()
                )

                counts = counts[
                    counts != ""
                ].value_counts(sort=False)

                for element, count in counts.items():
                    info.write(
                        f"{element:>4}: {int(count):,}\n"
                    )

                info.write("\n")

            # ==========================================================
            # CHAIN COMPOSITION
            # ==========================================================

            if chain_col:

                info.write("CHAIN COMPOSITION\n")
                info.write("-" * 72 + "\n")

                nonwater_balls = balls.loc[~water_mask]

                for chain_value, chain_atoms in nonwater_balls.groupby(
                    chain_col,
                    sort=False,
                    dropna=False,
                ):

                    chain_name = clean(chain_value)

                    line = (
                        f"Chain {chain_name}: "
                        f"{len(chain_atoms):,} atoms"
                    )

                    if res_name_col and res_seq_col:

                        residue_count = len(
                            chain_atoms[
                                [res_seq_col, res_name_col]
                            ].drop_duplicates()
                        )

                        line += (
                            f", {residue_count:,} residues"
                        )

                    info.write(line + "\n")

                info.write("\n")

            # ==========================================================
            # RESIDUE COMPOSITION
            # ==========================================================

            if res_name_col:

                info.write("RESIDUE COMPOSITION\n")
                info.write("-" * 72 + "\n")

                # Exclude bulk water here because it is summarized above.
                nonwater_names = res_names[~water_mask]
                counts = nonwater_names[
                    nonwater_names != ""
                ].value_counts(sort=False)

                for residue, count in counts.items():
                    info.write(
                        f"{residue:>6}: {int(count):,} atoms\n"
                    )

                info.write("\n")

        # ==============================================================
        # GROUPS
        # ==============================================================

        groups = getattr(sys, "groups", None)

        if groups:

            info.write("GROUPS\n")
            info.write("-" * 72 + "\n")

            setting_labels = [
                ("net_type", "Network Type"),
                ("surf_res", "Surface Resolution"),
                ("box_size", "Box Size"),
                ("max_vert", "Maximum Allowable Vertex"),
            ]

            for group in groups:

                group_name = getattr(
                    group,
                    "name",
                    "unnamed",
                )

                info.write(
                    f"Group: {group_name}\n"
                )

                info.write(
                    f"  Atoms:    "
                    f"{object_len(group, 'ball_ndxs', 'atms', 'atoms'):,}\n"
                )

                rsd_count = object_len(
                    group,
                    "rsds",
                    "residues",
                )

                chn_count = object_len(
                    group,
                    "chns",
                    "chains",
                )

                if rsd_count:
                    info.write(
                        f"  Residues: {rsd_count:,}\n"
                    )

                if chn_count:
                    info.write(
                        f"  Chains:   {chn_count:,}\n"
                    )

                group_dir = getattr(
                    group,
                    "dir",
                    None,
                )

                if group_dir:
                    info.write(
                        f"  Directory: {os.path.abspath(group_dir)}\n"
                    )

                net = getattr(
                    group,
                    "net",
                    None,
                )

                settings = getattr(
                    net,
                    "settings",
                    {},
                ) if net is not None else {}

                info.write(
                    "  Build Settings:\n"
                )

                printed_setting = False

                for key, label in setting_labels:

                    value = get_mapping_value(
                        settings,
                        key,
                        None,
                    )

                    if value is None:
                        continue

                    info.write(
                        f"    {label}: {value}\n"
                    )

                    printed_setting = True

                if not printed_setting:
                    info.write(
                        "    N/A\n"
                    )

                info.write("\n")

        # ==============================================================
        # INTERFACES
        # ==============================================================

        interfaces = get_interfaces(sys)

        if interfaces:

            info.write("INTERFACES\n")
            info.write("-" * 72 + "\n")

            for interface in interfaces:

                name = getattr(
                    interface,
                    "name",
                    "unnamed",
                )

                group1 = getattr(
                    getattr(interface, "group1", None),
                    "name",
                    None,
                )

                group2 = getattr(
                    getattr(interface, "group2", None),
                    "name",
                    None,
                )

                info.write(
                    f"Interface: {name}\n"
                )

                if group1 is not None:
                    info.write(
                        f"  Group 1: {group1}\n"
                    )

                if group2 is not None:
                    info.write(
                        f"  Group 2: {group2}\n"
                    )

                interface_dir = getattr(
                    interface,
                    "dir",
                    None,
                )

                if interface_dir:
                    info.write(
                        f"  Directory: "
                        f"{os.path.abspath(interface_dir)}\n"
                    )

                info.write("\n")

        else:

            info.write("INTERFACES\n")
            info.write("-" * 72 + "\n")
            info.write("No interfaces defined.\n\n")

        info.write("=" * 72 + "\n")
