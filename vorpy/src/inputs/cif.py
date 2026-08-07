import numpy as np
from pandas import DataFrame
from vorpy.src.objects import make_atom
from vorpy.src.objects import Sol, Chain
from vorpy.src.objects import Residue
from vorpy.src.chemistry import residue_names
from vorpy.src.chemistry import residue_atoms
from vorpy.src.inputs.fix_sol import fix_sol


def read_cif(sys, file=None):
    """
    Read a CIF/mmCIF molecular structure file into a VorPy System.

    The goal is to produce the same internal System structure as read_pdb()
    when the PDB and CIF files describe the same molecule.

    Unlike PDB files, CIF/mmCIF atom columns are not guaranteed to occur in
    a fixed order. This reader therefore finds the _atom_site field names
    first and uses them to determine the position of each value.

    Parameters
    ----------
    sys : System
        VorPy System object to populate.
    file : str, optional
        CIF file to load. If None, sys.files['base_file'] is used.

    Returns
    -------
    System
        The populated System object.
    """

    from vorpy.src.chemistry import my_masses

    # ----------------------------------------------------------------------
    # Resolve and read the input file
    # ----------------------------------------------------------------------

    if file is None:
        file = sys.files['base_file']

    with open(file, 'r') as rf:
        lines = rf.readlines()

    # Match read_pdb() behavior by resetting molecular collections before loading.
    sys.chains, sys.residues, sys.sol = [], [], None
    previous_res_seq, previous_chain = None, None

    atoms, data = [], []
    chains, resids = {}, {}
    atom_count, reset_checker = 0, 0
    printed_occ_warn = False

    # These names are treated as solvent throughout VorPy/read_pdb().
    solvent_names = {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'}

    # ----------------------------------------------------------------------
    # Find the _atom_site table layout
    # ----------------------------------------------------------------------
    #
    # CIF columns are described by lines such as:
    #
    # _atom_site.group_PDB
    # _atom_site.id
    # _atom_site.type_symbol
    # _atom_site.label_atom_id
    # ...
    #
    # We use these names rather than assuming a fixed column order.
    # ----------------------------------------------------------------------

    atom_site_fields = [line.strip() for line in lines if line.strip().startswith('_atom_site.')]

    if not atom_site_fields:
        raise ValueError(f'No _atom_site fields were found in CIF file: {file}')

    atom_site_index = {field: index for index, field in enumerate(atom_site_fields)}

    # ----------------------------------------------------------------------
    # Helper functions
    # ----------------------------------------------------------------------

    def get_value(values, *field_names, default=None):
        """
        Return the first existing, non-null value from the supplied CIF fields.

        CIF commonly uses "." and "?" to represent missing/not-applicable data,
        so those values are ignored here.
        """
        for field_name in field_names:
            if field_name not in atom_site_index:
                continue

            column = atom_site_index[field_name]

            if column >= len(values):
                continue

            value = values[column]

            if value not in {'.', '?', ''}:
                return value

        return default

    def get_raw_value(values, *field_names, default=None):
        """
        Return the first existing CIF value without treating "." or "?" as null.

        This is useful for fields such as alternate-location identifiers where
        "." specifically means that no alternate conformation is present.
        """
        for field_name in field_names:
            if field_name not in atom_site_index:
                continue

            column = atom_site_index[field_name]

            if column < len(values):
                return values[column]

        return default

    def safe_int(value, default=0):
        """Convert a CIF value to int while safely handling missing values."""
        if value in {None, '.', '?', ''}:
            return default

        try:
            return int(value)
        except ValueError:
            try:
                return int(float(value))
            except ValueError:
                return default

    def safe_float(value, default=None):
        """Convert a CIF value to float while safely handling missing values."""
        if value in {None, '.', '?', ''}:
            return default

        try:
            return float(value)
        except ValueError:
            return default

    # ----------------------------------------------------------------------
    # Parse the file
    # ----------------------------------------------------------------------

    for line in lines:
        stripped = line.strip()

        if not stripped:
            continue

        values = stripped.split()

        if not values:
            continue

        # Keep non-atom CIF information in sys.data.
        if values[0].upper() not in {'ATOM', 'HETATM'}:
            data.append(values)
            continue

        # An atom row should contain one value for every declared _atom_site field.
        if len(values) < len(atom_site_fields):
            raise ValueError(
                f"CIF atom row contains fewer values than its _atom_site header declares. "
                f"Expected {len(atom_site_fields)}, found {len(values)}.\nLine:\n{stripped}"
            )

        # ------------------------------------------------------------------
        # Alternate-location handling
        # ------------------------------------------------------------------

        # "." and "?" mean no alternate conformation. "A" is the primary alternate.
        alt_loc = get_raw_value(values, '_atom_site.label_alt_id', default='.')

        # Ignore alternate B/C/etc. conformations, matching the intended old behavior.
        if alt_loc not in {'A', '.', '?'}:
            if not printed_occ_warn:
                print('Warning! This molecule has multiple occupancy. Program will default to occupancy "A".')
                printed_occ_warn = True

            continue

        # ------------------------------------------------------------------
        # Coordinates
        # ------------------------------------------------------------------

        x = safe_float(get_value(values, '_atom_site.Cartn_x'))
        y = safe_float(get_value(values, '_atom_site.Cartn_y'))
        z = safe_float(get_value(values, '_atom_site.Cartn_z'))

        if x is None or y is None or z is None:
            raise ValueError(f"CIF atom is missing one or more Cartesian coordinates:\n{stripped}")

        loc = np.array([x, y, z])

        # ------------------------------------------------------------------
        # Atom identity
        # ------------------------------------------------------------------

        element = get_value(values, '_atom_site.type_symbol', default='')

        # Prefer author/PDB atom names when available.
        atom_name = get_value(values, '_atom_site.auth_atom_id', '_atom_site.label_atom_id', default=element)

        # Prefer author/PDB residue names when available.
        res_str = get_value(values, '_atom_site.auth_comp_id', '_atom_site.label_comp_id', default='UNK')

        # Prefer author/PDB residue numbering when available.
        res_seq = safe_int(get_value(values, '_atom_site.auth_seq_id', '_atom_site.label_seq_id', default=0), default=0)

        # Prefer auth_asym_id because this corresponds most directly to the PDB chain ID.
        chain_str = get_value(values, '_atom_site.auth_asym_id', '_atom_site.label_asym_id', default=' ')

        if chain_str in {None, '.', '?', ''}:
            chain_str = ' '

        # ------------------------------------------------------------------
        # Optional CIF metadata
        # ------------------------------------------------------------------

        occupancy = safe_float(get_value(values, '_atom_site.occupancy', default=1.0), default=1.0)
        b_factor = safe_float(get_value(values, '_atom_site.B_iso_or_equiv', default=0.0), default=0.0)
        charge = get_value(values, '_atom_site.pdbx_formal_charge', default='')

        pdb_ins_code = get_raw_value(values, '_atom_site.pdbx_PDB_ins_code', default='')
        if pdb_ins_code in {'.', '?'}:
            pdb_ins_code = ''

        cif_atom_id = get_value(values, '_atom_site.id', default=None)
        model_num = safe_int(get_value(values, '_atom_site.pdbx_PDB_model_num', default=1), default=1)

        # ------------------------------------------------------------------
        # Mass
        # ------------------------------------------------------------------

        element_lower = element.lower()

        if sys.type == 'mol' and element_lower in my_masses:
            mass = my_masses[element_lower]
        else:
            mass = 1

        # ------------------------------------------------------------------
        # Create atom
        # ------------------------------------------------------------------
        #
        # Use atom_count instead of _atom_site.id because VorPy expects the
        # same zero-based sequential indexing produced by read_pdb().
        # ------------------------------------------------------------------

        ball = make_atom(location=loc, system=sys, element=element, res_seq=res_seq, res_name=res_str,
                         chn_name=chain_str, name=atom_name, index=atom_count, mass=mass, radius=None,
                         occ_choice=alt_loc, occupancy=occupancy, b_factor=b_factor, charge=charge,
                         pdb_ins_code=pdb_ins_code, pdbx_PDB_model_num=model_num)

        # Keep the original CIF ID as metadata without using it as the VorPy index.
        ball['cif_id'] = cif_atom_id

        # ------------------------------------------------------------------
        # Normalize chain assignment
        # ------------------------------------------------------------------

        if chain_str == ' ':
            if res_str.lower() in solvent_names:
                chain_str = 'SOL'
            else:
                chain_str = 'A'

        elif sys.type == 'foam' and res_str.lower() != 'bub' and chain_str != '0':
            chain_str = 'SOL'
        if previous_chain == chain_str and previous_res_seq is not None and res_seq < previous_res_seq:
            reset_checker += 1

        previous_chain, previous_res_seq = chain_str, res_seq

        # CIF has no meaningful equivalent of PDB's line[17:20], so build the
        # residue identity directly from parsed chain/residue information.
        res_key = f'{chain_str}_{res_str}{res_seq}_{reset_checker}'
        chn_name = chain_str

        # ------------------------------------------------------------------
        # Create/retrieve Chain
        # ------------------------------------------------------------------

        if chn_name in chains:
            my_chn = chains[chn_name]
            my_chn.add_atom(ball['num'])
            ball['chn'] = my_chn

        else:
            if res_str.lower() in solvent_names or chn_name == 'SOL':
                my_chn = Sol(atoms=[ball['num']], residues=[], name=chn_name, sys=sys)
                sys.sol = my_chn
            else:
                my_chn = Chain(atoms=[ball['num']], residues=[], name=chn_name, sys=sys)
                sys.chains.append(my_chn)

            chains[chn_name] = my_chn
            ball['chn'] = my_chn

        # ------------------------------------------------------------------
        # Create/retrieve Residue
        # ------------------------------------------------------------------

        if res_key in resids:
            my_res = resids[res_key]
            my_res.atoms.append(ball['num'])

        else:
            my_res = Residue(sys=sys, atoms=[ball['num']], name=res_str, sequence=ball['res_seq'], chain=ball['chn'])
            resids[res_key] = my_res

            if res_str.lower() in solvent_names or chain_str == 'SOL':
                if sys.sol is None:
                    sys.sol = Sol(sys=sys, atoms=[], residues=[])

                sys.sol.residues.append(my_res)

            else:
                sys.residues.append(my_res)
                ball['chn'].residues.append(my_res)

        ball['res'] = my_res

        # Save the atom and move to the next zero-based VorPy index.
        atoms.append(ball)
        atom_count += 1


    # ----------------------------------------------------------------------
    # Final System cleanup
    # ----------------------------------------------------------------------

    # Ensure a solvent object exists even if the molecule contains no solvent.
    if sys.sol is None:
        sys.sol = Sol(sys=sys, atoms=[], residues=[])

    # Update custom residue chemistry information, matching read_pdb().
    for res in sys.residues:
        if res.name.lower() not in residue_names and res.chain.name != 'SOL':
            residue_names[res.name.lower()] = res.name.upper()
            residue_atoms[res.name.upper()] = {atoms[atom_ndx]['name'] for atom_ndx in res.atoms}

    # Convert the atom list into the normal VorPy balls DataFrame.
    sys.balls, sys.data = DataFrame(atoms), data

    # ----------------------------------------------------------------------
    # Solvent cleanup
    # ----------------------------------------------------------------------
    #
    # Match read_pdb(): solvent residues containing more than three atoms are
    # passed through fix_sol() so combined waters can be split appropriately.
    # ----------------------------------------------------------------------

    adjusted_residues = []

    for res in sys.sol.residues:
        if len(res.atoms) > 3:
            try:
                adjusted_residues += fix_sol(sys, res)
            except TypeError:
                print(res.atoms)
        else:
            adjusted_residues.append(res)

    sys.sol.residues = adjusted_residues

    return sys


if __name__ == '__main__':
    import tkinter as tk
    from tkinter import filedialog
    from vorpy.src.inputs.pdb import read_pdb


    def test_file(pdb, cif, verbose=False):
        """
        Load equivalent PDB and CIF files into separate System objects and
        compare the resulting molecular data.

        The PDB-loaded System is treated as the reference.
        """
        from vorpy.src.system.system import System
        import traceback

        print("\n" + "=" * 70)
        print("PDB / CIF LOADING TEST")
        print("=" * 70)

        # ------------------------------------------------------------------
        # Load PDB reference
        # ------------------------------------------------------------------
        try:
            print(f"\nLoading PDB: {pdb}")
            pdb_sys = System(file=pdb)
            print("\nPDB load: PASS")
        except Exception as exc:
            print("\nPDB load: FAIL")
            print(f"{type(exc).__name__}: {exc}")

            if verbose:
                traceback.print_exc()

            return False

        # ------------------------------------------------------------------
        # Load CIF
        # ------------------------------------------------------------------
        try:
            print(f"\nLoading CIF: {cif}")
            cif_sys = System(file=cif)
            print("\nCIF load: PASS")
        except Exception as exc:
            print("\nCIF load: FAIL")
            print(f"{type(exc).__name__}: {exc}")

            if verbose:
                traceback.print_exc()

            return False

        print("\n" + "-" * 70)
        print("SYSTEM COMPARISON")
        print("-" * 70)

        all_passed = True

        def check(label, pdb_value, cif_value):
            nonlocal all_passed

            passed = pdb_value == cif_value

            print(
                f"{'PASS' if passed else 'FAIL':4} | "
                f"{label:<25} "
                f"PDB={pdb_value!r}  CIF={cif_value!r}"
            )

            if not passed:
                all_passed = False

            return passed

        # ------------------------------------------------------------------
        # Basic system structure
        # ------------------------------------------------------------------
        check("number of atoms", len(pdb_sys.balls), len(cif_sys.balls))
        check("non-SOL residues", len(pdb_sys.residues), len(cif_sys.residues))
        check("chains", len(pdb_sys.chains), len(cif_sys.chains))
        check("SOL residues",
              len(pdb_sys.sol.residues),
              len(cif_sys.sol.residues))

        # ------------------------------------------------------------------
        # Compare SOL residue groupings
        # ------------------------------------------------------------------

        print("\n" + "-" * 70)
        print("SOL RESIDUE COMPARISON")
        print("-" * 70)

        # Compare residues by the actual VorPy atom indices assigned to them.
        # Since the PDB and CIF atom lists already align, matching atom sets
        # should represent the same physical solvent residue.
        pdb_sol_map = {tuple(sorted(res.atoms)): res for res in pdb_sys.sol.residues}
        cif_sol_map = {tuple(sorted(res.atoms)): res for res in cif_sys.sol.residues}

        pdb_sol_sets = set(pdb_sol_map)
        cif_sol_sets = set(cif_sol_map)

        pdb_only = sorted(pdb_sol_sets - cif_sol_sets)
        cif_only = sorted(cif_sol_sets - pdb_sol_sets)

        print(f"Matching SOL residues : {len(pdb_sol_sets & cif_sol_sets)}")
        print(f"PDB-only residues    : {len(pdb_only)}")
        print(f"CIF-only residues    : {len(cif_only)}")

        def print_sol_residue(label, atom_ndxs, sys_obj):
            """Print enough information to identify an unmatched solvent residue."""
            res = pdb_sol_map.get(atom_ndxs) if sys_obj is pdb_sys else cif_sol_map.get(atom_ndxs)

            print(f"\n{label}")
            print(f"  residue name     = {res.name}")
            print(f"  residue sequence = {res.seq}")
            print(f"  chain            = {res.chain.name}")
            print(f"  atom indices     = {list(atom_ndxs)}")

            for atom_ndx in atom_ndxs:
                atom = sys_obj.balls.iloc[atom_ndx]
                print(f"    {atom_ndx}: name={atom['name']}, element={atom['element']}, "
                      f"res={atom['res_name']} {atom['res_seq']}, loc={atom['loc']}")

        if pdb_only:
            print("\nPDB SOL RESIDUES NOT PRESENT AS IDENTICAL RESIDUES IN CIF:")

            for i, atom_ndxs in enumerate(pdb_only):
                print_sol_residue(f"PDB-only residue {i + 1}", atom_ndxs, pdb_sys)

        if cif_only:
            print("\nCIF SOL RESIDUES NOT PRESENT AS IDENTICAL RESIDUES IN PDB:")

            for i, atom_ndxs in enumerate(cif_only):
                print_sol_residue(f"CIF-only residue {i + 1}", atom_ndxs, cif_sys)

        # Also verify that both loaders classify exactly the same atoms as solvent.
        pdb_sol_atoms = set(atom for res in pdb_sys.sol.residues for atom in res.atoms)
        cif_sol_atoms = set(atom for res in cif_sys.sol.residues for atom in res.atoms)

        pdb_only_atoms = sorted(pdb_sol_atoms - cif_sol_atoms)
        cif_only_atoms = sorted(cif_sol_atoms - pdb_sol_atoms)

        print("\nSOL atom coverage:")
        print(f"  PDB SOL atoms     = {len(pdb_sol_atoms)}")
        print(f"  CIF SOL atoms     = {len(cif_sol_atoms)}")
        print(f"  PDB-only atoms    = {pdb_only_atoms}")
        print(f"  CIF-only atoms    = {cif_only_atoms}")

        # Can't safely do atom-by-atom comparison unless the counts match
        if len(pdb_sys.balls) != len(cif_sys.balls):
            print("\nAtom counts differ. Skipping atom-by-atom comparison.")
            return False

        print("\n" + "-" * 70)
        print("ATOM COMPARISON")
        print("-" * 70)

        pdb_balls = pdb_sys.balls
        cif_balls = cif_sys.balls

        # ------------------------------------------------------------------
        # Compare simple atom fields
        # ------------------------------------------------------------------
        compare_columns = [
            "name",
            "element",
            "res_name",
            "res_seq",
        ]

        for column in compare_columns:
            if column not in pdb_balls.columns:
                print(f"SKIP | {column:<25} missing from PDB balls")
                continue

            if column not in cif_balls.columns:
                print(f"FAIL | {column:<25} missing from CIF balls")
                all_passed = False
                continue

            pdb_values = pdb_balls[column].tolist()
            cif_values = cif_balls[column].tolist()

            mismatches = [
                i for i, (p, c) in enumerate(zip(pdb_values, cif_values))
                if p != c
            ]

            if not mismatches:
                print(f"PASS | {column:<25} all {len(pdb_values)} match")
            else:
                all_passed = False
                print(
                    f"FAIL | {column:<25} "
                    f"{len(mismatches)} / {len(pdb_values)} differ"
                )

                if verbose:
                    for i in mismatches[:10]:
                        print(
                            f"       atom {i}: "
                            f"PDB={pdb_values[i]!r}, "
                            f"CIF={cif_values[i]!r}"
                        )

        # ------------------------------------------------------------------
        # Compare coordinates
        # ------------------------------------------------------------------
        if "loc" in pdb_balls.columns and "loc" in cif_balls.columns:
            pdb_locs = np.vstack(pdb_balls["loc"].to_numpy())
            cif_locs = np.vstack(cif_balls["loc"].to_numpy())

            coordinate_diff = np.abs(pdb_locs - cif_locs)
            max_diff = np.max(coordinate_diff)

            if np.allclose(pdb_locs, cif_locs):
                print(
                    f"PASS | {'coordinates':<25} "
                    f"max difference = {max_diff}"
                )
            else:
                all_passed = False
                bad_atoms = np.where(
                    np.any(~np.isclose(pdb_locs, cif_locs), axis=1)
                )[0]

                print(
                    f"FAIL | {'coordinates':<25} "
                    f"{len(bad_atoms)} atoms differ; "
                    f"max difference = {max_diff}"
                )

                if verbose:
                    for i in bad_atoms[:10]:
                        print(
                            f"       atom {i}: "
                            f"PDB={pdb_locs[i]}, "
                            f"CIF={cif_locs[i]}"
                        )

        # ------------------------------------------------------------------
        # Compare assigned chain objects
        # ------------------------------------------------------------------
        if "chn" in pdb_balls.columns and "chn" in cif_balls.columns:
            pdb_chains = [
                getattr(chn, "name", None)
                for chn in pdb_balls["chn"]
            ]

            cif_chains = [
                getattr(chn, "name", None)
                for chn in cif_balls["chn"]
            ]

            mismatches = [
                i for i, (p, c) in enumerate(zip(pdb_chains, cif_chains))
                if p != c
            ]

            if not mismatches:
                print(
                    f"PASS | {'chain assignments':<25} "
                    f"all {len(pdb_chains)} match"
                )
            else:
                all_passed = False
                print(
                    f"FAIL | {'chain assignments':<25} "
                    f"{len(mismatches)} / {len(pdb_chains)} differ"
                )

                if verbose:
                    for i in mismatches[:10]:
                        print(
                            f"       atom {i}: "
                            f"PDB={pdb_chains[i]!r}, "
                            f"CIF={cif_chains[i]!r}"
                        )

        # ------------------------------------------------------------------
        # Final result
        # ------------------------------------------------------------------
        print("\n" + "=" * 70)

        if all_passed:
            print("OVERALL RESULT: PASS")
        else:
            print("OVERALL RESULT: FAIL")

        print("=" * 70)

        return all_passed

    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    pdb_filename = filedialog.askopenfilename(title='Get PDB File')
    cif_filename = filedialog.askopenfilename(title='Get CIF File')
    test_file(pdb=pdb_filename, cif=cif_filename, verbose=True)