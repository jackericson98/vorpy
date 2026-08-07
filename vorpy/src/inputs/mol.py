import re
import numpy as np

from pandas import DataFrame

from vorpy.src.objects.atom import make_atom
from vorpy.src.objects import Sol, Chain, Residue
from vorpy.src.chemistry import my_masses, residue_names, residue_atoms
from vorpy.src.inputs.fix_sol import fix_sol


def read_mol(sys, file=None):
    """
    Read a MOL or SDF file into a VorPy System.

    Supports:
        - MDL V2000 MOL/SDF
        - MDL V3000 MOL/SDF
        - Cartesian coordinates
        - Element symbols
        - Formal charges when available
        - Bond connectivity

    MOL/SDF generally does not contain reliable residue or chain information.
    Therefore this reader loads atoms and bonds but leaves sys.chains and
    sys.residues empty rather than inventing molecular hierarchy.

    Parameters
    ----------
    sys : System
        VorPy System object to populate.
    file : str, optional
        MOL/SDF file. If None, sys.files['base_file'] is used.

    Returns
    -------
    System
        Populated System.
    """

    if file is None:
        file = sys.files['base_file']

    with open(file, 'r') as rf:
        lines = rf.readlines()

    # MOL/SDF does not reliably contain PDB-style chain/residue information.
    sys.chains, sys.residues = [], []
    sys.sol = Sol(sys=sys, atoms=[], residues=[])

    atoms, bonds, data = [], [], []
    atom_id_to_index = {}

    # Determine whether this is V2000 or V3000.
    is_v3000 = any('V3000' in line for line in lines[:10])
    is_v2000 = any('V2000' in line for line in lines[:10])

    if is_v3000:
        atoms, bonds, data = _read_mol_v3000(sys, lines, atom_id_to_index)

    elif is_v2000:
        atoms, bonds, data = _read_mol_v2000(sys, lines, atom_id_to_index)

    else:
        raise ValueError(f'Could not determine MOL/SDF version for file: {file}')

    sys.balls = DataFrame(atoms)
    sys.bonds = bonds
    sys.data = data

    return sys


def _read_mol_v3000(sys, lines, atom_id_to_index):
    """
    Parse the CTAB portion of an MDL V3000 MOL/SDF file.

    V3000 atoms have the general form:

        M  V30 id element x y z atom_map [properties]

    V3000 bonds have the general form:

        M  V30 id type atom1 atom2 [properties]
    """

    atoms, raw_bonds, data = [], [], []
    section = None
    atom_count = 0

    # V3000 permits continued records ending in "-".
    # Join those before attempting to parse the file.
    processed_lines = []
    pending = ''

    for line in lines:
        stripped = line.rstrip('\n')

        if pending:
            stripped = pending + stripped.lstrip()

        if stripped.rstrip().endswith('-'):
            pending = stripped.rstrip()[:-1] + ' '
            continue

        processed_lines.append(stripped)
        pending = ''

    if pending:
        processed_lines.append(pending)

    for line in processed_lines:
        stripped = line.strip()

        if stripped == 'M  V30 BEGIN ATOM':
            section = 'ATOM'
            continue

        if stripped == 'M  V30 END ATOM':
            section = None
            continue

        if stripped == 'M  V30 BEGIN BOND':
            section = 'BOND'
            continue

        if stripped == 'M  V30 END BOND':
            section = None
            continue

        # --------------------------------------------------------------
        # Atom section
        # --------------------------------------------------------------

        if section == 'ATOM' and stripped.startswith('M  V30 '):
            info = stripped[7:].split()

            if len(info) < 6:
                raise ValueError(f'Malformed V3000 atom line:\n{line}')

            atom_id = int(info[0])
            element = info[1]
            location = np.array([float(info[2]), float(info[3]), float(info[4])])

            # Optional V3000 properties occur after the required fields.
            properties = _parse_mol_properties(info[6:])
            charge = properties.get('CHG', '')

            element_lower = element.lower()
            mass = my_masses[element_lower] if sys.type == 'mol' and element_lower in my_masses else 1

            # MOL does not contain a reliable atom name, so use the element.
            ball = make_atom(location=location, system=sys, element=element, name=element, index=atom_count,
                             mass=mass, radius=None, charge=charge)

            ball['mol_id'] = atom_id

            atoms.append(ball)
            atom_id_to_index[atom_id] = atom_count
            atom_count += 1

        # --------------------------------------------------------------
        # Bond section
        # --------------------------------------------------------------

        elif section == 'BOND' and stripped.startswith('M  V30 '):
            info = stripped[7:].split()

            if len(info) < 4:
                raise ValueError(f'Malformed V3000 bond line:\n{line}')

            bond_id = int(info[0])
            bond_type = info[1]
            atom1 = int(info[2])
            atom2 = int(info[3])

            raw_bonds.append((bond_id, atom1, atom2, bond_type))

        else:
            data.append(line)

    # Convert file atom IDs into VorPy's zero-based internal atom indices.
    bonds = []

    for bond_id, atom1, atom2, bond_type in raw_bonds:
        if atom1 not in atom_id_to_index or atom2 not in atom_id_to_index:
            raise ValueError(f'V3000 bond {bond_id} references an atom that was not loaded: {atom1}, {atom2}')

        ndx1, ndx2 = atom_id_to_index[atom1], atom_id_to_index[atom2]
        bonds.append((ndx1, ndx2, bond_type))

    return atoms, bonds, data


def _read_mol_v2000(sys, lines, atom_id_to_index):
    """
    Parse an MDL V2000 MOL/SDF file.

    V2000 uses fixed-width atom and bond records, so these should not be
    interpreted using arbitrary whitespace-field counts.
    """

    atoms, raw_bonds, data = [], [], []

    if len(lines) < 4:
        raise ValueError('V2000 MOL file is too short to contain a valid counts line.')

    counts_line = lines[3]

    try:
        num_atoms = int(counts_line[0:3])
        num_bonds = int(counts_line[3:6])
    except ValueError:
        raise ValueError(f'Could not parse V2000 counts line:\n{counts_line}')

    atom_start = 4
    bond_start = atom_start + num_atoms
    atom_count = 0

    # ------------------------------------------------------------------
    # Atoms
    # ------------------------------------------------------------------

    for line in lines[atom_start:bond_start]:
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            element = line[31:34].strip()
        except ValueError:
            raise ValueError(f'Could not parse V2000 atom line:\n{line}')

        location = np.array([x, y, z])
        element_lower = element.lower()
        mass = my_masses[element_lower] if sys.type == 'mol' and element_lower in my_masses else 1

        ball = make_atom(location=location, system=sys, element=element, name=element, index=atom_count,
                         mass=mass, radius=None)

        # V2000 atom numbering is implicitly one-based by atom-table position.
        mol_id = atom_count + 1
        ball['mol_id'] = mol_id

        atoms.append(ball)
        atom_id_to_index[mol_id] = atom_count
        atom_count += 1

    # ------------------------------------------------------------------
    # Bonds
    # ------------------------------------------------------------------

    for line in lines[bond_start:bond_start + num_bonds]:
        try:
            atom1 = int(line[0:3])
            atom2 = int(line[3:6])
            bond_type = line[6:9].strip()
        except ValueError:
            raise ValueError(f'Could not parse V2000 bond line:\n{line}')

        raw_bonds.append((atom1, atom2, bond_type))

    for atom1, atom2, bond_type in raw_bonds:
        if atom1 not in atom_id_to_index or atom2 not in atom_id_to_index:
            raise ValueError(f'V2000 bond references an atom that was not loaded: {atom1}, {atom2}')

        bonds.append((atom_id_to_index[atom1], atom_id_to_index[atom2], bond_type))

    # Everything after the bond table is retained as additional MOL/SDF data.
    data = lines[bond_start + num_bonds:]

    return atoms, bonds, data


def _parse_mol_properties(values):
    """Convert V3000 property tokens such as CHG=-1 into a dictionary."""

    properties = {}

    for value in values:
        if '=' not in value:
            continue

        key, val = value.split('=', 1)
        properties[key] = val

    return properties


def read_mol2(sys, file=None):
    """
    Read a Tripos MOL2 file into a VorPy System.

    MOL2 contains considerably more structural metadata than MOL/SDF.
    In addition to atoms and bonds, this reader uses the MOL2 SUBSTRUCTURE
    information to reconstruct residue and chain organization when available.

    Atom format:
        atom_id atom_name x y z atom_type subst_id subst_name charge ...

    Relevant SUBSTRUCTURE information:
        subst_id subst_name root_atom subst_type ... chain residue_name

    Parameters
    ----------
    sys : System
        VorPy System object to populate.
    file : str, optional
        MOL2 file. If None, sys.files['base_file'] is used.

    Returns
    -------
    System
        Populated System.
    """

    if file is None:
        file = sys.files['base_file']

    with open(file, 'r') as rf:
        lines = rf.readlines()

    sys.chains, sys.residues, sys.sol = [], [], None

    atoms, bonds, data = [], [], []
    chains, resids = {}, {}
    atom_id_to_index = {}
    substructures = {}

    solvent_names = {'sol', 'hoh', 'sod', 'out', 'cl', 'mg', 'na', 'k', 'ion', 'cla'}

    # ------------------------------------------------------------------
    # First pass: read SUBSTRUCTURE information
    # ------------------------------------------------------------------
    #
    # SUBSTRUCTURE occurs after the atoms in MOL2, but we want residue/chain
    # metadata available while constructing atoms. Reading the file twice is
    # simpler and safer than trying to retroactively reconstruct every atom.
    # ------------------------------------------------------------------

    section = None

    for line in lines:
        stripped = line.strip()

        if stripped.startswith('@<TRIPOS>'):
            section = stripped[9:]
            continue

        if section != 'SUBSTRUCTURE' or not stripped or stripped.startswith('#'):
            continue

        info = stripped.split()

        if len(info) < 2:
            continue

        subst_id = int(info[0])
        subst_name = info[1]
        subst_type = info[3] if len(info) > 3 else ''
        chain_name = info[5] if len(info) > 5 else ''
        residue_name = info[6] if len(info) > 6 else ''

        if chain_name in {'****', '***', '**', '*', '.', '?'}:
            chain_name = ''

        # If no explicit residue name exists, derive it from the substructure name.
        if residue_name in {'', '****', '.', '?'}:
            residue_name = re.sub(r'-?\d+$', '', subst_name)

        # Extract a trailing integer from names such as DA-73, SOL136 or NA3400.
        sequence_match = re.search(r'(-?\d+)$', subst_name)
        residue_sequence = int(sequence_match.group(1)) if sequence_match else subst_id

        substructures[subst_id] = {
            'name': subst_name,
            'type': subst_type,
            'chain': chain_name,
            'res_name': residue_name,
            'res_seq': residue_sequence
        }

    # ------------------------------------------------------------------
    # Second pass: atoms and bonds
    # ------------------------------------------------------------------

    section = None
    raw_bonds = []
    atom_count = 0

    for line in lines:
        stripped = line.strip()

        if stripped.startswith('@<TRIPOS>'):
            section = stripped[9:]
            continue

        if not stripped or stripped.startswith('#'):
            continue

        # --------------------------------------------------------------
        # Atoms
        # --------------------------------------------------------------

        if section == 'ATOM':
            info = stripped.split()

            if len(info) < 6:
                raise ValueError(f'Malformed MOL2 atom line:\n{line}')

            atom_id = int(info[0])
            atom_name = info[1]
            location = np.array([float(info[2]), float(info[3]), float(info[4])])
            atom_type = info[5]

            subst_id = int(info[6]) if len(info) > 6 else 0
            subst_name = info[7] if len(info) > 7 else ''
            charge = float(info[8]) if len(info) > 8 else 0.0

            # Tripos atom types include hybridization after a period:
            # O.3 -> O, C.ar -> C, N.pl3 -> N, etc.
            element = atom_type.split('.')[0]

            # Normalize common element capitalization.
            if len(element) == 1:
                element = element.upper()
            elif len(element) > 1:
                element = element[0].upper() + element[1:].lower()

            element_lower = element.lower()
            mass = my_masses[element_lower] if sys.type == 'mol' and element_lower in my_masses else 1

            # ----------------------------------------------------------
            # Obtain residue/chain information from SUBSTRUCTURE
            # ----------------------------------------------------------

            subst = substructures.get(subst_id, {})

            res_str = subst.get('res_name', '')
            res_seq = subst.get('res_seq', subst_id)
            chain_str = subst.get('chain', '')

            if not res_str:
                res_str = re.sub(r'-?\d+$', '', subst_name) or 'UNK'

            # MOL2 solvent/group records often use **** instead of a chain.
            if res_str.lower() in solvent_names:
                chain_str = 'SOL'
            elif not chain_str:
                chain_str = 'A'

            # ----------------------------------------------------------
            # Create VorPy atom
            # ----------------------------------------------------------

            ball = make_atom(location=location, system=sys, element=element, res_seq=res_seq, res_name=res_str,
                             chn_name=chain_str, name=atom_name, index=atom_count, mass=mass, radius=None, charge=charge)

            ball['mol2_id'] = atom_id
            ball['mol2_type'] = atom_type
            ball['mol2_subst_id'] = subst_id
            ball['mol2_subst_name'] = subst_name

            atom_id_to_index[atom_id] = atom_count

            # ----------------------------------------------------------
            # Chain construction
            # ----------------------------------------------------------

            if chain_str in chains:
                my_chn = chains[chain_str]
                my_chn.add_atom(ball['num'])
                ball['chn'] = my_chn

            else:
                if res_str.lower() in solvent_names or chain_str == 'SOL':
                    my_chn = Sol(atoms=[ball['num']], residues=[], name='SOL', sys=sys)
                    sys.sol = my_chn
                else:
                    my_chn = Chain(atoms=[ball['num']], residues=[], name=chain_str, sys=sys)
                    sys.chains.append(my_chn)

                chains[chain_str] = my_chn
                ball['chn'] = my_chn

            # ----------------------------------------------------------
            # Residue construction
            # ----------------------------------------------------------
            #
            # MOL2 provides a real substructure ID, so use it as the primary
            # residue identity instead of trying to infer boundaries only
            # from residue sequence numbers.
            # ----------------------------------------------------------

            res_key = subst_id

            if res_key in resids:
                my_res = resids[res_key]
                my_res.atoms.append(ball['num'])

            else:
                my_res = Residue(sys=sys, atoms=[ball['num']], name=res_str, sequence=res_seq, chain=ball['chn'])
                resids[res_key] = my_res

                if res_str.lower() in solvent_names or chain_str == 'SOL':
                    if sys.sol is None:
                        sys.sol = Sol(sys=sys, atoms=[], residues=[])

                    sys.sol.residues.append(my_res)

                else:
                    sys.residues.append(my_res)
                    ball['chn'].residues.append(my_res)

            ball['res'] = my_res

            atoms.append(ball)
            atom_count += 1

        # --------------------------------------------------------------
        # Bonds
        # --------------------------------------------------------------

        elif section == 'BOND':
            info = stripped.split()

            if len(info) < 4:
                raise ValueError(f'Malformed MOL2 bond line:\n{line}')

            bond_id = int(info[0])
            atom1 = int(info[1])
            atom2 = int(info[2])
            bond_type = info[3]

            raw_bonds.append((bond_id, atom1, atom2, bond_type))

        elif section not in {'ATOM', 'BOND'}:
            data.append(line)

    # ------------------------------------------------------------------
    # Convert MOL2 atom IDs in bonds to VorPy indices
    # ------------------------------------------------------------------

    for bond_id, atom1, atom2, bond_type in raw_bonds:
        if atom1 not in atom_id_to_index or atom2 not in atom_id_to_index:
            raise ValueError(f'MOL2 bond {bond_id} references an atom that was not loaded: {atom1}, {atom2}')

        bonds.append((atom_id_to_index[atom1], atom_id_to_index[atom2], bond_type))

    # Ensure SOL always exists.
    if sys.sol is None:
        sys.sol = Sol(sys=sys, atoms=[], residues=[])

    # ------------------------------------------------------------------
    # Store atoms before fix_sol(), since fix_sol() may reference sys.balls.
    # ------------------------------------------------------------------

    sys.balls = DataFrame(atoms)
    sys.bonds = bonds
    sys.data = data

    # Update VorPy chemistry dictionaries.
    for res in sys.residues:
        if res.name.lower() not in residue_names and res.chain.name != 'SOL':
            residue_names[res.name.lower()] = res.name.upper()
            residue_atoms[res.name.upper()] = {atoms[atom_ndx]['name'] for atom_ndx in res.atoms}

    # Keep solvent cleanup consistent with PDB/CIF.
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


def read_mol_old(sys, file=None):
    """
    Read and process a MOL format file into a system object.

    This function parses MOL format files, which contain molecular structure data including:
    - Atom coordinates
    - Element information
    - Bond connectivity

    Parameters:
    -----------
    sys : System
        The system object to populate with MOL data
    file : str, optional
        Path to the MOL file. If None, uses sys.files['base_file']

    Returns:
    --------
    None
        Modifies the system object in place by:
        - Creating atom objects from MOL data
        - Storing atoms in a pandas DataFrame
        - Storing bond information
        - Initializing empty lists for chains and residues
    """

    # Check the file variable and if it is none get the systems base file
    if file is None:
        file = sys.files['base_file']

    # Create the dictionary that holds the information from the
    file_dict = {'balls': [], 'Additional Lines': [], 'bonds': []}
    # Open the file
    with open(file, 'r') as rf:
        # Create the index for counting the atoms
        index = 0

        # Loop through the lines
        for line in rf.readlines():

            # Split the line
            line_info = line.split()

            # Check for if it is an atom dood
            if len(line_info) >= 10:

                # Pull the location
                location = np.array([float(_) for _ in line_info[:3]])

                # Create the ball
                ball = make_atom(sys, location=location, element=line[3], index=index)

                # Add the ball
                file_dict['balls'].append(ball)

                # Increment the index
                index += 1

            # If the length of the line is 4 it is the bonds
            elif len(line_info) == 4:

                # Add the bond to the
                file_dict['bonds'].append([int(_) for _ in line_info])

            # Otherwise add the
            else:
                # Add the line to the extra lines list
                file_dict['Additional Lines'].append(line)
    # Return the dataframe
    sys.balls = DataFrame(file_dict['balls'])
    sys.data = file_dict['Additional Lines']
    sys.chains, sys.residues = [], []


if __name__ == '__main__':
    import tkinter as tk
    from tkinter import filedialog
    from vorpy.src.system.system import System

    def test_file(pdb, mol, verbose=False):
        """
        Compare a MOL/SDF-loaded System against an equivalent PDB-loaded System.

        PDB is treated as the reference.

        MOL/SDF does not reliably contain atom names, residues, or chains, so
        this test focuses on information the formats actually share:
            - atom count
            - element identity
            - coordinates

        MOL/SDF-specific bond information is reported separately.
        """

        print("\n" + "=" * 70)
        print("PDB / MOL-SDF LOADING TEST")
        print("=" * 70)

        # --------------------------------------------------------------
        # Load PDB reference
        # --------------------------------------------------------------

        try:
            print(f"\nLoading PDB: {pdb}")
            pdb_sys = System(file=pdb)
            print("\nPDB load: PASS")
        except Exception as exc:
            print("\nPDB load: FAIL")
            print(f"{type(exc).__name__}: {exc}")

            if verbose:
                import traceback
                traceback.print_exc()

            return False

        # --------------------------------------------------------------
        # Load MOL/SDF
        # --------------------------------------------------------------

        try:
            print(f"\nLoading MOL/SDF: {mol}")
            mol_sys = System(file=mol)
            print("\nMOL/SDF load: PASS")
        except Exception as exc:
            print("\nMOL/SDF load: FAIL")
            print(f"{type(exc).__name__}: {exc}")

            if verbose:
                import traceback
                traceback.print_exc()

            return False

        print("\n" + "-" * 70)
        print("SYSTEM COMPARISON")
        print("-" * 70)

        all_passed = True

        def check(label, pdb_value, mol_value):
            nonlocal all_passed

            passed = pdb_value == mol_value
            print(f"{'PASS' if passed else 'FAIL':4} | {label:<25} PDB={pdb_value!r}  MOL/SDF={mol_value!r}")

            if not passed:
                all_passed = False

            return passed

        # --------------------------------------------------------------
        # Basic structure
        # --------------------------------------------------------------

        check("number of atoms", len(pdb_sys.balls), len(mol_sys.balls))

        print(f"INFO | PDB residues              {len(pdb_sys.residues)}")
        print(f"INFO | MOL/SDF residues          {len(mol_sys.residues)}")
        print(f"INFO | PDB chains                {len(pdb_sys.chains)}")
        print(f"INFO | MOL/SDF chains            {len(mol_sys.chains)}")
        print(f"INFO | MOL/SDF bonds             {len(getattr(mol_sys, 'bonds', []))}")

        if len(pdb_sys.balls) != len(mol_sys.balls):
            print("\nAtom counts differ. Skipping atom-by-atom comparison.")
            return False

        print("\n" + "-" * 70)
        print("ATOM COMPARISON")
        print("-" * 70)

        # --------------------------------------------------------------
        # Elements
        # --------------------------------------------------------------

        pdb_elements = pdb_sys.balls['element'].str.upper().tolist()
        mol_elements = mol_sys.balls['element'].str.upper().tolist()

        mismatches = [i for i, (pdb_element, mol_element) in enumerate(zip(pdb_elements, mol_elements))
                      if pdb_element != mol_element]

        if not mismatches:
            print(f"PASS | {'element':<25} all {len(pdb_elements)} match")
        else:
            all_passed = False
            print(f"FAIL | {'element':<25} {len(mismatches)} / {len(pdb_elements)} differ")

            if verbose:
                for i in mismatches[:20]:
                    print(f"       atom {i}: PDB={pdb_elements[i]!r}, MOL/SDF={mol_elements[i]!r}")

        # --------------------------------------------------------------
        # Coordinates
        # --------------------------------------------------------------

        pdb_locs = np.vstack(pdb_sys.balls['loc'].to_numpy())
        mol_locs = np.vstack(mol_sys.balls['loc'].to_numpy())

        coordinate_diff = np.abs(pdb_locs - mol_locs)
        max_diff = np.max(coordinate_diff)

        if np.allclose(pdb_locs, mol_locs):
            print(f"PASS | {'coordinates':<25} max difference = {max_diff}")
        else:
            all_passed = False
            bad_atoms = np.where(np.any(~np.isclose(pdb_locs, mol_locs), axis=1))[0]

            print(f"FAIL | {'coordinates':<25} {len(bad_atoms)} atoms differ; max difference = {max_diff}")

            if verbose:
                for i in bad_atoms[:20]:
                    print(f"       atom {i}: PDB={pdb_locs[i]}, MOL/SDF={mol_locs[i]}")

        # --------------------------------------------------------------
        # Final result
        # --------------------------------------------------------------

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

    pdb_filename = filedialog.askopenfilename(title='Get PDB File', filetypes=[('PDB files', '*.pdb'), ('All files', '*.*')])
    mol_filename = filedialog.askopenfilename(title='Get MOL/SDF File', filetypes=[('MOL/SDF files', '*.mol *.sdf'), ('All files', '*.*')])

    test_file(pdb=pdb_filename, mol=mol_filename, verbose=True)
