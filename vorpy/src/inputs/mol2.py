import re
import numpy as np

from pandas import DataFrame

from vorpy.src.objects.atom import make_atom
from vorpy.src.objects import Sol, Chain, Residue
from vorpy.src.chemistry import my_masses, residue_names, residue_atoms
from vorpy.src.inputs.fix_sol import fix_sol


def read_mol2(sys, file=None):
    """
    Read a Tripos MOL2 file into a VorPy System.

    MOL2 can contain atom names, atom types, bonds, substructures, residues
    and chain information, so this reader reconstructs the molecular hierarchy
    where that information is available.
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

    # ----------------------------------------------------------------------
    # First pass: read MOL2 substructures
    # ----------------------------------------------------------------------

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

        if residue_name in {'', '****', '.', '?'}:
            residue_name = re.sub(r'-?\d+$', '', subst_name)

        sequence_match = re.search(r'(-?\d+)$', subst_name)
        residue_sequence = int(sequence_match.group(1)) if sequence_match else subst_id

        substructures[subst_id] = {'name': subst_name, 'type': subst_type, 'chain': chain_name,
                                   'res_name': residue_name, 'res_seq': residue_sequence}

    # ----------------------------------------------------------------------
    # Second pass: atoms and bonds
    # ----------------------------------------------------------------------

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

            # Tripos types look like C.3, C.ar, N.pl3, O.2, etc.
            element = atom_type.split('.')[0]
            element = element.upper() if len(element) == 1 else element[0].upper() + element[1:].lower()

            element_lower = element.lower()
            mass = my_masses[element_lower] if sys.type == 'mol' and element_lower in my_masses else 1

            # Pull chain/residue information from the SUBSTRUCTURE section.
            subst = substructures.get(subst_id, {})
            res_str = subst.get('res_name', '')
            res_seq = subst.get('res_seq', subst_id)
            chain_str = subst.get('chain', '')

            if not res_str:
                res_str = re.sub(r'-?\d+$', '', subst_name) or 'UNK'

            if res_str.lower() in solvent_names:
                chain_str = 'SOL'
            elif not chain_str:
                chain_str = 'A'

            ball = make_atom(location=location, system=sys, element=element, res_seq=res_seq, res_name=res_str,
                             chn_name=chain_str, name=atom_name, index=atom_count, mass=mass, radius=None, charge=charge)

            ball['mol2_id'] = atom_id
            ball['mol2_type'] = atom_type
            ball['mol2_subst_id'] = subst_id
            ball['mol2_subst_name'] = subst_name

            atom_id_to_index[atom_id] = atom_count

            # --------------------------------------------------------------
            # Chain
            # --------------------------------------------------------------

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

            # --------------------------------------------------------------
            # Residue
            # --------------------------------------------------------------
            #
            # MOL2 has an explicit substructure ID, so use that directly as
            # the residue identity rather than inferring residue boundaries.
            # --------------------------------------------------------------

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

        elif section == 'BOND':
            info = stripped.split()

            if len(info) < 4:
                raise ValueError(f'Malformed MOL2 bond line:\n{line}')

            bond_id, atom1, atom2, bond_type = int(info[0]), int(info[1]), int(info[2]), info[3]
            raw_bonds.append((bond_id, atom1, atom2, bond_type))

        elif section not in {'ATOM', 'BOND'}:
            data.append(line)

    # Convert one-based MOL2 atom IDs to zero-based VorPy indices.
    for bond_id, atom1, atom2, bond_type in raw_bonds:
        if atom1 not in atom_id_to_index or atom2 not in atom_id_to_index:
            raise ValueError(f'MOL2 bond {bond_id} references an unknown atom: {atom1}, {atom2}')

        bonds.append((atom_id_to_index[atom1], atom_id_to_index[atom2], bond_type))

    if sys.sol is None:
        sys.sol = Sol(sys=sys, atoms=[], residues=[])

    sys.balls = DataFrame(atoms)
    sys.bonds = bonds
    sys.data = data

    # Update VorPy residue chemistry dictionaries.
    for res in sys.residues:
        if res.name.lower() not in residue_names and res.chain.name != 'SOL':
            residue_names[res.name.lower()] = res.name.upper()
            residue_atoms[res.name.upper()] = {atoms[atom_ndx]['name'] for atom_ndx in res.atoms}

    # Keep solvent handling consistent with PDB/CIF.
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
    from vorpy.src.system.system import System

    def test_file(pdb, mol2, verbose=False):
        """
        Compare a MOL2-loaded System against an equivalent PDB-loaded System.

        PDB is treated as the reference.

        MOL2 contains enough structural metadata that we can compare:
            - atom count
            - atom names
            - elements
            - coordinates
            - residue names
            - residue sequences
            - chain assignments
            - residue counts
            - chain counts
            - SOL residue counts
        """

        print("\n" + "=" * 70)
        print("PDB / MOL2 LOADING TEST")
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
        # Load MOL2
        # --------------------------------------------------------------

        try:
            print(f"\nLoading MOL2: {mol2}")
            mol2_sys = System(file=mol2)
            print("\nMOL2 load: PASS")
        except Exception as exc:
            print("\nMOL2 load: FAIL")
            print(f"{type(exc).__name__}: {exc}")

            if verbose:
                import traceback
                traceback.print_exc()

            return False

        print("\n" + "-" * 70)
        print("SYSTEM COMPARISON")
        print("-" * 70)

        all_passed = True

        def check(label, pdb_value, mol2_value):
            nonlocal all_passed

            passed = pdb_value == mol2_value
            print(f"{'PASS' if passed else 'FAIL':4} | {label:<25} PDB={pdb_value!r}  MOL2={mol2_value!r}")

            if not passed:
                all_passed = False

            return passed

        check("number of atoms", len(pdb_sys.balls), len(mol2_sys.balls))
        check("non-SOL residues", len(pdb_sys.residues), len(mol2_sys.residues))
        check("chains", len(pdb_sys.chains), len(mol2_sys.chains))
        check("SOL residues", len(pdb_sys.sol.residues), len(mol2_sys.sol.residues))

        if len(pdb_sys.balls) != len(mol2_sys.balls):
            print("\nAtom counts differ. Skipping atom-by-atom comparison.")
            return False

        print("\n" + "-" * 70)
        print("ATOM COMPARISON")
        print("-" * 70)

        # --------------------------------------------------------------
        # Shared atom fields
        # --------------------------------------------------------------

        compare_columns = ['name', 'element', 'res_name', 'res_seq']

        for column in compare_columns:
            if column not in pdb_sys.balls.columns or column not in mol2_sys.balls.columns:
                print(f"SKIP | {column:<25} missing from one loader")
                continue

            pdb_values = pdb_sys.balls[column].tolist()
            mol2_values = mol2_sys.balls[column].tolist()

            # Element capitalization should not create a false mismatch.
            if column == 'element':
                pdb_values = [str(value).upper() for value in pdb_values]
                mol2_values = [str(value).upper() for value in mol2_values]

            mismatches = [i for i, (pdb_value, mol2_value) in enumerate(zip(pdb_values, mol2_values))
                          if pdb_value != mol2_value]

            if not mismatches:
                print(f"PASS | {column:<25} all {len(pdb_values)} match")
            else:
                all_passed = False
                print(f"FAIL | {column:<25} {len(mismatches)} / {len(pdb_values)} differ")

                if verbose:
                    for i in mismatches[:20]:
                        print(f"       atom {i}: PDB={pdb_values[i]!r}, MOL2={mol2_values[i]!r}")

        # --------------------------------------------------------------
        # Coordinates
        # --------------------------------------------------------------

        pdb_locs = np.vstack(pdb_sys.balls['loc'].to_numpy())
        mol2_locs = np.vstack(mol2_sys.balls['loc'].to_numpy())

        coordinate_diff = np.abs(pdb_locs - mol2_locs)
        max_diff = np.max(coordinate_diff)

        if np.allclose(pdb_locs, mol2_locs):
            print(f"PASS | {'coordinates':<25} max difference = {max_diff}")
        else:
            all_passed = False
            bad_atoms = np.where(np.any(~np.isclose(pdb_locs, mol2_locs), axis=1))[0]

            print(f"FAIL | {'coordinates':<25} {len(bad_atoms)} atoms differ; max difference = {max_diff}")

            if verbose:
                for i in bad_atoms[:20]:
                    print(f"       atom {i}: PDB={pdb_locs[i]}, MOL2={mol2_locs[i]}")

        # --------------------------------------------------------------
        # Chain assignments
        # --------------------------------------------------------------

        pdb_chains = [getattr(chain, 'name', None) for chain in pdb_sys.balls['chn']]
        mol2_chains = [getattr(chain, 'name', None) for chain in mol2_sys.balls['chn']]

        mismatches = [i for i, (pdb_chain, mol2_chain) in enumerate(zip(pdb_chains, mol2_chains))
                      if pdb_chain != mol2_chain]

        if not mismatches:
            print(f"PASS | {'chain assignments':<25} all {len(pdb_chains)} match")
        else:
            all_passed = False
            print(f"FAIL | {'chain assignments':<25} {len(mismatches)} / {len(pdb_chains)} differ")

            if verbose:
                for i in mismatches[:20]:
                    print(f"       atom {i}: PDB={pdb_chains[i]!r}, MOL2={mol2_chains[i]!r}")

        # --------------------------------------------------------------
        # SOL residue grouping
        # --------------------------------------------------------------

        pdb_sol_sets = {tuple(sorted(res.atoms)) for res in pdb_sys.sol.residues}
        mol2_sol_sets = {tuple(sorted(res.atoms)) for res in mol2_sys.sol.residues}

        pdb_only = sorted(pdb_sol_sets - mol2_sol_sets)
        mol2_only = sorted(mol2_sol_sets - pdb_sol_sets)

        print("\n" + "-" * 70)
        print("SOL RESIDUE COMPARISON")
        print("-" * 70)

        print(f"Matching SOL residues : {len(pdb_sol_sets & mol2_sol_sets)}")
        print(f"PDB-only residues     : {len(pdb_only)}")
        print(f"MOL2-only residues    : {len(mol2_only)}")

        if pdb_only or mol2_only:
            all_passed = False

        if verbose and pdb_only:
            print("\nFirst PDB-only SOL residues:")
            for atoms in pdb_only[:10]:
                print(f"  {list(atoms)}")

        if verbose and mol2_only:
            print("\nFirst MOL2-only SOL residues:")
            for atoms in mol2_only[:10]:
                print(f"  {list(atoms)}")

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
    mol2_filename = filedialog.askopenfilename(title='Get MOL2 File', filetypes=[('MOL2 files', '*.mol2'), ('All files', '*.*')])

    test_file(pdb=pdb_filename, mol2=mol2_filename, verbose=True)
