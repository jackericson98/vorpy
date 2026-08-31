# Solvent

VorPy treats solvent as part of the molecular geometry rather than as a feature added after the Voronoi calculation.

Water molecules and ions remain part of the full System and can therefore influence cell boundaries, shared surfaces, contacts, molecular exposure, and interfaces even when the primary Group contains only solute atoms.

This page describes the current solvent representation, PDB solvent identification, water normalization, solvent-facing surface area, and interface-water analysis.

---

# Conceptual Model

For molecular systems, VorPy separates the primary molecular residue/chain collections from a dedicated solvent/environment collection.

Conceptually:

```text
System
│
├── balls
│   └── all atoms/balls in the loaded structure
│
├── residues
│   └── ordinary non-solvent residues
│
├── chains
│   └── ordinary non-solvent chains
│
└── sol
    ├── atoms
    └── residues
```

`System.balls` remains the complete atomic table.

The separation occurs at the molecular-object level: solvent/environment residues are stored under `System.sol` instead of being inserted into the ordinary `System.residues` collection.

This allows a Group to contain only the molecular atoms of interest while the surrounding solvent remains available to define the geometry of those atoms.

---

## Image Placeholder 1 — Solute and Solvent Representation

![Placeholder: solute and solvent organization](../../assets/docs/solvent/system_solvent_structure.png)

> **Suggested figure:** Show a protein or ligand surrounded by water and ions. On the right, show the VorPy object hierarchy: `System.balls` contains everything, ordinary solute residues/chains are stored separately, and waters/ions are collected in `System.sol`.

---

# The `Sol` Object

`Sol` is implemented as a subclass of `Chain`.

It stores:

- atom indices,
- solvent/environment residue objects,
- a name,
- volume and surface-area attributes.

The class currently defaults its name to `H2O`, and its docstring describes it as a solvent molecule.

In actual PDB parsing, however, `System.sol` is used more broadly as a **system-level solvent/environment container**. It can contain many water residues and recognized ion residues rather than representing one individual water molecule.

For scientific interpretation, the current object is therefore better understood as:

```text
System.sol
    =
recognized solvent/environment collection
```

rather than:

```text
one solvent molecule
```

Individual waters are represented by `Residue` objects inside `System.sol.residues`.

---

# PDB Solvent Identification

During PDB loading, VorPy reads both `ATOM` and `HETATM` records.

Solvent classification is not based solely on whether an atom came from `HETATM`. Instead, the current PDB parser checks residue names and chain assignment.

The parser currently recognizes the following residue names as part of the solvent/environment collection:

```text
SOL
HOH
SOD
OUT
CL
MG
NA
K
ION
CLA
```

If one of these residues has no explicit PDB chain identifier, VorPy assigns it to a synthetic chain named:

```text
SOL
```

An atom already associated with chain `SOL` is also treated as part of the solvent collection.

---

## Water and Ions Share the Solvent Container

The current PDB parser places both water-like residues and several ions in `System.sol`.

For example:

```text
SOL
HOH
```

and:

```text
NA
CL
MG
K
SOD
CLA
```

can all be represented within the same `Sol` object.

This is an implementation-level organization. It does **not** mean that VorPy treats an ion as a water molecule during later water-specific analyses.

Water-specific calculations should use a water-residue classification rather than simply assuming every residue in `System.sol` is water.

---

# Residue Organization

When a recognized solvent/environment residue is encountered, its `Residue` object is appended to:

```text
System.sol.residues
```

rather than:

```text
System.residues
```

Ordinary molecular residues are instead added to the normal residue collection and to their molecular chain.

Thus:

```text
ordinary residue
      ↓
System.residues
      ↓
ordinary Chain

recognized solvent/environment residue
      ↓
System.sol.residues
      ↓
Sol
```

This separation is why VorPy selection and analysis code may inspect both the ordinary residue collection and the solvent collection.

---

# Systems Without Solvent

If a PDB contains no recognized solvent, VorPy still initializes an empty `Sol` object.

Conceptually:

```text
System.sol ≠ None
System.sol.residues = []
System.sol.atoms = []
```

This provides a consistent downstream object interface while allowing solvent-free systems to remain valid.

---

# Normalizing Water Residues

Some PDB files do not present each water molecule as a clean three-atom residue.

For example, a parser may encounter one large residue containing atoms from many waters.

After PDB parsing, VorPy inspects every residue in `System.sol.residues`.

If a solvent residue contains more than three atoms, it is passed to:

```text
fix_sol(...)
```

Residues containing three or fewer atoms are retained directly.

---

## `fix_sol()` Water Reconstruction

`fix_sol()` attempts to reorganize a large solvent residue into individual water residues.

The current algorithm is geometric.

### Step 1 — Separate oxygen and hydrogen atoms

Each oxygen creates a new candidate `Residue` containing that oxygen.

Hydrogen atoms are collected separately.

Other element types are not used by the water-reconstruction step.

### Step 2 — Assign hydrogens by distance

For each hydrogen, VorPy searches the candidate oxygen residues and identifies the nearest oxygen.

The hydrogen is assigned to that oxygen when:

\[
d_{\mathrm{O-H}} < 1.5\ \text{Å}
\]

The current implementation therefore uses a **1.5 Å geometric cutoff**, rather than bond records, to reconstruct waters.

### Step 3 — Identify complete waters

A candidate oxygen residue containing exactly three atoms is considered complete:

```text
1 oxygen + 2 assigned atoms
```

Its atoms are then remapped to the newly created residue.

### Step 4 — Attempt incomplete-water recovery

Candidate residues containing fewer than three atoms receive another attempt to collect nearby remaining hydrogens within the same 1.5 Å cutoff.

The resulting residue collection is returned to the PDB reader and becomes the normalized `System.sol.residues` list.

---

## Image Placeholder 2 — Water Reconstruction

![Placeholder: distance-based water reconstruction](../../assets/docs/solvent/fix_sol_water_reconstruction.png)

> **Suggested figure:** Start with a single oversized `SOL` residue containing several oxygen and hydrogen atoms. Draw a 1.5 Å neighborhood around each oxygen, then show the atoms divided into individual three-atom water residues.

---

# Important `fix_sol()` Limitations

The current water-normalization code should be interpreted as a pragmatic repair step rather than a chemical bond-perception algorithm.

It assumes:

- water oxygens can be identified from the atom element,
- water hydrogens can be identified from the atom element,
- the nearest appropriate O-H pair lies within 1.5 Å,
- a normal explicit-water representation should contain approximately one oxygen and two hydrogens.

The current implementation also contains edge-case behavior that is still under validation. In particular, the hydrogen list is modified while it is being iterated, and leftover/incomplete hydrogen handling contains paths that can associate the same remaining hydrogen list with multiple generated residues.

For well-formed modern PDB/GROMACS waters, `fix_sol()` may never be needed because individual three-atom solvent residues are already present.

---

# Solvent in Group Geometry

A VorPy Group does not have to contain solvent atoms for solvent to influence its geometry.

When a Group Network is constructed, the complete System coordinates and radii remain available while the Group's ball indices determine which cells belong to the Group.

Conceptually:

```text
full system
├── selected Group atoms
└── surrounding atoms, including solvent
          ↓
      define geometry
          ↓
Group-owned Voronoi cells
```

Therefore a protein-only Group can still have boundaries defined against explicit water molecules.

This distinction is important:

> **Group membership** determines which cells are analyzed as part of the Group; it does not imply that surrounding solvent has been removed from the geometric problem.

---

# Solvent-Facing Surface Area

VorPy can classify Group boundary surfaces according to whether the opposing atom is recognized as water.

For a Group atom \(i\), a shared surface with neighboring atom \(j\) contributes to solvent-interfacial surface area when:

```text
i belongs to the analyzed Group
and
j is recognized as water
```

The same logic is applied symmetrically when the Group atom is the second defining atom.

At the residue level:

\[
A_{\mathrm{solv}}(R)
=
\sum_{\substack{s\in\partial R\\
\text{opposing atom is water}}}
A_s
\]

and similarly for chains.

This quantity is a subset of the object's total boundary surface area.

---

## Water Names Used by Surface-Area Summaries

The current Group information/export logic recognizes these water residue names for solvent-interfacial surface area:

```text
SOL
HOH
WAT
H2O
TIP3
TIP3P
TIP4
TIP4P
SPC
SPCE
```

Ions are not included in this water-specific surface-area classification.

Therefore:

```text
Solvent-Interfacial SA
```

in the current Group summary is more precisely a **water-interfacial surface area** under this residue-name definition.

---

## Residue-Level Interpretation

For a residue, VorPy can distinguish:

```text
Total Boundary SA
Inter-Residue SA
Solvent-Interfacial SA
```

The solvent-interfacial component is the area of residue boundary surfaces whose opposing atom is a recognized water atom.

Internal surfaces between atoms belonging to the same residue do not define the residue boundary.

---

## Chain-Level Interpretation

The analogous chain quantities distinguish:

```text
Total Boundary SA
Inter-Chain SA
Solvent-Interfacial SA
```

A surface internal to the same selected chain is excluded from the chain boundary.

A surface between a selected chain atom and recognized water contributes to the chain's solvent-interfacial area.

---

## Solvent Exposure

A useful normalized geometric exposure measure is:

\[
f_{\mathrm{solv}}
=
\frac{A_{\mathrm{solv}}}{A_{\mathrm{boundary}}}
\]

or as a percentage:

\[
\%_{\mathrm{solv}}
=
100
\frac{A_{\mathrm{solv}}}{A_{\mathrm{boundary}}}.
\]

This measures what fraction of an atom, residue, chain, or other molecular boundary is geometrically adjacent to recognized water in the VorPy partition.

It should be described as a **Voronoi surface-based solvent exposure measure**, not automatically equated with conventional SASA.

---

## Image Placeholder 3 — Solvent-Interfacial Surface Area

![Placeholder: solvent-facing Voronoi surfaces](../../assets/docs/solvent/solvent_interfacial_surface.png)

> **Suggested figure:** Show a residue with several shared Voronoi surfaces. Color surfaces against the same residue separately from surfaces against another biomolecular residue and surfaces against water. Label the summed water-facing surfaces as **Solvent-Interfacial SA**.

---

# Solvent-Interfacial Area Is Not SASA

VorPy's solvent-interfacial surface area is defined from the spatial-partition network.

It is based on **shared Voronoi surfaces between molecular atoms and explicit water atoms**.

This differs conceptually from a probe-based solvent-accessible surface area calculation.

A VorPy water-facing surface answers:

```text
Which part of this cell boundary is shared with explicit solvent?
```

rather than:

```text
Where could the center of an idealized solvent probe travel?
```

The two quantities can be scientifically related, but they are not interchangeable definitions.

---

# System Solvent and Ion Summaries

The current System information workflow independently classifies water and ions from the atomic table.

Its recognized water names include:

```text
SOL
HOH
WAT
H2O
TIP3
TIP3P
TIP4
TIP4P
SPC
SPCE
```

The ion list is broader and includes common names for sodium, potassium, chloride, magnesium, calcium, zinc, iron, manganese, copper, and several other ions.

The System information file can report:

```text
Molecular Atoms
Non-Water Residues
Waters
Water Atoms
Ions
Non-Water Chains

SYSTEM TOTALS
Total Atoms
Total Residues

SOLVENT / ION COMPOSITION
Water Molecules
Water Atoms
individual ion counts
```

This composition summary is derived from the loaded atomic table rather than from solved Voronoi geometry.

---

# Current Classification Consistency

The project currently contains more than one solvent/water classification list.

For example:

- the PDB parser uses a relatively small hard-coded solvent/environment set,
- the System information exporter recognizes a broader set of water and ion names,
- the Group solvent-interfacial surface classifier uses its own water-name set.

These definitions overlap but are not currently identical.

That means an unusual water or ion residue could potentially be:

```text
recognized by one subsystem
but not another
```

until the classification tables are centralized.

For reproducible analyses, report the residue naming convention used in the input structure, particularly when using nonstandard water or ion names.

---

# Solvent Selections

The CLI selection system can create Groups from solvent-related categories including:

```text
water
sol
ions
```

These selections use VorPy's molecular classification logic rather than requiring the user to enumerate atom indices manually.

This allows solvent itself to become the analyzed Group when desired.

However, the ordinary default molecular Group behavior generally focuses on non-solvent molecular residues unless the user explicitly requests solvent.

See [Selections](selections.md) for the selection grammar and [Groups](groups.md) for how selected residues are converted into Group-owned cells.

---

# Solvent Around Interfaces

Interface analysis has an additional solvent workflow.

This workflow is distinct from general `System.sol` organization and from global solvent-interfacial surface area.

VorPy can identify **interface-associated waters** from the retained Interface Network.

A water residue is considered interface-associated when at least one of its atoms appears among the defining balls of a retained Interface Network surface.

Conceptually:

```text
retained Interface Network surfaces
             ↓
collect defining ball indices
             ↓
compare with water-residue atoms
             ↓
water has at least one matching atom
             ↓
interface-associated water
```

The current definition uses the complete retained Interface Network surface set, not only direct Group-1/Group-2 surfaces.

---

# Interface-Local Water Geometry

For each identified interface water, VorPy gathers the retained Interface Network geometry associated with that water.

This includes:

- water atom indices,
- associated retained surfaces,
- associated edges,
- associated vertices,
- summed retained surface area,
- summed Contact Area.

This describes the portion of the **Interface Network** involving that water.

It is not necessarily the complete Voronoi geometry of the water molecule.

---

# Complete Water Cells

To obtain the complete cell geometry of an interface-associated water, VorPy can build a normal Group containing that water residue.

Conceptually:

```text
interface-associated water
          ↓
Group(residues=[water])
          ↓
normal VorPy Network
          ↓
complete water cell(s)
```

This creates a second, distinct geometric description:

1. **Interface-local water geometry** — geometry retained within the Interface Network.
2. **Complete water-cell geometry** — the full cell geometry obtained by solving the water as its own Group in the complete System environment.

The distinction should be preserved in scientific interpretation.

---

## Image Placeholder 4 — Interface Water

![Placeholder: interface-local water versus complete water cell](../../assets/docs/solvent/interface_water_geometry.png)

> **Suggested figure:** Place one water between two molecular Groups. Panel A shows only Interface Network surfaces/edges involving the water. Panel B shows the same water rebuilt independently with its complete Voronoi cell geometry.

---

# Current Interface-Water Status

The interface-water workflow is still under validation.

The current Interface build can identify multiple associated waters, but complete water-Group construction is temporarily limited to the first detected water during the current test implementation.

The water-cell volume summary also still requires final validation of its authoritative volume source.

These are implementation limitations, not changes to the underlying distinction between interface-local and complete water-cell geometry.

See [Interfaces](interfaces.md) for the full interface workflow.

---

# Recommended Reporting

For a solvent-related VorPy analysis, report:

- input structure and solvent model,
- water residue names used in the structure,
- ion residue names where relevant,
- VorPy network type,
- atomic radius model if modified,
- surface resolution,
- whether solvent exposure refers to absolute water-interfacial area or normalized exposure,
- whether a reported water is based on general solvent adjacency or the Interface Network definition,
- whether water geometry is interface-local or a complete independently rebuilt water Group,
- VorPy version.

When comparing different structures or simulation packages, confirm that solvent and ion residue naming is interpreted consistently.

---

# Summary

VorPy retains explicit solvent as part of the spatial-partition geometry.

At the System level:

```text
System.balls
```

contains the complete atomic system, while recognized solvent/environment residues are organized under:

```text
System.sol
```

rather than the ordinary residue and chain collections.

The current PDB reader recognizes a hard-coded set of water and ion residue names and can use a geometric `fix_sol()` procedure to split oversized water residues into candidate individual waters.

At the geometric-analysis level, VorPy can measure surfaces shared between selected molecular objects and recognized water atoms, providing a Voronoi surface-based solvent exposure descriptor.

At molecular interfaces, VorPy additionally identifies interface-associated waters and can distinguish their Interface Network geometry from independently rebuilt complete water-cell geometry.

The primary development need is to centralize solvent/water/ion classification and complete validation of the water-reconstruction and interface-water workflows.
