# Groups

Groups define the atoms whose Voronoi cells VorPy builds, analyzes, summarizes, and exports as a single molecular unit.

A Group is therefore more than a named atom selection. It connects a molecular selection from the parent System to a dedicated VorPy Network and provides the context for group-level geometry, surrounding layers, exports, and interfaces.

---

## What Is a Group?

Every VorPy Group belongs to a parent System.

The System contains the complete molecular structure and its molecular classifications. A Group selects a subset of that System and converts the selected atoms into a list of system-level ball indices. Those indices define the cells owned by the Group's Network.

Conceptually:

```text
System
│
├── all atoms
├── residues
├── chains
├── molecular classifications
│
└── Group selection
      │
      ├── selected atoms
      └── Group Network
            ├── cells
            ├── surfaces
            ├── edges
            ├── vertices
            ├── geometry
            └── surrounding layers
```

The Group keeps references to the parent System, its selected molecular objects, its selected ball indices, its Network, build settings, analysis results, layer information, exports, and interface relationships.

---

## Image Placeholder 1 — System, Group, and Network

![Placeholder: System to Group to Network](../../assets/docs/groups/system_group_network.png)

> **Suggested figure:** Show a complete molecular System on the left, highlight a subset of atoms in the center as the Group, and show the corresponding solved Voronoi cells on the right. Use the same highlighted atoms through all three stages. A small callout should distinguish **selected/group atoms** from **surrounding atoms needed to define their cell boundaries**.

---

## How Groups Are Defined

VorPy can currently construct Groups from several kinds of molecular selections.

### Direct structural selections

The Group object accepts selections based on:

- atoms,
- residues,
- chains,
- molecules.

Atom selections ultimately become system ball indices. Residue and chain selections are expanded to the atoms belonging to those molecular objects.

### Broad molecular categories

The command-selection layer also supports broad molecular categories such as:

```text
protein
dna
rna
ligand
water
sol
ions
others
```

These category selectors are resolved through VorPy's molecular classification system and converted into matching residues before Group construction.

### Specific residue identity

Residues can be selected by molecular identity using a residue name and sequence number, for example conceptually:

```text
THY 124
```

Specific atoms within a residue can also be selected when an atom name is supplied.

### Index-based selections

Atoms, residues, and chains can also be selected using internal indices or index ranges.

Because molecular file numbering and VorPy's internal indices are not necessarily the same concept, documentation and command examples should distinguish **internal index selection** from **structural identity selection**.

> **Note:** Exact command-line spellings and aliases are documented separately in the CLI selection reference. This page describes the Group model rather than serving as the complete command grammar.

---

## Default Group Behavior

When no explicit grouping command is supplied for a molecular system, VorPy creates one standard Group from the System's non-solvent residue collection.

This default is useful for conventional molecular analyses because the molecular solute can be analyzed as the Group while solvent atoms remain available in the parent System and network geometry.

VorPy also supports an explicit full-system grouping path in which all system balls are selected.

These two cases should not be confused:

```text
Default molecular Group
    selected non-solvent residues

Full-system Group
    selected all system balls
```

The distinction becomes especially important for solvent exposure and interface analysis.

---

## Selected Atoms vs. Surrounding Atoms

A Group owns the cells corresponding to its selected `ball_ndxs`, but the Network retains the geometry of the complete parent System.

This is an important feature of molecular Voronoi analysis.

A cell belonging to a selected atom may be bounded by atoms that are not themselves members of the Group. Those surrounding atoms are still necessary to determine the selected atom's Voronoi boundary correctly.

Therefore:

> **Group membership determines which cells are analyzed as belonging to the Group; it does not imply that surrounding System atoms are removed from the geometric problem.**

This distinction is particularly important for molecular solutes surrounded by solvent.

---

## Group Networks

Each Group can own a dedicated Network.

The Network receives:

- the complete System ball locations,
- the complete System ball radii,
- the complete System ball masses,
- the Group's selected ball indices,
- the Group's build settings,
- a reference to the parent System.

The selected ball indices determine the cells assigned to the Group, while complete-system geometry remains available during network construction.

A Group may be initialized without immediately solving the Network, allowing Group definition and network construction to remain separate operations.

---

## Group Build Settings

Groups store the settings used to construct and analyze their Network.

Current defaults include:

| Setting | Default | Purpose |
|---|---:|---|
| Surface resolution | `0.2` | Controls surface triangulation resolution |
| Box size | `1.25` | Spatial-search/build setting |
| Maximum allowable vertex | `40` | Vertex-search/build limit |
| Network type | `aw` | Spatial partitioning scheme |
| Build type | `all` | Network construction scope |

Additional visualization, splitting, rounding, and output settings are also stored with the Group.

Because Groups may be analyzed independently, different Groups can in principle carry different network settings. When comparing Groups quantitatively, the relevant settings should therefore be reported.

---

## Group-Level Analysis

After the Network is built and analyzed, VorPy can summarize geometry across the Group.

Current Group-level quantities include:

- Voronoi volume,
- van der Waals volume,
- surface area,
- density,
- mass,
- center-related quantities,
- moment-of-inertia information,
- integrated mean curvature,
- integrated mean-curvature squared,
- integrated Gaussian curvature,
- area-weighted mean curvature,
- area-weighted Gaussian curvature.

Volume-based Group quantities are accumulated from complete cells.

The Group surface area and integrated curvature quantities are calculated over the Group's **first exposed surface layer**, rather than by summing every surface attached to every atom. This prevents internal Group boundaries from being treated as part of the exposed Group boundary.

---

## Group Boundary

VorPy identifies the Group boundary from surfaces that cross out of the current atom layer.

For a surface shared by two atoms:

```text
both atoms in current layer
    -> internal surface
    -> not part of that layer's outer boundary

one atom in current layer
one atom outside
    -> boundary surface
```

The first such surface layer defines the exposed Group boundary used for Group surface-area and curvature summaries.

This makes the Group boundary a direct consequence of the solved Voronoi adjacency network.

---

## Image Placeholder 2 — Internal vs. Boundary Surfaces

![Placeholder: Internal and boundary group surfaces](../../assets/docs/groups/internal_vs_boundary_surfaces.png)

> **Suggested figure:** Highlight a small cluster of selected atoms. Draw their Voronoi cells and distinguish internal Group–Group surfaces from surfaces crossing from a Group atom to a surrounding atom. Label the latter **Group boundary / layer 1 surfaces**. This figure can also reinforce why simply summing every atom surface would double-count or include internal boundaries.

---

## Surrounding Layers

Groups can be expanded into successive topological layers using the solved Voronoi network.

Layer 0 is the Group body:

```text
Layer 0
    selected Group atoms
```

The next atom layer is found by following Voronoi surfaces that cross outward from the current layer. The atom on the opposite side of each boundary surface becomes a member of the next layer.

This process can be repeated outward:

```text
Group / Layer 0
      ↓
Layer 1 neighbors
      ↓
Layer 2 neighbors
      ↓
Layer 3 neighbors
      ↓
...
```

VorPy stores layer-specific atoms, surfaces, edges, vertices, and summary information.

### Keeping residues together

By default, surrounding-layer construction uses `group_resids=True`.

When an atom is added to a surrounding layer and belongs to a residue, the remaining atoms from that residue are also added to that atom layer.

This allows molecular layers to preserve residue identity rather than fragmenting residues solely according to atom-level Voronoi adjacency.

---

## Image Placeholder 3 — Topological Group Layers

![Placeholder: Voronoi neighbor layers around a group](../../assets/docs/groups/group_layers.png)

> **Suggested figure:** Show a central molecular Group surrounded by two or three shells of atoms. Label **Layer 0: Group**, **Layer 1: direct Voronoi neighbors**, and **Layer 2: next topological layer**. If possible, show one residue whose atoms are kept together when residue grouping is enabled.

---

## Groups and Interfaces

Groups provide the molecular selections used to define interfaces.

A complete Group may participate in multiple interfaces. VorPy stores interface metadata separately from the Group's primary Network so that the same Group can be related to multiple opposing selections.

Conceptually:

```text
Group A ───── Interface A:B ───── Group B

Group A ───── Interface A:C ───── Group C
```

An interface therefore does not replace either Group. It describes a geometric relationship between Group selections.

See [Interfaces](interfaces.md) for the interface-specific workflow and calculations.

---

## Groups and Contacts

Contacts are defined at the solved-cell level, while Groups determine how those contacts are interpreted at larger molecular scales.

For example, a cell boundary can represent:

```text
Group atom ↔ Group atom
    internal Group relationship

Group atom ↔ surrounding atom
    Group boundary

Group A atom ↔ Group B atom
    potential contribution to an interface
```

This is why Group definitions are foundational for contact, interface, and solvent analyses.

See [Contacts](contacts.md) for the contact definitions.

---

## Group Outputs

Group information is available through several output types.

### `info.txt`

The Group information file summarizes the Group and its solved Network, including:

- molecular composition,
- build settings,
- timing,
- network topology,
- volume and surface-area quantities,
- curvature quantities,
- representative surface energy,
- molecular classification summaries.

### Network logs

Network logs contain atom-, surface-, edge-, and vertex-level information associated with the Group Network.

### Geometry and visualization exports

Groups can export combinations of:

- atoms,
- atom surfaces,
- atom edges,
- atom vertices,
- complete surfaces,
- separate surfaces,
- shell surfaces,
- edges,
- shell edges,
- vertices,
- shell vertices,
- surrounding atoms,
- layer outputs,
- logs,
- information summaries.

The export reference documents the exact presets and command syntax.

---

## Example: EDTA

For the default EDTA example, VorPy creates an additively weighted Group containing 32 selected molecular atoms.

The resulting Group Network contains:

```text
32 selected atoms
500 vertices
892 edges
424 surfaces
```

The Group information file then summarizes its geometry, including volume, exposed surface area, curvature quantities, and representative surface energy.

The same System also contains surrounding water and other non-Group atoms. These atoms provide environmental geometry without becoming part of the 32-atom Group itself.

This provides a useful example of the distinction between:

```text
System
    complete molecular environment

Group
    selected EDTA atoms

Group Network
    Voronoi cells owned by those selected atoms,
    solved in the context of the complete System
```

---

## Interpretation Guidelines

When reporting or comparing Group results, record at minimum:

- input structure,
- Group definition,
- network type,
- atomic radius model where relevant,
- surface resolution,
- other non-default build settings,
- VorPy version.

For comparisons between Groups, use equivalent selection logic and compatible network settings unless the settings themselves are the variable being studied.

Group definitions should also be described in molecular terms whenever possible. For example, reporting a residue, chain, protein, DNA, or ligand selection is usually more reproducible than reporting only an undocumented list of internal indices.

---

## Summary

A VorPy Group connects a molecular selection to a solved spatial-partitioning Network.

The Group determines **which cells belong to the analysis**, while the parent System continues to provide the surrounding molecular geometry required to construct those cells correctly.

Groups provide the basis for:

- group-level volume and exposed surface geometry,
- curvature summaries,
- topological surrounding layers,
- molecular contacts,
- solvent relationships,
- interfaces,
- exports and visualization.

They are therefore the central organizational unit between the complete molecular System and VorPy's cell-level geometric analysis.
