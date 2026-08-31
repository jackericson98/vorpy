# Summaries

VorPy produces information at several levels of the molecular and geometric hierarchy.

A **summary** is not a replacement for the detailed network data. Instead, summaries provide progressively higher-level interpretations of the same calculation:

```text
System
  ↓
Group / Interface
  ↓
Residue / Chain
  ↓
Atom / Cell
  ↓
Surface / Edge / Vertex
```

The appropriate output depends on the scientific question.

Use the System information file to understand the molecular model and configured analyses, the Group or Interface information file to understand a solved geometric object, and the network logs when atom-, surface-, edge-, or vertex-level data are required.

---

# Summary Hierarchy

VorPy separates human-readable summaries from detailed machine-readable network records.

```text
Input molecular structure
        │
        ▼
      System
        │
        ├── System info
        │     molecular composition
        │     solvent / ions
        │     groups
        │     interfaces
        │
        ├── Group
        │     ├── Group info
        │     │     geometry
        │     │     build settings
        │     │     timing
        │     │     residue / chain summaries
        │     │
        │     └── Network logs
        │           atoms
        │           surfaces
        │           edges
        │           vertices
        │
        └── Interface
              ├── Interface info
              └── Interface-compatible logs
```

The distinction is important because the same numerical quantity may appear at more than one level for different purposes.

---

## Image Placeholder 1 — VorPy Output Hierarchy

![Placeholder: VorPy summary hierarchy](../../assets/docs/summaries/summary_hierarchy.png)

> **Suggested figure:** Show the System at the top, branching to Groups and Interfaces. Under a Group, show `info.txt` and the network log. Under the log, show Atom, Surface, Edge, and Vertex tables. Add Residue and Chain summaries as aggregated interpretations of the Group Network.

---

# System Summary

The System information file is the top-level manifest for a molecular calculation.

It is intended to describe the **input molecular system and configured VorPy analyses**, rather than to duplicate every solved geometric quantity from every Group or Interface.

Current System summaries can include:

```text
FILE LOCATIONS
MOLECULAR COMPOSITION
SYSTEM TOTALS
SOLVENT / ION COMPOSITION
ELEMENT COMPOSITION
CHAIN COMPOSITION
RESIDUE COMPOSITION
GROUPS
INTERFACES
```

Solved cell geometry, curvature, surface energy, and detailed network topology belong primarily to the Group or Interface output directories.

---

## Molecular Composition

The System summary derives molecular composition from the loaded atomic table.

Depending on the available metadata, it can report:

```text
Molecular Atoms
Non-Water Residues
Waters
Water Atoms
Ions
Non-Water Chains
Total Atoms
Total Residues
```

Water and ions are classified separately in the current information exporter.

See [Solvent](solvent.md) for the exact classification behavior and current consistency limitations.

---

## Element, Chain, and Residue Composition

The System summary can also report elemental and molecular composition directly from the input structure.

Element composition is an atom count by element.

Chain composition reports non-water chains and, when residue identifiers are available, the number of residues represented in each chain.

Residue composition reports non-water residue names and atom counts.

These quantities describe the **input molecular model**. They should not be confused with the solved geometric residue and chain summaries produced from a Group Network.

---

## Configured Groups and Interfaces

The System information file records which Groups and Interfaces belong to the System.

For Groups, the manifest can include the Group name, selected atom/residue counts, and principal build settings such as:

```text
Network Type
Surface Resolution
Box Size
Maximum Allowable Vertex
```

The Interface section provides system-level context for the requested or solved molecular relationships.

The System information file is therefore the best starting point when opening an unfamiliar VorPy result directory.

---

# Group Summary

The Group `info.txt` is the principal human-readable summary of a solved Group Network.

A current Group summary is organized into sections such as:

```text
COMPOSITION
BUILD INFORMATION
BUILD TIMING
VORONOI NETWORK
GROUP GEOMETRY
SURFACE CURVATURE
SURFACE ENERGY ESTIMATE
SURFACE CLASSIFICATION
CHAIN COMPOSITION
RESIDUE COMPOSITION
```

Not every section necessarily contains meaningful values for every molecular system.

---

## Group Composition

The Group composition section reports the molecular selection represented by the Group, including values such as:

```text
Atoms
Residues
Chains
Mass
```

These values describe the selected Group, not the entire parent System.

A Group may therefore contain 32 analyzed atoms while the parent System contains hundreds or thousands of surrounding atoms that remain available to define the cells.

See [Groups](groups.md).

---

# Build Information and Timing

The Group summary records the network settings that produced the geometry.

Important settings include:

```text
Network Type
Surface Resolution
Box Size
Maximum Allowable Vertex
```

Build timing currently separates major stages such as:

```text
Vertex Time
Connection Time
Surface Building Time
Analysis Time
Total Time
```

The analysis implementation also contains a more detailed profiler that can separate setup/cache construction, surface gathering, completeness checks, curvature, neighbors, contacts, point properties, center-of-mass calculations, and assignment overhead.

The detailed profiler is primarily a performance/development diagnostic; the stable summary output retains the higher-level timing categories.

---

# Network Topology Summary

The Group information file records the size of the solved Voronoi Network:

```text
Vertices
Edges
Surfaces
```

These are topology counts for the Group's solved Network, not molecular atom counts.

They provide a useful compact description of geometric complexity and are especially valuable when comparing network settings or partitioning schemes.

---

# Group Geometry

The Group summary currently reports quantities including:

\[
V_{\mathrm{group}}
\]

for Group Voronoi volume,

\[
V_{\mathrm{vdW}}
\]

for reconstructed van der Waals / non-overlap volume,

and:

\[
A_{\mathrm{group}}
\]

for the exposed Group boundary surface area.

It can also report:

```text
Surface Area / Volume
Density
```

The meaning of `Density` requires special care.

---

## Current `Density` Definition

The current Group implementation defines:

\[
\mathrm{Density}
=
\frac{V_{\mathrm{vdW}}}{V_{\mathrm{group}}}.
\]

Despite the output label, this is **not a mass density**.

It is a dimensionless occupied-volume fraction / geometric packing ratio comparing the reconstructed van der Waals volume with the Voronoi volume assigned to the Group.

For scientific writing, a term such as **van der Waals volume fraction**, **occupied-volume fraction**, or **packing fraction** is more descriptive than simply `Density`.

Until the output field is renamed, publications and downstream analyses should explicitly define this quantity.

---

# Complete Cells and Group Aggregation

Cell-level analysis marks whether a solved cell is geometrically complete.

Group volume, van der Waals volume, mass, and related aggregate quantities are currently accumulated from cells marked complete.

This means incomplete boundary/pathological cells can be excluded from some Group-level aggregates.

The completeness criterion should therefore be considered when interpreting Groups containing truncated or otherwise incomplete geometry.

See [Groups](groups.md) for the Group aggregation model.

---

# Group Boundary Surface Area

Group surface area is not simply the sum of the total surface area of every selected atom.

VorPy distinguishes:

```text
selected atom ↔ selected atom
    internal Group surface

selected atom ↔ outside atom
    Group boundary surface
```

The exposed Group boundary is built from surfaces crossing out of the Group's first atom layer.

This avoids treating internal Group-Group cell boundaries as externally exposed molecular surface.

---

# Curvature Summary

The Group summary aggregates curvature over the exposed Group boundary.

Current fields include:

```text
Integrated Mean Curvature
Integrated Mean Curvature Squared
Integrated Gaussian Curvature
Area-Weighted Mean Curvature
Area-Weighted Gaussian Curvature
```

For a Group boundary with total area \(A\),

\[
\bar{H}_A
=
\frac{\int H\,dA}{A}
\]

and:

\[
\bar{K}_A
=
\frac{\int K\,dA}{A}.
\]

Integrated quantities preserve the geometry of the entire Group boundary, while area-weighted quantities provide normalized descriptors useful for comparing differently sized Groups.

See [Curvature](../theory/curvature.md).

---

# Representative Surface Energy

The Group information file can report:

```text
Representative Surface Energy
Representative Energy / Area
```

The current representative model is:

\[
\frac{E_{\mathrm{rep}}}{k_B T}
=
2\int H^2\,dA.
\]

This is a curvature-dependent reference bending metric.

It is **not** a calibrated molecular free energy and should not be interpreted as including solvent thermodynamics, electrostatics, dispersion, or conformational entropy.

See [Surface Energy](../theory/surface_energy.md).

---

# Residue Summaries

VorPy can aggregate solved Group geometry by molecular residue.

A current residue row can contain:

```text
Residue identity
Residue sequence
Atom count
Volume
Total Boundary SA
Inter-Residue SA
Solvent-Interfacial SA
```

These quantities are geometric aggregates derived from the solved Group Network.

They are different from the System-level residue composition, which is primarily an input-structure count.

---

## Residue Volume

Residue volume is obtained from the Voronoi cell volumes associated with atoms assigned to that residue.

Conceptually:

\[
V_R
=
\sum_{i\in R} V_i.
\]

As with other Group volume aggregates, incomplete-cell handling can affect the available result.

---

## Total Residue Boundary Surface Area

For a residue, a surface shared by two atoms belonging to the same residue is internal and is not part of the residue boundary.

A surface contributes to the residue boundary when one defining atom belongs to the residue and the opposing atom does not.

Thus:

```text
atom in residue R ↔ atom in residue R
    internal residue surface

atom in residue R ↔ atom outside residue R
    residue boundary surface
```

`Total Boundary SA` therefore describes the complete geometric boundary of the residue within the solved environment.

---

## Inter-Residue Surface Area

`Inter-Residue SA` is more restrictive.

It represents surfaces between atoms belonging to different residues represented within the analyzed molecular Group.

This helps distinguish residue-residue molecular contacts from other parts of the residue boundary.

---

## Solvent-Interfacial Surface Area

`Solvent-Interfacial SA` currently represents surfaces whose opposing atom is recognized as water by the Group output classifier.

As discussed in [Solvent](solvent.md), the current implementation is more precisely **water-interfacial surface area** than an all-solvent metric because ions are not included in that water-specific classification.

---

## Image Placeholder 2 — Residue Summary Decomposition

![Placeholder: residue surface-area decomposition](../../assets/docs/summaries/residue_summary_decomposition.png)

> **Suggested figure:** Highlight one residue and divide its boundary into three categories: surfaces against the same residue (internal and excluded), surfaces against other molecular residues, and surfaces against water. Show how the latter two contribute to `Total Boundary SA`, while the water portion is additionally classified as solvent-interfacial.

---

# Chain Summaries

Chain summaries use the same general logic at a larger molecular scale.

A current chain row can report:

```text
Chain identity
Atom count
Residue count
Volume
Total Boundary SA
Inter-Chain SA
Solvent-Interfacial SA
```

For chain boundaries:

```text
same selected chain ↔ same selected chain
    internal

selected chain ↔ different selected chain
    inter-chain

selected chain ↔ recognized water
    solvent-interfacial
```

The same surface may belong to more than one descriptive category.

For example, water-interfacial area is part of total chain boundary area rather than a quantity that should be added to it again.

---

# Atom / Cell Summaries

The most detailed per-cell measurements are stored on `net.balls` and exported through the network log.

The analysis stage currently calculates or stores atom/cell fields spanning several categories:

```text
Geometry
    Volume
    Surface Area
    van der Waals / non-overlap volume
    overlap volume

Curvature
    maximum mean curvature
    maximum Gaussian curvature
    integrated mean curvature
    integrated mean-curvature squared
    integrated Gaussian curvature

Shape
    sphericity
    isoperimetric quotient

Topology
    complete cell
    inner ball
    number of neighbors
    nearest neighbor
    neighbor distances

Contact / overlap
    number of overlaps
    contact area
    surface overlap

Spatial descriptors
    minimum / maximum point distance
    center-related quantities
    moment of inertia
    bounding box
```

The exact availability of expensive fields can depend on analysis settings.

---

# Optional Analysis Features

The current analyzer preserves a master `complicated` switch while also allowing individual expensive analyses to be enabled or disabled.

Feature switches include:

```text
spikes
contacts
second_neighbors
com
moi
bounding_box
```

When an option is `None`, it inherits the value of `complicated`.

This means two otherwise similar network builds can contain different detailed summary fields if expensive analyses were disabled.

For reproducibility, record non-default analysis settings when comparing these descriptors.

---

# Network Logs

The standard network log is the principal machine-readable quantitative record.

It contains multiple CSV sections rather than one flat table:

```text
build information
group information
Atoms
Surfaces
Edges
Vertices
```

The log preserves both molecular identity and geometric/topological quantities.

---

## Build Information in Logs

The build section currently records fields such as:

```text
Name
Location
Completion Date
Network Type
Surface Resolution
Box Size
Maximum Allowable Vertex
Total Time
Vertex Time
Connect Time
Surface Building Time
Analysis Time
Maximum Found Vertex
VorPy Version
```

These fields provide provenance for the detailed numerical data that follow.

---

# Group Information in Logs

Before writing detailed topology tables, the Group log calls the Group summary calculation and writes a compact compatibility row containing:

```text
Name
Volume
Surface Area
Mass
Density
Center of Mass
VDW Volume
VDW Center of Mass
Moment of Inertia
Spatial Moment of Inertia
```

Because this CSV section is parsed positionally by the existing log reader, its schema is part of the current backward-compatibility contract.

Changes to this row should therefore be versioned and tested carefully.

---

# Atom Rows in Logs

Current atom log rows include molecular identifiers, coordinates, radius, geometric quantities, curvature, representative surface energy, shape descriptors, neighbor information, overlap/contact information, center-related quantities, bounding boxes, and explicit neighbor lists.

Importantly, the standard Group log writer currently skips an atom when:

```text
Surface Area == 0
```

or when:

```text
Complete Cell? == False
```

Therefore the `Atoms` section of a Group log should **not automatically be interpreted as a complete copy of `System.balls` or even every requested Group atom**.

It is a record of solved atom/cell rows that pass the writer's current output criteria.

---

# Surface Rows in Logs

Each surface row identifies its two defining balls and records quantities including:

```text
Surface Area
Mean Curvature
Average Mean Curvature
Gaussian Curvature
Average Gaussian Curvature
Integrated Mean Curvature
Integrated Mean Curvature Squared
Integrated Gaussian Curvature
Representative Surface Energy
Ball 1 Volume Contribution
Ball 2 Volume Contribution
Contact Area
Overlap
```

These rows are the most direct quantitative representation of shared cell boundaries.

They are particularly useful for reconstructing atom-atom contacts and molecular interfaces.

---

# Edge and Vertex Rows

Edge rows record the defining three balls and edge length.

Vertex rows record the defining four balls, vertex location, and vertex radius.

Together with the Surface and Atom sections, these tables expose the topology and geometry underlying the higher-level summaries.

---

# Interface Summaries

Interfaces require their own summary semantics because an Interface Network is not necessarily a collection of complete closed cells.

The Interface information output can report:

```text
Interface identity
Group 1 / Group 2 definitions
Atom membership
Network topology
Direct inter-group surfaces
Group 1 internal surfaces
Group 2 internal surfaces
Supporting surfaces
Surface-area statistics
Curvature statistics
Contact / overlap statistics
Interface-associated water information
```

See [Interfaces](interfaces.md) for the exact surface classifications.

---

## Interface Logs and Compatibility Fields

Interface logs intentionally preserve the standard Group-log section ordering and schemas so the existing log reader can parse them.

This requires a `group information` compatibility row even though several ordinary Group quantities are not physically defined for an Interface Network.

For the current Interface writer:

```text
Surface Area
Mass
Center of Mass
```

are populated where meaningful, while closed-volume-dependent fields such as:

```text
Volume
Density
VDW Volume
```

are written as undefined (`NaN`) because an Interface-only Network does not necessarily enclose a closed volume.

This distinction is scientifically important: compatibility with a file schema does not imply that every Group metric has a meaningful Interface analogue.

---

# Network Checkpoint vs. Analysis Summary

VorPy also contains a network checkpoint writer.

A network checkpoint is conceptually different from a scientific summary.

The checkpoint records enough network construction data to help restore or inspect network state, including:

```text
network settings
vertex topology and coordinates
edge topology
surface topology
surface points
surface triangles
```

It is oriented toward reconstruction and geometry storage.

The standard analysis log, by contrast, stores calculated atom/surface metrics and molecular metadata.

Conceptually:

```text
Network checkpoint
    "How was the geometric network represented?"

Analysis log
    "What quantities were calculated from that network?"

Info summary
    "What are the main scientifically interpretable results?"
```

These files should not be described as interchangeable.

---

## Image Placeholder 3 — Summary vs. Log vs. Checkpoint

![Placeholder: summary, analysis log, and checkpoint roles](../../assets/docs/summaries/output_roles.png)

> **Suggested figure:** Three columns. `info.txt` = compact scientific interpretation; `*_logs.csv` = detailed atom/surface/edge/vertex analysis tables; `*_net.csv` = network checkpoint/geometry state. Show arrows from the same solved Network into all three outputs.

---

# Important Current Validation Issues

Several summary/output fields require care in the current implementation.

### `Density`

`Density` is currently \(V_{\mathrm{vdW}} / V_{\mathrm{Voronoi}}\), not mass per volume.

### `VDW Center of Mass`

The Group implementation currently accumulates a mass-weighted position numerator but has an inconsistent denominator in the `vdw_com` calculation. This field remains under validation and should not yet be treated as an authoritative physical center of mass.

### `Non-Overlap Volume` and `Overlap Volume` in Group logs

The current standard Group log header orders the fields as:

```text
Non-Overlap Volume
Overlap Volume
```

but the writer currently outputs:

```text
atom.olap_vol
atom.vdw_vol
```

in those positions.

The values are therefore reversed relative to the current column labels and should be corrected before those two exported columns are used as authoritative labeled values.

### `Isometric Quotient`

The analysis calculation calls an isoperimetric-quotient function but stores/exports the historical field name:

```text
Isometric Quotient
```

The terminology should be standardized so the mathematical descriptor and the output label agree.

### Incomplete cells

Incomplete cells can be excluded from Group aggregate calculations and from standard Group atom-log rows.

When counts differ between input atoms, Group membership, and exported atom rows, cell completeness is one possible explanation.

---

# Which Summary Should I Use?

For interpreting an analysis:

| Question | Primary output |
|---|---|
| What molecular system was loaded? | System `info.txt` |
| How many waters/ions/chains/residues exist? | System `info.txt` |
| What Groups and Interfaces were configured? | System `info.txt` |
| What settings produced this Group geometry? | Group `info.txt` |
| What is the Group volume or exposed boundary area? | Group `info.txt` |
| What are the residue/chain geometric summaries? | Group `info.txt` |
| What are the integrated curvature values? | Group `info.txt` |
| Which atom is adjacent to which atom? | Network log |
| What is one atom's volume/contact/curvature? | Network log |
| What is the area of one shared surface? | Network log |
| Which balls define one edge or vertex? | Network log |
| What is the direct Group 1–Group 2 interface? | Interface `info.txt` |
| What topology/mesh state was stored for reconstruction? | Network checkpoint |

---

# Recommended Scientific Workflow

A practical interpretation workflow is:

```text
1. System info
      ↓
   verify molecular composition and analysis definitions

2. Group / Interface info
      ↓
   inspect settings, topology, and aggregate geometry

3. Residue / chain summaries
      ↓
   identify molecular-scale patterns

4. Network logs
      ↓
   investigate specific atoms and shared surfaces

5. Geometry exports
      ↓
   inspect the result visually in three dimensions
```

This approach keeps detailed cell-level results connected to the molecular and geometric context that produced them.

---

# Reproducibility

When reporting a VorPy summary, preserve enough context to reproduce the geometry.

At minimum, retain:

```text
Input structure
VorPy version
Group / Interface definition
Network type
Atomic radius model when modified
Surface resolution
Maximum vertex setting
Other non-default analysis settings
```

For quantitative comparisons between structures or partitioning schemes, use equivalent molecular selections and compatible build/analysis settings unless those settings are intentionally part of the comparison.

---

# Summary

VorPy's output hierarchy is designed to move from detailed geometric data to molecular interpretation.

The **System summary** describes the molecular model and configured analyses.

The **Group summary** describes one solved molecular selection and aggregates its boundary geometry, curvature, energy reference metric, residue geometry, and chain geometry.

The **Interface summary** describes the geometry connecting two molecular selections and distinguishes direct inter-group surfaces from supporting network geometry.

The **network log** exposes detailed atom, surface, edge, and vertex records.

The **network checkpoint** stores network state for reconstruction rather than serving as a scientific summary.

These levels should be interpreted together: summaries provide context and aggregation, while the network tables preserve the detailed geometry from which those summaries are derived.
