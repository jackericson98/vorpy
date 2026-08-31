# Contacts

VorPy defines molecular contacts from the topology and geometry of the solved Voronoi network.

Two atoms are considered **Voronoi neighbors** when their Voronoi cells share a surface. In this sense, a contact is not defined by a user-selected radial cutoff. Instead, it is determined by the local spatial partition produced by the selected Voronoi scheme.

This distinction is important because VorPy can describe both:

- **contact topology** — which atoms are neighbors, and
- **contact geometry** — the geometry associated with the shared boundary between neighboring cells.

---

## What Is a Voronoi Contact?

For two atoms \(i\) and \(j\), a Voronoi contact exists when their solved cells share a Voronoi surface.

Conceptually:

```text
Atom i  <---- shared Voronoi surface ---->  Atom j
```

The corresponding atoms are therefore direct neighbors in the Voronoi network.

This gives VorPy a topological definition of molecular adjacency:

```text
shared surface
      |
      v
cell i ----- cell j
  |             |
atom i        atom j
```

A nearby atom that does not share a Voronoi surface is **not** a Voronoi neighbor, even if its center lies within a conventional radial contact cutoff.

> **Important:** A Voronoi contact represents geometric adjacency. It does not, by itself, imply a hydrogen bond, salt bridge, favorable van der Waals interaction, electrostatic interaction, or other specific chemical interaction.

---

## Image Placeholder 1 — Definition of a Contact

> **Suggested figure:** Show 3–5 atoms as spheres together with their Voronoi cells. Highlight one shared surface between two atoms. On that surface, distinguish the **complete Voronoi surface** from the smaller **contact-area region lying inside the overlapping balls**. Include at least one nearby atom that does **not** share a Voronoi surface.
>
> Suggested filename:
>
> `assets/docs/contacts/voronoi_contact_definition.png`

![Placeholder: Voronoi contact definition](../../assets/docs/contacts/voronoi_contact_definition.png)

**Figure 1. Voronoi definition of an atomic contact.** Two atoms are Voronoi neighbors when their cells share a surface. Spatial proximity alone does not guarantee that two atoms share a Voronoi boundary.

---

## Voronoi Contacts vs. Distance-Cutoff Contacts

Many molecular-analysis methods define contacts using an interatomic distance threshold:

```text
distance(i, j) <= cutoff
```

This is useful for many purposes, but the resulting contact network depends directly on the selected cutoff.

VorPy instead determines adjacency from the spatial decomposition:

```text
cell i shares a surface with cell j
```

This means the local molecular environment influences which atoms become neighbors.

A Voronoi contact therefore answers a different question:

> **Which atoms directly share the partitioned molecular space?**

rather than:

> **Which atoms lie within a selected distance of one another?**

Neither definition is universally superior. They describe different geometric relationships and should be interpreted according to the scientific question being asked.

---

## Image Placeholder 2 — Distance Cutoff vs. Voronoi Contact

> **Suggested figure:** Use the same small molecular configuration twice. On the left, show a circular or spherical radial cutoff around one atom. On the right, show the Voronoi cells and shared surfaces. Highlight atoms that are included by one definition but excluded by the other.
>
> Suggested filename:
>
> `assets/docs/contacts/cutoff_vs_voronoi_contacts.png`

![Placeholder: Distance cutoff versus Voronoi contacts](../../assets/docs/contacts/cutoff_vs_voronoi_contacts.png)

**Figure 2. Distance-based and Voronoi definitions of molecular contact.** A distance cutoff identifies atoms within a selected radius, whereas a Voronoi contact identifies atoms whose cells share a boundary.

---

## Contact Topology

The solved network stores the neighbor relationships associated with each atom.

At the atom level, the current network log contains fields including:

| Field | Meaning |
|---|---|
| `Number of Neighbors` | Number of Voronoi-neighbor cells associated with the atom |
| `neighbors` | Indices of the neighboring atoms/cells |
| `Closest Neighbor` | Index of the nearest neighbor according to the current VorPy neighbor-distance calculation |
| `Closest Neighbor Distance` | Distance value associated with the closest-neighbor calculation |
| `Layer Distance Average` | Summary of the current neighbor-layer distance calculation |
| `Layer Distance RMSD` | Spread associated with the neighbor-layer distance calculation |
| `Number of Overlaps` | Number of overlaps reported for the atom |
| `Contact Area` | Atom-level contact-area quantity calculated by VorPy |

The `neighbors` field provides the explicit connectivity of the Voronoi network.

For example, an atom with:

```text
Number of Neighbors = 12
neighbors = [1, 3, 4, 5, ...]
```

has twelve cells sharing surfaces with its own cell.

### Neighbor Count

The number of neighbors provides a simple measure of local network connectivity.

High neighbor counts generally indicate that a cell shares boundaries with many surrounding cells, while lower values indicate fewer direct Voronoi connections.

Neighbor count should not automatically be interpreted as coordination number in a specific chemical sense. The value is defined by the selected spatial partition.

---

## Contact Geometry

A Voronoi contact is associated with a specific shared surface.

The surface-level portion of the network log identifies the two cells associated with every surface using:

```text
Ball 1
Ball 2
```

and records geometric quantities including:

| Field | Meaning |
|---|---|
| `Surface Area` | Total area of the solved Voronoi surface |
| `Contact Area` | Contact-area quantity associated with the surface |
| `Overlap` | Overlap quantity associated with the two balls |
| `Ball 1 Volume Contribution` | Surface-associated volume contribution for Ball 1 |
| `Ball 2 Volume Contribution` | Surface-associated volume contribution for Ball 2 |

Curvature and representative surface-energy quantities are also stored on each surface, allowing the geometry of a contact boundary to be analyzed beyond its area alone.

### Surface Area and Contact Area Are Distinct

The current VorPy log reports both:

```text
Surface Area
Contact Area
```

for a shared surface.

These are intentionally different quantities.

- **Surface Area** is the area of the complete solved Voronoi boundary between two neighboring cells.
- **Contact Area** is a mesh-based estimate of the portion of that boundary lying within the atomic sphere associated with the cell being analyzed.

For weighted Voronoi surfaces, points on the ideal shared boundary that lie inside one neighboring sphere are also expected to lie inside the other under the corresponding equal-distance condition. Contact area therefore serves as a practical surface-based representation of the overlap region between neighboring balls.

It should nevertheless be interpreted as a **geometric overlap descriptor**, not as an exact analytic sphere-intersection area.

---

## How VorPy Calculates Contact Area

VorPy calculates contact area directly from the triangulated Voronoi surface.

For each surface attached to an atom, every mesh point is tested against that atom's sphere:

```text
||p - c|| <= r
```

where:

- `p` is a surface mesh point,
- `c` is the ball center,
- `r` is the ball radius.

The implementation uses the equivalent squared-distance comparison:

```text
||p - c||² <= r²
```

Each triangle is then classified according to how many of its three vertices lie inside the sphere.

### Triangle Classification

```text
3 vertices inside
    -> full triangle area contributes to Contact Area

0 vertices inside
    -> triangle contributes no Contact Area

1 or 2 vertices inside
    -> full triangle area contributes to Contact Area
```

This last case is important.

VorPy currently does **not** analytically clip a triangle where it crosses the sphere boundary. Instead, a mixed triangle is counted in full if at least one of its vertices lies inside the sphere.

Contact area is therefore a discretized approximation whose precision depends partly on the surface triangulation resolution.

As surface resolution is refined, the boundary triangles become smaller and the approximation should more closely represent the portion of the Voronoi surface contained within the overlapping balls.

### Practical Interpretation

A useful interpretation is:

> **Contact Area is the triangulated Voronoi surface area associated with the region where neighboring atomic spheres overlap.**

More precisely, it is the sum of surface triangles classified as at least partially contained within the analyzed ball.

Because the calculation is performed on the Voronoi surface rather than directly on the sphere-sphere intersection geometry, contact area is complementary to overlap depth and overlap volume rather than equivalent to either quantity.

---

## Overlap Information

VorPy also records the degree to which neighboring atomic spheres overlap.

For two neighboring balls \(i\) and \(j\), VorPy first calculates the surface-to-surface separation:

```text
d_surface = ||c_j - c_i|| - r_i - r_j
```

where `c_i` and `c_j` are the ball centers and `r_i` and `r_j` are their radii.

This gives:

```text
d_surface > 0    separated spheres
d_surface = 0    tangent spheres
d_surface < 0    overlapping spheres
```

The surface-level `Overlap` value is then:

```text
Overlap = max(-d_surface, 0)
```

and therefore has units of length.

At the atom level, `Number of Overlaps` counts first-layer Voronoi neighbors whose surface-to-surface separation is negative.

### Closest Neighbor Distance

The logged `Closest Neighbor Distance` is therefore not a center-to-center distance.

It is the minimum **surface-to-surface separation** among the atom's Voronoi neighbors.

Negative values are expected and indicate sphere overlap.

### Three Complementary Overlap Descriptors

VorPy currently provides three related geometric descriptions:

| Quantity | Units | Interpretation |
|---|---:|---|
| `Overlap` | Å | Linear penetration depth between neighboring spheres |
| `Contact Area` | Å² | Voronoi boundary area associated with the overlap region |
| `Overlap Volume` | Å³ | Sphere volume excluded by the reconstructed non-overlapping volume |

These quantities describe different aspects of the same local geometry and should not be treated as interchangeable.

---

## Contacts Depend on the Partitioning Scheme

VorPy supports three primary spatial partitioning schemes:

- additively weighted (`aw`)
- power / Laguerre (`pow`)
- primitive (`prm`)

Because each scheme partitions space differently, the resulting contact network can also differ.

Changing the partitioning scheme may affect:

- which atoms are neighbors,
- the number of neighbors,
- the geometry of shared surfaces,
- contact area,
- interfacial area,
- higher-level residue or group contacts.

This is especially important for systems containing atoms or particles with different radii.

A contact should therefore always be interpreted together with the network type used to generate it.

---

## Image Placeholder 3 — Contact Networks Across AW, POW, and PRM

> **Suggested figure:** Use the same small molecule or local molecular environment solved with AW, POW, and PRM. Highlight contacts that are common to all three schemes and contacts that appear only in one or two schemes.
>
> Suggested filename:
>
> `assets/docs/contacts/contact_scheme_comparison.png`

![Placeholder: AW, POW, and PRM contact comparison](../../assets/docs/contacts/contact_scheme_comparison.png)

**Figure 3. Partitioning-scheme dependence of Voronoi contacts.** The same molecular coordinates can produce different neighbor networks under additively weighted, power, and primitive spatial decompositions.

---

## Molecular Hierarchy

Raw Voronoi contacts are defined between cells, but the corresponding atoms carry molecular identity.

This allows contacts to be interpreted at several structural levels, including:

```text
atom <-> atom
residue <-> residue
chain <-> chain
molecule <-> molecule
group <-> group
solute <-> solvent
```

The atom rows in the network logs include molecular identifiers such as:

- atom name,
- residue name,
- residue sequence,
- chain,
- atom index.

This makes it possible to map cell-level neighbor relationships back onto the molecular structure.

For example, an atom-level contact can be classified as:

```text
same residue
different residue
same chain
different chain
protein-solvent
protein-DNA
protein-ligand
```

depending on the molecular identities associated with the two cells.

---

## Contacts and Interfaces

Contacts and interfaces are closely related but describe different levels of organization.

A **contact** describes the relationship between two neighboring cells.

An **interface** aggregates shared boundaries connecting two molecular groups.

Conceptually:

```text
Atom contacts
     |
     v
Shared Voronoi surfaces
     |
     v
Group-group interface
```

For example, a protein-DNA interface can contain many individual atom-atom Voronoi contacts.

This allows VorPy to move naturally between:

- individual atomic neighbors,
- residue-level contact patterns,
- complete intermolecular interfaces.

See [Interfaces](interfaces.md).

---

## Solvent Contacts

When solvent is included in the solved system, solvent cells can appear as neighbors of molecular cells.

These contacts allow VorPy to distinguish between surfaces associated with:

- internal molecular neighbors,
- neighboring molecular groups,
- surrounding solvent.

This relationship contributes to analyses of solvent-facing surface area and molecular exposure.

See [Solvent Analysis](solvent.md).

---

## Where Contact Information Is Stored

VorPy distributes contact-related information across the normal output hierarchy rather than producing a single dedicated contact file.

### Network Logs

The network log contains the most detailed contact information.

Atom records include neighbor counts, neighbor indices, overlap information, and contact-related quantities.

Surface records identify the two balls connected by each surface and contain surface-level geometric information.

### Group `info.txt`

The group information file provides the context for the network from which the contacts were generated, including:

- group composition,
- network type,
- surface resolution,
- build settings,
- network size,
- group geometry,
- curvature and energy summaries,
- molecular composition summaries.

This is important because contact relationships should be interpreted together with the network type and group definition that produced them.

### System `info.txt`

The system information file describes the molecular system and configured groups/interfaces.

It therefore provides the system-level context necessary to interpret network contacts across:

- molecular composition,
- groups,
- interfaces,
- solvent and ion composition.

---

## EDTA Example

A default EDTA calculation provides a useful small example of the contact representation.

The group contains 32 atoms and was solved using the additively weighted (`aw`) network with a surface resolution of `0.2`. The resulting network contains 500 vertices, 892 edges, and 424 surfaces.

The atom-level network log then records the neighbor list for each of those atoms, while the surface section records the two balls associated with each of the 424 shared surfaces.

This is useful for illustrating the distinction between:

```text
network topology
    -> which cells are connected

surface geometry
    -> how the connected cells share space

molecular identity
    -> which atoms/residues/groups those cells represent
```

---

## Image Placeholder 4 — Real VorPy EDTA Contact Example

> **Suggested figure:** Generate this one directly from the EDTA example used in the documentation. Pick one central EDTA atom and show:
>
> 1. the atom as a sphere,
> 2. its solved Voronoi cell,
> 3. all Voronoi-neighbor atoms,
> 4. one or two highlighted shared surfaces,
> 5. labels containing the corresponding atom indices/names.
>
> A second panel could show the same local neighborhood as a simple contact graph.
>
> Suggested filename:
>
> `assets/docs/contacts/edta_contact_example.png`

![Placeholder: EDTA VorPy contact example](../../assets/docs/contacts/edta_contact_example.png)

**Figure 4. Example VorPy contact neighborhood for EDTA.** A local Voronoi neighborhood can be visualized both as shared three-dimensional cell boundaries and as a network of atom-atom connections.

---

## Interpretation Guidelines

When reporting VorPy contacts, record enough information for the calculation to be reproduced.

At minimum, report:

- molecular structure or input model,
- VorPy version,
- partitioning scheme,
- atomic radii or radius model,
- relevant group selection,
- relevant build settings.

For quantitative surface comparisons, also record the surface-resolution settings used to construct the geometry.

### Recommended Terminology

Use:

- **Voronoi neighbor**
- **Voronoi contact**
- **shared Voronoi surface**
- **contact area**
- **contact count**
- **neighbor count**

Avoid treating a Voronoi contact as synonymous with:

- chemical bond,
- hydrogen bond,
- salt bridge,
- energetic interaction,
- experimentally observed binding contact.

Those interpretations require additional chemical or energetic criteria.

---

## Summary

VorPy contacts describe the topology and geometry of molecular spatial adjacency.

The key ideas are:

1. two cells are Voronoi neighbors when they share a surface,
2. the corresponding atoms form a Voronoi contact,
3. the network log stores explicit atom neighbor lists,
4. shared surfaces provide geometric information about each connection,
5. molecular metadata maps those contacts onto residues, chains, molecules, and groups,
6. contact networks can differ between AW, POW, and PRM decompositions,
7. collections of contacts form the basis for higher-level interface analysis.

This makes contacts a bridge between VorPy's low-level geometric network and its higher-level molecular analysis.
