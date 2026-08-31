# Interfaces

Interfaces describe the geometric boundary between two molecular selections in VorPy.

An Interface is constructed from two Groups and owns a dedicated Network containing the topology relevant to their interaction. This allows VorPy to isolate inter-group surfaces, contacts, curvature, overlap, supporting geometry, and interface-associated solvent without requiring the complete Group networks to be treated as the interface itself.

---

## What Is an Interface?

Conceptually, an interface begins with two molecular selections:

```text
Group 1                     Group 2
   │                           │
   └──────── Interface ────────┘
                 │
                 ▼
       dedicated interface Network
                 │
        ├── vertices
        ├── edges
        ├── surfaces
        ├── direct inter-group surfaces
        ├── internal/supporting surfaces
        └── interface-associated waters
```

The two Groups remain independent molecular selections. The Interface stores references to them, builds its own Network, and records the relationship back on the participating Groups.

A single Group may therefore participate in multiple Interfaces.

---

## Image Placeholder 1 — Two Groups and Their Interface

![Placeholder: Group 1 and Group 2 interface](../../assets/docs/interfaces/group1_group2_interface.png)

> **Suggested figure:** Show two molecular Groups with distinct highlighted atoms. Overlay their Voronoi cells and emphasize the shared Group 1–Group 2 surfaces. Include a few nearby cells that help close the geometry but do not represent direct inter-group contact.

---

## Interface Definitions

An Interface is normally defined from:

- `Group 1`
- `Group 2`

If a second Group is not supplied, VorPy constructs a concrete **surrounding Group** from nearby system atoms.

The surrounding selection is spatially restricted rather than containing every atom outside Group 1. This keeps the dedicated interface calculation focused on nearby geometry.

The resulting Interface stores:

- the two source Groups,
- their system-level ball-index sets,
- a stable interface identifier,
- independent network settings copied from Group 1,
- its dedicated Network,
- interface metadata,
- interface-associated water geometry,
- full water Groups where requested.

---

## Interface Identity

For an explicit pair of Groups, VorPy constructs a stable identifier from their Group identifiers.

The two identifiers are sorted before being combined, so:

```text
Group A – Group B
```

and

```text
Group B – Group A
```

resolve to the same pairwise interface identifier.

This is useful when the same molecular relationship is encountered from different analysis paths.

---

## Dedicated Interface Network

The Interface owns a Network separate from the complete Group Networks.

The Network receives:

- the complete System coordinates,
- the complete System radii,
- the complete System masses,
- the combined atom indices from the two interface sides,
- the original Group 1 and Group 2 index sets,
- a copy of Group 1's network settings.

The complete System geometry remains available so that surrounding atoms can participate in geometric validity checks.

The two Group index sets are additionally supplied to the Network as interface-side membership. Interface vertex construction can then restrict retained topology to geometry involving both interface sides.

> **Important:** An Interface Network is not simply the union of two already-built Group Networks. It is its own network calculation, constructed specifically for the pairwise relationship.

---

## Image Placeholder 2 — Interface Network Filtering

![Placeholder: retained and rejected interface topology](../../assets/docs/interfaces/interface_network_filtering.png)

> **Suggested figure:** Show several candidate Voronoi vertices/surfaces around two Groups. Highlight retained geometry that involves both interface sides and gray out geometry belonging only to unrelated surrounding atoms. The figure should make clear that complete-system geometry can participate in validity checks even when the retained topology is interface-focused.

---

## Direct and Supporting Surfaces

VorPy classifies surfaces in the completed Interface Network according to the defining balls on each side of the surface.

### Direct inter-group surfaces

A surface is a **direct interface surface** when:

```text
one defining ball belongs to Group 1
and
the other defining ball belongs to Group 2
```

These surfaces are the clearest geometric representation of direct Group 1–Group 2 adjacency.

### Group 1 internal surfaces

Both defining balls belong to Group 1.

### Group 2 internal surfaces

Both defining balls belong to Group 2.

### Supporting surfaces

At least one defining ball falls outside the two explicit interface Groups.

Supporting surfaces may be necessary to construct or close the dedicated Interface Network but do not represent direct Group 1–Group 2 contact.

This classification is important because the complete Interface Network contains more geometric information than the direct molecular interface alone.

---

## Interface Surface Area

The most direct inter-group surface-area quantity is obtained by summing the areas of the **direct Group 1–Group 2 surfaces**.

This should be distinguished from:

```text
Full interface-network surface area
    all retained surfaces in the Interface Network

Direct inter-group surface area
    only surfaces with one defining ball in each Group
```

For molecular interpretation, the direct inter-group set is generally the appropriate set when asking how much Voronoi boundary is shared by the two Groups.

---

## Interface Contacts and Overlap

Each retained surface can also carry the same surface-level contact and overlap descriptors used elsewhere in VorPy.

For a collection of interface surfaces, VorPy can summarize:

- total surface area,
- mean/minimum/maximum surface area,
- curvature statistics,
- integrated curvature,
- total Contact Area,
- number of surfaces with positive Contact Area,
- total overlap measure,
- number of surfaces with positive overlap.

The **direct inter-group surface statistics** are particularly useful because they restrict those quantities to Group 1–Group 2 boundaries.

Contact Area should not be confused with Voronoi surface area. Contact Area is the mesh-based estimate of the portion of a shared Voronoi surface associated with overlap between the neighboring balls, while surface area refers to the complete retained Voronoi boundary.

See [Contacts](contacts.md) for the exact Contact Area and overlap definitions.

---

## Image Placeholder 3 — Interface Surface Area vs. Contact Area

![Placeholder: interface surface and contact area](../../assets/docs/interfaces/interface_surface_contact_area.png)

> **Suggested figure:** Show two atoms from opposing Groups sharing a Voronoi surface. Highlight the complete shared surface as **Direct Interface Surface Area** and the smaller sphere-overlap-associated portion as **Contact Area**. Place this within a larger Group 1–Group 2 molecular interface to connect atom-level and group-level interpretation.

---

## Interface Curvature

Interface surface collections can also be summarized using local and integrated curvature quantities.

Current interface information output includes, where available:

- surface-average mean curvature,
- area-weighted surface-average mean curvature,
- minimum/maximum surface-average mean curvature,
- maximum local mean curvature,
- surface-average Gaussian curvature,
- area-weighted surface-average Gaussian curvature,
- minimum/maximum surface-average Gaussian curvature,
- maximum local Gaussian curvature,
- integrated mean curvature,
- integrated squared mean curvature,
- integrated Gaussian curvature,
- area-normalized integrated curvature.

These quantities can be reported separately for:

- the complete Interface Network,
- direct Group 1–Group 2 surfaces,
- Group 1 internal surfaces,
- Group 2 internal surfaces,
- supporting surfaces.

This separation is important when curvature is intended to describe the molecular interface rather than all geometry needed to construct the Interface Network.

---

# Interface-Associated Water

VorPy can identify water molecules associated with a completed Interface Network and can build complete water-cell Groups for those waters.

This functionality provides a way to move from:

```text
Which waters participate in the retained interface geometry?
```

to:

```text
What is the complete Voronoi geometry of each participating water molecule?
```

This is useful for studying solvent molecules occupying or contacting molecular interfaces.

---

## Identifying Interface Waters

VorPy first gathers every ball index appearing on a retained Interface Network surface.

It then searches the System's solvent-water residues.

A water residue is identified as interface-associated when **at least one atom from that water appears among the defining balls of a retained Interface Network surface**.

Conceptually:

```text
retained interface surfaces
          │
          ▼
collect defining ball indices
          │
          ▼
compare against water-residue atoms
          │
          ▼
water has at least one matching atom
          │
          ▼
interface-associated water
```

This is a topology-based definition tied to the retained Interface Network.

> **Current implementation note:** This test uses the complete retained Interface Network surface set, not only the direct Group 1–Group 2 surface subset. The exact scientific interpretation of "interface-associated water" should therefore be reported with this definition until a narrower criterion is deliberately implemented.

---

## Extracting Water-Associated Interface Geometry

For each identified water residue, VorPy extracts retained Interface Network geometry associated with that water.

The stored geometry includes:

- the water's atom indices,
- retained surfaces defined by at least one water atom,
- associated edges,
- associated vertices,
- summed retained surface area,
- summed retained Contact Area.

Additional retained edges or vertices directly defined by a water atom are also included even if they are not reached through the initially collected surface list.

This produces an interface-local geometric description of each participating water.

---

## Building Complete Water Cells

The interface-local geometry does not by itself represent the complete Voronoi cells of the water molecule.

To obtain complete cell geometry, VorPy can construct a normal Group for each identified water residue:

```text
interface-associated water
          │
          ▼
create Group(residues=[water])
          │
          ▼
build normal VorPy Network
          │
          ▼
complete water-cell geometry
```

The water Group uses the parent System and copies the primary interface Group's settings. Because it is built as a normal Group, the complete System geometry remains available to define the water's cells.

The resulting water summary can contain:

- Group name,
- water atom indices,
- vertex count,
- edge count,
- surface count,
- total surface area,
- volume,
- full vertex table,
- full edge table,
- full surface table.

This provides both:

1. **interface-local water geometry** from the Interface Network, and
2. **complete water-cell geometry** from a dedicated water Group.

That distinction should be preserved in analysis and reporting.

---

## Image Placeholder 4 — Interface Water and Complete Water Cell

![Placeholder: interface-associated water and complete cell](../../assets/docs/interfaces/interface_water_cell.png)

> **Suggested figure:** Show two molecular Groups forming an interface with one water molecule between them. Panel A highlights only the retained Interface Network geometry involving the water. Panel B shows the same water rebuilt as its own Group with its complete Voronoi cell(s). Label these **Interface-local water geometry** and **Complete water-cell geometry**.

---

## Current Water-Analysis Status

The water-cell workflow is currently under active validation.

The Interface build identifies all touching waters, but the present implementation intentionally passes only the **first identified water** into the complete water-Group builder as a temporary test.

Therefore:

> Interface-water detection can currently identify multiple waters, but complete water-cell construction during `Interface.build()` is temporarily limited to the first identified water.

This limitation should remain visible in development documentation and testing until the temporary restriction is removed.

The water-volume summary also contains an implementation note indicating that its authoritative volume source still needs to be finalized.

These limitations do not change the conceptual workflow, but they should be resolved before the water-cell feature is presented as a fully validated production analysis.

---

## Interface Outputs

VorPy can export interface-specific:

- atoms,
- surfaces,
- edges,
- vertices,
- logs,
- `info.txt`,
- participating Group information.

The interface `info.txt` includes:

### Interface definition

- interface identifier,
- Group 1 name,
- Group 2 name.

### Atom membership

- Group 1 atom count,
- Group 2 atom count,
- unique Interface Network atom count,
- atoms shared by both definitions.

### Network topology

- vertices,
- edges,
- surfaces,
- edges per vertex,
- surfaces per vertex.

### Surface classification

- direct Group 1–Group 2 surfaces,
- Group 1 internal surfaces,
- Group 2 internal surfaces,
- supporting surfaces.

### Geometric statistics

Surface area, curvature, contact, and overlap statistics are written separately for the different surface classes.

This makes the information file useful both for understanding the complete dedicated network and for extracting the direct molecular interface.

---

## Interfaces Without Complete Group Networks

An Interface can be defined and analyzed without requiring both source Groups to have already built complete Networks.

Interface-specific Group information can therefore record:

- the molecular selection,
- the number of selected atoms,
- whether a complete Group Network exists,
- selection counts,
- interface metadata.

This helps separate **molecular Group definition** from **whether a complete standalone Group Network has been solved**.

---

## Interpretation Guidelines

When reporting a VorPy Interface, specify:

- the two Group definitions,
- whether the second side was explicit or automatically surrounding,
- network type,
- atomic radius model where relevant,
- surface resolution,
- other non-default build settings,
- whether reported surface area refers to the full Interface Network or direct inter-group surfaces,
- whether water analysis used interface-local geometry or complete water-cell Groups,
- VorPy version.

For solvent analyses, also report the operational definition used for an interface-associated water.

---

## Summary

A VorPy Interface is a dedicated spatial-partitioning Network constructed for the relationship between two molecular selections.

It separates:

- direct inter-group surfaces,
- Group-internal surfaces,
- supporting geometry,
- surface/contact/overlap metrics,
- curvature,
- interface-associated solvent.

The water-analysis extension additionally allows participating waters to be identified from retained interface topology and rebuilt as normal Groups so that their complete Voronoi cells can be analyzed.

Together, these features provide a geometric description of both the molecular boundary itself and the solvent structures associated with that boundary.
