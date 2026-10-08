# Representation cache contract

This is a persistence extension point, not a geometry implementation.

The current reader rejects unknown object kinds and fields. Therefore v2 is not
an open extension format. New representation-cache fields are reserved by
archive format version 3; the reader still accepts old v1 and v2 archives.
An archive labelled v2 that already contains the reserved representation fields
is rejected rather than silently treated as a v2 document.
Future adapter kinds must be registered before archives containing them are
read.

## Cache layout

An interface may later expose two optional mappings:

```text
representation_caches[key] -> registered cached adapter
representation_links[key]  -> IDs and provenance only
```

The key is a stable representation name, for example:

```text
physical_interface
dual_contact_complex
generator_complex
molecular_contact_surface:A
molecular_contact_surface:B
dual_separator
```

The codec stores registered adapter objects. It does not require every adapter
to contain surfaces, area, or curvature. Link records must use archive IDs, not
Python object IDs.

## Adapter registration

Agent #3's final dataclasses should be added to
`vorpy/src/io/scientific_adapters.py`'s `registry()` with a unique adapter kind.
Nested custom records need their own registry entries. Adapter fields must be
explicit and codec-supported; arbitrary constructors, imports, and pickle are
not allowed. Include an `archive_id` and `archive_provenance` when the record
has independent scientific identity. Network/interface ownership is checked
from those stable provenance IDs.

## Molecular contact surface envelope

Required future metadata:

- definition and schema version
- side (`A` or `B`)
- molecular surface convention
- system and frame identity
- source network and interface IDs
- selected-contact query ID
- contact-selection provenance
- radius model and probe radius
- canonical environment (`dry` or `solvent_competing`)
- source atom stable IDs
- support/status and provenance

Pending geometry fields, exact shape owned by Agent #3:

- spherical patch records
- circular boundary arcs
- junction vertices
- connected-component IDs
- atom/contact provenance on geometric records
- area and topology summaries
- mean/Gaussian curvature summaries
- unresolved/support diagnostics

## Capability and dependencies

Capability reads `representation_links` and cache metadata only. It never
constructs geometry. A link may carry dependency states such as:

```text
atom_centers_radii
stable_atom_ids
selected_ab_contacts
molecular_boundary_convention
spherical_clipping
```

This permits messages such as “selected contacts available” or “clipping
geometry unavailable” while the representation itself is `NOT_CALCULATED`.

## Round trip

After a real molecular contact surface is cached:

```text
save -> fresh process -> load -> inspect/export
```

must use only the cached adapter and links. It must not solve a network, select
contacts, clip spheres, rebuild arrangements, or integrate curvature.

Old archives without these optional mappings load normally. They report the
requested molecular-contact representation as `NOT_CALCULATED`; no cache is
synthesized during load.

## Cache-only visualization

`export_molecular_contact_surface()` and
`export_molecular_contact_bundle()` consume loaded canonical records only:

- canonical: `SphericalPatch`, `CircularArc`, `JunctionVertex`, and available
  `RefinedSeamPiece` records
- derived: rendering triangles, OBJ boundary/junction files, CSV mappings,
  and visualization manifests

The renderer does not reload OBJ/OFF output as scientific state. Missing
representations are reported as `NOT_CALCULATED`; partial topology does not
prevent rendering known patches and arcs. Bundle paths are centralized under
`interfaces/<interaction>/` with independent A/B representation directories
and a `visualization/` manifest directory.
