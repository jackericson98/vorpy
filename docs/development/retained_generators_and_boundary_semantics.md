# Retained molecular generators and boundary semantics

FORGE keeps bulk solvent and unselected solvent ions out of the default
`MOLECULE_ONLY` generator universe. A caller may explicitly retain
System-global rows when a residue classified as solvent/ion is a molecular
component, for example:

```python
system = System(file="EDTA.pdb", retained_system_indices=(32,))
geometry = prepare_network_geometry(system, SolventMode.MOLECULE_ONLY)
```

The selection is deterministic, deduplicated, and validated against the
unchanged System row numbering. Retained rows are recorded separately from
`solvent_local_indices`; therefore EDTA Mg (system ID 32) can participate as a
molecular competitor without re-admitting the 579 water rows. The resulting
Network stores `retained_system_indices`, `retained_generator_ids`, and a
retention digest. These fields are included in solve-context membership
identity and version-3 network archives. Legacy archives load with an empty
retention selection.

## Supported entry point

The supported F14 entry point is the public `System` constructor. It accepts
System-global row IDs before any geometry is built:

```python
from vorpy.src.system import System

system = System(
    file="EDTA.pdb",
    make_dir=False,
    retained_system_indices=(32,),  # EDTA Mg; System-global row ID
)
```

The default remains `retained_system_indices=()`: molecule-only preparation
excludes unselected solvent and solvent ions. Retained IDs are explicit and
are not inferred from atom names, residue names, or serial numbers. The
current command-line parser has solvent-mode flags but no retained-ID flag;
the constructor/API is therefore the stable entry point until CLI ownership
and frame propagation are coordinated. Workbench controls are intentionally
out of scope for this producer milestone.

## Boundary semantics

`BoundaryMode.SHELL` adds deterministic virtual generators with stable
`boundary:<ordinal>` IDs. Boundary identity is independent of solvent
classification: virtual generators are never physical solvent, are not added
to the System atom table, and remain separately represented by
`boundary_indices`, `boundary_generators`, and `is_boundary_generator`.

Faces, areas, and volumes incident to a virtual shell are computational-domain
quantities. They depend on shell spacing, offset, and placement and must not be
reported as intrinsic physical observables without an independent boundary
convergence or physical-boundary argument. A shell can make an otherwise open
molecular cell complete, but it does not certify that the shell-truncated
volume is a physical molecular volume.

`BoundaryMode.NONE` leaves genuinely unbounded cells open. Explicit physical
solvent is controlled separately by `SOLVENT_CONTEXT` or `FULL_SOLVENT`; these
modes do not convert virtual shell geometry into physical solvent.

## Boundary-layer identity contract

The layer export has separate identity namespaces:

* `Group.ball_ndxs` and `Group.layer_atoms` contain non-negative System-global
  row IDs and are valid inputs to System/PDB access.
* `Network.balls["num"]`, `layer_net_atoms`, and topology incidence contain
  Network-local topology IDs. They may refer to physical or virtual boundary
  generators.
* Virtual generators have `generator_system_indices == None`, stable IDs such
  as `boundary:0`, and are identified separately by `boundary_indices` and
  `is_boundary_generator`.

`Network.balls["system_num"] == -1` is only the internal DataFrame encoding of
that missing System mapping. It is not a physical System ID and must never be
copied into `Group.layer_atoms`. A missing/invalid mapping is omitted from the
System layer and remains available only in the topology layer; archive
validation should therefore continue to reject `-1` in `group.layer_atoms`
while accepting valid topology IDs in `group.layer_net_atoms`.

### Consumer validation

VECTOR may validate a loaded group layer by requiring every `layer_atoms`
value to be in `range(len(network.sys.balls))`, every `layer_net_atoms` value
to be present in `network.balls["num"]`, and every boundary topology value to
be marked by `boundary_indices`/`is_boundary_generator`. A `-1` in
`layer_atoms` is malformed metadata, not a virtual-boundary reference; a
boundary reference is valid only in the Network topology namespace.
