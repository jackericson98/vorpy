# Selections

Selections define which molecular objects become members of a VorPy Group.

They provide the connection between a molecular structure—atoms, residues, chains, molecular classifications, and structural identifiers—and the system-level ball indices used by VorPy's geometric network calculations.

Selections can be created through either the command-line interface (CLI) or the graphical user interface (GUI). The two interfaces currently expose different selection capabilities.

---

## Selection and Group Are Different Concepts

A **Selection** identifies molecular objects.

A **Group** converts those selected objects into the set of atoms whose Voronoi cells belong to a network analysis.

Conceptually:

```text
molecular structure
      │
      ▼
   Selection
      │
      ├── atoms
      ├── residues
      ├── chains
      └── molecular categories
      │
      ▼
     Group
      │
      ▼
system ball indices
      │
      ▼
VorPy Network
```

Several different selection descriptions can therefore resolve to the same final Group atoms.

---

# Selection Concepts

VorPy currently uses three useful selection concepts.

## 1. Positional / index selection

Objects can be selected by their internal VorPy index.

Examples include:

```text
atom index 15
atoms 15-25

residue index 4
residues 4-8

chain index 0
```

Internal indices begin at **0**.

## 2. Molecular classification

The CLI can select broad molecular categories using VorPy's chemistry classification system.

Current category selectors include:

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

These selectors identify matching residues and use those residues to define Group membership.

## 3. Structural identity

The CLI can also select molecular objects by structural identity.

Examples include:

```text
THY 124
ALA 25 CA
```

These refer to residue/atom identity rather than to internal list position.

This distinction is important:

```text
residue index 124
```

is not necessarily the same object as:

```text
THY 124
```

The former is positional; the latter is structural.

---

# Command-Line Selections

CLI Group selection begins with the exact top-level flag:

```text
-g
```

Each new `-g` starts a new Group definition.

For example:

```text
python vorpy structure.pdb -g protein
```

defines one Group.

A second `-g` begins another Group:

```text
python vorpy structure.pdb -g protein -g dna
```

Conceptually:

```text
Group 0
    protein

Group 1
    dna
```

This distinction is especially important when two Groups are used to define an Interface.

---

## Combining Selections Within One Group

Multiple selection clauses can be combined within the same `-g` definition using a conjunction.

The canonical documented form is:

```text
and
```

For example:

```text
-g THY 124 and ALA 25
```

Conceptually:

```text
Group 0
    THY 124
    ALA 25
```

The command module currently recognizes several conjunction aliases internally:

```text
&
and
nd
also
+
&&
```

For documentation and reproducible command examples, `and` should be preferred.

---

## Top-Level Flags vs. Selection Aliases

VorPy contains many historical aliases for molecular-object names, but the current high-level command parser recognizes the top-level Group flag specifically as:

```text
-g
```

The object aliases apply **inside** the Group selection command.

For this reason, documentation should distinguish:

```text
-g
```

from selection keywords such as:

```text
atom
residue
chain
```

and should use a small set of canonical spellings rather than reproducing every historical typo/novelty alias.

---

# Canonical CLI Object Selectors

For scientific documentation, the following canonical names are recommended.

| Object | Canonical selector | Current alias family |
|---|---|---|
| Atom | `atom` | `a`, `as`, `atom`, `atoms`, `at`, `ats`, `am`, `ams` |
| Residue | `residue` | `r`, `rs`, `residue`, `residues`, `resid`, `resids`, `res`, `ress`, `reses`, `rdue`, `rdues` |
| Chain | `chain` | implemented through the combined chain/molecule alias family |
| Loaded index/group | `index` | `i`, `is`, `index`, `indexs`, `indexes`, `indices`, `ndx`, `ndxs`, `ndex`, `group`, `g`, `grp`, `n` |

The code also contains a combined chain/molecule alias family:

```text
m
ms
molecule
molecules
mols
ml
mls
c
cs
chain
```

However, in the currently reviewed Group-selection implementation, this entire alias family resolves through the **chain-selection path**.

Therefore `molecule` should **not** currently be documented as an independent working CLI molecular-object selector.

---

## Atom Selections

The CLI atom-selection helper supports several identifier forms.

### Atom internal index

Conceptually:

```text
-g atom 5
```

### Atom index range

Conceptually:

```text
-g atom 1-10
```

The range is inclusive.

### Atom name

An atom name can identify all matching atom names in the System.

Conceptually:

```text
-g atom CA
```

### Element

An element identifier can identify atoms of that element.

Conceptually:

```text
-g atom C
```

### Residue + sequence + atom

A specific structural atom can also be selected through residue identity:

```text
-g ALA 25 CA
```

This route identifies residue `ALA 25` and then the requested atom within it.

> The legacy helper for one residue/sequence/atom lookup contains a known implementation typo and is tracked for correction. A second active structural-selection path performs the same conceptual operation.

---

## Residue Selections

Residues can be selected in several ways.

### Residue internal index

Conceptually:

```text
-g residue 4
```

### Residue index range

Conceptually:

```text
-g residue 4-8
```

### Residue name + sequence number

Conceptually:

```text
-g THY 124
```

This searches for a residue whose name and sequence identifier match the request.

### Residue name

A residue name can also be used to collect matching residues of that type.

This is distinct from a broad molecular category such as `protein` or `dna`.

---

## Chain Selections

The current chain helper supports:

```text
single internal chain index
chain-index range
```

For example, conceptually:

```text
-g chain 0
```

or:

```text
-g chain 0-2
```

The reviewed helper does **not** currently resolve a structural chain identifier such as:

```text
chain A
```

That capability remains a planned improvement.

---

## Molecular Categories

Broad molecular-category selections are resolved separately from the atom/residue/chain object helpers.

Current categories include:

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

Matching residues are gathered from the System's ordinary and solvent residue collections as appropriate.

Examples:

```text
-g protein
```

```text
-g dna
```

```text
-g ligand
```

These are among the most convenient selectors for whole molecular classes.

---

## Additive Group Membership

Selection clauses accumulated within one Group are additive.

For example:

```text
-g protein and THY 124
```

conceptually combines both selection results.

During Group construction, selected atoms are converted to system ball indices and duplicate ball indices are suppressed.

This allows an atom reached through more than one selection route to remain a single Group cell.

---

# GUI Selections

The current GUI Selection panel exposes four selectable object classes:

```text
Atoms/Balls
Residues
Chains
Molecules
```

The GUI panel is explicitly index-based and displays:

```text
Indices start at 0
```

For each object class, the user can enter either:

- one index, or
- a start and end index.

---

## Image Placeholder 1 — GUI Selection Panel

![Placeholder: VorPy GUI selection panel](../../assets/docs/selections/gui_selection_panel.png)

> **Suggested figure:** Screenshot the current Group Selection panel showing the object dropdown, Index and Range entries, Add/Remove/Clear controls, and Selection Tracker. Highlight the **Indices start at 0** note.

---

## Inclusive GUI Ranges

GUI ranges include both endpoints.

For example:

```text
start = 4
end = 8
```

produces:

```text
4, 5, 6, 7, 8
```

The same range-building logic is used for all four displayed object classes.

Negative values are rejected by the Add path because the entry must contain digits and the converted value must not be negative.

---

## Selection Tracker

The GUI stores current selection indices separately for:

- balls,
- residues,
- chains,
- molecules.

The Selection Tracker compresses consecutive indices into human-readable ranges.

For example:

```text
0, 1, 2, 3, 7, 8, 12
```

is displayed as:

```text
0-3, 7-8, 12
```

The tracker is a display representation of the stored index lists rather than a separate molecular-selection grammar.

---

## Adding, Removing, and Clearing

The Selection panel provides:

- **Add**
- **Remove**
- **Clear**

Adding duplicate indices does not intentionally create duplicate stored selections.

Removing a range removes each index in that inclusive range.

Clear resets all four displayed selection collections.

An internal undo implementation exists, but the Undo button is currently commented out and is therefore not a documented user-facing feature.

---

# How GUI Selections Leave the Panel

Each Group tab owns its own `SelectionFrame`.

The parent Groups panel stores that frame alongside the Group's:

- build settings,
- export settings,
- name,
- optional logs file.

When Group settings are requested, the parent panel returns:

```text
name
build_settings
export_settings
selections
```

where `selections` is the SelectionFrame's current selection dictionary.

Therefore the reviewed GUI path confirms:

```text
SelectionFrame
      │
      ▼
GroupsFrame group settings
      │
      ▼
GUI run_group(...)
```

The final conversion is performed directly by the GUI's `run_group()` method. The four SelectionFrame collections are passed to the `Group` constructor as:

```text
balls      → atoms
residues   → residues
chains     → chains
molecules  → molecules
```

The GUI therefore preserves the object class selected in the Selection panel when it constructs a Group.

The build settings are passed at the same time, including surface resolution, box size, maximum vertex size, network type, surface coloring, and vertex/edge colors. The resulting Group is then appended to `System.groups` if it is not already present.

This closes the GUI selection path:

```text
SelectionFrame
      ↓
GroupsFrame
      ↓
VorPyGUI.run_group()
      ↓
Group(...)
      ↓
Group.process_inputs()
      ↓
ball_ndxs
```

For atoms, residues, and chains, this path is now source-validated. For molecules, the GUI handoff reaches `Group(molecules=...)` correctly, but the downstream Group expansion bug remains.


---

## Loaded Logs and GUI Selections

When a Group is reconstructed from logs, the Groups panel repopulates the Selection Tracker using the Group's logged ball indices.

The reconstructed GUI selection state becomes:

```text
balls      = logged Group ball indices
residues   = None
chains     = None
molecules  = None
```

This means log-loaded Group membership is represented in the GUI primarily as explicit ball-index membership rather than attempting to reconstruct the original higher-level residue/chain selection expression.

That is an important distinction between **recovered Group membership** and **recovered original selection intent**.

---

## Image Placeholder 2 — Selection Flow Through the GUI

![Placeholder: GUI selection data flow](../../assets/docs/selections/gui_selection_flow.png)

> **Suggested figure:** Show `SelectionFrame` on the left, the stored dictionary in the middle (`balls`, `residues`, `chains`, `molecules`), and a Group/Network on the right. Add a separate arrow from a logs file back into explicit ball indices.

---

# Current CLI and GUI Capability Matrix

| Selection concept | CLI | GUI Selection panel |
|---|---|---|
| Atom/ball internal index | Supported | Supported |
| Atom/ball index range | Supported | Supported |
| Atom name | Supported | Not exposed |
| Element-based atom selection | Supported | Not exposed |
| Residue internal index | Supported | Supported |
| Residue index range | Supported | Supported |
| Residue name | Supported | Not exposed |
| Residue name + sequence | Supported | Not exposed |
| Residue + sequence + atom name | Supported, with one legacy helper bug tracked | Not exposed |
| Chain internal index | Supported | Supported |
| Chain index range | Supported | Supported |
| Structural chain ID (`A`) | Not currently implemented | Not exposed |
| Independent molecule-object selection | Current CLI alias maps to chain path | GUI passes molecule indices to `Group(molecules=...)`, but Group molecule expansion currently needs correction |
| Protein category | Supported | Not exposed |
| DNA category | Supported | Not exposed |
| RNA category | Supported | Not exposed |
| Ligand category | Supported | Not exposed |
| Water/solvent category | Supported | Not exposed |
| Ion category | Supported | Not exposed |

The GUI and CLI should therefore not currently be described as equivalent selection interfaces.

---

# Internal Index vs. Structural Identifier

The most important selection distinction for reproducibility is:

```text
internal index
```

versus:

```text
identifier from the molecular structure
```

For example:

```text
residue index 124
```

means the residue at internal position 124.

By contrast:

```text
THY 124
```

means a thymine residue whose sequence identifier is 124.

The two may refer to completely different objects.

Similarly, a future:

```text
chain A
```

selector should refer to molecular chain identity rather than merely internal chain-list position.

---

## Selection Reproducibility

For scientific reporting, structural selections should generally be preferred when they are available and unambiguous.

Examples:

```text
THY 124
protein
DNA
ligand
```

communicate molecular meaning directly.

Internal indices remain valuable for:

- scripting,
- debugging,
- generated selections,
- exact access to VorPy's object collections,
- log-reconstructed Group membership.

When reporting an index-based selection, retain the input structure and software version so the selection can be reconstructed.

---

# Known Limitations and Active Improvements

## Structural chain identifiers

The current chain-selection helper supports chain indices and ranges but not chain IDs such as `A`.

## Independent molecule selection

The GUI exposes a Molecules category, and `Group` contains a molecule-selection field.

However:

1. the current CLI `molecule` aliases resolve through the chain path, and
2. the reviewed Group input-processing code does not yet expand `self.mols` into Group atoms.

Molecule selection should therefore be considered incomplete until both paths are reconciled and tested.

## GUI molecular identity

The current Selection panel does not expose protein/DNA/RNA/ligand categories, named residues, atom names, elements, or chain IDs.

## Final GUI Group conversion

The complete GUI handoff is now validated: `SelectionFrame` stores the selected indices, `GroupsFrame` retains the SelectionFrame for that Group tab, and `VorPyGUI.run_group()` passes each selection collection directly into the corresponding `Group` constructor argument.

---

# Recommended Long-Term Selection Model

The user-facing selection grammar should make three concepts explicit.

### Positional

```text
atom 20
residue 15
chain 0
```

### Molecular classification

```text
protein
dna
rna
ligand
water
ions
```

### Structural identity

```text
chain A
THY 124
ALA 25 CA
```

The same conceptual model should eventually be available from both CLI and GUI.

---

## Image Placeholder 3 — Selection Types

![Placeholder: positional, classification, and structural selections](../../assets/docs/selections/selection_types.png)

> **Suggested figure:** Three columns labeled **Positional / Index**, **Molecular Classification**, and **Structural Identity**. Show examples such as `atom 25`, `protein`, and `THY 124`, with arrows resolving each selection to system atom indices and then a VorPy Group.

---

# Relationship to Groups and Interfaces

Selections define molecular membership; they do not perform geometric calculations.

```text
Selection
    ↓
Group
    ↓
Network
```

For interfaces:

```text
Selection A       Selection B
    ↓                 ↓
 Group A           Group B
      \             /
       \           /
         Interface
```

See [Groups](groups.md) for Group construction and [Interfaces](interfaces.md) for pairwise molecular interfaces.

---

# Summary

VorPy selections translate molecular intent into Group membership.

The current CLI supports:

- atom indices and ranges,
- atom names and elements,
- residue indices and ranges,
- residue identities,
- residue-specific atom identities,
- chain indices and ranges,
- broad molecular categories,
- additive clauses within one Group using `and`,
- separate Groups using repeated `-g`.

The current GUI Selection panel supports zero-based index/range selections for:

- Atoms/Balls,
- Residues,
- Chains,
- Molecules.

The two interfaces are therefore not yet feature-equivalent.

For future development, VorPy should preserve a clear distinction between **positional selection**, **molecular classification**, and **structural identity**, while making all three available consistently across CLI and GUI.
