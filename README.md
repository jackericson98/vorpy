# VorPy

![VorPy Logo](assets/VorpyLogo.svg)

**Quantitative molecular geometry through three-dimensional spatial partitioning**

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![License](https://img.shields.io/badge/License-MIT-blue)
![Status](https://img.shields.io/badge/Status-Active-success)

---

## Overview

**VorPy** is a Python framework for constructing and analyzing three-dimensional Voronoi partitions of molecular and other spherical systems.

VorPy represents atoms or other objects as spheres and partitions space using three complementary geometric schemes:

- **Additively weighted (AW)** Voronoi diagrams
- **Power (Laguerre)** diagrams
- **Primitive (unweighted)** Voronoi diagrams

From these decompositions, VorPy can quantify **volume, surface area, neighbors, contacts, interfaces, exposure, surface geometry, curvature, and related molecular properties**. Analyses can be performed at the level of atoms, residues, chains, molecules, user-defined groups, and interfaces between groups.

VorPy provides both a **graphical user interface (GUI)** for interactive analysis and a **command-line interface (CLI)** for reproducible, automated, and large-scale workflows.

> Molecular structure is fundamentally geometric. Changes in packing, exposure, local contacts, and intermolecular interfaces alter how atoms occupy and share space. VorPy provides a common spatial framework for quantifying these relationships.

---

## Installation

### PyPI

```bash
pip install vorpy3
```

### Conda

```bash
conda install vorpy3
```

VorPy currently targets **Python 3.10+**. Explicit multi-version support will be listed here once the test suite and end-to-end regression suite have been validated across supported interpreters and platforms.

---

## Quick Start

Launch the graphical interface:

```bash
python vorpy
```

Open the file browser:

```bash
python vorpy browse
```

Analyze a structure directly:

```bash
python vorpy path/to/structure.pdb
```

VorPy currently supports `.pdb`, `.cif`, `.mol`, `.mol2`, and `.gro` files.

---

## Spatial Partitioning Schemes

### Additively weighted Voronoi (`aw`)

Accounts directly for sphere radii through distance to sphere surfaces, making it well suited to radius-aware molecular geometry.

### Power / Laguerre (`pow`)

Incorporates sphere size through power distance while retaining planar cell boundaries.

### Primitive (`prm`)

Depends only on generator positions. Sphere radii do not affect boundary placement, making it useful as an unweighted geometric reference.

VorPy supports direct comparison between these schemes so the effects of the chosen spatial model can be evaluated explicitly.

---

## Molecular Geometry

Depending on the requested analysis and exports, VorPy can provide information about:

- Voronoi cell volume
- molecular and cell surface area
- neighboring atoms and cells
- contact topology
- interfacial surface area
- solvent-facing exposure
- cell shape descriptors
- local and integrated curvature
- triangulated surfaces
- edges and vertices
- atom, residue, chain, group, and system summaries

---

## Interface Analysis

Interface analysis is one of VorPy's central molecular-analysis capabilities.

Groups of atoms can be defined independently and their shared Voronoi boundaries analyzed as explicit geometric interfaces. This supports protein-protein, protein-DNA, protein-RNA, protein-ligand, solute-solvent, chain-chain, residue-residue, and arbitrary user-defined interfaces.

Because interfaces are derived from shared Voronoi boundaries, VorPy can describe not only **which atoms are neighbors**, but also the **geometry and area of the boundary shared between them**.

See [`docs/molecular_analysis/interfaces.md`](docs/molecular_analysis/interfaces.md).

---

## Groups and Selections

VorPy analyses can operate on the full structure or on selected subsets.

Selections can be based on:

- atom indices
- residue indices
- chain indices
- molecule indices
- molecular classes such as **protein**, **ligand**, **DNA**, and **RNA**
- structural identities such as a specific chain or residue

Examples of increasingly specific selections include:

```text
protein
dna
chain A
THY 124
```

The exact CLI syntax for identity-based selectors is being expanded and will be documented in the full selection reference.

See [`docs/molecular_analysis/selections.md`](docs/molecular_analysis/selections.md).

---

## Curvature and Representative Surface Energy

VorPy can calculate local and integrated curvature quantities, including `int_mean_curv`, `int_mean_curv_sq`, `int_gauss_curv`, and `surf_energy`.

The current `surf_energy` quantity is a **representative curvature-dependent bending metric**, not a calibrated molecular free energy. It does not by itself represent solvent thermodynamics, electrostatics, dispersion, conformational entropy, or a system-specific fitted bending modulus.

See [`docs/theory/curvature.md`](docs/theory/curvature.md) and [`docs/theory/surface_energy.md`](docs/theory/surface_energy.md).

---

## Exports

VorPy supports both **export presets** and **individual export components**.

Preset bundles include:

```text
tiny
small
medium
large
all
```

Additional output types can be requested independently or appended to a preset:

```bash
python vorpy example.pdb -e small
python vorpy example.pdb -e small and shell
```

See [`docs/cli/exports.md`](docs/cli/exports.md).

---

## Graphical User Interface

The GUI provides an interactive environment for loading structures, defining groups, changing radii and masses, configuring build settings, selecting exports, running individual groups or complete systems, and working with interfaces and visualization output.

The detailed GUI guide is maintained in [`docs/gui/`](docs/gui/).

---

## Command-Line Interface

Basic syntax:

```bash
python vorpy <input_file> [options]
```

Examples:

```bash
python vorpy example.pdb
python vorpy example.pdb -s nt pow
python vorpy example.pdb -g a 0-100
python vorpy example.pdb -e small and shell
```

Full CLI documentation is available in [`docs/cli/`](docs/cli/).

---

## Documentation

The full documentation is organized by topic:

- [Installation](docs/installation.md)
- [Quick Start](docs/quick_start.md)
- [Spatial Partitioning Theory](docs/theory/)
- [Molecular Analysis](docs/molecular_analysis/)
- [Command-Line Interface](docs/cli/)
- [Graphical User Interface](docs/gui/)
- [Input and Output Formats](docs/formats/)
- [Development and Testing](docs/development/)
- [Terminology and References](docs/reference/)

The documentation is intentionally stored as ordinary Markdown so it can be read directly in GitHub, PyCharm, and other editors. It can later serve as the source for a dedicated documentation website.

---

## Validation and Python Compatibility

VorPy currently targets Python 3.10+. Before individual Python versions are advertised as explicitly supported, they should pass unit tests, numerical regression tests, parser/export integration tests, and representative end-to-end molecular calculations.

See [`docs/development/testing.md`](docs/development/testing.md).

---

## Scientific Reference

**Ericson, J. M.; Wolpert, N.; Poon, G. M. K.**  
*Evaluation of weighted Voronoi decompositions of physicochemical ensembles.*  
**Physical Chemistry Chemical Physics** 2025, **27**, 16204-16218.  
DOI: `10.1039/D5CP00763A`

---

## Citation

```bibtex
@article{Ericson2025WeightedVoronoi,
  author  = {Ericson, John M. and Wolpert, Nicola and Poon, Gregory M. K.},
  title   = {Evaluation of weighted Voronoi decompositions of physicochemical ensembles},
  journal = {Physical Chemistry Chemical Physics},
  year    = {2025},
  volume  = {27},
  pages   = {16204--16218},
  doi     = {10.1039/D5CP00763A}
}
```

---

## License

VorPy is released under the **MIT License**. See [`LICENSE`](LICENSE).

---

## Project Status

VorPy is under active development. For publication-grade analyses, record the VorPy version, Python version, partitioning scheme, atomic radii, relevant build settings, and export settings used for each calculation.
