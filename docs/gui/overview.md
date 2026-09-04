# Graphical User Interface

VorPy's official GUI is the PySide6/PyVista Analysis Studio in `vorpy.workbench`.

A typical workflow is:

1. Launch `vorpy` with no command-line arguments.
2. Create a project or open a PDB structure and inspect it in the molecular viewer.
3. Build named groups from the Structure selection browser.
4. Save the project as a versioned `.vpyworkbench.json` file. The project references the molecular file rather than embedding it and restores groups through stable atom identities.
5. Choose Atomic Voronoi, Power, or Primitive geometry and set the maximum vertex distance.
6. Run the analysis. The solve stays outside the Qt GUI thread and reports native VorPy progress.
7. Inspect Voronoi vertex and edge layers, atom or residue selections, and result metrics.
8. Adjust layer visibility, opacity, and color or save a high-resolution screenshot.

When a selection is active, Run VorPy snapshots its atom indices before starting the background job and passes them to VorPy’s concrete `Group`; with no selection, the run uses the complete structure.

Project commands are available from the File menu: New Project, Open Project, Save Project, and Save Project As. Unsaved projects are marked with an asterisk in the window title, and the Workbench prompts before discarding changes. Structure paths inside the project directory are stored relatively so a project folder can be moved as a unit.

The left workflow panel follows the Structure → Selection → Groups → Interfaces sequence. Solve, Analysis, and Results occupy the enlarged bottom tray as tabs. Small molecules (up to about 1,500 non-water/non-ion atoms) automatically enable sticks alongside the cartoon representation. Structure starts with a centered Browse action; additional structures appear in its loaded-structures list and can be switched without losing each molecule’s selections, groups, interfaces, rendering, or statistics. Visual and Layers remain on the right of the viewer for presentation and result inspection.

The Studio also opens completed VorPy output directories containing a PDB plus OFF, PLY, or VTP geometry. Command-line invocations with arguments retain the existing VorPy CLI for batch and scripted analysis.
