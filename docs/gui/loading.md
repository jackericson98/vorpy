# Loading Structures in the GUI

The GUI supports PDB, CIF, MOL, MOL2, and GRO input structures. After loading, the interface displays the selected structure name and molecular metadata derived from the file.

System-level controls include input structure information, file locations, system exports, radii/mass adjustments, and reset/reload behavior.

**TODO:** Add screenshots and exact field descriptions from the current GUI.

## New workbench: molecule-first PDB loading

The workbench displays the main molecule before preparing water and ions. Once the
molecule appears, the loading dialog closes and solvent loading continues in the
background, with progress in the status bar. You can rotate, zoom, select atoms,
and create groups during this stage. Solving and saving become available when the
complete system is ready.

Solvent completion preserves the camera, display settings, selections, and groups.
The final atom indices follow the original file order, including interleaved solvent.
If background loading fails, the previously loaded structure is restored.

The selection browser uses a virtual list so large solvent populations do not require
creating a separate widget item for every atom. All atoms remain searchable once loaded.
