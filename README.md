<h1 align="center" style="font-size: 4em; line-height: 1.2;">
  <img src="assets/VorpyIcon.svg" alt="VorPy Logo" width="160" height="160"
       style="vertical-align: -36px; margin-right: 10px;"/>
  VorPy
</h1>

<p align="center">
  <strong>Comprehensive Voronoi Analysis for Molecular and Geometric Systems</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Python-3.10%2B-blue"/>
  <img src="https://img.shields.io/badge/License-MIT-blue"/>
  <img src="https://img.shields.io/badge/Status-Active-success"/>
</p>

---

VorPy is a **comprehensive Voronoi-diagram analysis toolkit** for **scientific and molecular applications**.  
It supports **additively weighted**, **power**, and **primitive** Voronoi constructions from multiple input formats, and provides both a **command-line interface** and an intuitive **graphical UI** for interactive exploration.  
Built for performance and clarity, VorPy can produce **lightweight analytical outputs** or generate **rich 3D visualizations** for geometric inspection and publication-quality figures.

---

## 🚀 Quick Usage

Install from PyPI:
```bash
   pip install vorpy3
```

### Launch the GUI
Run VorPy with no arguments to open the graphical interface:
```bash
   python vorpy
```
or from within Python:
```python
import vorpy as vp
vp.run()
```

### Launch File Browser
Open a file selection dialog:
```bash
   python vorpy browse
```

### Run from Command Line
To solve directly from the command line, provide a file path:
```bash
   python vorpy path/to/file.pdb
```

VorPy automatically switches to **command-line mode** when a third argument is detected.  
Any additional arguments after the file adjust the default settings (see below).

---

## 🖥️ VorPy GUI

<p align="center">
  <img src="https://github.com/user-attachments/assets/43f33001-3ef3-427e-9293-eef0382d754c"
       alt="VorPy GUI" width="80%"/>
</p>

## 🖥️ VorPy GUI Overview

The VorPy graphical interface provides an intuitive workspace for configuring, solving, and visualizing Voronoi analyses.  
Each interface element serves a clear purpose, from managing input files to customizing groups and exports.

---

### 🔹 System-Level Controls

1. **Input File Name**  
   Displays the name of the selected input file (with folder path and extension removed).  
   This name becomes the root directory for all generated outputs, including subfolders for each defined group.

2. **Input File Information**  
   Shows general details about the input file, such as atom count, residue grouping, and other molecular metadata.  
   Helps confirm that the correct structure has been loaded before analysis.

3. **Input Locations**  
   Lists the current paths for input files, output directories, and any additional data sources (e.g., radii or log files).  
   Use this section to verify or modify file paths as needed.

4. **Radii and Mass Adjustments**  
   Opens a dedicated window for changing **atomic radii** and **masses**, either globally by element or selectively for individual atoms within specific residues.  
   (See [Radii/Mass Adjustments](#radii-mass-adjustments) for details.)

5. **System Exports**  
   Provides access to system-level export settings, such as **logs**, **ball files**, and **set-atom configurations**.  
   Launches a separate dialog where you can choose export types and formats.  
   (See [System](#system) for configuration options.)

6. **Reset**  
   Clears all loaded data, including input files, group settings, and system parameters.  
   Use this option to start fresh without restarting the application.

---

### 🔹 Group Management and Analysis

7. **Groups Section**  
   Displays all currently defined groups and their associated parameters.  
   Each group functions as an independent analysis unit with its own selection, build settings, and export options.  
   (See [Groups](#groups-gui) for more details.)

8. **Group Name**  
   Allows you to rename the active group.  
   The group name determines the name of its output subdirectory within the main results folder.

9. **Group Selections**  
   Lets you add or remove **atoms**, **balls**, **residues**, **chains**, or **molecules** by specifying indices or ranges.  
   If the range is left blank, only the single object corresponding to the entered index is selected.  
   Indices begin at **0** and correspond to the order of objects as they appear in the input file.

10. **Group Selection Tracker**  
    Displays a live record of all atoms or structures currently assigned to the active group.  
    Useful for confirming selections and avoiding duplicates.

11. **Group Build Settings**  
    Shows all computation parameters applied to the active group (e.g., network type, vertex count, box multiplier, surface resolution).  
    Settings can be modified interactively before running the solver.

12. **Group Exports**  
    Lists available export configurations for the current group.  
    Includes three standard presets — **small**, **medium**, and **large** — as well as a **custom exports** button that provides access to every available export type.

13. **Run Group**  
    Executes the Voronoi solver for the **active group only**, using its assigned selections and settings.  
    Exports the selected data to that group’s subdirectory.

14. **Add/Delete Groups**  
    Allows dynamic creation and removal of groups within the session.  
    The delete function applies to the current group and requires confirmation before removal to prevent accidental loss.

---

### 🔹 Global and Support Tools

15. **Help**  
    Opens the in-application help window, which explains each GUI element, available settings, and analysis options.  
    Acts as a quick reference for both new and experienced users.  
    (See [Help Window](#help-window).)

16. **Run All Groups**  
    Executes all defined groups sequentially, along with system-level exports.  
    This serves as the **main run function** for VorPy, automating the complete analysis workflow from input to export.

---

### 💡 Tip
You can run individual groups to test parameters or use **Run All Groups** for batch processing.  
System and group-level settings can be mixed freely — VorPy applies the most recent changes immediately and reflects them in the corresponding export directories.
---

## ⚙️ Command Line Interface

VorPy’s **command-line mode** uses the same executable as the GUI.  
It activates automatically when a file path or the keyword `browse` is supplied.

### Examples
```bash
   python vorpy                     # Launch GUI
   python vorpy browse              # Open file selection window
   python vorpy example.pdb         # Run CLI mode on example.pdb
   python vorpy example.pdb -s sr 0.1 and mv 80 -e small  # Override default parameters. 
```

### File Argument
The first argument after `vorpy` is the **path to the input file** (`.pdb`, `.mol`, `.gro`, `.cif`).

### Additional Arguments
Subsequent arguments modify the default settings.

#### Load Options (`-l`)
Load supplemental data files (vertices, logs, or index files):
```bash
   python vorpy example.pdb -l example_logs.csv
```

#### Setting Options (`-s`)
Adjust simulation or visualization parameters:
```bash
   python vorpy example.pdb -s nt pow and sr 0.05 and mv 80.0
```
- `nt` — Network Type: `aw`, `pow`, `prm`, or `com 'type1' 'type2'`
- `sr` — Surface Resolution (default: 0.2)
- `mv` — Maximum Vertex Count (default: 40)
- `bm` — Box Multiplier (default: 1.25)
- `ss` — Surface Scheme: curvature (`curv`), inside/outside (`nout`), or distance (`dist`)
- `sc` — Colormap: any [matplotlib colormap](https://matplotlib.org/stable/gallery/color/colormap_reference.html)

#### Group Options (`-g`)
Select atoms, residues, chains, or balls for analysis:
```bash
   python vorpy example.pdb -g a 0-100
```

#### Export Options (`-e`)
Choose export type and scope:
```bash
   python vorpy -e large and pdb and logs
```

---

## 🧪 Examples

The following examples demonstrate how to run VorPy from the command line.
To execute these commands, make sure your terminal’s working directory is set to the root directory for vorpy.
Each command can include multiple settings or flags in any order or combination. If two consecutive options contradict each other, VorPy will apply the last-specified command, ensuring that the most recent parameters take precedence.

---

### **Example 1 – Simple Molecular Visualization**
```bash
   python vorpy EDTA_Mg
```
This runs VorPy using the **built-in EDTA_Mg molecule**, which models the chelating agent EDTA bound to a magnesium ion.  
By default, VorPy computes the **additively weighted Voronoi diagram**, where each atom’s radius influences its spatial region.  
The results include geometric metrics (volume, surface area, sphericity) and a default visualization of the structure.

🧩 *Purpose:* Demonstrates the simplest possible use — one command, no extra flags, producing Voronoi partitions and visual outputs.

---

### **Example 2 – Change the Settings**
```bash
   python vorpy cambrin -s sr 0.05 and mv 80.0
```

This runs the **built-in Cambrin molecule**, but with modified computation settings:
- `sr 0.05` sets the **surface resolution** to 0.05 Å — very fine, meaning the surface will be calculated at high detail.  
- `mv 80` sets the **maximum vertex count** per Voronoi cell to 80 — an extremely high cap that ensures **all potential vertices** are evaluated (for most cases, values around 5–10 are sufficient).

🧩 *Purpose:* Demonstrates how to adjust **surface detail** and **computation depth** using the `-s` (settings) flag.

---

### **Example 3 – Compare Network Types**
```bash
   python vorpy EDTA -s nt compare prm pow -g mg
```
This compares two different Voronoi construction types — **primitive** and **power** — for the **Mg atom** in the built-in EDTA molecule:
- `nt compare prm pow` tells VorPy to compute **two networks**: one using the *primitive* Voronoi diagram and another using the *power* (weighted) version.  
- `-g mg` limits the analysis to the magnesium atom group.  

VorPy will generate a **comparison document** summarizing differences in geometric and topological metrics between the two network types.

🧩 *Purpose:* Highlights how VorPy can directly compare different geometric decomposition schemes on the same atomic subset.

---

### **Example 4 – Shell Visualization**
```bash
   python vorpy hairpin -s ss nout and sr 0.01 -e shell and pdb
```
This runs VorPy on the **hairpin molecule**, producing a **shell visualization**:
- `ss nout` colors the surface based on **inside vs. outside** curvature (a spatial separation scheme).  
- `sr 0.01` sets a **very high surface resolution**, ideal for fine structural mapping.  
- `-e shell and pdb` exports both a **shell visualization** and the corresponding **PDB structure file**.

The result is a highly detailed visual representation showing the molecular boundary and how it interacts with surrounding atoms.

🧩 *Purpose:* Demonstrates how to generate and export high-resolution, color-mapped surface models.

---

### ✅ **Summary of What These Teach**
| Example | Focus | Key Feature Demonstrated |
|----------|--------|---------------------------|
| 1 | Simple run | Default additively weighted Voronoi computation |
| 2 | Custom settings | Adjusting surface resolution and vertex limits |
| 3 | Network comparison | Comparing primitive vs. power diagrams |
| 4 | Visualization exports | Generating detailed color-mapped shell surfaces |

---

## 🧩 Requirements and Dependencies

- **Python 3.8+**
- **Dependencies:**
  - numpy  
  - scipy  
  - matplotlib  
  - pandas  
  - Pillow  
  - tkinter (for GUI)  
  - pytest (optional for tests)

Install dependencies via:
```bash
   pip install -r requirements.txt
```

---

## 🪪 License

This project is licensed under the [MIT License](LICENSE).

---

## 🧾 Citation

If you use VorPy in your research, please cite:

```bibtex
@software{vorpy2024,
  author = {John Ericson},
  title = {VORPY: A Python package for Voronoi analysis of molecular structures},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/jackericson98/vorpy}
}
```

---

## 📬 Contact

- **Email:** [jericson1@gsu.edu](mailto:jericson1@gsu.edu)  
- **Site:** [ericsonlabs.com](https://ericsonlabs.com)  
- **Phone:** +1 (404)-413-5491
