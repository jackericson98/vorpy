"""Build and export the dry-2KAI Power molecular contact-surface prototype."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vorpy.src.geometry.molecular_contact import (  # noqa: E402
    AtomBall,
    SelectedContact,
    build_molecular_contact_surface_pair,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "dual_side_incidence_audit" / "2kai_power" / "molecular_contact_surface"


def stable_id(row: dict[str, str]) -> str:
    return (
        f"{row['index']}|pdb={row['pdb_serial']}|group={row['group']}|"
        f"{row['chain']}:{row['residue_name']}{row['residue_number']}:{row['atom_name']}"
    )


def main() -> None:
    with (ROOT / "2kai_cazals_atoms.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    atoms = tuple(
        AtomBall(
            int(row["index"]),
            stable_id(row),
            row["group"],
            (float(row["x"]), float(row["y"]), float(row["z"])),
            float(row["expanded_radius_A"]),
        )
        for row in rows
    )
    with (ROOT / "2kai_alpha0_edges.csv").open(newline="", encoding="utf-8") as stream:
        selected_rows = list(csv.DictReader(stream))
    selected_pairs = {
        tuple(sorted((int(row["facet_1_i"]), int(row["facet_1_j"]))))
        for row in selected_rows
    }
    selected_pairs |= {
        tuple(sorted((int(row["facet_2_i"]), int(row["facet_2_j"]))))
        for row in selected_rows
    }
    contacts = tuple(
        SelectedContact(
            f"power-alpha0:{a}:{b}",
            a if next(atom for atom in atoms if atom.generator_index == a).group == "A" else b,
            b if next(atom for atom in atoms if atom.generator_index == a).group == "A" else a,
            "dry_2kai_power",
            "c_0_and_c_1_c_2_interface",
        )
        for a, b in sorted(selected_pairs)
    )
    side_a, side_b = build_molecular_contact_surface_pair(
        atoms,
        contacts,
        source_network_id="dry_2kai_power",
        source_interface_id="c_0_and_c_1_c_2_interface",
        contact_selection_source="power_alpha0_interface",
        contact_selection_id="dry_2kai_power:selected_power_alpha0_contacts",
    )
    OUT.mkdir(parents=True, exist_ok=True)
    for label, surface in (("A", side_a), ("B", side_b)):
        (OUT / f"surface_{label}.json").write_text(
            json.dumps(
                {
                    "summary": surface.summary(),
                    "patches": [asdict(value) for value in surface.patches],
                    "arcs": [asdict(value) for value in surface.arcs],
                    "junctions": [asdict(value) for value in surface.junctions],
                    "refined_seams": [asdict(value) for value in surface.refined_seams],
                    "topological_cuts": [asdict(value) for value in surface.topological_cuts],
                    "curvature": asdict(surface.curvature) if surface.curvature else None,
                },
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
    summary = {
        "input_atoms": {"A": sum(atom.group == "A" for atom in atoms), "B": sum(atom.group == "B" for atom in atoms)},
        "selected_contacts": len(contacts),
        "direct_contact_atoms": {
            "A": len({contact.atom_a_index for contact in contacts}),
            "B": len({contact.atom_b_index for contact in contacts}),
        },
        "A": side_a.summary(),
        "B": side_b.summary(),
        "canonical_geometry": "spherical patches with circular-arc boundaries; no mesh is authoritative",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
