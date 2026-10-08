"""Export 2KAI molecular-contact visualization from a saved .vpy archive.

The OBJ files are render meshes only.  The JSON spherical patches, arcs, and
junctions remain the scientific source geometry.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vorpy.src.geometry.visualization import export_molecular_contact_bundle  # noqa: E402
from vorpy.src.io import load_session  # noqa: E402


def load_cached_surfaces(archive: Path) -> dict[str, object]:
    """Find the cached A/B surfaces without loading any raw molecular input."""
    session = load_session(archive)
    for system in getattr(session, "systems", ()):
        for interface in getattr(system, "ifaces", ()) or ():
            caches = getattr(interface, "representation_caches", {}) or {}
            keys = ("molecular_contact_surface:A", "molecular_contact_surface:B")
            if all(key in caches for key in keys):
                return {key: caches[key] for key in keys}
    raise ValueError("Archive has no cached molecular_contact_surface:A/B pair")


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="saved .vpy session")
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    parser.add_argument("--interaction", default="2KAI_AB")
    args = parser.parse_args(argv)
    manifest = export_molecular_contact_bundle(
        load_cached_surfaces(args.archive),
        args.output_dir,
        interaction=args.interaction,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


if __name__ == "__main__":
    main()
