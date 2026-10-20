"""Small, deterministic helpers for the on-disk output layout."""

import hashlib
import re
import shutil
from pathlib import Path


def slug(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._") or "network"


def interface_folder_name(group1, group2):
    first = slug(getattr(group1, "name", group1) or "group_1")
    second = slug(getattr(group2, "name", group2) or "surrounding")
    return f"interface_{first}_{second}"


def unique_output_directory(root, name, occupied=()):
    """Return a non-colliding child path without creating it."""
    root = Path(root)
    occupied = {Path(path).resolve() for path in occupied}
    candidate = root / slug(name)
    suffix = 2
    while candidate.resolve() in occupied or candidate.exists():
        candidate = root / f"{slug(name)}_{suffix}"
        suffix += 1
    return candidate


def canonical_layout_enabled(system):
    """Return whether the normalized export path owns the canonical layout."""
    return bool(
        getattr(system, "_canonical_output_layout", False)
        or (
            getattr(system, "_export_plan", None) is not None
            and getattr(system, "balls", None) is not None
        )
    )


def canonical_root(system):
    return Path(system.files["dir"])


def canonical_system_directory(system):
    return canonical_root(system) / "system"


def canonical_groups_directory(system):
    return canonical_root(system) / "groups"


def canonical_interfaces_directory(system):
    return canonical_root(system) / "interfaces"


def canonical_group_directory(system, group):
    group_id = getattr(group, "group_id", None) or getattr(group, "name", "group")
    return canonical_groups_directory(system) / slug(group_id)


def canonical_interface_directory(system, interface):
    interface_id = getattr(interface, "interface_id", None) or getattr(interface, "name", "interface")
    return canonical_interfaces_directory(system) / slug(interface_id)


def copy_source_once(system):
    """Copy the original input bytes once, refusing a different collision."""
    source = (getattr(system, "files", {}) or {}).get("base_file")
    if not source:
        return None
    source = Path(source)
    if not source.is_file():
        return None
    destination = canonical_root(system) / source.name
    if destination.resolve() == source.resolve():
        return destination
    if destination.exists():
        if destination.stat().st_size != source.stat().st_size:
            raise FileExistsError(f"Refusing to overwrite unrelated source copy: {destination}")
        if hashlib.sha256(destination.read_bytes()).digest() != hashlib.sha256(source.read_bytes()).digest():
            raise FileExistsError(f"Refusing to overwrite unrelated source copy: {destination}")
        return destination
    shutil.copyfile(source, destination)
    return destination


def prepare_canonical_root(system):
    """Create the shared root and source copy used by normalized exports."""
    root = canonical_root(system)
    root.mkdir(parents=True, exist_ok=True)
    system._canonical_output_layout = True
    copy_source_once(system)
    canonical_system_directory(system).mkdir(parents=True, exist_ok=True)
    canonical_groups_directory(system).mkdir(parents=True, exist_ok=True)
    canonical_interfaces_directory(system).mkdir(parents=True, exist_ok=True)
    for interface in getattr(system, "ifaces", ()) or ():
        directory = canonical_interface_directory(system, interface)
        directory.mkdir(parents=True, exist_ok=True)
        interface.dir = str(directory)
    return root
