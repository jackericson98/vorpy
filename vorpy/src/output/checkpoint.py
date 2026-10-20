"""VECTOR filesystem adapter for the ATLAS checkpoint contract."""

from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import tempfile
import zipfile

from vorpy.src.io import save_network
from vorpy.src.results.checkpoint import (
    CheckpointAtom,
    CheckpointChain,
    CheckpointContainer,
    CheckpointDocument,
    CheckpointGroup,
    CheckpointInterfaceIdentity,
    CheckpointNetworkState,
    CheckpointPayload,
    CheckpointResidue,
    CheckpointResultState,
    CheckpointSource,
    load_checkpoint_container,
    make_checkpoint_dependency,
    prepare_checkpoint_update,
)
from vorpy.src.results.model import ProvenanceCatalog, Status


class CheckpointSourceUnavailable(RuntimeError):
    """The archive is usable, but exact original source bytes are unavailable."""

    def __init__(self, *, filename=None, source_path=None, source_sha256=None):
        self.filename = filename
        self.source_path = source_path
        self.source_sha256 = source_sha256
        label = filename or source_path or "the original source"
        super().__init__(
            "Exact original source bytes are unavailable for canonical checkpoint "
            f"emission: {label}"
        )


def _json_safe(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    return str(value)


def _row_value(row, *names, default=""):
    for name in names:
        try:
            value = row.get(name)
        except AttributeError:
            value = None
        if value is not None:
            try:
                if value != value:  # NaN
                    continue
            except Exception:
                pass
            return value
    return default


def _int_value(value, default):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _network_ids(network):
    values = []
    balls = getattr(network, "balls", None)
    if balls is None or not hasattr(balls, "iterrows"):
        return ()
    for position, (_, row) in enumerate(balls.iterrows()):
        value = _row_value(row, "stable_id", "generator_stable_id", "num", default=None)
        values.append(str(value) if value is not None else f"generator:{position}")
    return tuple(dict.fromkeys(values))


def _network_payload(network):
    """Reuse the existing v3 serializer as one opaque, lossless payload."""
    with tempfile.TemporaryDirectory(prefix="vorpy-checkpoint-") as directory:
        path = Path(directory) / "network.vpy"
        try:
            save_network(network, path, write_timing_manifest=False)
            return "network-v3", path.read_bytes()
        except Exception:
            # Synthetic/source-view networks can intentionally have system
            # indices that are not valid v3 archive references. Preserve the
            # available tables as an explicit partial payload instead of
            # silently dropping checkpoint geometry.
            tables = {}
            for name in ("balls", "edges", "verts", "surfs"):
                table = getattr(network, name, None)
                tables[name] = table.to_json(orient="split", default_handler=str) if table is not None else None
            payload = json.dumps({
                "format": "vorpy.checkpoint.network-geometry",
                "settings": _json_safe(getattr(network, "settings", {}) or {}),
                "tables": tables,
            }, sort_keys=True, separators=(",", ":")).encode("utf-8")
            return "network-geometry", payload


def _network_status(network, payload):
    return Status.PARTIAL if payload else Status.UNRESOLVED


def _archive_source_metadata(system):
    metadata = getattr(system, "archive_metadata", None)
    if not isinstance(metadata, dict):
        metadata = {}
    source_path = metadata.get("original_structure_path")
    filename = metadata.get("original_structure_filename")
    if not filename and source_path:
        filename = Path(source_path).name
    return {
        "path": str(source_path) if source_path else None,
        "filename": str(filename) if filename else None,
        "sha256": metadata.get("source_sha256"),
    }


def describe_checkpoint_source(system):
    """Describe exact-source availability without reconstructing source bytes."""
    metadata = _archive_source_metadata(system)
    source_path = (getattr(system, "files", {}) or {}).get("base_file")
    if source_path and Path(source_path).is_file():
        digest = hashlib.sha256(Path(source_path).read_bytes()).hexdigest()
        return {
            "available": True,
            "origin": "filesystem",
            "path": str(Path(source_path).resolve()),
            "filename": Path(source_path).name,
            "sha256": digest,
            "recorded_sha256": metadata["sha256"],
        }
    document = getattr(system, "_checkpoint_document", None)
    source = getattr(document, "source", None)
    if source is not None:
        return {
            "available": True,
            "origin": "checkpoint",
            "path": None,
            "filename": source.filename,
            "sha256": source.sha256,
            "recorded_sha256": metadata["sha256"],
        }
    return {
        "available": False,
        "origin": None,
        "path": metadata["path"],
        "filename": metadata["filename"],
        "sha256": None,
        "recorded_sha256": metadata["sha256"],
    }


def write_checkpoint_status(system, *, reason, error=None):
    """Record an honest source/checkpoint outcome beside canonical outputs."""
    from vorpy.src.output.layout import canonical_system_directory

    directory = canonical_system_directory(system)
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / "checkpoint_status.json"
    source = describe_checkpoint_source(system)
    payload = {
        "format": "vorpy.checkpoint_export_status",
        "version": 1,
        "status": "NOT_WRITTEN",
        "reason": reason,
        "source": {
            "exact_bytes_available": source["available"],
            "origin": source["origin"],
            "filename": source["filename"],
            "path": source["path"],
            "sha256": source["sha256"],
            "original_recorded_sha256": source["recorded_sha256"],
        },
        "scientific_archive_reloadable": bool(
            getattr(system, "_loaded_from_archive", False)
        ),
        "checkpoint": {
            "written": False,
            "path": str(
                (Path(system.files["dir"]) / f"{system.name}_checkpoint.vpy").resolve()
            ),
            "existing_preserved": bool(
                (Path(system.files["dir"]) / f"{system.name}_checkpoint.vpy").is_file()
            ),
        },
    }
    if error is not None:
        payload["error"] = str(error)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=directory,
            mode="w",
            encoding="utf-8",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    system._checkpoint_export_status = payload
    return destination
def _checkpoint_source(system):
    metadata = _archive_source_metadata(system)
    source_path = (getattr(system, "files", {}) or {}).get("base_file")
    if source_path and Path(source_path).is_file():
        source_path = Path(source_path)
        source_bytes = source_path.read_bytes()
        source = CheckpointSource(source_path.name, source_bytes)
        recorded = metadata["sha256"]
        if recorded and source.sha256 != recorded:
            raise ValueError(
                "Original source hash does not match the archived identity: "
                f"recorded={recorded!r}, actual={source.sha256!r}"
            )
        return source
    document = getattr(system, "_checkpoint_document", None)
    source = getattr(document, "source", None)
    if source is not None:
        recorded = metadata["sha256"]
        if recorded and source.sha256 != recorded:
            raise ValueError(
                "Embedded checkpoint source hash does not match the archived identity: "
                f"recorded={recorded!r}, actual={source.sha256!r}"
            )
        return source
    raise CheckpointSourceUnavailable(
        filename=metadata["filename"],
        source_path=metadata["path"],
        source_sha256=metadata["sha256"],
    )


def _owner_networks(system):
    owners = []
    for group in getattr(system, "groups", ()) or ():
        network = getattr(group, "net", None) or getattr(group, "source_network", None)
        if network is not None:
            owners.append(("group", getattr(group, "group_id", None) or getattr(group, "name", "group"), network))
    for interface in getattr(system, "ifaces", ()) or ():
        network = getattr(interface, "net", None) or getattr(interface, "source_network", None)
        if network is not None:
            owners.append(("interface", getattr(interface, "interface_id", None) or getattr(interface, "name", "interface"), network))
    return owners


def _table_len(network, name):
    table = getattr(network, name, None)
    return len(table) if table is not None else 0


def build_checkpoint_container(system):
    source = _checkpoint_source(system)
    balls = getattr(system, "balls", None)
    if balls is None:
        raise ValueError("A canonical checkpoint requires source atom metadata")

    group_objects = list(getattr(system, "groups", ()) or ())
    group_by_id = {}
    for group in group_objects:
        group_id = str(getattr(group, "group_id", None) or getattr(group, "name", "group"))
        group_by_id[group_id] = group
    for interface in getattr(system, "ifaces", ()) or ():
        for group in (getattr(interface, "group1", None), getattr(interface, "group2", None)):
            if group is not None:
                group_id = str(getattr(group, "group_id", None) or getattr(group, "name", "group"))
                group_by_id.setdefault(group_id, group)

    atom_group_ids = {}
    for group_id, group in group_by_id.items():
        for index in getattr(group, "ball_ndxs", ()) or ():
            atom_group_ids.setdefault(int(index), []).append(group_id)

    atoms = []
    residues = {}
    chains = {}
    for position, (_, row) in enumerate(balls.iterrows()):
        atom_id = f"atom:{position}"
        chain_name = str(_row_value(row, "chain_name", "chain", "Chain", default="_") or "_")
        residue_name = str(_row_value(row, "res_name", "residue", "res", "Residue", default="UNK") or "UNK")
        residue_seq = _int_value(_row_value(row, "res_seq", "residue_sequence", "Residue Sequence", default=position + 1), position + 1)
        chain_id = f"chain:{chain_name}"
        residue_id = f"residue:{chain_name}:{residue_seq}:{residue_name}"
        location = _row_value(row, "loc", "coordinates", default=(0.0, 0.0, 0.0))
        coordinates = tuple(float(value) for value in location)
        if len(coordinates) != 3:
            coordinates = (0.0, 0.0, 0.0)
        atoms.append(CheckpointAtom(
            atom_id=atom_id,
            source_index=position,
            atom_name=str(_row_value(row, "name", "atom_name", default="ATOM") or "ATOM"),
            element=str(_row_value(row, "element", "Element", default="X") or "X"),
            residue_id=residue_id,
            residue_name=residue_name,
            residue_seq=residue_seq,
            chain_id=chain_id,
            coordinates=coordinates,
            radius=max(0.0, float(_row_value(row, "rad", "radius", default=0.0) or 0.0)),
            group_ids=tuple(atom_group_ids.get(position, ())),
        ))
        residues.setdefault(residue_id, {"chain_id": chain_id, "name": residue_name, "sequence": residue_seq, "atom_ids": []})["atom_ids"].append(atom_id)
        chain = chains.setdefault(chain_id, {"name": chain_name, "residue_ids": [], "atom_ids": []})
        if residue_id not in chain["residue_ids"]:
            chain["residue_ids"].append(residue_id)
        chain["atom_ids"].append(atom_id)

    residue_records = tuple(
        CheckpointResidue(residue_id=key, **value)
        for key, value in residues.items()
    )
    chain_records = tuple(
        CheckpointChain(chain_id=key, **value)
        for key, value in chains.items()
    )

    payloads = []
    dependencies = []
    networks = []
    network_ids = {}
    for scope, owner_id, network in _owner_networks(system):
        key = id(network)
        if key in network_ids:
            continue
        network_id = f"network:{scope}:{owner_id}"
        network_ids[key] = network_id
        generator_ids = _network_ids(network)
        topology = {
            "balls": _table_len(network, "balls"),
            "edges": _table_len(network, "edges"),
            "verts": _table_len(network, "verts"),
            "surfs": _table_len(network, "surfs"),
            "generator_ids": generator_ids,
        }
        coordinate_rows = []
        radii = []
        network_balls = getattr(network, "balls", None)
        for _, row in (network_balls.iterrows() if network_balls is not None else ()):
            coordinate_rows.append(tuple(float(value) for value in _row_value(row, "loc", default=(0.0, 0.0, 0.0))))
            radii.append(float(_row_value(row, "rad", default=0.0) or 0.0))
        dependency = make_checkpoint_dependency(
            source_network_id=str(getattr(network, "source_network_id", None) or f"producer:{system.name}"),
            solve_universe_id=str(getattr(network, "solve_universe_id", None) or f"solve:{system.name}"),
            source_sha256=source.sha256,
            coordinates=coordinate_rows,
            radii=radii,
            solver_settings=_json_safe(getattr(network, "settings", {}) or {}),
            generator_ids=generator_ids,
            topology=topology,
        )
        dependencies.append(dependency)
        try:
            payload_kind, content = _network_payload(network)
        except Exception:
            payload_kind, content = None, b""
        status = _network_status(network, content)
        payload = None
        if content:
            payload = CheckpointPayload(
                kind=payload_kind,
                content=content,
                dependency_digest=dependency.content_digest,
                status=status,
            )
            payloads.append(payload)
        networks.append(CheckpointNetworkState(
            network_id=network_id,
            source_network_id=dependency.source_network_id,
            scope=scope,
            scheme=str((getattr(network, "settings", {}) or {}).get("net_type", "aw")).lower(),
            dependency_digest=dependency.content_digest,
            generator_ids=generator_ids,
            status=status,
            completeness={
                "balls": topology["balls"], "edges": topology["edges"],
                "verts": topology["verts"], "surfs": topology["surfs"],
            },
            payload_ref=payload.payload_id if payload is not None else None,
            settings=_json_safe(getattr(network, "settings", {}) or {}),
        ))

    groups = []
    for group_id, group in group_by_id.items():
        network = getattr(group, "net", None) or getattr(group, "source_network", None)
        groups.append(CheckpointGroup(
            group_id=group_id,
            name=str(getattr(group, "name", group_id)),
            atom_ids=tuple(f"atom:{int(index)}" for index in (getattr(group, "ball_ndxs", ()) or ()) if 0 <= int(index) < len(atoms)),
            generator_ids=_network_ids(network) if network is not None else (),
            network_id=network_ids.get(id(network)),
            status=Status.CERTIFIED,
        ))
    interfaces = []
    for interface in getattr(system, "ifaces", ()) or ():
        first = getattr(interface, "group1", None)
        second = getattr(interface, "group2", None)
        if first is None or second is None:
            continue
        first_id = str(getattr(first, "group_id", None) or getattr(first, "name", "group_a"))
        second_id = str(getattr(second, "group_id", None) or getattr(second, "name", "group_b"))
        if first_id != second_id:
            interfaces.append(CheckpointInterfaceIdentity(
                system_id=f"system:{system.name}",
                group_a_id=first_id,
                group_b_id=second_id,
                display_group_a_id=first_id,
                display_group_b_id=second_id,
            ))

    network_by_owner = {}
    for group in group_objects:
        network = getattr(group, "net", None) or getattr(group, "source_network", None)
        if network is not None:
            network_by_owner[("group", str(getattr(group, "group_id", None) or getattr(group, "name", "group")))] = network_ids.get(id(network))
    for interface in getattr(system, "ifaces", ()) or ():
        network = getattr(interface, "net", None) or getattr(interface, "source_network", None)
        if network is not None:
            network_by_owner[("interface", str(getattr(interface, "interface_id", None) or getattr(interface, "name", "interface")))] = network_ids.get(id(network))

    result_records = []
    provenance_values = []
    candidates = list(getattr(system, "cli_interface_results", ()) or ())
    candidates.extend(getattr(system, "results", ()) or ())
    for owner in list(group_objects) + list(getattr(system, "ifaces", ()) or ()):
        candidates.extend(getattr(owner, "typed_results", ()) or ())
        candidates.extend(getattr(owner, "results", ()) or ())
    seen_result_ids = set()
    dependency_by_network = {item.network_id: item.dependency_digest for item in networks}
    for result in candidates:
        identity = getattr(result, "identity", None)
        provenance = getattr(result, "provenance", None)
        if identity is None or provenance is None or not hasattr(result, "to_dict"):
            continue
        result_id = str(getattr(result, "result_id", None) or getattr(identity, "result_id", ""))
        if not result_id or result_id in seen_result_ids:
            continue
        network_ref = None
        interface_id = getattr(identity, "interface_id", None)
        group_id = getattr(identity, "group_id", None)
        if interface_id:
            network_ref = network_by_owner.get(("interface", str(interface_id)))
        if network_ref is None and group_id:
            network_ref = network_by_owner.get(("group", str(group_id)))
        if network_ref is None:
            continue
        status = Status(getattr(result, "status", Status.UNRESOLVED))
        if status is Status.CERTIFIED:
            network = next(item for item in networks if item.network_id == network_ref)
            if network.scheme == "prm":
                status = Status.UNRESOLVED
        try:
            provenance_values.append(provenance)
            result_records.append({
                "result_id": result_id,
                "representation": str(getattr(identity, "representation", "unknown")),
                "status": status,
                "dependency_digest": dependency_by_network[network_ref],
                "typed_result": _json_safe(result.to_dict()),
                "provenance_ref": provenance.content_digest,
                "network_id": network_ref,
            })
            seen_result_ids.add(result_id)
        except (AttributeError, TypeError, ValueError):
            continue

    document = CheckpointDocument(
        system_id=f"system:{system.name}",
        source=source,
        atoms=tuple(atoms),
        residues=residue_records,
        chains=chain_records,
        groups=tuple(groups),
        interfaces=tuple(interfaces),
        dependencies=tuple(dependencies),
        networks=tuple(networks),
        results=tuple(
            CheckpointResultState(**item)
            for item in result_records
        ),
        provenance=ProvenanceCatalog.from_provenances(provenance_values),
    )
    return CheckpointContainer(document=document, payloads=tuple(payloads))


def _load_checkpoint(path):
    with zipfile.ZipFile(path, "r") as archive:
        manifest = archive.read("checkpoint/manifest.json")
        members = {
            name: archive.read(name)
            for name in archive.namelist()
            if name.startswith("checkpoint/payloads/")
        }
    return load_checkpoint_container(manifest, members)


def write_checkpoint(system, destination=None):
    """Write one validated ATLAS container with same-directory atomic replace."""
    root = Path(system.files["dir"])
    destination = Path(destination or root / f"{system.name}_checkpoint.vpy")
    if destination.exists() and getattr(system, "_fail_if_exists", False):
        raise FileExistsError(f"Refusing to overwrite checkpoint: {destination}")
    replacement = build_checkpoint_container(system)
    current = _load_checkpoint(destination) if destination.exists() else None
    plan = prepare_checkpoint_update(replacement, current=current)
    bundle = plan.save_bundle()
    plan.validate_written(bundle)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix=f".{destination.name}.", suffix=".tmp", dir=root, delete=False) as handle:
            temporary = Path(handle.name)
            with zipfile.ZipFile(handle, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("checkpoint/manifest.json", bundle.manifest_bytes)
                for member, content in sorted(bundle.payload_members.items()):
                    archive.writestr(member, content)
            handle.flush()
            os.fsync(handle.fileno())
        _load_checkpoint(temporary)
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    system._last_checkpoint_path = str(destination.resolve())
    system._last_checkpoint_bytes = destination.stat().st_size
    system._last_checkpoint_digest = replacement.content_digest
    return destination
