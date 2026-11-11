"""ATLAS-owned checkpoint identity and validation contracts.

This module defines the payload that VECTOR may encode in the one shared
``<system>_checkpoint.vpy`` archive.  It deliberately performs no file I/O,
geometry solve, topology merge, or archive migration.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
from collections.abc import Mapping

from .model import ProvenanceCatalog, Status, stable_id


CHECKPOINT_FORMAT = "vorpy.system.checkpoint"
CHECKPOINT_SCHEMA_VERSION = 1
CHECKPOINT_CONTAINER_FORMAT = "vorpy.system.checkpoint.container"
CHECKPOINT_CONTAINER_VERSION = 1
CHECKPOINT_CONTAINER_PREFIX = "checkpoint-container-v1"
CHECKPOINT_PAYLOAD_PREFIX = "checkpoint-payload-v1"
CHECKPOINT_MANIFEST_MEMBER = "checkpoint/manifest.json"
CHECKPOINT_PAYLOAD_MEMBER_PREFIX = "checkpoint/payloads/"
CHECKPOINT_FILENAME_SUFFIX = "_checkpoint.vpy"
CHECKPOINT_DEPENDENCY_PREFIX = "checkpoint-dependency-v1"
CHECKPOINT_FEATURE_PREFIX = "checkpoint-feature-v1"
GROUP_C_EXTENSION_SCHEMA_VERSION = 1
_SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
_SAFE_SYSTEM_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class CheckpointContractError(ValueError):
    """Raised when checkpoint identity or cross-reference evidence is invalid."""


class CheckpointPersistenceError(CheckpointContractError):
    """Base class for explicit checkpoint serialization/persistence failures."""


class CheckpointLoadError(CheckpointPersistenceError):
    """Raised when a checkpoint manifest or referenced payload cannot be loaded."""


class CheckpointSaveError(CheckpointPersistenceError):
    """Raised when a checkpoint cannot produce a validated write bundle."""


class CheckpointUpdateError(CheckpointPersistenceError):
    """Raised when an update precondition or post-write validation fails."""


def _sha256(value, label):
    if not isinstance(value, str) or not _SHA256_RE.fullmatch(value):
        raise CheckpointContractError(f"{label} must be a SHA-256 hex digest")
    return value.lower()


def _payload_digest(value):
    try:
        encoded = json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise CheckpointContractError("Checkpoint digest payload is not JSON-safe") from error
    return hashlib.sha256(encoded).hexdigest()


def _prefixed(value, prefix, label):
    expected = f"{prefix}:"
    if not isinstance(value, str) or not value.startswith(expected):
        raise CheckpointContractError(f"{label} must use the {prefix}:<sha256> form")
    _sha256(value.rsplit(":", 1)[-1], label)
    if len(value) != len(prefix) + 65:
        raise CheckpointContractError(f"{label} must use the {prefix}:<sha256> form")
    return value


def _text(value, label):
    if not isinstance(value, str) or not value:
        raise CheckpointContractError(f"{label} must be a non-empty string")
    return value


def _unique_strings(values, label):
    values = tuple(values)
    if any(not isinstance(value, str) or not value for value in values):
        raise CheckpointContractError(f"{label} must contain non-empty strings")
    if len(set(values)) != len(values):
        raise CheckpointContractError(f"{label} must not contain duplicates")
    return values


@dataclass(frozen=True)
class CheckpointPolicy:
    """One-root checkpoint policy; VECTOR owns the eventual filesystem write."""

    enabled: bool = True
    atomic_replace: bool = True
    preserve_previous_valid: bool = True
    same_directory_temporary: bool = True
    filename_suffix: str = CHECKPOINT_FILENAME_SUFFIX

    def __post_init__(self):
        if type(self.enabled) is not bool:
            raise CheckpointContractError("checkpoint enabled must be boolean")
        for name in (
            "atomic_replace",
            "preserve_previous_valid",
            "same_directory_temporary",
        ):
            if type(getattr(self, name)) is not bool:
                raise CheckpointContractError(f"{name} must be boolean")
        if self.filename_suffix != CHECKPOINT_FILENAME_SUFFIX:
            raise CheckpointContractError("checkpoint filename suffix is fixed")

    @staticmethod
    def _validate_system_name(system_name):
        if not isinstance(system_name, str) or not _SAFE_SYSTEM_NAME_RE.fullmatch(system_name):
            raise CheckpointContractError(
                "system name must contain only letters, digits, '.', '_' or '-'"
            )
        return system_name

    def filename(self, system_name):
        return self._validate_system_name(system_name) + self.filename_suffix

    def path(self, output_root, system_name):
        root = Path(output_root)
        if root.name in {"", "."}:
            raise CheckpointContractError("checkpoint output root must be explicit")
        return root / self.filename(system_name)


@dataclass(frozen=True)
class CheckpointSource:
    """Byte-preserved molecular input embedded in a new checkpoint."""

    filename: str
    content: bytes

    def __post_init__(self):
        if not isinstance(self.filename, str) or not self.filename:
            raise CheckpointContractError("source filename must be non-empty")
        if Path(self.filename).name != self.filename or self.filename in {".", ".."}:
            raise CheckpointContractError("source filename must be a basename")
        if not isinstance(self.content, (bytes, bytearray)):
            raise CheckpointContractError("source content must be bytes")
        object.__setattr__(self, "content", bytes(self.content))

    @property
    def size_bytes(self):
        return len(self.content)

    @property
    def sha256(self):
        return hashlib.sha256(self.content).hexdigest()

    @property
    def source_id(self):
        return f"source-v1:{self.sha256}"

    def as_dict(self):
        return {
            "filename": self.filename,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "source_id": self.source_id,
            "content_base64": base64.b64encode(self.content).decode("ascii"),
        }

    @classmethod
    def from_mapping(cls, value):
        if not isinstance(value, Mapping):
            raise CheckpointContractError("checkpoint source must be a mapping")
        try:
            content = base64.b64decode(value["content_base64"], validate=True)
        except (KeyError, TypeError, ValueError) as error:
            raise CheckpointContractError("checkpoint source content is invalid") from error
        result = cls(value["filename"], content)
        if (
            value.get("size_bytes") != result.size_bytes
            or value.get("sha256") != result.sha256
            or value.get("source_id", result.source_id) != result.source_id
        ):
            raise CheckpointContractError("checkpoint source digest or size mismatch")
        return result


@dataclass(frozen=True)
class CheckpointPayload:
    """Opaque, lossless bytes for network geometry, analysis, or cache state."""

    kind: str
    content: bytes
    dependency_digest: str | None = None
    status: Status = Status.CERTIFIED
    payload_id: str = field(init=False)

    def __post_init__(self):
        _text(self.kind, "payload kind")
        if not isinstance(self.content, (bytes, bytearray)):
            raise CheckpointContractError("payload content must be bytes")
        object.__setattr__(self, "content", bytes(self.content))
        if self.dependency_digest is not None:
            _prefixed(self.dependency_digest, CHECKPOINT_DEPENDENCY_PREFIX, "payload dependency_digest")
        object.__setattr__(self, "status", Status(self.status))
        object.__setattr__(
            self,
            "payload_id",
            stable_id(
                CHECKPOINT_PAYLOAD_PREFIX,
                {"kind": self.kind, "sha256": self.sha256},
            ),
        )

    @property
    def size_bytes(self):
        return len(self.content)

    @property
    def sha256(self):
        return hashlib.sha256(self.content).hexdigest()

    @property
    def member_name(self):
        return CHECKPOINT_PAYLOAD_MEMBER_PREFIX + self.payload_id.rsplit(":", 1)[-1]

    def as_dict(self):
        return {
            "payload_id": self.payload_id,
            "kind": self.kind,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "dependency_digest": self.dependency_digest,
            "status": self.status.value,
            "member": self.member_name,
        }

    @classmethod
    def from_mapping(cls, value, content):
        if not isinstance(value, Mapping):
            raise CheckpointLoadError("checkpoint payload metadata must be a mapping")
        try:
            result = cls(
                kind=value["kind"],
                content=content,
                dependency_digest=value.get("dependency_digest"),
                status=value.get("status", Status.CERTIFIED.value),
            )
        except (KeyError, TypeError, ValueError) as error:
            raise CheckpointLoadError("checkpoint payload metadata is invalid") from error
        if value.get("payload_id") != result.payload_id:
            raise CheckpointLoadError("checkpoint payload identity mismatch")
        if value.get("size_bytes") != result.size_bytes or value.get("sha256") != result.sha256:
            raise CheckpointLoadError("checkpoint payload digest or size mismatch")
        if value.get("member") != result.member_name:
            raise CheckpointLoadError("checkpoint payload member mismatch")
        return result


@dataclass(frozen=True)
class CheckpointAtom:
    """Stable source atom identity and the geometry needed for source-free reopen."""

    atom_id: str
    source_index: int
    atom_name: str
    element: str
    residue_id: str
    residue_name: str
    residue_seq: int
    chain_id: str
    coordinates: tuple[float, float, float]
    radius: float
    group_ids: tuple[str, ...] = ()

    def __post_init__(self):
        _text(self.atom_id, "atom_id")
        if isinstance(self.source_index, bool) or not isinstance(self.source_index, int) or self.source_index < 0:
            raise CheckpointContractError("source_index must be a nonnegative integer")
        for name in ("atom_name", "element", "residue_id", "residue_name", "chain_id"):
            _text(getattr(self, name), name)
        if isinstance(self.residue_seq, bool) or not isinstance(self.residue_seq, int):
            raise CheckpointContractError("residue_seq must be an integer")
        coordinates = tuple(float(value) for value in self.coordinates)
        if len(coordinates) != 3 or not all(math.isfinite(value) for value in coordinates):
            raise CheckpointContractError("atom coordinates must be three finite numbers")
        if not isinstance(self.radius, (int, float)) or isinstance(self.radius, bool):
            raise CheckpointContractError("atom radius must be numeric")
        if not math.isfinite(float(self.radius)) or float(self.radius) < 0:
            raise CheckpointContractError("atom radius must be finite and nonnegative")
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "radius", float(self.radius))
        object.__setattr__(self, "group_ids", _unique_strings(self.group_ids, "atom group_ids"))

    def as_dict(self):
        return {
            "atom_id": self.atom_id,
            "source_index": self.source_index,
            "atom_name": self.atom_name,
            "element": self.element,
            "residue_id": self.residue_id,
            "residue_name": self.residue_name,
            "residue_seq": self.residue_seq,
            "chain_id": self.chain_id,
            "coordinates": self.coordinates,
            "radius": self.radius,
            "group_ids": self.group_ids,
        }


@dataclass(frozen=True)
class CheckpointResidue:
    residue_id: str
    chain_id: str
    name: str
    sequence: int
    atom_ids: tuple[str, ...]

    def __post_init__(self):
        for name in ("residue_id", "chain_id", "name"):
            _text(getattr(self, name), name)
        if isinstance(self.sequence, bool) or not isinstance(self.sequence, int):
            raise CheckpointContractError("residue sequence must be an integer")
        object.__setattr__(self, "atom_ids", _unique_strings(self.atom_ids, "residue atom_ids"))

    def as_dict(self):
        return {
            "residue_id": self.residue_id,
            "chain_id": self.chain_id,
            "name": self.name,
            "sequence": self.sequence,
            "atom_ids": self.atom_ids,
        }


@dataclass(frozen=True)
class CheckpointChain:
    chain_id: str
    name: str
    residue_ids: tuple[str, ...]
    atom_ids: tuple[str, ...]

    def __post_init__(self):
        _text(self.chain_id, "chain_id")
        _text(self.name, "chain name")
        object.__setattr__(self, "residue_ids", _unique_strings(self.residue_ids, "chain residue_ids"))
        object.__setattr__(self, "atom_ids", _unique_strings(self.atom_ids, "chain atom_ids"))

    def as_dict(self):
        return {
            "chain_id": self.chain_id,
            "name": self.name,
            "residue_ids": self.residue_ids,
            "atom_ids": self.atom_ids,
        }


@dataclass(frozen=True)
class CheckpointGroup:
    group_id: str
    name: str
    atom_ids: tuple[str, ...]
    generator_ids: tuple[str, ...] = ()
    network_id: str | None = None
    status: Status = Status.CERTIFIED

    def __post_init__(self):
        _text(self.group_id, "group_id")
        _text(self.name, "group name")
        object.__setattr__(self, "atom_ids", _unique_strings(self.atom_ids, "group atom_ids"))
        object.__setattr__(self, "generator_ids", _unique_strings(self.generator_ids, "group generator_ids"))
        if self.network_id is not None:
            _text(self.network_id, "group network_id")
        object.__setattr__(self, "status", Status(self.status))

    def as_dict(self):
        return {
            "group_id": self.group_id,
            "name": self.name,
            "atom_ids": self.atom_ids,
            "generator_ids": self.generator_ids,
            "network_id": self.network_id,
            "status": self.status.value,
        }


@dataclass(frozen=True)
class CheckpointInterfaceIdentity:
    """One unordered physical interface with a deterministic reference direction."""

    system_id: str
    group_a_id: str
    group_b_id: str
    display_group_a_id: str | None = None
    display_group_b_id: str | None = None
    interface_id: str = field(init=False)

    def __post_init__(self):
        _text(self.system_id, "interface system_id")
        first, second = sorted((_text(self.group_a_id, "group_a_id"), _text(self.group_b_id, "group_b_id")))
        if first == second:
            raise CheckpointContractError("interface groups must be distinct")
        display_a = self.group_a_id if self.display_group_a_id is None else _text(
            self.display_group_a_id, "display_group_a_id"
        )
        display_b = self.group_b_id if self.display_group_b_id is None else _text(
            self.display_group_b_id, "display_group_b_id"
        )
        if {display_a, display_b} != {first, second}:
            raise CheckpointContractError("display orientation must use the interface groups")
        object.__setattr__(self, "group_a_id", first)
        object.__setattr__(self, "group_b_id", second)
        object.__setattr__(self, "display_group_a_id", display_a)
        object.__setattr__(self, "display_group_b_id", display_b)
        object.__setattr__(
            self,
            "interface_id",
            stable_id("interface-v1", {"system_id": self.system_id, "group_ids": (first, second)}),
        )

    @property
    def reference_orientation(self):
        return (self.group_a_id, self.group_b_id)

    @property
    def display_orientation(self):
        return (self.display_group_a_id, self.display_group_b_id)

    def as_dict(self):
        return {
            "system_id": self.system_id,
            "interface_id": self.interface_id,
            "group_ids": (self.group_a_id, self.group_b_id),
            "reference_orientation": self.reference_orientation,
            "display_orientation": self.display_orientation,
        }


@dataclass(frozen=True)
class CheckpointDependency:
    """Identity of every input that can invalidate a cached network/result."""

    source_network_id: str
    solve_universe_id: str
    source_sha256: str
    coordinates_digest: str
    radii_digest: str
    solver_settings_digest: str
    generator_universe_digest: str
    topology_digest: str

    def __post_init__(self):
        for name in ("source_network_id", "solve_universe_id"):
            _text(getattr(self, name), name)
        for name in (
            "source_sha256",
            "coordinates_digest",
            "radii_digest",
            "solver_settings_digest",
            "generator_universe_digest",
            "topology_digest",
        ):
            object.__setattr__(self, name, _sha256(getattr(self, name), name))

    @property
    def content_digest(self):
        return stable_id(f"{CHECKPOINT_DEPENDENCY_PREFIX}", self.as_dict())

    def as_dict(self):
        return {
            "source_network_id": self.source_network_id,
            "solve_universe_id": self.solve_universe_id,
            "source_sha256": self.source_sha256,
            "coordinates_digest": self.coordinates_digest,
            "radii_digest": self.radii_digest,
            "solver_settings_digest": self.solver_settings_digest,
            "generator_universe_digest": self.generator_universe_digest,
            "topology_digest": self.topology_digest,
        }


def make_checkpoint_dependency(
    *,
    source_network_id,
    solve_universe_id,
    source_sha256,
    coordinates,
    radii,
    solver_settings,
    generator_ids,
    topology,
):
    """Build dependency digests from explicit source/geometry inputs."""
    return CheckpointDependency(
        source_network_id=source_network_id,
        solve_universe_id=solve_universe_id,
        source_sha256=source_sha256,
        coordinates_digest=_payload_digest(coordinates),
        radii_digest=_payload_digest(radii),
        solver_settings_digest=_payload_digest(solver_settings),
        generator_universe_digest=_payload_digest(generator_ids),
        topology_digest=_payload_digest(topology),
    )


def validate_checkpoint_dependency(expected, cached):
    """Return a status without rebuilding a stale network or analysis cache."""
    if cached is None:
        return Status.NOT_CALCULATED
    if not isinstance(expected, CheckpointDependency) or not isinstance(cached, CheckpointDependency):
        return Status.UNRESOLVED
    return Status.CERTIFIED if expected.content_digest == cached.content_digest else Status.STALE


@dataclass(frozen=True)
class CheckpointNetworkState:
    network_id: str
    source_network_id: str
    scope: str
    scheme: str
    dependency_digest: str
    generator_ids: tuple[str, ...]
    feature_ids: tuple[str, ...] = ()
    status: Status = Status.NOT_CALCULATED
    completeness: Mapping = field(default_factory=dict)
    payload_ref: str | None = None
    settings: Mapping = field(default_factory=dict)
    solvent_policy: Mapping = field(default_factory=dict)
    competitor_ids: tuple[str, ...] = ()

    def __post_init__(self):
        for name in ("network_id", "source_network_id", "dependency_digest"):
            _text(getattr(self, name), name)
        if self.scope not in {"system", "group", "interface"}:
            raise CheckpointContractError("network scope must be system, group, or interface")
        if self.scheme not in {"aw", "pow", "prm"}:
            raise CheckpointContractError("network scheme must be aw, pow, or prm")
        _prefixed(self.dependency_digest, CHECKPOINT_DEPENDENCY_PREFIX, "network dependency_digest")
        object.__setattr__(self, "generator_ids", _unique_strings(self.generator_ids, "network generator_ids"))
        object.__setattr__(self, "feature_ids", _unique_strings(self.feature_ids, "network feature_ids"))
        if not isinstance(self.completeness, Mapping):
            raise CheckpointContractError("network completeness must be a mapping")
        for name in ("settings", "solvent_policy"):
            if not isinstance(getattr(self, name), Mapping):
                raise CheckpointContractError(f"network {name} must be a mapping")
        object.__setattr__(
            self,
            "competitor_ids",
            _unique_strings(self.competitor_ids, "network competitor_ids"),
        )
        object.__setattr__(self, "status", Status(self.status))
        if self.payload_ref is None:
            object.__setattr__(self, "payload_ref", self.network_id)

    def as_dict(self):
        return {
            "network_id": self.network_id,
            "source_network_id": self.source_network_id,
            "scope": self.scope,
            "scheme": self.scheme,
            "dependency_digest": self.dependency_digest,
            "generator_ids": self.generator_ids,
            "feature_ids": self.feature_ids,
            "status": self.status.value,
            "completeness": dict(self.completeness),
            "payload_ref": self.payload_ref,
            "settings": dict(self.settings),
            "solvent_policy": dict(self.solvent_policy),
            "competitor_ids": self.competitor_ids,
        }


@dataclass(frozen=True)
class CheckpointCacheState:
    cache_id: str
    kind: str
    dependency_digest: str
    status: Status
    payload_ref: str
    provenance_ref: str | None = None

    def __post_init__(self):
        for name in ("cache_id", "kind", "payload_ref"):
            _text(getattr(self, name), name)
        _prefixed(self.dependency_digest, CHECKPOINT_DEPENDENCY_PREFIX, "cache dependency_digest")
        if self.provenance_ref is not None:
            _prefixed(self.provenance_ref, "provenance-v1", "cache provenance_ref")
        object.__setattr__(self, "status", Status(self.status))

    def as_dict(self):
        return {
            "cache_id": self.cache_id,
            "kind": self.kind,
            "dependency_digest": self.dependency_digest,
            "status": self.status.value,
            "payload_ref": self.payload_ref,
            "provenance_ref": self.provenance_ref,
        }


@dataclass(frozen=True)
class CheckpointResultState:
    result_id: str
    representation: str
    status: Status
    dependency_digest: str
    typed_result: Mapping
    provenance_ref: str | None = None
    network_id: str | None = None

    def __post_init__(self):
        _text(self.result_id, "result_id")
        _text(self.representation, "result representation")
        _prefixed(self.dependency_digest, CHECKPOINT_DEPENDENCY_PREFIX, "result dependency_digest")
        if not isinstance(self.typed_result, Mapping):
            raise CheckpointContractError("typed_result must be a JSON mapping")
        if self.provenance_ref is not None:
            _prefixed(self.provenance_ref, "provenance-v1", "result provenance_ref")
        if self.network_id is not None:
            _text(self.network_id, "result network_id")
        object.__setattr__(self, "status", Status(self.status))

    def as_dict(self):
        return {
            "result_id": self.result_id,
            "representation": self.representation,
            "status": self.status.value,
            "dependency_digest": self.dependency_digest,
            "typed_result": dict(self.typed_result),
            "provenance_ref": self.provenance_ref,
            "network_id": self.network_id,
        }


@dataclass(frozen=True)
class CheckpointFeatureIdentity:
    """Stable source feature identity used by the deferred Group C contract."""

    kind: str
    source_network_id: str
    source_id: str
    generator_ids: tuple[str, ...]
    topology_digest: str
    feature_id: str = field(init=False)

    def __post_init__(self):
        _text(self.kind, "feature kind")
        _text(self.source_network_id, "feature source_network_id")
        _text(self.source_id, "feature source_id")
        object.__setattr__(self, "generator_ids", _unique_strings(self.generator_ids, "feature generator_ids"))
        object.__setattr__(self, "topology_digest", _sha256(self.topology_digest, "feature topology_digest"))
        object.__setattr__(
            self,
            "feature_id",
            stable_id(
                CHECKPOINT_FEATURE_PREFIX,
                {
                    "kind": self.kind,
                    "source_network_id": self.source_network_id,
                    "source_id": self.source_id,
                    "generator_ids": self.generator_ids,
                },
            ),
        )

    def as_dict(self):
        return {
            "feature_id": self.feature_id,
            "kind": self.kind,
            "source_network_id": self.source_network_id,
            "source_id": self.source_id,
            "generator_ids": self.generator_ids,
            "topology_digest": self.topology_digest,
        }


@dataclass(frozen=True)
class GroupCExtensionContract:
    """Append-only Group C identity contract; topology merge is intentionally deferred."""

    system_id: str
    base_source_network_id: str
    group_id: str
    atom_ids: tuple[str, ...]
    generator_ids: tuple[str, ...]
    feature_identities: tuple[CheckpointFeatureIdentity, ...]
    dependency_digest: str
    status: Status = Status.NOT_CALCULATED
    schema_version: int = GROUP_C_EXTENSION_SCHEMA_VERSION

    @property
    def merge_supported(self):
        return False

    def __post_init__(self):
        for name in ("system_id", "base_source_network_id", "group_id"):
            _text(getattr(self, name), name)
        object.__setattr__(self, "atom_ids", _unique_strings(self.atom_ids, "Group C atom_ids"))
        object.__setattr__(self, "generator_ids", _unique_strings(self.generator_ids, "Group C generator_ids"))
        if any(not isinstance(item, CheckpointFeatureIdentity) for item in self.feature_identities):
            raise CheckpointContractError("Group C feature_identities must be typed identities")
        feature_ids = tuple(item.feature_id for item in self.feature_identities)
        _unique_strings(feature_ids, "Group C feature identities")
        _prefixed(self.dependency_digest, CHECKPOINT_DEPENDENCY_PREFIX, "Group C dependency_digest")
        object.__setattr__(self, "status", Status(self.status))
        if self.schema_version != GROUP_C_EXTENSION_SCHEMA_VERSION:
            raise CheckpointContractError("unsupported Group C extension schema")

    def as_dict(self):
        return {
            "schema_version": self.schema_version,
            "system_id": self.system_id,
            "base_source_network_id": self.base_source_network_id,
            "group_id": self.group_id,
            "atom_ids": self.atom_ids,
            "generator_ids": self.generator_ids,
            "feature_identities": [item.as_dict() for item in self.feature_identities],
            "dependency_digest": self.dependency_digest,
            "status": self.status.value,
            "merge_supported": False,
        }


@dataclass(frozen=True)
class CheckpointDocument:
    """Self-contained checkpoint manifest independent of the archive byte writer."""

    system_id: str
    source: CheckpointSource
    atoms: tuple[CheckpointAtom, ...]
    residues: tuple[CheckpointResidue, ...]
    chains: tuple[CheckpointChain, ...]
    groups: tuple[CheckpointGroup, ...]
    interfaces: tuple[CheckpointInterfaceIdentity, ...]
    dependencies: tuple[CheckpointDependency, ...]
    networks: tuple[CheckpointNetworkState, ...] = ()
    caches: tuple[CheckpointCacheState, ...] = ()
    results: tuple[CheckpointResultState, ...] = ()
    provenance: ProvenanceCatalog = field(default_factory=ProvenanceCatalog)
    group_extensions: tuple[GroupCExtensionContract, ...] = ()
    schema_version: int = CHECKPOINT_SCHEMA_VERSION

    def __post_init__(self):
        _text(self.system_id, "checkpoint system_id")
        if not isinstance(self.source, CheckpointSource):
            raise CheckpointContractError("checkpoint source must be typed")
        if not isinstance(self.provenance, ProvenanceCatalog):
            raise CheckpointContractError("checkpoint provenance must be a ProvenanceCatalog")
        if self.schema_version != CHECKPOINT_SCHEMA_VERSION:
            raise CheckpointContractError("unsupported checkpoint schema version")
        self._validate_unique(self.atoms, "atom_id", "atoms")
        self._validate_unique(self.residues, "residue_id", "residues")
        self._validate_unique(self.chains, "chain_id", "chains")
        self._validate_unique(self.groups, "group_id", "groups")
        self._validate_unique(self.interfaces, "interface_id", "interfaces")
        self._validate_unique(self.networks, "network_id", "networks")
        self._validate_unique(self.caches, "cache_id", "caches")
        self._validate_unique(self.results, "result_id", "results")
        if any(not isinstance(item, CheckpointDependency) for item in self.dependencies):
            raise CheckpointContractError("dependencies contain an invalid typed record")
        dependency_ids = tuple(item.content_digest for item in self.dependencies)
        if len(set(dependency_ids)) != len(dependency_ids):
            raise CheckpointContractError("dependencies contain duplicate identities")
        atom_ids = {item.atom_id for item in self.atoms}
        residue_ids = {item.residue_id for item in self.residues}
        chain_ids = {item.chain_id for item in self.chains}
        group_ids = {item.group_id for item in self.groups}
        dependency_ids = set(dependency_ids)
        network_by_id = {network.network_id: network for network in self.networks}
        for atom in self.atoms:
            if atom.residue_id not in residue_ids or atom.chain_id not in chain_ids:
                raise CheckpointContractError("atom membership references an unknown residue or chain")
            if any(group_id not in group_ids for group_id in atom.group_ids):
                raise CheckpointContractError("atom membership references an unknown group")
        for residue in self.residues:
            if residue.chain_id not in chain_ids or any(atom_id not in atom_ids for atom_id in residue.atom_ids):
                raise CheckpointContractError("residue membership is not closed")
        for chain in self.chains:
            if any(residue_id not in residue_ids for residue_id in chain.residue_ids):
                raise CheckpointContractError("chain references an unknown residue")
            if any(atom_id not in atom_ids for atom_id in chain.atom_ids):
                raise CheckpointContractError("chain references an unknown atom")
        for group in self.groups:
            if any(atom_id not in atom_ids for atom_id in group.atom_ids):
                raise CheckpointContractError("group references an unknown atom")
        for interface in self.interfaces:
            if interface.system_id != self.system_id or any(
                group_id not in group_ids for group_id in interface.reference_orientation
            ):
                raise CheckpointContractError("interface identity is outside this checkpoint")
        for network in self.networks:
            if network.dependency_digest not in dependency_ids:
                raise CheckpointContractError("network references a missing dependency")
            if any(competitor_id not in atom_ids for competitor_id in network.competitor_ids):
                raise CheckpointContractError("network competitor universe references an unknown atom")
        for cache in self.caches:
            if cache.dependency_digest not in dependency_ids:
                raise CheckpointContractError("cache references a missing dependency")
            if cache.provenance_ref is not None:
                try:
                    self.provenance.resolve(cache.provenance_ref)
                except (KeyError, ValueError) as error:
                    raise CheckpointContractError("cache references missing provenance") from error
        for result in self.results:
            if result.dependency_digest not in dependency_ids:
                raise CheckpointContractError("result references a missing dependency")
            if result.network_id is None:
                if result.status is Status.CERTIFIED:
                    raise CheckpointContractError(
                        "certified Result must identify its source network"
                    )
            else:
                network = network_by_id.get(result.network_id)
                if network is None:
                    raise CheckpointContractError("result references an unknown network")
                if result.dependency_digest != network.dependency_digest:
                    raise CheckpointContractError("result dependency differs from its network")
                if network.scheme == "prm" and result.status is Status.CERTIFIED:
                    raise CheckpointContractError(
                        "PRM typed Results remain unresolved until scientific support is approved"
                    )
            if result.provenance_ref is not None:
                try:
                    self.provenance.resolve(result.provenance_ref)
                except (KeyError, ValueError) as error:
                    raise CheckpointContractError("result references missing provenance") from error
        for extension in self.group_extensions:
            if validate_group_c_extension(self, extension) is not Status.CERTIFIED:
                raise CheckpointContractError("invalid Group C extension identity")

    @staticmethod
    def _validate_unique(values, field_name, label):
        values = tuple(values)
        if any(not hasattr(value, field_name) for value in values):
            raise CheckpointContractError(f"{label} contain an invalid typed record")
        identifiers = tuple(getattr(value, field_name) for value in values)
        if len(set(identifiers)) != len(identifiers):
            raise CheckpointContractError(f"{label} contain duplicate identities")

    def as_dict(self):
        return {
            "format": CHECKPOINT_FORMAT,
            "schema_version": self.schema_version,
            "system_id": self.system_id,
            "source": self.source.as_dict(),
            "atoms": [item.as_dict() for item in self.atoms],
            "residues": [item.as_dict() for item in self.residues],
            "chains": [item.as_dict() for item in self.chains],
            "groups": [item.as_dict() for item in self.groups],
            "interfaces": [item.as_dict() for item in self.interfaces],
            "dependencies": [item.as_dict() | {"dependency_digest": item.content_digest} for item in self.dependencies],
            "networks": [item.as_dict() for item in self.networks],
            "caches": [item.as_dict() for item in self.caches],
            "results": [item.as_dict() for item in self.results],
            "provenance": self.provenance.as_dict(),
            "group_extensions": [item.as_dict() for item in self.group_extensions],
        }

    def serialize(self):
        return json.dumps(
            self.as_dict(),
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )

    @classmethod
    def from_mapping(cls, value):
        if not isinstance(value, Mapping) or value.get("format") != CHECKPOINT_FORMAT:
            raise CheckpointContractError("not a VorPy checkpoint contract")
        if value.get("schema_version") != CHECKPOINT_SCHEMA_VERSION:
            raise CheckpointContractError("unsupported checkpoint schema version")
        dependencies = []
        for item in value.get("dependencies", ()):
            item = dict(item)
            item.pop("dependency_digest", None)
            dependencies.append(CheckpointDependency(**item))
        interfaces = []
        for item in value.get("interfaces", ()):
            groups = tuple(item["group_ids"])
            display = tuple(item.get("display_orientation", groups))
            interfaces.append(
                CheckpointInterfaceIdentity(
                    item["system_id"], groups[0], groups[1], display[0], display[1]
                )
            )
        features = []
        for item in value.get("group_extensions", ()):
            features.append(
                GroupCExtensionContract(
                    schema_version=item.get("schema_version", GROUP_C_EXTENSION_SCHEMA_VERSION),
                    system_id=item["system_id"],
                    base_source_network_id=item["base_source_network_id"],
                    group_id=item["group_id"],
                    atom_ids=tuple(item.get("atom_ids", ())),
                    generator_ids=tuple(item.get("generator_ids", ())),
                    feature_identities=tuple(
                        CheckpointFeatureIdentity(
                            kind=feature["kind"],
                            source_network_id=feature["source_network_id"],
                            source_id=feature["source_id"],
                            generator_ids=tuple(feature.get("generator_ids", ())),
                            topology_digest=feature["topology_digest"],
                        )
                        for feature in item.get("feature_identities", ())
                    ),
                    dependency_digest=item["dependency_digest"],
                    status=item.get("status", Status.NOT_CALCULATED.value),
                )
            )
        return cls(
            system_id=value["system_id"],
            source=CheckpointSource.from_mapping(value["source"]),
            atoms=tuple(CheckpointAtom(**item) for item in value.get("atoms", ())),
            residues=tuple(CheckpointResidue(**item) for item in value.get("residues", ())),
            chains=tuple(CheckpointChain(**item) for item in value.get("chains", ())),
            groups=tuple(CheckpointGroup(**item) for item in value.get("groups", ())),
            interfaces=tuple(interfaces),
            dependencies=tuple(dependencies),
            networks=tuple(CheckpointNetworkState(**item) for item in value.get("networks", ())),
            caches=tuple(CheckpointCacheState(**item) for item in value.get("caches", ())),
            results=tuple(CheckpointResultState(**item) for item in value.get("results", ())),
            provenance=ProvenanceCatalog.from_mapping(value.get("provenance", {})),
            group_extensions=tuple(features),
            schema_version=value.get("schema_version"),
        )

    @classmethod
    def deserialize(cls, payload):
        try:
            value = json.loads(payload)
        except (TypeError, json.JSONDecodeError) as error:
            raise CheckpointContractError("invalid checkpoint contract JSON") from error
        return cls.from_mapping(value)

    @property
    def content_digest(self):
        return stable_id("checkpoint-v1", self.as_dict())


@dataclass(frozen=True)
class CheckpointContainer:
    """Validated manifest plus opaque payloads for one shared checkpoint."""

    document: CheckpointDocument
    payloads: tuple[CheckpointPayload, ...] = ()
    container_version: int = CHECKPOINT_CONTAINER_VERSION

    def __post_init__(self):
        if not isinstance(self.document, CheckpointDocument):
            raise CheckpointContractError("checkpoint container document must be typed")
        if self.container_version != CHECKPOINT_CONTAINER_VERSION:
            raise CheckpointContractError("unsupported checkpoint container version")
        if any(not isinstance(item, CheckpointPayload) for item in self.payloads):
            raise CheckpointContractError("checkpoint payloads must be typed")
        payload_ids = tuple(item.payload_id for item in self.payloads)
        if len(set(payload_ids)) != len(payload_ids):
            raise CheckpointContractError("checkpoint payloads contain duplicate identities")
        payload_map = set(payload_ids)
        dependency_ids = {item.content_digest for item in self.document.dependencies}
        for payload in self.payloads:
            if payload.dependency_digest is not None and payload.dependency_digest not in dependency_ids:
                raise CheckpointContractError("payload references a missing dependency")
        for network in self.document.networks:
            if network.status in {Status.CERTIFIED, Status.PARTIAL} and network.payload_ref not in payload_map:
                raise CheckpointContractError("calculated network references a missing payload")
        for cache in self.document.caches:
            if cache.status in {Status.CERTIFIED, Status.PARTIAL} and cache.payload_ref not in payload_map:
                raise CheckpointContractError("calculated cache references a missing payload")

    @property
    def system_id(self):
        return self.document.system_id

    @property
    def content_digest(self):
        return stable_id(CHECKPOINT_CONTAINER_PREFIX, self.as_manifest())

    def as_manifest(self):
        return {
            "format": CHECKPOINT_CONTAINER_FORMAT,
            "container_version": self.container_version,
            "system_id": self.system_id,
            "document": self.document.as_dict(),
            "payloads": [item.as_dict() for item in self.payloads],
        }

    def serialize(self):
        return json.dumps(
            self.as_manifest(),
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )

    def payload_members(self):
        return {item.member_name: item.content for item in self.payloads}

    def save_bundle(self):
        try:
            manifest_bytes = self.serialize().encode("utf-8")
            return CheckpointSaveBundle(
                manifest_bytes=manifest_bytes,
                payload_members=self.payload_members(),
                container_digest=self.content_digest,
            )
        except CheckpointPersistenceError:
            raise
        except (TypeError, ValueError, OSError) as error:
            raise CheckpointSaveError("checkpoint manifest could not be serialized") from error

    @classmethod
    def from_mapping(cls, value, payload_members=None):
        if not isinstance(value, Mapping) or value.get("format") != CHECKPOINT_CONTAINER_FORMAT:
            raise CheckpointLoadError("not a VorPy checkpoint container")
        if value.get("container_version") != CHECKPOINT_CONTAINER_VERSION:
            raise CheckpointLoadError("unsupported checkpoint container version")
        members = dict(payload_members or {})
        try:
            document = CheckpointDocument.from_mapping(value["document"])
        except (KeyError, TypeError, ValueError) as error:
            raise CheckpointLoadError("checkpoint document is invalid") from error
        payloads = []
        try:
            for metadata in value.get("payloads", ()):
                member = metadata["member"]
                if member not in members:
                    raise CheckpointLoadError(f"missing checkpoint payload member: {member}")
                payloads.append(CheckpointPayload.from_mapping(metadata, members[member]))
            result = cls(
                document=document,
                payloads=tuple(payloads),
                container_version=value.get("container_version"),
            )
        except CheckpointLoadError:
            raise
        except (KeyError, TypeError, ValueError) as error:
            raise CheckpointLoadError("checkpoint payload set is invalid") from error
        if value.get("system_id") != result.system_id:
            raise CheckpointLoadError("checkpoint system identity mismatch")
        return result

    @classmethod
    def deserialize(cls, payload, payload_members=None):
        try:
            value = json.loads(payload)
        except (TypeError, json.JSONDecodeError) as error:
            raise CheckpointLoadError("invalid checkpoint container JSON") from error
        return cls.from_mapping(value, payload_members)


@dataclass(frozen=True)
class CheckpointSaveBundle:
    """Writer-neutral bytes and members that VECTOR persists atomically."""

    manifest_bytes: bytes
    payload_members: Mapping[str, bytes]
    container_digest: str

    def __post_init__(self):
        if not isinstance(self.manifest_bytes, (bytes, bytearray)):
            raise CheckpointSaveError("checkpoint manifest bytes are required")
        object.__setattr__(self, "manifest_bytes", bytes(self.manifest_bytes))
        try:
            _prefixed(self.container_digest, CHECKPOINT_CONTAINER_PREFIX, "container_digest")
        except CheckpointContractError as error:
            raise CheckpointSaveError("container_digest is invalid") from error
        members = dict(self.payload_members)
        for name, content in members.items():
            if not isinstance(name, str) or not name.startswith(CHECKPOINT_PAYLOAD_MEMBER_PREFIX):
                raise CheckpointSaveError("checkpoint payload member name is invalid")
            if not isinstance(content, (bytes, bytearray)):
                raise CheckpointSaveError("checkpoint payload member must be bytes")
            members[name] = bytes(content)
        object.__setattr__(self, "payload_members", members)

    def load(self):
        try:
            container = load_checkpoint_container(self.manifest_bytes, self.payload_members)
        except CheckpointLoadError as error:
            raise CheckpointSaveError("written checkpoint bundle does not load") from error
        if container.content_digest != self.container_digest:
            raise CheckpointSaveError("written checkpoint bundle digest mismatch")
        return container


@dataclass(frozen=True)
class CheckpointUpdatePlan:
    """Precondition and post-write checks for VECTOR's atomic replace step."""

    replacement: CheckpointContainer
    expected_current_digest: str | None = None

    def __post_init__(self):
        if not isinstance(self.replacement, CheckpointContainer):
            raise CheckpointUpdateError("checkpoint replacement must be a container")
        if self.expected_current_digest is not None:
            try:
                _prefixed(
                    self.expected_current_digest,
                    CHECKPOINT_CONTAINER_PREFIX,
                    "expected_current_digest",
                )
            except CheckpointContractError as error:
                raise CheckpointUpdateError("expected_current_digest is invalid") from error

    def validate_current(self, current):
        if current is not None and not isinstance(current, CheckpointContainer):
            raise CheckpointUpdateError("current checkpoint is not a container")
        if self.expected_current_digest is not None:
            if current is None:
                raise CheckpointUpdateError("expected an existing checkpoint")
            if current.content_digest != self.expected_current_digest:
                raise CheckpointUpdateError("current checkpoint changed before update")
        return current

    def save_bundle(self):
        try:
            return self.replacement.save_bundle()
        except CheckpointSaveError as error:
            raise CheckpointUpdateError("checkpoint replacement cannot be saved") from error

    def validate_written(self, bundle):
        if not isinstance(bundle, CheckpointSaveBundle):
            raise CheckpointUpdateError("checkpoint writer returned an invalid bundle")
        try:
            loaded = bundle.load()
        except CheckpointSaveError as error:
            raise CheckpointUpdateError("checkpoint replacement failed validation") from error
        if loaded.content_digest != self.replacement.content_digest:
            raise CheckpointUpdateError("checkpoint replacement identity changed")
        return loaded


def load_checkpoint_container(payload, payload_members=None):
    """Load and validate a manifest plus its opaque payload members."""
    try:
        return CheckpointContainer.deserialize(payload, payload_members)
    except CheckpointLoadError:
        raise
    except (CheckpointContractError, TypeError, ValueError, KeyError) as error:
        raise CheckpointLoadError("checkpoint container failed validation") from error


def save_checkpoint_container(container):
    """Return a validated writer bundle; VECTOR owns filesystem persistence."""
    if not isinstance(container, CheckpointContainer):
        raise CheckpointSaveError("checkpoint save requires a typed container")
    try:
        return container.save_bundle()
    except CheckpointSaveError:
        raise
    except CheckpointContractError as error:
        raise CheckpointSaveError("checkpoint save failed validation") from error


def prepare_checkpoint_update(replacement, *, current=None, expected_current_digest=None):
    """Build an atomic-update plan without touching the destination path."""
    if not isinstance(replacement, CheckpointContainer):
        raise CheckpointUpdateError("checkpoint replacement must be a container")
    if current is not None and not isinstance(current, CheckpointContainer):
        raise CheckpointUpdateError("current checkpoint is not a container")
    if current is not None:
        current_digest = current.content_digest
        if expected_current_digest is not None and expected_current_digest != current_digest:
            raise CheckpointUpdateError("provided current checkpoint digest is stale")
        expected_current_digest = current_digest
    return CheckpointUpdatePlan(replacement, expected_current_digest)


def validate_group_c_extension(document, extension):
    """Validate append-only Group C identity without claiming a topology merge."""
    if not isinstance(document, CheckpointDocument) or not isinstance(extension, GroupCExtensionContract):
        return Status.UNRESOLVED
    if extension.system_id != document.system_id:
        return Status.UNRESOLVED
    if extension.group_id in {group.group_id for group in document.groups}:
        return Status.UNRESOLVED
    dependencies = {item.content_digest: item for item in document.dependencies}
    dependency = dependencies.get(extension.dependency_digest)
    if dependency is None or dependency.source_network_id != extension.base_source_network_id:
        return Status.STALE if dependency is None else Status.UNRESOLVED
    if any(item.source_network_id != extension.base_source_network_id for item in extension.feature_identities):
        return Status.UNRESOLVED
    existing_features = {
        feature_id
        for network in document.networks
        for feature_id in network.feature_ids
    }
    if any(item.feature_id in existing_features for item in extension.feature_identities):
        return Status.UNRESOLVED
    if extension.status is Status.CERTIFIED:
        return Status.NOT_SUPPORTED
    return Status.CERTIFIED


__all__ = [
    "CHECKPOINT_DEPENDENCY_PREFIX",
    "CHECKPOINT_CONTAINER_FORMAT",
    "CHECKPOINT_CONTAINER_PREFIX",
    "CHECKPOINT_CONTAINER_VERSION",
    "CHECKPOINT_FEATURE_PREFIX",
    "CHECKPOINT_FILENAME_SUFFIX",
    "CHECKPOINT_FORMAT",
    "CHECKPOINT_MANIFEST_MEMBER",
    "CHECKPOINT_PAYLOAD_MEMBER_PREFIX",
    "CHECKPOINT_PAYLOAD_PREFIX",
    "CHECKPOINT_SCHEMA_VERSION",
    "CheckpointContainer",
    "CheckpointAtom",
    "CheckpointCacheState",
    "CheckpointChain",
    "CheckpointContractError",
    "CheckpointDependency",
    "CheckpointDocument",
    "CheckpointFeatureIdentity",
    "CheckpointGroup",
    "CheckpointInterfaceIdentity",
    "CheckpointLoadError",
    "CheckpointNetworkState",
    "CheckpointPayload",
    "CheckpointPolicy",
    "CheckpointResidue",
    "CheckpointResultState",
    "CheckpointPersistenceError",
    "CheckpointSaveBundle",
    "CheckpointSaveError",
    "CheckpointSource",
    "CheckpointUpdateError",
    "CheckpointUpdatePlan",
    "GroupCExtensionContract",
    "GROUP_C_EXTENSION_SCHEMA_VERSION",
    "load_checkpoint_container",
    "make_checkpoint_dependency",
    "prepare_checkpoint_update",
    "save_checkpoint_container",
    "validate_checkpoint_dependency",
    "validate_group_c_extension",
]
