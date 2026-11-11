"""Results-owned identity checks for shared solved-network dependencies.

This module consumes producer identity; it never creates or rebuilds geometry.
The solve and view IDs are deliberately kept separate: a view selection may
change a dependency digest without changing the solve-universe identity.
"""

from __future__ import annotations

import json
from hashlib import sha256
from collections.abc import Mapping
from dataclasses import dataclass

from .model import ResultModel, Status, stable_id


ANALYSIS_CONTEXT_KINDS = frozenset({"GLOBAL_NETWORK", "ISOLATED_GROUP"})
SOURCE_IDENTITY_STATUSES = frozenset({
    "AUTHORITATIVE",
    "PROVISIONAL",
    "UNRESOLVED",
})
SOLVENT_MODES = frozenset({
    "MOLECULE_ONLY",
    "SOLVENT_CONTEXT",
    "FULL_SOLVENT",
})


def canonical_retained_generator_digest(
    retained_system_indices, retained_generator_ids
) -> str:
    """Return F12's raw digest for ordered retained-generator membership."""
    payload = json.dumps(
        {
            "system_indices": tuple(retained_system_indices),
            "generator_ids": tuple(retained_generator_ids),
        },
        separators=(",", ":"),
        ensure_ascii=True,
        sort_keys=True,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _validated_retained_values(system_indices, generator_ids, digest):
    if isinstance(system_indices, (str, bytes)) or isinstance(
        generator_ids, (str, bytes)
    ):
        return None
    try:
        indices = tuple(system_indices)
        generator_ids = tuple(generator_ids)
    except TypeError:
        return None
    if (
        any(isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in indices)
        or tuple(sorted(set(indices))) != indices
        or any(not isinstance(value, str) or not value for value in generator_ids)
        or len(indices) != len(generator_ids)
        or len(set(generator_ids)) != len(generator_ids)
        or not isinstance(digest, str)
        or not digest
        or digest != canonical_retained_generator_digest(indices, generator_ids)
    ):
        return None
    return indices, generator_ids, digest


def is_provisional_source_network_id(
    source_network_id: str,
    solve_universe_id: str,
) -> bool:
    """Return whether F7's temporary solve-derived namespace is being used."""
    return source_network_id == f"network:{solve_universe_id.removeprefix('solve:')}"


@dataclass(frozen=True)
class SharedNetworkIdentity:
    """Canonical identity projected into Results provenance.

    ``solve_universe_id`` is the F6 ``SolveContextKey.solve_universe_id``;
    ``selection_id`` is the F7 ``NetworkView.selection_id``.  Selection fields
    are intentionally absent from the solve identity payload.  A nonempty
    ``source_network_id`` is only a namespace; ``source_identity_status``
    carries the producer authority needed for certification.
    """

    source_network_id: str
    solve_universe_id: str
    frame_id: str
    analysis_context_kind: str
    geometry_scheme: str
    partition_scheme: str
    scope_mode: str
    source_identity_status: str = "UNRESOLVED"
    solvent_mode: str = "MOLECULE_ONLY"
    solvent_generator_digest: str | None = None
    solvent_cells_complete: bool | None = None
    retained_system_indices: tuple[int, ...] = ()
    retained_generator_ids: tuple[str, ...] = ()
    retained_generator_digest: str | None = None
    selection_id: str | None = None
    group_A_selection_id: str | None = None
    group_B_selection_id: str | None = None
    orientation: str | None = None

    def __post_init__(self):
        for name in (
            "source_network_id",
            "solve_universe_id",
            "frame_id",
            "geometry_scheme",
            "partition_scheme",
            "scope_mode",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a non-empty stable identifier")
        if self.analysis_context_kind not in ANALYSIS_CONTEXT_KINDS:
            raise ValueError("analysis_context_kind must be global or isolated")
        if self.source_identity_status not in SOURCE_IDENTITY_STATUSES:
            raise ValueError("Unknown source identity status")
        if (
            self.source_identity_status == "AUTHORITATIVE"
            and is_provisional_source_network_id(
                self.source_network_id, self.solve_universe_id
            )
        ):
            raise ValueError(
                "A provisional source namespace cannot be authoritative"
            )
        if not isinstance(self.solvent_mode, str) or self.solvent_mode not in SOLVENT_MODES:
            raise ValueError("Unknown solvent mode")
        if self.solvent_generator_digest is not None and (
            not isinstance(self.solvent_generator_digest, str)
            or not self.solvent_generator_digest
        ):
            raise ValueError("solvent_generator_digest must be a stable identifier")
        if self.solvent_cells_complete is not None and not isinstance(
            self.solvent_cells_complete, bool
        ):
            raise ValueError("solvent_cells_complete must be boolean or null")
        if self.solvent_mode == "MOLECULE_ONLY":
            if self.solvent_generator_digest or self.solvent_cells_complete:
                raise ValueError(
                    "MOLECULE_ONLY cannot contain participating solvent evidence"
                )
        elif not self.solvent_generator_digest:
            raise ValueError("Solvent-participating modes require generator digest")
        if self.solvent_mode == "SOLVENT_CONTEXT" and self.solvent_cells_complete:
            raise ValueError("SOLVENT_CONTEXT does not certify complete solvent cells")
        if self.solvent_mode == "FULL_SOLVENT" and self.solvent_cells_complete is not True:
            raise ValueError("FULL_SOLVENT requires complete solvent cells")
        if self.retained_system_indices is None or self.retained_generator_ids is None:
            raise ValueError("retained-generator evidence must be complete when supplied")
        retained_indices = tuple(self.retained_system_indices)
        retained_ids = tuple(self.retained_generator_ids)
        retention_supplied = bool(
            retained_indices or retained_ids or self.retained_generator_digest is not None
        )
        if retention_supplied:
            if _validated_retained_values(
                retained_indices, retained_ids, self.retained_generator_digest
            ) is None:
                raise ValueError("inconsistent retained-generator evidence")
        elif self.retained_generator_digest is not None:
            raise ValueError("retained-generator digest requires membership evidence")
        object.__setattr__(self, "retained_system_indices", retained_indices)
        object.__setattr__(self, "retained_generator_ids", retained_ids)
        if self.selection_id is not None and (
            not isinstance(self.selection_id, str)
            or not self.selection_id.startswith("view:")
        ):
            raise ValueError("selection_id must be the canonical view:<digest> ID")
        ordered = (
            self.group_A_selection_id,
            self.group_B_selection_id,
            self.orientation,
        )
        if any(value is not None for value in ordered) and (
            not self.selection_id or not all(ordered)
        ):
            raise ValueError(
                "A/B interface identity requires both ordered selections and orientation"
            )

    @property
    def provisional_source_network(self) -> bool:
        return is_provisional_source_network_id(
            self.source_network_id,
            self.solve_universe_id,
        )

    @property
    def authoritative_source_network(self) -> bool:
        """Whether producer authority, not namespace presence, is attested."""
        return (
            self.source_identity_status == "AUTHORITATIVE"
            and not self.provisional_source_network
        )

    def context_payload(self) -> dict[str, object]:
        """Return identity data that can describe a solve context.

        View/group selection IDs are excluded by construction.
        """
        value = {
            "source_network_id": self.source_network_id,
            "source_identity_status": self.source_identity_status,
            "solve_universe_id": self.solve_universe_id,
            "frame_id": self.frame_id,
            "analysis_context_kind": self.analysis_context_kind,
            "geometry_scheme": self.geometry_scheme,
            "partition_scheme": self.partition_scheme,
            "scope_mode": self.scope_mode,
            "solvent_mode": self.solvent_mode,
            "solvent_generator_digest": self.solvent_generator_digest,
            "solvent_cells_complete": self.solvent_cells_complete,
        }
        # Omit absent F12 evidence so legacy shared-network dependency IDs
        # remain stable; explicit retention is part of the solve payload.
        if self.retained_generator_digest is not None:
            value.update({
                "retained_system_indices": self.retained_system_indices,
                "retained_generator_ids": self.retained_generator_ids,
                "retained_generator_digest": self.retained_generator_digest,
            })
        return value

    def payload(self) -> dict[str, object]:
        value = self.context_payload()
        value.update({
            "selection_id": self.selection_id,
            "group_A_selection_id": self.group_A_selection_id,
            "group_B_selection_id": self.group_B_selection_id,
            "orientation": self.orientation,
        })
        return value


@dataclass(frozen=True)
class ResultsDependency:
    """Immutable dependency identity for a cached Results representation."""

    identity: SharedNetworkIdentity
    source_id: str
    source_revision: str | None = None
    source_digest: str | None = None
    upstream_dependency_digests: tuple[str, ...] = ()
    required: bool = True
    status: Status = Status.CERTIFIED

    def __post_init__(self):
        if not isinstance(self.source_id, str) or not self.source_id:
            raise ValueError("source_id must be a non-empty stable identifier")
        if not self.source_revision and not self.source_digest:
            raise ValueError("dependency requires a source revision or digest")
        if type(self.required) is not bool:
            raise ValueError("required must be boolean")
        if not isinstance(self.status, Status):
            raise ValueError("status must be a Results Status")
        if any(
            not isinstance(value, str) or not value
            for value in self.upstream_dependency_digests
        ):
            raise ValueError("upstream dependency digests must be non-empty strings")
        object.__setattr__(
            self,
            "upstream_dependency_digests",
            tuple(sorted(self.upstream_dependency_digests)),
        )

    @property
    def dependency_digest(self) -> str:
        """Canonical digest; status is evidence, not identity."""
        return stable_id(
            "dependency-v1",
            {
                "identity": self.identity.payload(),
                "source_id": self.source_id,
                "source_revision": self.source_revision,
                "source_digest": self.source_digest,
                "upstream_dependency_digests": self.upstream_dependency_digests,
                "required": self.required,
            },
        )


@dataclass(frozen=True)
class SharedResultsCollection:
    """Immutable handoff of Results sharing one solved source context."""

    results: tuple[ResultModel, ...]
    source_network_id: str
    solve_universe_id: str

    def __post_init__(self):
        if not isinstance(self.source_network_id, str) or not self.source_network_id:
            raise ValueError("source_network_id must be a non-empty stable identifier")
        if not isinstance(self.solve_universe_id, str) or not self.solve_universe_id:
            raise ValueError("solve_universe_id must be a non-empty stable identifier")
        results = tuple(self.results)
        if any(not isinstance(result, ResultModel) for result in results):
            raise TypeError("Results collection entries must be typed Results")
        selections = []
        for result in results:
            identity = result.identity
            if (
                identity.source_network_id != self.source_network_id
                or identity.solve_universe_id != self.solve_universe_id
                or not identity.selection_id
            ):
                raise ValueError("Results collection entries have incompatible lineage")
            selections.append(identity.selection_id)
        if len(set(selections)) != len(selections):
            raise ValueError("Results collection requires distinct selection identities")
        object.__setattr__(self, "results", results)

    def __iter__(self):
        return iter(self.results)

    def __len__(self):
        return len(self.results)

    def __getitem__(self, item):
        return self.results[item]

    @property
    def selection_ids(self) -> tuple[str, ...]:
        return tuple(result.identity.selection_id for result in self.results)

    def by_selection(self, selection_id: str) -> ResultModel:
        for result in self.results:
            if result.identity.selection_id == selection_id:
                return result
        raise KeyError(selection_id)


def collect_shared_results(
    results, *, source_network_id: str, solve_universe_id: str
) -> SharedResultsCollection:
    """Validate a producer tuple for NEXUS/Workbench handoff."""
    return SharedResultsCollection(
        tuple(results),
        source_network_id=source_network_id,
        solve_universe_id=solve_universe_id,
    )


def validate_results_dependency(
    expected: ResultsDependency,
    cached: ResultsDependency | None,
) -> Status:
    """Validate a cached dependency without refreshing or rebuilding it."""
    if cached is None:
        return Status.NOT_CALCULATED
    if not isinstance(cached, ResultsDependency):
        return Status.UNRESOLVED
    if not expected.identity.authoritative_source_network:
        return Status.UNRESOLVED
    if not cached.identity.authoritative_source_network:
        return Status.UNRESOLVED
    if (
        expected.identity.source_identity_status
        != cached.identity.source_identity_status
    ):
        return Status.UNRESOLVED
    if cached.status is Status.STALE or expected.status is Status.STALE:
        return Status.STALE
    if expected.dependency_digest != cached.dependency_digest:
        return Status.STALE
    if expected.status is not Status.CERTIFIED:
        return expected.status
    return cached.status


def validate_solve_context_identity(
    identity: SharedNetworkIdentity,
    solve_context,
) -> Status:
    """Validate Results identity against an explicit producer solve key.

    The producer owns the solve-key object.  This duck-typed boundary keeps
    Results independent of FORGE while requiring the canonical F6 fields.
    """
    if not isinstance(identity, SharedNetworkIdentity):
        return Status.UNRESOLVED
    if not identity.authoritative_source_network or solve_context is None:
        return Status.UNRESOLVED
    fields = (
        ("solve_universe_id", "solve_universe_id"),
        ("frame_id", "frame_id"),
        ("geometry_scheme", "scheme"),
        ("partition_scheme", "partition_scheme"),
        ("scope_mode", "scope_mode"),
    )
    try:
        if not all(
            getattr(identity, left) == getattr(solve_context, right)
            for left, right in fields
        ):
            return Status.UNRESOLVED
        solvent_fields = (
            "solvent_mode",
            "solvent_generator_digest",
            "solvent_cells_complete",
        )
        solvent_present = tuple(hasattr(solve_context, name) for name in solvent_fields)
        if any(solvent_present):
            if not all(solvent_present):
                return Status.UNRESOLVED
            if any(
                getattr(identity, name) != getattr(solve_context, name)
                for name in solvent_fields
            ):
                return Status.UNRESOLVED
        retained_provenance = getattr(
            solve_context, "retained_generator_provenance", None
        )
        if retained_provenance is not None and validate_retained_generator_provenance(
            identity, retained_provenance
        ) is not Status.CERTIFIED:
            return Status.UNRESOLVED
        return Status.CERTIFIED
    except AttributeError:
        return Status.UNRESOLVED


def validate_solvent_provenance(
    identity: SharedNetworkIdentity,
    provenance,
) -> Status:
    """Validate F10's explicit solvent provenance against Results identity."""
    if not isinstance(identity, SharedNetworkIdentity):
        return Status.UNRESOLVED
    if not isinstance(provenance, Mapping):
        return Status.UNRESOLVED
    fields = (
        "solvent_mode",
        "solvent_generator_digest",
        "solvent_cells_complete",
    )
    if any(field not in provenance for field in fields):
        return Status.UNRESOLVED
    return (
        Status.CERTIFIED
        if all(getattr(identity, field) == provenance[field] for field in fields)
        else Status.UNRESOLVED
    )


def validate_retained_generator_provenance(
    identity: SharedNetworkIdentity,
    provenance,
) -> Status:
    """Validate F12 retained-generator evidence without inspecting geometry."""
    if not isinstance(identity, SharedNetworkIdentity) or not isinstance(
        provenance, Mapping
    ):
        return Status.UNRESOLVED
    required = {
        "policy",
        "retained_system_indices",
        "retained_generator_ids",
        "retained_generator_digest",
        "boundary_independent",
    }
    if not required.issubset(provenance):
        return Status.UNRESOLVED
    if provenance["policy"] != "explicit_system_indices" or provenance[
        "boundary_independent"
    ] is not True:
        return Status.UNRESOLVED
    try:
        default_retention = (
            tuple(provenance["retained_system_indices"]) == ()
            and tuple(provenance["retained_generator_ids"]) == ()
            and provenance["retained_generator_digest"] is None
        )
    except TypeError:
        default_retention = False
    if default_retention:
        return (
            Status.CERTIFIED
            if not identity.retained_system_indices
            and not identity.retained_generator_ids
            and identity.retained_generator_digest is None
            else Status.UNRESOLVED
        )
    values = _validated_retained_values(
        provenance["retained_system_indices"],
        provenance["retained_generator_ids"],
        provenance["retained_generator_digest"],
    )
    if values is None:
        return Status.UNRESOLVED
    indices, generator_ids, digest = values
    return (
        Status.CERTIFIED
        if (
            identity.retained_system_indices == indices
            and identity.retained_generator_ids == generator_ids
            and identity.retained_generator_digest == digest
        )
        else Status.UNRESOLVED
    )


__all__ = [
    "ANALYSIS_CONTEXT_KINDS",
    "SOURCE_IDENTITY_STATUSES",
    "SOLVENT_MODES",
    "ResultsDependency",
    "SharedResultsCollection",
    "SharedNetworkIdentity",
    "canonical_retained_generator_digest",
    "collect_shared_results",
    "is_provisional_source_network_id",
    "validate_solve_context_identity",
    "validate_solvent_provenance",
    "validate_retained_generator_provenance",
    "validate_results_dependency",
]
