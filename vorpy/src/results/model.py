"""Immutable scientific result contracts, independent of kernels and exporters.

No model constructor computes geometry or aggregates scientific measurements.
IDs hash identity only, never values or file locations. Mapping payloads are
snapshotted recursively so later edits to a cache cannot mutate a result.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass, replace
from enum import Enum
from types import MappingProxyType
from typing import ClassVar

SCHEMA_VERSION = "1.0"


class Status(str, Enum):
    CERTIFIED = "CERTIFIED"
    PARTIAL = "PARTIAL"
    UNRESOLVED = "UNRESOLVED"
    NOT_CALCULATED = "NOT_CALCULATED"
    NOT_SUPPORTED = "NOT_SUPPORTED"
    STALE = "STALE"


class GeometryLineageOrigin(str, Enum):
    """How the geometry backing a result was obtained.

    This is provenance, not scientific certification.  In particular,
    ``EXACT_REUSE`` never promotes a quantity's :class:`Status`.
    """

    DIRECT_SOLVE = "DIRECT_SOLVE"
    EXACT_REUSE = "EXACT_REUSE"
    MATERIALIZED_VIEW = "MATERIALIZED_VIEW"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True)
class GeometryLineage:
    """Validated, optional lineage for geometry-backed scientific results.

    The record is intentionally separate from ``ResultIdentity`` and
    ``Status``. A legacy result with no ``geometry_lineage`` detail remains
    a legacy result; it is not silently labelled as a direct solve.
    """

    origin: GeometryLineageOrigin | str
    schema_version: int = 1
    solve_universe_id: str | None = None
    geometry_source_id: str | None = None
    geometry_version: str | None = None
    solver_settings_digest: str | None = None
    solver_settings: Mapping | None = None
    solver_settings_ref: str | None = None
    environment_digest: str | None = None
    environment: Mapping | None = None
    environment_ref: str | None = None
    source_network_id: str | None = None
    source_interface_id: str | None = None
    source_view_id: str | None = None
    selection_scope: Mapping | None = None
    orientation: Mapping | None = None
    compatibility: Mapping | None = None
    completeness: Mapping | None = None
    reason: str | None = None

    _IDENTIFIERS = (
        "solve_universe_id",
        "geometry_source_id",
        "geometry_version",
        "solver_settings_digest",
        "solver_settings_ref",
        "environment_digest",
        "environment_ref",
        "source_network_id",
        "source_interface_id",
        "source_view_id",
        "reason",
    )
    _MAPPINGS = (
        "solver_settings",
        "environment",
        "selection_scope",
        "orientation",
        "compatibility",
        "completeness",
    )

    def __post_init__(self):
        try:
            origin = GeometryLineageOrigin(self.origin)
        except (TypeError, ValueError) as exc:
            raise ValueError("Unknown geometry lineage origin") from exc
        object.__setattr__(self, "origin", origin)
        if type(self.schema_version) is not int or self.schema_version != 1:
            raise ValueError("Unsupported geometry lineage schema version")

        for name in self._IDENTIFIERS:
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value):
                raise ValueError(f"{name} must be null or a non-empty string")

        for name in self._MAPPINGS:
            value = getattr(self, name)
            if value is not None:
                if not isinstance(value, Mapping):
                    raise TypeError(f"{name} must be a JSON mapping or null")
                if name in {"selection_scope", "orientation", "completeness"} and not value:
                    raise ValueError(f"{name} must be explicit when supplied")
                object.__setattr__(self, name, freeze(value))

        if self.origin is GeometryLineageOrigin.UNRESOLVED:
            if not self.reason:
                raise ValueError("UNRESOLVED geometry lineage requires a reason")
            return

        required = {
            "solve_universe_id": self.solve_universe_id,
            "geometry_source_id": self.geometry_source_id,
            "geometry_version": self.geometry_version,
            "solver_settings_digest": self.solver_settings_digest,
            "environment_digest": self.environment_digest,
            "source_network_id": self.source_network_id,
            "selection_scope": self.selection_scope,
            "orientation": self.orientation,
            "completeness": self.completeness,
        }
        missing = tuple(name for name, value in required.items() if value is None)
        if missing:
            raise ValueError(
                "Resolved geometry lineage is incomplete: " + ", ".join(missing)
            )
        if self.solver_settings is None and self.solver_settings_ref is None:
            raise ValueError(
                "Resolved geometry lineage requires solver settings or a resolvable reference"
            )
        if self.environment is None and self.environment_ref is None:
            raise ValueError(
                "Resolved geometry lineage requires environment metadata or a resolvable reference"
            )
        if not isinstance(self.completeness.get("completeness_proven"), bool):
            raise ValueError("completeness requires a boolean completeness_proven field")

        if self.origin is GeometryLineageOrigin.EXACT_REUSE:
            if self.source_view_id is None:
                raise ValueError("EXACT_REUSE requires source_view_id")
            if not isinstance(self.compatibility, Mapping):
                raise ValueError("EXACT_REUSE requires compatibility evidence")
            required_compatibility = {
                "decision",
                "generator_compatibility",
                "competitor_compatibility",
                "solve_context",
            }
            missing = required_compatibility.difference(self.compatibility)
            if missing:
                raise ValueError(
                    "EXACT_REUSE compatibility is incomplete: "
                    + ", ".join(sorted(missing))
                )
            if self.compatibility["decision"] != GeometryLineageOrigin.EXACT_REUSE.value:
                raise ValueError("EXACT_REUSE compatibility has the wrong decision")
            for name in required_compatibility - {"decision"}:
                evidence = self.compatibility[name]
                if not isinstance(evidence, Mapping) or evidence.get("match") is not True:
                    raise ValueError(
                        f"EXACT_REUSE requires positive {name} evidence"
                    )
            if self.completeness["completeness_proven"] is not True:
                raise ValueError("EXACT_REUSE requires proven completeness")

    @classmethod
    def from_mapping(cls, value):
        """Construct lineage from its JSON-safe mapping representation."""
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError("geometry_lineage must be a mapping")
        names = {item.name for item in fields(cls) if not item.name.startswith("_")}
        unknown = set(value).difference(names)
        if unknown:
            raise ValueError(
                "Unknown geometry lineage fields: " + ", ".join(sorted(map(str, unknown)))
            )
        return cls(**{name: value[name] for name in names if name in value})

    def as_dict(self):
        """Return detached JSON-safe lineage state, including explicit nulls."""
        return dictionary(self)


class MolecularSurfaceConvention(str, Enum):
    EXPANDED_UNION_OF_BALLS = "expanded_union_of_balls"
    SOLVENT_ACCESSIBLE = "solvent_accessible"
    SOLVENT_EXCLUDED = "solvent_excluded"


class OrientationConvention(str, Enum):
    MOLECULAR_OUTWARD_NORMAL = "outward_molecular_normal"
    PHYSICAL_INTERFACE_A_TO_B = "physical_interface_A_to_B"


INTERFACE_LOG_BASE_COLUMNS = (
    "interface_id",
    "representation",
    "side",
    "quantity",
    "value",
    "units",
    "status",
)
"""Seven-column payload prefix; not a valid ``interfaces/logs.csv`` header."""

INTERFACE_LOG_COLUMNS = INTERFACE_LOG_BASE_COLUMNS + ("provenance",)
"""Authoritative eight-column CSV header; provenance is nullable but present."""

INTERFACE_LOG_STATUS_VALUES = tuple(status.value for status in Status)

INTERFACE_LOG_REPRESENTATIONS = frozenset(
    {
        "physical_partition_interface",
        "molecular_contact_surface",
        "dual_contact_complex",
        "dual_generator_complex",
        "dual_separator",
        "alpha_selection",
        "water_analysis",
        "timing",
    }
)

RESULT_REPRESENTATION_TO_INTERFACE_LOG = MappingProxyType(
    {
        "voronoi": "physical_partition_interface",
        "dual": "dual_contact_complex",
        "dual_a": "dual_generator_complex",
        "dual_b": "dual_generator_complex",
        "dual_generator_complex": "dual_generator_complex",
        "dual_separator": "dual_separator",
        "molecular_contact_surface": "molecular_contact_surface",
    }
)

INTERFACE_LOG_QUANTITY_UNITS = MappingProxyType(
    {
        "area": "Å²",
        "contact_count": "1",
        "residue_contact_count": "1",
        "selected_ab_contact_count": "1",
        "direct_contact_atom_count": "1",
        "contributing_surface_atom_count": "1",
        "occluding_same_group_atom_count": "1",
        "partner_clipping_atom_count": "1",
        "unresolved_contact_count": "1",
        "carrier_atom_count": "1",
        "geometric_edge_count": "1",
        "spherical_patch_count": "1",
        "circular_arc_count": "1",
        "junction_record_count": "1",
        "refined_seam_count": "1",
        "boundary_arc_count": "1",
        "vertex_count": "1",
        "edge_count": "1",
        "topological_edge_count": "1",
        "topological_cut_count": "1",
        "patch_attachment_record_count": "1",
        "face_count": "1",
        "connected_component_count": "1",
        "boundary_component_count": "1",
        "euler_characteristic": "1",
        "orientable": "1",
        "genus": "1",
        "nonorientable_genus": "1",
        "manifold_vertex_count": "1",
        "singular_vertex_count": "1",
        "chain_complex_valid": "1",
        "boundary_of_boundary_zero": "1",
        "beta0": "1",
        "beta1": "1",
        "beta2": "1",
        "mean_curvature_per_area": "Å^-1",
        "gaussian_curvature_per_area": "Å^-2",
        "smooth_integrated_mean_curvature": "Å",
        "physical_seam_integrated_mean_curvature": "Å",
        "total_integrated_mean_curvature": "Å",
        "complete_integrated_mean_curvature": "Å",
        "boundary_scalar_mean_curvature": "Å",
        "intrinsic_integrated_gaussian_curvature": "1",
        "smooth_integrated_gaussian_curvature": "1",
        "intrinsic_discrete_gaussian_curvature": "1",
        "intrinsic_seam_gaussian_curvature": "1",
        "exterior_boundary_geodesic_curvature": "1",
        "interior_junction_gaussian_defect": "1",
        "boundary_corner_turning": "1",
        "completed_gauss_bonnet": "1",
        "gauss_bonnet_expected": "1",
        "gauss_bonnet_residual": "1",
        "alpha_candidate_count": "1",
        "alpha_selected_count": "1",
        "alpha_selected_generator_pairs": "1",
        "alpha_selected_system_pairs": "1",
        "alpha_selected_physical_system_pairs": "1",
    }
)

INTERFACE_LOG_NATIVE_UNIT_QUANTITIES = frozenset({"alpha_value"})

RESULT_CURVATURE_QUANTITY_TO_INTERFACE_LOG = MappingProxyType(
    {
        "mean_curvature.integrated_mean_curvature": "total_integrated_mean_curvature",
        "mean_curvature.complete_total": "complete_integrated_mean_curvature",
        "mean_curvature.smooth_face_contribution": "smooth_integrated_mean_curvature",
        "mean_curvature.edge_contribution": "physical_seam_integrated_mean_curvature",
        "mean_curvature.boundary_contribution": "boundary_scalar_mean_curvature",
        "gaussian_curvature.integrated_gaussian_curvature": "intrinsic_integrated_gaussian_curvature",
        "gaussian_curvature.smooth_face_contribution": "smooth_integrated_gaussian_curvature",
        "gaussian_curvature.intrinsic_discrete_contribution": "intrinsic_discrete_gaussian_curvature",
        "gaussian_curvature.intrinsic_seam_contribution": "intrinsic_seam_gaussian_curvature",
        "gaussian_curvature.interior_vertex_defects": "interior_junction_gaussian_defect",
        "gaussian_curvature.boundary_geodesic_contribution": "exterior_boundary_geodesic_curvature",
        "gaussian_curvature.boundary_corner_contribution": "boundary_corner_turning",
        "gaussian_curvature.gauss_bonnet_completion": "completed_gauss_bonnet",
        "gaussian_curvature.gauss_bonnet_expected": "gauss_bonnet_expected",
        "gaussian_curvature.gauss_bonnet_residual": "gauss_bonnet_residual",
    }
)

RESULT_MOLECULAR_QUANTITY_TO_INTERFACE_LOG = MappingProxyType(
    {
        "molecular_selection.selected_contacts": "selected_ab_contact_count",
        "molecular_selection.unresolved_contacts": "unresolved_contact_count",
    }
)

RESULT_METRIC_QUANTITY_TO_INTERFACE_LOG = MappingProxyType(
    {
        "metrics.area": "area",
        "metrics.n_contacts": "contact_count",
        "metrics.n_residue_contacts": "residue_contact_count",
    }
)


INTERACTION_COMPARISON_QUANTITIES = (
    "area_A_contact",
    "area_partition_interface",
    "area_B_contact",
    "area_ratio_A_to_partition",
    "area_ratio_B_to_partition",
    "area_ratio_A_to_B",
    "integrated_mean_curvature_A",
    "integrated_mean_curvature_partition",
    "integrated_mean_curvature_B",
    "integrated_gaussian_curvature_A",
    "integrated_gaussian_curvature_partition",
    "integrated_gaussian_curvature_B",
)


_ENVIRONMENT_ALIASES = {
    "dry": "dry",
    "solvent_competing": "solvent_competing",
}
_NETWORK_SCOPE_ALIASES = {
    "system": "system",
    "interface": "interface",
    "system network": "system",
    "interface/dedicated": "interface",
}
GENERATOR_COMPLEX_SIDES = frozenset({"A", "B", "AB/contact"})


def normalize_environment(value):
    try:
        return _ENVIRONMENT_ALIASES[value]
    except (KeyError, TypeError):
        raise ValueError(
            "Environment must be dry or solvent_competing"
        ) from None


def normalize_network_scope(value):
    try:
        return _NETWORK_SCOPE_ALIASES[value]
    except (KeyError, TypeError):
        raise ValueError("Network scope must be system or interface") from None


def freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({str(k): freeze(v) for k, v in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(freeze(v) for v in value)
    if value is None or isinstance(value, (str, bool, int, float, Enum)):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError("Nonfinite values cannot enter a scientific result")
        return value
    if is_dataclass(value) and value.__dataclass_params__.frozen:
        return value
    raise TypeError(
        f"Result payload is not immutable JSON state: {type(value).__name__}"
    )


def dictionary(value):
    """Return detached, JSON-safe Python state; performs no file IO."""
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        omitted = {
            "contact_selection_source",
            "contact_selection_id",
            "molecular_surface_convention",
            "definition_version",
            "partner_id",
            "source_network_id",
            "source_interface_id",
        "interaction_id",
        "orientation_convention",
        "orientation_convention",
        }
        return {
            f.name: dictionary(getattr(value, f.name))
            for f in fields(value)
            if not (
                type(value).__name__ == "ResultIdentity"
                and f.name in omitted
                and getattr(value, f.name) is None
            )
        }
    if isinstance(value, Mapping):
        return {str(k): dictionary(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [dictionary(v) for v in value]
    return value


def stable_id(prefix, payload):
    encoded = json.dumps(
        dictionary(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return prefix + ":" + hashlib.sha256(encoded).hexdigest()


def _stale(value):
    """Invalidate certification without modifying retained cached measurements."""
    if isinstance(value, Quantity):
        return replace(value, status=Status.STALE)
    if type(value).__name__ in {
        "MolecularContactSelection",
        "MolecularContactSurfaceGeometry",
        "MolecularPatchProvenance",
        "InteractionLinks",
        "MolecularGeometryCounts",
        "TopologicalCellComplex",
    }:
        updates = {f.name: _stale(getattr(value, f.name)) for f in fields(value)}
        if "status" in updates:
            updates["status"] = Status.STALE
        return replace(value, **updates)
    if isinstance(value, (MeanCurvature, GaussianCurvature, Topology, ContactRecord)):
        updates = {f.name: _stale(getattr(value, f.name)) for f in fields(value)}
        if isinstance(value, Topology):
            updates["orientation_status"] = Status.STALE
        if isinstance(value, ContactRecord):
            updates["selection_status"] = "STALE"
        return replace(value, **updates)
    if isinstance(value, Mapping):
        return {k: _stale(v) for k, v in value.items()}
    if isinstance(value, tuple):
        return tuple(_stale(v) for v in value)
    return value


@dataclass(frozen=True)
class Provenance:
    source: str
    source_revision: str | None = None
    scientific_state_id: str | None = None
    details: Mapping = field(default_factory=dict)
    diagnostics: tuple[str, ...] = ()

    def __post_init__(self):
        if not self.source:
            raise ValueError("Provenance requires a source")
        if not isinstance(self.details, Mapping):
            raise TypeError("Provenance details must be a mapping")
        details = dict(self.details)
        if "geometry_lineage" in details and details["geometry_lineage"] is not None:
            details["geometry_lineage"] = GeometryLineage.from_mapping(
                details["geometry_lineage"]
            ).as_dict()
        object.__setattr__(self, "details", freeze(details))
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))

    @property
    def geometry_lineage(self):
        """Return validated lineage, or ``None`` when legacy data omitted it."""
        value = self.details.get("geometry_lineage")
        return None if value is None else GeometryLineage.from_mapping(value)


@dataclass(frozen=True)
class ResultIdentity:
    result_kind: str
    system: str
    frame: str | int | None
    network_id: str
    network_scope: str
    environment: str
    partition: str
    representation: str
    radius_configuration_id: str
    group_id: str | None = None
    interface_id: str | None = None
    group_A: str | None = None
    group_B: str | None = None
    side: str | None = None
    orientation: str | None = None
    alpha_method: str | None = None
    alpha_value: float | None = None
    alpha_units: str | None = None
    probe_radius: float | None = None
    probe_units: str | None = None
    system_id: str | None = None
    contact_selection_source: str | None = None
    contact_selection_id: str | None = None
    molecular_surface_convention: str | None = None
    definition_version: str | None = None
    partner_id: str | None = None
    source_network_id: str | None = None
    source_interface_id: str | None = None
    interaction_id: str | None = None
    orientation_convention: str | None = None

    def __post_init__(self):
        if self.result_kind not in {
            "system",
            "network",
            "group",
            "interface",
            "molecular_contact_surface",
        }:
            raise ValueError("Unknown result kind")
        for name in ("system", "network_id", "radius_configuration_id"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name):
                raise ValueError(f"{name} must be an explicit stable identifier")
        if self.frame is not None and (
            isinstance(self.frame, bool) or not isinstance(self.frame, (str, int))
        ):
            raise ValueError("Frame must be a string, integer, or null")
        for name in ("alpha_value", "probe_radius"):
            value = getattr(self, name)
            if value is not None:
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(value)
                ):
                    raise ValueError(f"{name} must be a finite number")
                object.__setattr__(self, name, float(value) if value != 0 else 0.0)
        object.__setattr__(self, "network_scope", normalize_network_scope(self.network_scope))
        object.__setattr__(self, "environment", normalize_environment(self.environment))
        if self.partition not in {"aw", "power"}:
            raise ValueError("Use canonical partition names")
        if self.representation not in {
            "voronoi",
            "dual_a",
            "dual_b",
            "dual_generator_complex",
            "dual_separator",
            "alpha_selection",
            "molecular_contact_surface",
            "dual",
            "molecular",
        }:
            raise ValueError("Unknown representation")
        if self.side not in {None, "shared", "A", "B", "AB/contact"}:
            raise ValueError("Unknown side")
        if self.system_id is not None and not self.system_id:
            raise ValueError("Stable system namespace cannot be empty")
        if self.result_kind == "group" and not self.group_id:
            raise ValueError("Group results require group_id")
        if self.result_kind == "interface":
            if not all(
                (self.interface_id, self.group_A, self.group_B, self.orientation)
            ):
                raise ValueError(
                    "Interface identity requires groups, ID, and orientation"
                )
            expected = {
                "voronoi": {"shared"},
                "dual": {"shared"},
                "dual_a": {"A"},
                "dual_b": {"B"},
                "alpha_selection": {"shared"},
            }
            if self.representation in expected:
                valid_sides = expected[self.representation]
            elif self.representation == "dual_generator_complex":
                valid_sides = GENERATOR_COMPLEX_SIDES
            elif self.representation == "dual_separator":
                valid_sides = {"shared", "AB/contact"}
            else:
                valid_sides = set()
            if self.side not in valid_sides:
                raise ValueError("Interface representation and side disagree")
        if self.molecular_surface_convention is not None:
            object.__setattr__(
                self,
                "molecular_surface_convention",
                MolecularSurfaceConvention(self.molecular_surface_convention).value,
            )
        if self.orientation_convention is not None:
            object.__setattr__(
                self,
                "orientation_convention",
                OrientationConvention(self.orientation_convention).value,
            )
        if self.result_kind == "molecular_contact_surface" and self.orientation_convention != OrientationConvention.MOLECULAR_OUTWARD_NORMAL.value:
            raise ValueError("molecular contact surfaces require outward molecular normals")
        if self.result_kind == "molecular_contact_surface":
            required = {
                "system_id": self.system_id,
                "interface_id": self.interface_id,
                "group_A": self.group_A,
                "group_B": self.group_B,
                "orientation": self.orientation,
                "partner_id": self.partner_id,
                "contact_selection_source": self.contact_selection_source,
                "contact_selection_id": self.contact_selection_id,
                "molecular_surface_convention": self.molecular_surface_convention,
                "definition_version": self.definition_version,
                "source_network_id": self.source_network_id,
                "source_interface_id": self.source_interface_id,
                "interaction_id": self.interaction_id,
            }
            if any(not isinstance(value, str) or not value for value in required.values()):
                raise ValueError(
                    "Molecular contact-surface identity is incomplete"
                )
            if self.representation != "molecular_contact_surface" or self.side not in {"A", "B"}:
                raise ValueError("Molecular contact surfaces require side A or B")
        if self.alpha_value is not None and (
            not self.alpha_method
            or not self.alpha_units
            or not math.isfinite(self.alpha_value)
        ):
            raise ValueError("Alpha value requires finite value, convention, and units")
        if self.probe_radius is not None and (
            not self.probe_units
            or not math.isfinite(self.probe_radius)
            or self.probe_radius < 0
        ):
            raise ValueError("Probe requires finite nonnegative radius and units")

    @property
    def result_id(self):
        return stable_id("result-v1", self)


@dataclass(frozen=True)
class Quantity:
    name: str
    value: float | int | bool | None
    units: str
    status: Status
    method: str
    scope: str
    support_count: int | None = None
    total_count: int | None = None
    provenance: Provenance | None = None

    def __post_init__(self):
        object.__setattr__(self, "status", Status(self.status))
        if not all((self.name, self.units, self.method, self.scope)):
            raise ValueError("Quantities require name, units, method, and scope")
        if self.value is not None and (
            not isinstance(self.value, (int, float, bool))
            or not math.isfinite(self.value)
        ):
            raise ValueError(
                "Quantity values must be finite scalar numbers or booleans"
            )
        if self.status in {Status.CERTIFIED, Status.PARTIAL} and self.value is None:
            raise ValueError("Certified/partial quantities require a value")
        if (
            self.status
            in {Status.UNRESOLVED, Status.NOT_CALCULATED, Status.NOT_SUPPORTED}
            and self.value is not None
        ):
            raise ValueError(
                "Unavailable quantities must use null, never zero or a partial sum"
            )
        for count in (self.support_count, self.total_count):
            if count is not None and (
                isinstance(count, bool) or not isinstance(count, int) or count < 0
            ):
                raise ValueError("Coverage counts must be nonnegative integers")
        if (
            self.support_count is not None
            and self.total_count is not None
            and self.support_count > self.total_count
        ):
            raise ValueError("Support count exceeds total count")


def canonical_interface_log_representation(representation):
    """Return the stable long-form log name for a Results representation."""
    try:
        return RESULT_REPRESENTATION_TO_INTERFACE_LOG[representation]
    except KeyError:
        if representation in INTERFACE_LOG_REPRESENTATIONS:
            return representation
        raise ValueError(f"Unsupported interface-log representation: {representation}") from None


def _interface_log_json_value(value):
    value = dictionary(value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("Interface-log values must be finite")
        return value
    if isinstance(value, (list, dict)):
        json.dumps(value, ensure_ascii=False, allow_nan=False)
        return value
    raise TypeError(f"Unsupported interface-log value type: {type(value).__name__}")


def serialize_interface_log_value(value):
    """Serialize a scalar or structured log value deterministically for CSV."""
    value = _interface_log_json_value(value)
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)


def serialize_interface_log_provenance(provenance):
    """Serialize optional provenance as canonical JSON, or an empty CSV field."""
    if provenance is None:
        return ""
    if not isinstance(provenance, (Provenance, Mapping)):
        raise TypeError("Interface-log provenance must be Provenance or a mapping")
    value = _interface_log_json_value(provenance)
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"), sort_keys=True)


@dataclass(frozen=True)
class InterfaceLogRow:
    """Validated long-form row consumed by the consolidated interface log writer."""

    interface_id: str
    representation: str
    side: str
    quantity: str
    value: object | None
    units: str
    status: Status
    provenance: Provenance | Mapping | None = None

    def __post_init__(self):
        if not self.interface_id or not self.quantity or not self.units:
            raise ValueError("Interface-log rows require interface, quantity, and units")
        object.__setattr__(
            self,
            "representation",
            canonical_interface_log_representation(self.representation),
        )
        if self.side not in {"shared", "A", "B", "AB/contact"}:
            raise ValueError("Interface-log rows require a canonical side")
        object.__setattr__(self, "status", Status(self.status))
        if self.status in {Status.CERTIFIED, Status.PARTIAL} and self.value is None:
            raise ValueError("Certified/partial log rows require a value")
        if self.status in {
            Status.UNRESOLVED,
            Status.NOT_CALCULATED,
            Status.NOT_SUPPORTED,
        } and self.value is not None:
            raise ValueError("Unavailable log rows require an empty value")
        expected_units = INTERFACE_LOG_QUANTITY_UNITS.get(self.quantity)
        if expected_units is not None and self.units != expected_units:
            raise ValueError(f"Invalid units for interface-log quantity: {self.quantity}")
        if self.quantity == "alpha_value":
            pass  # Alpha retains the producer's explicit native units.
        elif self.quantity.startswith("alpha_"):
            if self.units != "1":
                raise ValueError("Alpha counts and selections are dimensionless")
        elif self.quantity.startswith("timing_"):
            if self.units != "s":
                raise ValueError("Timing quantities require seconds")
        elif self.quantity.startswith("water_"):
            if self.units != "1":
                raise ValueError("Water-analysis quantities require dimensionless units")
        elif expected_units is None:
            raise ValueError(f"Unsupported interface-log quantity: {self.quantity}")
        _interface_log_json_value(self.value)
        if self.provenance is not None and not isinstance(self.provenance, (Provenance, Mapping)):
            raise TypeError("Interface-log provenance must be Provenance or a mapping")

    def as_csv_row(self):
        return {
            "interface_id": self.interface_id,
            "representation": self.representation,
            "side": self.side,
            "quantity": self.quantity,
            "value": serialize_interface_log_value(self.value),
            "units": self.units,
            "status": self.status.value,
            "provenance": serialize_interface_log_provenance(self.provenance),
        }


@dataclass(frozen=True)
class Environment:
    environment: str
    explicit_solvent_present: bool | None = None
    solvent_generator_count: int | None = None
    solvent_interface_competitors: int | None = None
    provenance: Provenance | None = None

    def __post_init__(self):
        object.__setattr__(self, "environment", normalize_environment(self.environment))
        for count in (self.solvent_generator_count, self.solvent_interface_competitors):
            if count is not None and (
                isinstance(count, bool) or not isinstance(count, int) or count < 0
            ):
                raise ValueError("Solvent counts must be nonnegative integers")
        if self.solvent_generator_count is not None and (
            self.solvent_generator_count > 0
        ) != (self.environment == "solvent_competing"):
            raise ValueError(
                "Environment contradicts participating solvent generator count"
            )
        if self.environment == "dry" and self.solvent_interface_competitors:
            raise ValueError("Dry solve cannot have solvent interface competitors")
        if self.explicit_solvent_present is False and (
            self.solvent_generator_count or self.solvent_interface_competitors
        ):
            raise ValueError(
                "Competing explicit solvent contradicts recorded solvent absence"
            )
        if (
            self.solvent_interface_competitors is not None
            and self.solvent_generator_count is not None
        ) and self.solvent_interface_competitors > self.solvent_generator_count:
            raise ValueError("Interface competitors exceed solvent generators")


@dataclass(frozen=True)
class MeanCurvature:
    """Signed integrated mean-curvature accounting.

    ``integrated_mean_curvature`` is the total H integral when the requested
    surface is complete; ``complete_total`` is the separately named complete
    total and must stay null/status-bearing when completeness is not certified.
    ``smooth_face_contribution`` is the smooth face term, ``edge_contribution``
    is the signed physical-seam term (one half of the seam dihedral integral),
    and ``boundary_contribution`` is a scalar boundary term only when a
    producer actually calculates one. Unsupported free-boundary omission is
    ``None/NOT_SUPPORTED``, never ``0.0/CERTIFIED``.
    """

    integrated_mean_curvature: Quantity
    complete_total: Quantity
    smooth_face_contribution: Quantity
    edge_contribution: Quantity
    boundary_contribution: Quantity
    convention: str

    def __post_init__(self):
        expected = {
            "integrated_mean_curvature": "integrated_mean_curvature",
            "complete_total": "complete_mean_curvature",
            "smooth_face_contribution": "smooth_face_contribution",
            "edge_contribution": "edge_contribution",
            "boundary_contribution": "boundary_contribution",
        }
        for attr, name in expected.items():
            q = getattr(self, attr)
            if not isinstance(q, Quantity) or q.name != name or q.units != "Å":
                raise ValueError("Invalid mean-curvature quantity contract")
        if not self.convention:
            raise ValueError("Mean curvature requires a convention")


@dataclass(frozen=True)
class GaussianCurvature:
    """Intrinsic Gaussian accounting before and after boundary completion.

    ``integrated_gaussian_curvature`` is the intrinsic surface total before
    exterior boundary geodesic and corner terms. It is not the smooth face
    term and it is not the completed Gauss--Bonnet total.
    ``smooth_face_contribution`` is smooth integrated K; intrinsic seam and
    interior-junction terms are retained separately.
    ``gauss_bonnet_completion``, ``gauss_bonnet_expected``, and
    ``gauss_bonnet_residual`` are independent audit quantities.
    """

    integrated_gaussian_curvature: Quantity
    smooth_face_contribution: Quantity
    intrinsic_discrete_contribution: Quantity
    intrinsic_seam_contribution: Quantity
    interior_vertex_defects: Quantity
    boundary_geodesic_contribution: Quantity
    boundary_corner_contribution: Quantity
    gauss_bonnet_completion: Quantity
    gauss_bonnet_expected: Quantity
    gauss_bonnet_residual: Quantity
    convention: str

    def __post_init__(self):
        for f in fields(self):
            if f.name != "convention":
                q = getattr(self, f.name)
                if not isinstance(q, Quantity) or q.name != f.name or q.units != "1":
                    raise ValueError("Invalid Gaussian-curvature quantity contract")
        if not self.convention:
            raise ValueError("Gaussian curvature requires a convention")

    @property
    def intrinsic_total_gaussian_curvature(self):
        """Explicit semantic alias for the pre-boundary intrinsic total."""
        return self.integrated_gaussian_curvature

    @property
    def smooth_integrated_gaussian_curvature(self):
        """Explicit semantic alias for smooth integrated K."""
        return self.smooth_face_contribution


@dataclass(frozen=True)
class Topology:
    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    orientation_status: Status = Status.NOT_CALCULATED
    genus_by_component: tuple[int, ...] | None = None
    genus_method: str | None = None

    def __post_init__(self):
        object.__setattr__(self, "quantities", freeze(self.quantities))
        object.__setattr__(self, "orientation_status", Status(self.orientation_status))
        _quantity_mapping(self.quantities)
        if self.genus_by_component is not None:
            if not self.genus_method:
                raise ValueError(
                    "Genus requires an explicit scientific justification/method"
                )
            if any(
                isinstance(g, bool) or not isinstance(g, int) or g < 0
                for g in self.genus_by_component
            ):
                raise ValueError("Invalid component genus")
            object.__setattr__(
                self, "genus_by_component", tuple(self.genus_by_component)
            )


@dataclass(frozen=True)
class GeneratorComplex:
    """Dimension-counted generator complex; not assumed to be a surface."""

    side: str
    simplex_counts: Mapping[int, Quantity] = field(default_factory=dict)
    topology: Mapping[str, Quantity] = field(default_factory=dict)
    status: Status = Status.NOT_CALCULATED
    provenance: Provenance | None = None

    def __post_init__(self):
        if self.side not in GENERATOR_COMPLEX_SIDES:
            raise ValueError("Unknown generator-complex side")
        object.__setattr__(self, "status", Status(self.status))
        for dimension, quantity in self.simplex_counts.items():
            if isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 0:
                raise ValueError("Simplex dimensions must be nonnegative integers")
            if not isinstance(quantity, Quantity) or quantity.units != "1":
                raise ValueError("Simplex counts require dimensionless quantities")
            if quantity.name != f"simplex_count_dim_{dimension}":
                raise ValueError("Simplex counts must use dimension-specific names")
            if quantity.value is not None and (
                isinstance(quantity.value, bool)
                or not isinstance(quantity.value, int)
                or quantity.value < 0
            ):
                raise ValueError("Simplex counts must be nonnegative integers")
        # Dimensions are typed integer keys in the in-memory contract.  The
        # generic JSON freezer stringifies mapping keys, so use a read-only
        # mapping here and let ``dictionary()`` perform serialization later.
        object.__setattr__(self, "simplex_counts", MappingProxyType(dict(self.simplex_counts)))
        object.__setattr__(self, "topology", freeze(self.topology))
        _quantity_mapping(self.topology)


def _string_tuple(values, label):
    values = tuple(values)
    if any(not isinstance(value, str) or not value for value in values):
        raise ValueError(f"{label} must contain nonempty strings")
    return values


@dataclass(frozen=True)
class MolecularContactSelection:
    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    selected_ab_contact_ids: tuple[str, ...] = ()
    direct_contact_atom_ids: tuple[str, ...] = ()
    contributing_surface_atom_ids: tuple[str, ...] = ()
    occluding_same_group_atom_ids: tuple[str, ...] = ()
    partner_clipping_atom_ids: tuple[str, ...] = ()
    unresolved_contact_ids: tuple[str, ...] = ()
    status: Status = Status.NOT_CALCULATED

    def __post_init__(self):
        object.__setattr__(self, "status", Status(self.status))
        object.__setattr__(self, "quantities", freeze(self.quantities))
        _quantity_mapping(self.quantities)
        for name in (
            "selected_ab_contact_ids",
            "direct_contact_atom_ids",
            "contributing_surface_atom_ids",
            "occluding_same_group_atom_ids",
            "partner_clipping_atom_ids",
            "unresolved_contact_ids",
        ):
            object.__setattr__(self, name, _string_tuple(getattr(self, name), name))


@dataclass(frozen=True)
class MolecularContactSurfaceGeometry:
    """Optional geometry quantities; no surface field is mandatory."""

    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    status: Status = Status.NOT_CALCULATED

    def __post_init__(self):
        object.__setattr__(self, "status", Status(self.status))
        object.__setattr__(self, "quantities", freeze(self.quantities))
        _quantity_mapping(self.quantities)


MOLECULAR_GEOMETRY_QUANTITIES = frozenset(
    {
        "carrier_atom_count",
        "geometric_edge_count",
        "spherical_patch_count",
        "circular_arc_count",
        "junction_record_count",
        "refined_seam_count",
        "boundary_arc_count",
    }
)

TOPOLOGY_CELL_QUANTITIES = frozenset(
    {
        "vertex_count",
        "edge_count",
        "topological_edge_count",
        "topological_cut_count",
        "patch_attachment_record_count",
        "face_count",
        "connected_component_count",
        "boundary_component_count",
        "euler_characteristic",
        "orientable",
        "genus",
        "nonorientable_genus",
        "manifold_vertex_count",
        "singular_vertex_count",
        "chain_complex_valid",
        "boundary_of_boundary_zero",
        "beta0",
        "beta1",
        "beta2",
    }
)


@dataclass(frozen=True)
class MolecularGeometryCounts:
    """Counts of geometric records, never implicit topological cells."""

    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    status: Status = Status.NOT_CALCULATED

    def __post_init__(self):
        quantities = _quantity_mapping(self.quantities)
        for name, quantity in quantities.items():
            if name not in MOLECULAR_GEOMETRY_QUANTITIES:
                raise ValueError(f"unsupported molecular geometry quantity: {name}")
            if quantity.name != name or quantity.units != "1":
                raise ValueError(f"invalid schema for molecular geometry quantity: {name}")
            if quantity.value is not None and (
                not isinstance(quantity.value, int) or isinstance(quantity.value, bool) or quantity.value < 0
            ):
                raise ValueError(f"molecular geometry counts must be non-negative integers: {name}")
        object.__setattr__(self, "quantities", quantities)


@dataclass(frozen=True)
class TopologicalCellComplex:
    """A certified or explicitly incomplete cell complex for topology only."""

    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    status: Status = Status.NOT_CALCULATED
    provenance: Provenance | None = None

    def __post_init__(self):
        quantities = _quantity_mapping(self.quantities)
        boolean_names = {
            "orientable",
            "chain_complex_valid",
            "boundary_of_boundary_zero",
        }
        count_names = {
            "vertex_count",
            "edge_count",
            "topological_edge_count",
            "topological_cut_count",
            "patch_attachment_record_count",
            "face_count",
            "connected_component_count",
            "boundary_component_count",
            "manifold_vertex_count",
            "singular_vertex_count",
            "genus",
            "nonorientable_genus",
            "beta0",
            "beta1",
            "beta2",
        }
        integer_names = count_names | {
            "euler_characteristic",
            "genus",
            "nonorientable_genus",
        }
        for name, quantity in quantities.items():
            if name not in TOPOLOGY_CELL_QUANTITIES:
                raise ValueError(f"unsupported topological cell quantity: {name}")
            if quantity.name != name or quantity.units != "1":
                raise ValueError(f"invalid schema for topological cell quantity: {name}")
            if quantity.value is not None and name in boolean_names and not isinstance(quantity.value, bool):
                raise ValueError(f"topological boolean must be true, false, or null: {name}")
            if quantity.value is not None and name in count_names and (
                not isinstance(quantity.value, int) or isinstance(quantity.value, bool) or quantity.value < 0
            ):
                raise ValueError(f"topological counts must be non-negative integers: {name}")
            if quantity.value is not None and name in integer_names and (
                not isinstance(quantity.value, int) or isinstance(quantity.value, bool)
            ):
                raise ValueError(f"topological integer must be null or an integer: {name}")
        object.__setattr__(self, "quantities", quantities)


@dataclass(frozen=True)
class MolecularPatchProvenance:
    patch_id: str
    source_atom_stable_id: str | None = None
    source_sphere_center: tuple[float, float, float] | None = None
    source_radius: float | None = None
    source_radius_units: str | None = None
    partner_contact_atom_ids: tuple[str, ...] = ()
    parent_selected_contact_ids: tuple[str, ...] = ()
    clipping_occlusion_atom_ids: tuple[str, ...] = ()
    component_id: str | int | None = None
    boundary_arc_ids: tuple[str, ...] = ()
    junction_ids: tuple[str, ...] = ()
    status: Status = Status.NOT_CALCULATED

    def __post_init__(self):
        if not self.patch_id:
            raise ValueError("Molecular patches require patch_id")
        object.__setattr__(self, "status", Status(self.status))
        if self.source_sphere_center is not None:
            center = tuple(float(value) for value in self.source_sphere_center)
            if len(center) != 3 or any(not math.isfinite(value) for value in center):
                raise ValueError("Source sphere center must contain three finite values")
            object.__setattr__(self, "source_sphere_center", center)
        if self.source_radius is not None:
            if not math.isfinite(self.source_radius) or self.source_radius < 0:
                raise ValueError("Source sphere radius must be finite and nonnegative")
            if not self.source_radius_units:
                raise ValueError("Source sphere radius requires units")
            object.__setattr__(self, "source_radius", float(self.source_radius))
        for name in (
            "partner_contact_atom_ids",
            "parent_selected_contact_ids",
            "clipping_occlusion_atom_ids",
            "boundary_arc_ids",
            "junction_ids",
        ):
            object.__setattr__(self, name, _string_tuple(getattr(self, name), name))


_INTERACTION_COMPARISON_UNITS = {
    name: ("\u00c5\u00b2" if name.startswith("area_") and "ratio" not in name
           else "\u00c5" if name.startswith("integrated_mean_") else "1")
    for name in INTERACTION_COMPARISON_QUANTITIES
}


@dataclass(frozen=True)
class InteractionLinks:
    interaction_id: str
    result_ids: Mapping[str, str] = field(default_factory=dict)
    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    status: Status = Status.NOT_CALCULATED
    contact_selection_source: str | None = None
    contact_selection_id: str | None = None
    orientation_conventions: Mapping[str, str] = field(default_factory=dict)
    comparison_quantities: Mapping[str, Quantity] = field(default_factory=dict)

    def __post_init__(self):
        if not self.interaction_id:
            raise ValueError("Interaction links require interaction_id")
        object.__setattr__(self, "status", Status(self.status))
        result_ids = dict(self.result_ids)
        if any(
            not isinstance(key, str)
            or not key
            or not isinstance(value, str)
            or not value
            for key, value in result_ids.items()
        ):
            raise ValueError("Interaction result IDs require nonempty strings")
        object.__setattr__(self, "result_ids", freeze(result_ids))
        object.__setattr__(self, "quantities", freeze(self.quantities))
        _quantity_mapping(self.quantities)


def generator_complex_extension(
    complex_, *, namespace="vorpy.dual_generator_complex", schema_version="1"
):
    """Wrap a generator complex in the existing extension channel."""
    if not isinstance(complex_, GeneratorComplex):
        raise TypeError("Expected GeneratorComplex")
    return RepresentationExtension(
        namespace,
        schema_version,
        complex_.status,
        {
            "kind": "generator_complex",
            "side": complex_.side,
            "simplex_counts": {
                str(dimension): dictionary(quantity)
                for dimension, quantity in complex_.simplex_counts.items()
            },
            "topology": dictionary(complex_.topology),
        },
        complex_.provenance,
    )


@dataclass(frozen=True)
class AtomIdentity:
    system_id: str
    chain: str
    residue_number: str
    insertion_code: str
    residue_name: str
    atom_name: str
    alternate_location: str = ""
    system_atom_id: str | None = None
    established_stable_id: str | None = None
    generator_index: int | None = None
    network_id: str | None = None

    def __post_init__(self):
        for name in (
            "system_id",
            "chain",
            "residue_number",
            "insertion_code",
            "residue_name",
            "atom_name",
            "alternate_location",
        ):
            if not isinstance(getattr(self, name), str):
                raise TypeError(f"Atom {name} must be a string")
        if self.system_atom_id is not None and (
            not isinstance(self.system_atom_id, str) or not self.system_atom_id
        ):
            raise ValueError("System atom ID must be a nonempty string")
        if not self.system_id:
            raise ValueError(
                "Atom identity requires a stable molecular system namespace"
            )
        if not self.system_atom_id and not all(
            (self.residue_number, self.residue_name, self.atom_name)
        ):
            raise ValueError(
                "Atom requires stable system atom ID or molecular identity; generator index is insufficient"
            )
        if self.generator_index is not None and not self.network_id:
            raise ValueError("Network-local generator index requires network ID")
        if self.generator_index is not None and (
            isinstance(self.generator_index, bool)
            or not isinstance(self.generator_index, int)
        ):
            raise ValueError("Generator index must be an integer")

    @property
    def stable_atom_id(self):
        key = (
            {"system_id": self.system_id, "system_atom_id": self.system_atom_id}
            if self.system_atom_id
            else {
                "system_id": self.system_id,
                "chain": self.chain,
                "residue_number": self.residue_number,
                "insertion_code": self.insertion_code,
                "residue_name": self.residue_name,
                "atom_name": self.atom_name,
                "alternate_location": self.alternate_location,
            }
        )
        return stable_id("atom-v1", key)

    def to_dict(self):
        return {"stable_atom_id": self.stable_atom_id, **dictionary(self)}


@dataclass(frozen=True)
class ContactRecord:
    result_id: str
    identity: ResultIdentity
    atom_A: AtomIdentity
    atom_B: AtomIdentity
    selection_status: str
    dual_feature_id: str | None = None
    physical_feature_id: str | int | None = None
    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    provenance: Provenance | None = None

    def __post_init__(self):
        if (
            self.identity.result_kind != "interface"
            or self.result_id != self.identity.result_id
        ):
            raise ValueError("Contacts require their owning interface result identity")
        if self.selection_status not in {
            "SELECTED",
            "REJECTED",
            "UNRESOLVED",
            "NOT_CALCULATED",
            "NOT_SUPPORTED",
            "STALE",
        }:
            raise ValueError("Unknown contact selection status")
        if (
            self.atom_A.system_id != (self.identity.system_id or self.identity.system)
            or self.atom_B.system_id != self.atom_A.system_id
        ):
            raise ValueError("Contact atoms belong to a different system namespace")
        object.__setattr__(self, "quantities", freeze(self.quantities))
        _quantity_mapping(self.quantities)

    @property
    def molecular_pair_id(self):
        """Frame/representation-independent ordered chemical pair identifier."""
        return stable_id(
            "pair-v1",
            {
                "atom_A": self.atom_A.stable_atom_id,
                "atom_B": self.atom_B.stable_atom_id,
            },
        )

    @property
    def contact_id(self):
        return stable_id(
            "contact-v1",
            {
                "result_id": self.result_id,
                "atom_A": self.atom_A.stable_atom_id,
                "atom_B": self.atom_B.stable_atom_id,
                "dual_feature_id": self.dual_feature_id,
                "physical_feature_id": self.physical_feature_id,
            },
        )

    def to_dict(self):
        return {
            "contact_id": self.contact_id,
            "molecular_pair_id": self.molecular_pair_id,
            **dictionary(self),
            "atom_A": self.atom_A.to_dict(),
            "atom_B": self.atom_B.to_dict(),
        }


@dataclass(frozen=True)
class RepresentationExtension:
    namespace: str
    schema_version: str
    status: Status
    payload: Mapping = field(default_factory=dict)
    provenance: Provenance | None = None

    def __post_init__(self):
        if not self.namespace or not self.schema_version:
            raise ValueError("Extensions require a namespace and version")
        object.__setattr__(self, "status", Status(self.status))
        object.__setattr__(self, "payload", freeze(self.payload))


def _quantity_mapping(mapping):
    if any(not isinstance(q, Quantity) or key != q.name for key, q in mapping.items()):
        raise ValueError("Quantity mappings must be keyed by the quantity name")
    return mapping


def _quantities(value):
    if isinstance(value, Quantity):
        yield value
    elif is_dataclass(value):
        for f in fields(value):
            yield from _quantities(getattr(value, f.name))
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _quantities(item)
    elif isinstance(value, tuple):
        for item in value:
            yield from _quantities(item)


@dataclass(frozen=True)
class ResultModel:
    identity: ResultIdentity
    provenance: Provenance
    status: Status = Status.NOT_CALCULATED
    environment: Environment | None = None
    metrics: Mapping[str, Quantity] = field(default_factory=dict)
    topology: Topology | None = None
    mean_curvature: MeanCurvature | None = None
    gaussian_curvature: GaussianCurvature | None = None
    contacts: tuple[ContactRecord, ...] = ()
    selection: Mapping[str, Quantity] = field(default_factory=dict)
    extensions: tuple[RepresentationExtension, ...] = ()
    schema_version: str = SCHEMA_VERSION
    KIND: ClassVar[str | None] = None

    def __post_init__(self):
        if self.KIND and self.identity.result_kind != self.KIND:
            raise ValueError("Result type and identity kind disagree")
        object.__setattr__(self, "status", Status(self.status))
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError("Unsupported result schema version")
        if not isinstance(self.identity, ResultIdentity) or not isinstance(
            self.provenance, Provenance
        ):
            raise TypeError("Results require typed identity and provenance")
        topology_type = (
            (Topology, TopologicalCellComplex)
            if self.KIND == "molecular_contact_surface"
            else Topology
        )
        for name, cls in (
            ("environment", Environment),
            ("topology", topology_type),
            ("mean_curvature", MeanCurvature),
            ("gaussian_curvature", GaussianCurvature),
        ):
            if getattr(self, name) is not None and not isinstance(
                getattr(self, name), cls
            ):
                raise TypeError(f"{name} requires {cls.__name__}")
        for name in ("metrics", "selection"):
            object.__setattr__(self, name, freeze(getattr(self, name)))
            _quantity_mapping(getattr(self, name))
        object.__setattr__(self, "contacts", tuple(self.contacts))
        object.__setattr__(self, "extensions", tuple(self.extensions))
        if any(not isinstance(c, ContactRecord) for c in self.contacts) or any(
            not isinstance(e, RepresentationExtension) for e in self.extensions
        ):
            raise TypeError("Contacts and extensions require typed records")
        if (
            self.environment
            and self.environment.environment != self.identity.environment
        ):
            raise ValueError("Result identity and environment disagree")
        if any(c.result_id != self.result_id for c in self.contacts):
            raise ValueError("Contact belongs to a different result")
        if len({c.contact_id for c in self.contacts}) != len(self.contacts):
            raise ValueError("Duplicate contact IDs")
        if len({e.namespace for e in self.extensions}) != len(self.extensions):
            raise ValueError("Duplicate extension namespace")
        if self.status in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED} and any(
            q.value is not None
            for name in (
                "metrics",
                "selection",
                "topology",
                "mean_curvature",
                "gaussian_curvature",
                "contacts",
            )
            for q in _quantities(getattr(self, name))
        ):
            raise ValueError("Unavailable result cannot contain calculated quantities")
        if self.status == Status.STALE:
            for name in (
                "metrics",
                "selection",
                "topology",
                "mean_curvature",
                "gaussian_curvature",
                "contacts",
            ):
                object.__setattr__(self, name, freeze(_stale(getattr(self, name))))
            object.__setattr__(
                self,
                "extensions",
                tuple(replace(e, status=Status.STALE) for e in self.extensions),
            )

    @property
    def result_id(self):
        return self.identity.result_id

    def to_dict(self):
        return {
            "result_id": self.result_id,
            **dictionary(self),
            "contacts": [c.to_dict() for c in self.contacts],
        }


@dataclass(frozen=True)
class SystemResult(ResultModel):
    KIND = "system"


@dataclass(frozen=True)
class NetworkResult(ResultModel):
    KIND = "network"


@dataclass(frozen=True)
class GroupResult(ResultModel):
    KIND = "group"


@dataclass(frozen=True)
class InterfaceResult(ResultModel):
    KIND = "interface"


@dataclass(frozen=True)
class DualSimplexRecord:
    """An incidence record from a dual complex, with no geometric measure.

    ``generator_ids`` are network-local generator indices, as declared by the
    source dual cache.  ``primal_feature_ids`` are stable serialized feature
    references (for example ``"surf:17"``), not physical measurements.
    """

    simplex_id: str
    dimension: int
    generator_ids: tuple[int, ...]
    primal_feature_ids: tuple[str, ...] = ()
    bounded: bool | None = None
    complete: bool | None = None
    supported: bool | None = None
    status: Status = Status.NOT_CALCULATED
    provenance: Provenance | None = None

    def __post_init__(self):
        if isinstance(self.dimension, bool) or not isinstance(self.dimension, int) or self.dimension < 0:
            raise ValueError("Dual-simplex dimensions must be nonnegative integers")
        generator_ids = tuple(self.generator_ids)
        if len(generator_ids) != self.dimension + 1:
            raise ValueError("Dual-simplex dimension and generator count disagree")
        if (
            any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in generator_ids)
            or tuple(sorted(generator_ids)) != generator_ids
            or len(set(generator_ids)) != len(generator_ids)
        ):
            raise ValueError("Dual-simplex generator IDs must be sorted unique nonnegative integers")
        expected_id = f"{self.dimension}:" + ",".join(map(str, generator_ids))
        if self.simplex_id != expected_id:
            raise ValueError("Dual-simplex ID must match its dimension and generator IDs")
        primal_feature_ids = tuple(self.primal_feature_ids)
        if any(not isinstance(value, str) or not value for value in primal_feature_ids):
            raise ValueError("Dual primal-feature IDs must be non-empty strings")
        for name in ("bounded", "complete", "supported"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, bool):
                raise TypeError(f"Dual-simplex {name} must be a boolean or null")
        if self.provenance is not None and not isinstance(self.provenance, Provenance):
            raise TypeError("Dual-simplex provenance must be typed")
        object.__setattr__(self, "generator_ids", generator_ids)
        object.__setattr__(self, "primal_feature_ids", primal_feature_ids)
        object.__setattr__(self, "status", Status(self.status))


def _dual_generator_pairs(values, name):
    pairs = tuple(tuple(pair) for pair in values)
    if any(
        len(pair) != 2
        or any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in pair)
        or pair != tuple(sorted(pair))
        or pair[0] == pair[1]
        for pair in pairs
    ):
        raise ValueError(f"{name} must contain sorted pairs of distinct nonnegative generator IDs")
    if len(set(pairs)) != len(pairs):
        raise ValueError(f"{name} contains duplicate generator pairs")
    return pairs


def _selection_pairs(values, name):
    pairs = tuple(tuple(pair) for pair in values)
    if any(len(pair) != 2 or any(not isinstance(value, str) or not value for value in pair) for pair in pairs):
        raise ValueError(f"{name} must contain pairs of non-empty stable identifiers")
    return pairs


@dataclass(frozen=True)
class DualContactComplexResult(ResultModel):
    """Dual contact incidences; deliberately not a physical surface result."""

    KIND = "interface"

    contact_features: tuple[DualSimplexRecord, ...] = ()

    def __post_init__(self):
        super().__post_init__()
        if self.identity.representation != "dual":
            raise ValueError("Dual-contact results require the dual representation")
        object.__setattr__(self, "contact_features", tuple(self.contact_features))
        if any(
            not isinstance(feature, DualSimplexRecord) or feature.dimension != 1
            for feature in self.contact_features
        ):
            raise TypeError("Dual-contact results require typed dimension-one features")
        if len({feature.simplex_id for feature in self.contact_features}) != len(self.contact_features):
            raise ValueError("Dual-contact results contain duplicate features")
        if self.status in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED} and self.contact_features:
            raise ValueError("Unavailable dual-contact results cannot contain incidences")
        if self.status == Status.STALE:
            object.__setattr__(
                self,
                "contact_features",
                tuple(_stale(feature) for feature in self.contact_features),
            )


@dataclass(frozen=True)
class DualGeneratorComplexResult(ResultModel):
    """Traversed generator complex with explicit simplices and side identity."""

    KIND = "interface"

    generator_complex: GeneratorComplex | None = None
    simplex_records: tuple[DualSimplexRecord, ...] = ()

    def __post_init__(self):
        super().__post_init__()
        if self.identity.representation not in {"dual_a", "dual_b", "dual_generator_complex"}:
            raise ValueError("Dual-generator results require a generator-complex representation")
        if self.generator_complex is not None and not isinstance(self.generator_complex, GeneratorComplex):
            raise TypeError("generator_complex must be GeneratorComplex")
        if self.generator_complex is not None and self.generator_complex.side != self.identity.side:
            raise ValueError("Generator-complex side disagrees with result identity")
        object.__setattr__(self, "simplex_records", tuple(self.simplex_records))
        if any(not isinstance(record, DualSimplexRecord) for record in self.simplex_records):
            raise TypeError("Dual-generator results require typed simplex records")
        if len({record.simplex_id for record in self.simplex_records}) != len(self.simplex_records):
            raise ValueError("Dual-generator results contain duplicate simplex records")
        if self.status in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED} and (
            self.generator_complex is not None or self.simplex_records
        ):
            raise ValueError("Unavailable dual-generator results cannot contain incidences")
        if self.status == Status.STALE:
            if self.generator_complex is not None:
                object.__setattr__(self, "generator_complex", _stale(self.generator_complex))
            object.__setattr__(
                self,
                "simplex_records",
                tuple(_stale(record) for record in self.simplex_records),
            )


@dataclass(frozen=True)
class AlphaSelectionResult(ResultModel):
    """A standalone alpha-selected dual subset, separate from primal geometry."""

    KIND = "interface"

    selected_dual_feature_ids: tuple[str, ...] = ()
    selected_generator_pairs: tuple[tuple[int, int], ...] = ()
    selected_system_pairs: tuple[tuple[str, str], ...] = ()
    selected_physical_system_pairs: tuple[tuple[str, str], ...] = ()
    selected_dual_features: tuple[DualSimplexRecord, ...] = ()

    def __post_init__(self):
        super().__post_init__()
        if self.identity.representation != "alpha_selection":
            raise ValueError("Alpha-selection results require the alpha_selection representation")
        if self.identity.alpha_value is None:
            raise ValueError("Alpha-selection results require an explicit alpha identity")
        object.__setattr__(
            self,
            "selected_dual_feature_ids",
            _string_tuple(self.selected_dual_feature_ids, "selected_dual_feature_ids"),
        )
        if len(set(self.selected_dual_feature_ids)) != len(self.selected_dual_feature_ids):
            raise ValueError("Alpha-selection results contain duplicate dual feature IDs")
        object.__setattr__(
            self,
            "selected_generator_pairs",
            _dual_generator_pairs(self.selected_generator_pairs, "selected_generator_pairs"),
        )
        object.__setattr__(
            self,
            "selected_system_pairs",
            _selection_pairs(self.selected_system_pairs, "selected_system_pairs"),
        )
        object.__setattr__(
            self,
            "selected_physical_system_pairs",
            _selection_pairs(self.selected_physical_system_pairs, "selected_physical_system_pairs"),
        )
        object.__setattr__(self, "selected_dual_features", tuple(self.selected_dual_features))
        if any(
            not isinstance(feature, DualSimplexRecord) or feature.dimension != 1
            for feature in self.selected_dual_features
        ):
            raise TypeError("Alpha-selection results require typed dimension-one dual features")
        feature_ids = {feature.simplex_id for feature in self.selected_dual_features}
        if not feature_ids.issubset(set(self.selected_dual_feature_ids)):
            raise ValueError("Selected dual feature records must be named by the alpha selection")
        if self.status in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED} and any(
            (
                self.selected_dual_feature_ids,
                self.selected_generator_pairs,
                self.selected_system_pairs,
                self.selected_physical_system_pairs,
                self.selected_dual_features,
            )
        ):
            raise ValueError("Unavailable alpha selections cannot contain selected identities")
        if self.status == Status.STALE:
            object.__setattr__(
                self,
                "selected_dual_features",
                tuple(_stale(feature) for feature in self.selected_dual_features),
            )


@dataclass(frozen=True)
class InteractionLinks:
    """Links the three interaction representations without owning their state."""

    interaction_id: str
    result_ids: Mapping[str, str] = field(default_factory=dict)
    quantities: Mapping[str, Quantity] = field(default_factory=dict)
    status: Status = Status.NOT_CALCULATED
    contact_selection_source: str | None = None
    contact_selection_id: str | None = None
    orientation_conventions: Mapping[str, str] = field(default_factory=dict)
    comparison_quantities: Mapping[str, Quantity] = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.interaction_id, str) or not self.interaction_id:
            raise ValueError("interaction_id must be a non-empty string")
        if any(
            not isinstance(role, str) or not role or not isinstance(result_id, str) or not result_id
            for role, result_id in self.result_ids.items()
        ):
            raise ValueError("interaction result IDs must be non-empty strings")
        object.__setattr__(self, "result_ids", MappingProxyType(dict(self.result_ids)))
        object.__setattr__(self, "quantities", _quantity_mapping(self.quantities))
        if (self.contact_selection_source is None) != (self.contact_selection_id is None):
            raise ValueError("contact selection source and ID must be supplied together")
        if self.contact_selection_source is not None:
            if not self.contact_selection_source or not self.contact_selection_id:
                raise ValueError("contact selection provenance must be non-empty")

        normalized_orientations = {}
        for role, value in self.orientation_conventions.items():
            if not isinstance(role, str) or not role:
                raise ValueError("orientation roles must be non-empty strings")
            normalized_orientations[role] = OrientationConvention(value).value
        expected_orientations = {
            "molecular_contact_surface:A": OrientationConvention.MOLECULAR_OUTWARD_NORMAL.value,
            "molecular_contact_surface:B": OrientationConvention.MOLECULAR_OUTWARD_NORMAL.value,
            "physical_partition_interface": OrientationConvention.PHYSICAL_INTERFACE_A_TO_B.value,
        }
        for role, expected in expected_orientations.items():
            if role in normalized_orientations and normalized_orientations[role] != expected:
                raise ValueError(f"invalid orientation convention for {role}")
        if any(
            name.startswith("integrated_mean_curvature_")
            or name.startswith("integrated_gaussian_curvature_")
            for name in self.comparison_quantities
        ) and any(role not in normalized_orientations for role in expected_orientations):
            raise ValueError("signed curvature comparisons require all interaction orientations")
        object.__setattr__(self, "orientation_conventions", MappingProxyType(normalized_orientations))

        comparison = _quantity_mapping(self.comparison_quantities)
        for name, quantity in comparison.items():
            if name not in INTERACTION_COMPARISON_QUANTITIES:
                raise ValueError(f"unsupported interaction comparison quantity: {name}")
            if quantity.name != name or quantity.units != _INTERACTION_COMPARISON_UNITS[name]:
                raise ValueError(f"invalid schema for interaction comparison quantity: {name}")
        if self.status in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED} and any(
            quantity.value is not None for quantity in comparison.values()
        ):
            raise ValueError("unavailable interaction links cannot contain calculated comparisons")
        object.__setattr__(self, "comparison_quantities", comparison)


@dataclass(frozen=True)
class MolecularContactSurfaceResult(ResultModel):
    KIND = "molecular_contact_surface"

    molecular_selection: MolecularContactSelection | None = None
    geometry: MolecularContactSurfaceGeometry | None = None
    geometric_counts: MolecularGeometryCounts | None = None
    topology: TopologicalCellComplex | None = None
    patch_provenance: tuple[MolecularPatchProvenance, ...] = ()
    singular_junction_contribution: Quantity | None = None
    interaction_links: InteractionLinks | None = None

    def __post_init__(self):
        super().__post_init__()
        if self.geometric_counts is not None and not isinstance(self.geometric_counts, MolecularGeometryCounts):
            raise TypeError("geometric_counts must be MolecularGeometryCounts")
        if self.topology is not None and not isinstance(self.topology, TopologicalCellComplex):
            raise TypeError("topology must be TopologicalCellComplex")
        if self.status == Status.STALE:
            if self.geometric_counts is not None:
                object.__setattr__(self, "geometric_counts", _stale(self.geometric_counts))
            if self.topology is not None:
                object.__setattr__(self, "topology", _stale(self.topology))
        if self.molecular_selection is not None and not isinstance(
            self.molecular_selection, MolecularContactSelection
        ):
            raise TypeError("molecular_selection requires MolecularContactSelection")
        if self.geometry is not None and not isinstance(
            self.geometry, MolecularContactSurfaceGeometry
        ):
            raise TypeError("geometry requires MolecularContactSurfaceGeometry")
        object.__setattr__(self, "patch_provenance", tuple(self.patch_provenance))
        if any(not isinstance(patch, MolecularPatchProvenance) for patch in self.patch_provenance):
            raise TypeError("patch_provenance requires typed patch records")
        if self.singular_junction_contribution is not None and (
            not isinstance(self.singular_junction_contribution, Quantity)
            or self.singular_junction_contribution.name
            != "singular_junction_contribution"
        ):
            raise ValueError("Invalid singular junction contribution")
        if self.interaction_links is not None and not isinstance(
            self.interaction_links, InteractionLinks
        ):
            raise TypeError("interaction_links requires InteractionLinks")
        if self.interaction_links is not None and (
            self.interaction_links.interaction_id != self.identity.interaction_id
        ):
            raise ValueError("Interaction link identity disagrees with result identity")
        quantities = []
        if self.molecular_selection is not None:
            quantities.extend(self.molecular_selection.quantities.values())
        if self.geometry is not None:
            quantities.extend(self.geometry.quantities.values())
        if self.singular_junction_contribution is not None:
            quantities.append(self.singular_junction_contribution)
        if self.status in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED} and any(
            quantity.value is not None for quantity in quantities
        ):
            raise ValueError("Unavailable molecular surface cannot contain values")
        if self.status == Status.STALE:
            object.__setattr__(self, "molecular_selection", _stale(self.molecular_selection))
            object.__setattr__(self, "geometry", _stale(self.geometry))
            object.__setattr__(
                self,
                "patch_provenance",
                tuple(_stale(patch) for patch in self.patch_provenance),
            )
            object.__setattr__(
                self,
                "singular_junction_contribution",
                _stale(self.singular_junction_contribution),
            )
            object.__setattr__(self, "interaction_links", _stale(self.interaction_links))


# Future table contract only. No CSV writer or scientific derivations live here.
METRICS_COLUMNS = (
    "result_id",
    "system",
    "frame",
    "network_id",
    "network_scope",
    "group_id",
    "interface_id",
    "group_A",
    "group_B",
    "environment",
    "partition",
    "representation",
    "side",
    "orientation",
    "orientation_convention",
    "interaction_id",
    "contact_selection_source",
    "contact_selection_id",
    "molecular_surface_convention",
    "definition_version",
    "partner_id",
    "source_network_id",
    "source_interface_id",
    "alpha_method",
    "alpha_value",
    "alpha_units",
    "probe_radius",
    "probe_units",
    "radius_configuration_id",
    "n_atoms_A",
    "n_atoms_B",
    "n_interface_atoms_A",
    "n_interface_atoms_B",
    "n_interface_residues_A",
    "n_interface_residues_B",
    "n_contacts",
    "n_residue_contacts",
    "area",
    "area_status",
    "area_scope",
    "integrated_mean_curvature",
    "mean_curvature_per_area",
    "mean_curvature_status",
    "mean_curvature_scope",
    "complete_mean_curvature",
    "complete_mean_curvature_status",
    "integrated_gaussian_curvature",
    "gaussian_curvature_per_area",
    "gaussian_curvature_status",
    "gaussian_curvature_scope",
    "gauss_bonnet_completion",
    "gauss_bonnet_status",
    "vertices",
    "edges",
    "faces",
    "euler_characteristic",
    "components",
    "boundary_edges",
    "boundary_loops",
    "nonmanifold_edges",
    "topology_status",
    "selected_features",
    "rejected_features",
    "unresolved_features",
    "selection_status",
    "explicit_solvent_present",
    "solvent_generator_count",
    "solvent_interface_competitors",
)

HUMAN_SUMMARY_LEVELS = MappingProxyType(
    {
        "headline": (
            "system",
            "group_id",
            "interface_id",
            "group_A",
            "group_B",
            "environment",
            "partition",
            "representation",
            "area",
            "integrated_mean_curvature",
            "integrated_gaussian_curvature",
            "euler_characteristic",
            "n_contacts",
        ),
        "support": (
            "topology",
            "mean_curvature",
            "gaussian_curvature",
            "selection",
            "quantity_status_and_scope",
        ),
        "provenance": (
            "network_id",
            "network_scope",
            "frame",
            "alpha_method",
            "alpha_value",
            "alpha_units",
            "radius_configuration_id",
            "probe_radius",
            "orientation",
            "provenance",
            "diagnostics",
        ),
    }
)
