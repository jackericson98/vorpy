"""Small opt-in group-analysis prototype for a solved :class:`NetworkView`.

This module deliberately does not materialize a ``Network`` and is not wired
into ``Group.build``.  It proves the ownership and area/volume part of shared
group analysis while the legacy curvature pipeline remains Network-shaped.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral, Real

from vorpy.src.calculations import calc_surf_sa
from vorpy.src.network.views import NetworkViewError


@dataclass(frozen=True)
class GroupViewAnalysis:
    """Read-only group metrics and identity evidence from one view.

    ``area`` includes every exposed selected-cell surface.  The two optional
    area components make the shell-boundary scope explicit: boundary area is
    not a physical molecular interface measurement.
    """

    source_network_id: str
    solve_universe_id: str
    selection_id: str
    system_ids: tuple[int, ...]
    boundary_surface_ids: tuple[object, ...]
    competitor_system_ids: tuple[int, ...]
    area: float
    volume: float
    complete: bool
    boundary_dependent_area: float = 0.0
    physical_interface_area: float = 0.0


def system_ids_for_group(group, source_network) -> tuple[int, ...]:
    """Validate a Group's System-level IDs against one solved source network.

    ``Group.ball_ndxs`` is a System selection namespace.  It must not be used
    as a compact ``Network.balls`` row position when solvent rows are absent
    from the solve universe.

    shortcut: this accepts only already-global ``ball_ndxs``; upgrade when
    FORGE publishes a typed producer selection map for ambiguous legacy groups.
    """
    values = getattr(group, "ball_ndxs", None)
    if values is None:
        raise NetworkViewError("Group selection has no System-level ball IDs.")
    try:
        values = tuple(values)
    except TypeError as error:
        raise NetworkViewError("Group selection must contain System-level IDs.") from error
    if not values or any(not isinstance(value, Integral) or isinstance(value, bool) for value in values):
        raise NetworkViewError("Group selection must contain integer System-level IDs.")
    selected = tuple(sorted(set(int(value) for value in values)))

    try:
        context = source_network.get_solve_context_key()
    except (AttributeError, TypeError, ValueError) as error:
        raise NetworkViewError(f"Unsupported source solve context: {error}") from error
    available = {
        system_id
        for _, system_id in context.local_to_system
        if system_id is not None
    }
    missing = set(selected).difference(available)
    if missing:
        raise NetworkViewError(
            "Group selection is outside the solved generator universe: "
            + ", ".join(str(value) for value in sorted(missing))
        )
    return selected


def make_group_view(group, source_network, *, source_network_id=None):
    """Create one read-only view from a Group without invoking a build."""
    if source_network_id is not None:
        if not isinstance(source_network_id, str) or not source_network_id.strip():
            raise NetworkViewError("source_network_id must be a non-empty string.")
        existing_source_id = getattr(source_network, "source_network_id", None)
        if existing_source_id and existing_source_id != source_network_id:
            raise NetworkViewError(
                "Explicit source_network_id conflicts with the source Network."
            )
    system_ids = system_ids_for_group(group, source_network)
    view = source_network.make_view(
        group_ids=system_ids,
        source_network_id=source_network_id,
    )
    if source_network_id is not None:
        source_network.source_network_id = source_network_id
    return view


def shared_network_identity(view, source_network):
    """Project one F10/F12-backed view into the Results shared-network identity."""
    try:
        context = source_network.get_solve_context_key()
    except (AttributeError, TypeError, ValueError) as error:
        raise NetworkViewError(f"Unsupported source solve context: {error}") from error
    if view.solve_universe_id != context.solve_universe_id:
        raise NetworkViewError("View and source network solve identities disagree.")

    provenance = getattr(source_network, "solvent_provenance", None)
    if not isinstance(provenance, dict) or not provenance.get("solvent_mode"):
        raise NetworkViewError("Results identity requires authoritative F10 solvent provenance.")
    retention = getattr(source_network, "retained_generator_provenance", None)
    if not isinstance(retention, dict):
        raise NetworkViewError("Results identity requires authoritative F12 retention provenance.")

    from vorpy.src.results.shared_network import SharedNetworkIdentity

    try:
        return SharedNetworkIdentity(
            source_network_id=view.source_network_id,
            source_identity_status=(
                "AUTHORITATIVE"
                if view.source_network_id_is_authoritative
                else "PROVISIONAL"
            ),
            solve_universe_id=view.solve_universe_id,
            frame_id=view.frame_id,
            analysis_context_kind="GLOBAL_NETWORK",
            geometry_scheme=view.scheme,
            partition_scheme=context.partition_scheme,
            scope_mode=context.scope_mode,
            solvent_mode=provenance["solvent_mode"],
            solvent_generator_digest=provenance.get("solvent_generator_digest"),
            solvent_cells_complete=provenance.get("solvent_cells_complete"),
            retained_system_indices=retention.get("retained_system_indices", ()),
            retained_generator_ids=retention.get("retained_generator_ids", ()),
            retained_generator_digest=retention.get("retained_generator_digest"),
            selection_id=view.selection_id,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise NetworkViewError(f"Invalid F10 solvent provenance: {error}") from error


def attach_group_view_results(
    analysis_result,
    source_network,
    groups,
    *,
    system: str,
    radius_configuration_id: str,
    source_revision: str | None = None,
    source_network_id: str | None = None,
    cached_dependencies: Mapping[str, object] | None = None,
):
    """Attach typed Results for several selections of one solved network.

    This is an explicit integration seam for PRISM/ATLAS consumers. It only
    creates immutable views and never calls ``Network.build``. Producer source
    IDs are preferred. Without a producer source ID, NetworkView may retain a
    solve-scoped provisional bookkeeping ID, but Results never promotes it (or
    an archive ID) to authoritative source identity for the analysis.
    """
    from vorpy.src.results.adapters import adapt_group_view_analysis
    from vorpy.src.results.shared_network import ResultsDependency
    from vorpy.src.network.aw_geometry_reuse import build_geometry_lineage_metadata
    from vorpy.src.results import make_geometry_lineage

    context = source_network.get_solve_context_key()
    provenance = getattr(source_network, "solvent_provenance", None)
    if not isinstance(provenance, dict) or not provenance.get("solvent_mode"):
        raise NetworkViewError(
            "Typed group Results require authoritative F10 solvent provenance."
        )
    retention = getattr(source_network, "retained_generator_provenance", None)
    if not isinstance(retention, dict):
        raise NetworkViewError(
            "Typed group Results require authoritative F12 retention provenance."
        )
    source_digest = getattr(context, "digest", None) or context.solve_universe_id
    settings = getattr(source_network, "settings", {})
    result_partition = settings.get("net_type") if isinstance(settings, dict) else None
    typed_results = []
    for group in groups:
        view = make_group_view(
            group, source_network, source_network_id=source_network_id
        )
        geometry_lineage = make_geometry_lineage(
            build_geometry_lineage_metadata(
                source_network,
                orientation={"convention": "source_view"},
            ),
            view_evidence=view,
        ).as_dict()
        analysis = analyze_group_view(view)
        identity = shared_network_identity(view, source_network)
        dependency = ResultsDependency(
            identity=identity,
            source_id=f"{view.source_network_id}:{view.selection_id}",
            source_revision=source_revision,
            source_digest=source_digest,
        )
        group_id = (
            getattr(group, "group_id", None)
            or getattr(group, "name", None)
            or view.selection_id
        )
        typed_results.append(
            adapt_group_view_analysis(
                analysis,
                identity,
                dependency,
                group_id=str(group_id),
                system=system,
                radius_configuration_id=radius_configuration_id,
                results_partition=result_partition,
                cached_dependency=(cached_dependencies or {}).get(view.selection_id),
                source_solve_context=context,
                source_solvent_provenance=provenance,
                source_retained_generator_provenance=retention,
                geometry_lineage=geometry_lineage,
            )
        )
    analysis_result.typed_results = tuple(typed_results)
    return analysis_result.typed_results


def analyze_group_view(view) -> GroupViewAnalysis:
    """Analyze group volume and exposed area without rebuilding or copying a network.

    The view owns the selection and preserves full-source incidence.  A surface
    is exposed when it touches a selected generator and either an unselected
    physical generator or a boundary generator.
    """
    if getattr(view, "selection_kind", None) != "group":
        raise NetworkViewError("Group analysis requires a group NetworkView.")
    if not getattr(view.completeness, "proven", False):
        raise NetworkViewError("Group analysis requires proven view completeness.")

    selected_ids = set(view.system_ids)
    volume = 0.0
    for cell in view.cells:
        if not bool(cell.data.get("complete", False)):
            raise NetworkViewError(
                f"Group cell {cell.source_id!r} is not complete."
            )
        try:
            volume += float(cell.data["vol"])
        except (KeyError, TypeError, ValueError) as error:
            raise NetworkViewError(
                f"Group cell {cell.source_id!r} has no usable volume."
            ) from error

    area = 0.0
    boundary_dependent_area = 0.0
    physical_interface_area = 0.0
    boundary_surface_ids = []
    competitor_ids = set()
    for surface in view.surfaces:
        # Legacy solved surface tables do not carry a per-row complete flag;
        # NetworkView completeness is the authoritative feature evidence.
        if "complete" in surface.data and not bool(surface.data["complete"]):
            raise NetworkViewError(
                f"Group surface {surface.source_id!r} is not complete."
            )
        generator_ids = tuple(surface.generator_system_ids)
        physical_ids = {value for value in generator_ids if value is not None}
        if not physical_ids.intersection(selected_ids):
            continue

        competitors = physical_ids.difference(selected_ids)
        exposed = bool(competitors) or any(value is None for value in generator_ids)
        if not exposed:
            continue

        raw_area = surface.data.get("sa")
        if raw_area is None:
            try:
                raw_area = calc_surf_sa(
                    tris=surface.data["tris"], points=surface.data["points"]
                )
            except (KeyError, TypeError, ValueError) as error:
                raise NetworkViewError(
                    f"Group surface {surface.source_id!r} has no usable area."
                ) from error
        if not isinstance(raw_area, Real):
            raise NetworkViewError(
                f"Group surface {surface.source_id!r} has no usable area."
            )
        area += float(raw_area)
        if any(value is None for value in generator_ids):
            boundary_dependent_area += float(raw_area)
        elif competitors:
            physical_interface_area += float(raw_area)
        boundary_surface_ids.append(surface.source_id)
        competitor_ids.update(competitors)

    return GroupViewAnalysis(
        source_network_id=view.source_network_id,
        solve_universe_id=view.solve_universe_id,
        selection_id=view.selection_id,
        system_ids=tuple(view.system_ids),
        boundary_surface_ids=tuple(boundary_surface_ids),
        competitor_system_ids=tuple(sorted(competitor_ids)),
        area=area,
        volume=volume,
        complete=True,
        boundary_dependent_area=boundary_dependent_area,
        physical_interface_area=physical_interface_area,
    )


__all__ = [
    "GroupViewAnalysis",
    "analyze_group_view",
    "attach_group_view_results",
    "make_group_view",
    "shared_network_identity",
    "system_ids_for_group",
]
