"""Selected source-network rows used by physical exporters."""

from __future__ import annotations

from copy import copy
from dataclasses import dataclass


def _complete_value(value):
    if value is None:
        return False
    try:
        if value != value:  # NaN without importing NumPy into this small module.
            return False
    except (TypeError, ValueError):
        return False
    try:
        return bool(value)
    except (TypeError, ValueError):
        return False


def _validate_materialized_interface_subset(view, network, positions):
    """Permit only the locally proven portion of a partial interface view."""
    if getattr(view, "selection_kind", None) != "interface":
        raise ValueError("Cannot export an incomplete materialized geometry view.")

    completeness = getattr(view, "completeness", None)
    unresolved = tuple(getattr(completeness, "unresolved_surface_ids", ()) or ())
    if not unresolved or not positions["surface"]:
        raise ValueError(
            "Cannot export an incomplete materialized interface without a "
            "supported physical subset."
        )

    for kind, table_name in (
        ("surface", "surfs"),
        ("edge", "edges"),
        ("vertex", "verts"),
        ("cell", "balls"),
    ):
        table = getattr(network, table_name, None)
        if table is None or "complete" not in table.columns:
            continue
        for position in positions[kind]:
            if not _complete_value(table.iloc[int(position)].get("complete")):
                raise ValueError(
                    f"Cannot export locally incomplete materialized {kind} "
                    f"row {position}."
                )


@dataclass(frozen=True)
class ExportScope:
    """A validated selection over one solved Network."""

    network: object
    view: object | None
    cells: tuple[int, ...]
    surfaces: tuple[int, ...]
    edges: tuple[int, ...]
    vertices: tuple[int, ...]
    support_generators: tuple[int, ...]

    @property
    def is_materialized_view(self) -> bool:
        return self.view is not None

    def source_ids(self, kind: str) -> tuple[object, ...]:
        if self.view is not None:
            return tuple(feature.source_id for feature in self.view.features(kind))
        table_name = {
            "cell": "balls", "surface": "surfs", "edge": "edges", "vertex": "verts",
        }[kind]
        table = getattr(self.network, table_name, None)
        positions = {
            "cell": self.cells, "surface": self.surfaces,
            "edge": self.edges, "vertex": self.vertices,
        }[kind]
        return tuple(table.index[int(position)] for position in positions)

    def scoped_network(self):
        """Return a small table-only proxy for CSV/info export.

        Incidence values remain in the source topology namespace. Mesh writers
        use the source Network and the original row positions directly.
        """
        if self.view is None and all(
            positions == tuple(range(len(getattr(self.network, name))))
            for name, positions in (
                ("balls", self.cells), ("surfs", self.surfaces),
                ("edges", self.edges), ("verts", self.vertices),
            )
        ):
            return self.network

        scoped = copy(self.network)
        for name, positions in (
            ("balls", self.support_generators),
            ("surfs", self.surfaces),
            ("edges", self.edges),
            ("verts", self.vertices),
        ):
            table = getattr(self.network, name, None)
            setattr(scoped, name, table.iloc[list(positions)].copy())

        boundary = set(int(value) for value in (getattr(self.network, "boundary_indices", ()) or ()))
        scoped.boundary_indices = tuple(sorted(boundary.intersection(self.support_generators)))
        return scoped

    def mesh_network(self):
        """Return the source network with only harmless display defaults filled."""
        settings = getattr(self.network, "settings", {}) or {}
        required = {"surf_scheme", "surf_col", "scheme_factor", "net_type"}
        if required.issubset(settings):
            return self.network
        mesh_network = copy(self.network)
        mesh_network.settings = dict(settings)
        mesh_network.settings.setdefault("surf_scheme", "int_mean_curv")
        mesh_network.settings.setdefault("surf_col", "coolwarm")
        mesh_network.settings.setdefault("scheme_factor", "log")
        mesh_network.settings.setdefault("net_type", "aw")
        return mesh_network


def export_scope(owner) -> ExportScope:
    """Resolve a Group/Interface to its explicit source selection."""
    view = getattr(owner, "source_view", None)
    network = (
        getattr(owner, "source_network", None)
        if view is not None else getattr(owner, "net", None)
    )
    if network is None:
        network = getattr(owner, "net", None)
    if network is None:
        raise ValueError("Export requires a solved source Network.")

    if view is None:
        def count(name):
            table = getattr(network, name, None)
            return 0 if table is None else len(table)

        return ExportScope(
            network=network,
            view=None,
            cells=tuple(range(count("balls"))),
            surfaces=tuple(range(count("surfs"))),
            edges=tuple(range(count("edges"))),
            vertices=tuple(range(count("verts"))),
            support_generators=tuple(range(count("balls"))),
        )

    positions = {}
    for kind, table_name in (
        ("cell", "balls"),
        ("surface", "surfs"),
        ("edge", "edges"),
        ("vertex", "verts"),
    ):
        table = getattr(network, table_name, None)
        values = tuple(int(feature.topology_index) for feature in view.features(kind))
        if any(value < 0 or value >= len(table) for value in values):
            raise ValueError(f"Materialized {kind} selection is outside the source table.")
        source_ids = tuple(table.index[value] for value in values)
        selected_ids = tuple(feature.source_id for feature in view.features(kind))
        if source_ids != selected_ids:
            raise ValueError(f"Materialized {kind} source IDs do not match the source table.")
        positions[kind] = values

    if not getattr(getattr(view, "completeness", None), "proven", False):
        _validate_materialized_interface_subset(view, network, positions)

    support = tuple(sorted({int(value) for value in view.support_generator_local_indices}))
    if any(value < 0 or value >= len(network.balls) for value in support):
        raise ValueError("Materialized generator support is outside the source table.")

    return ExportScope(
        network=network,
        view=view,
        cells=positions["cell"],
        surfaces=positions["surface"],
        edges=positions["edge"],
        vertices=positions["vertex"],
        support_generators=support,
    )


__all__ = ["ExportScope", "export_scope"]
