"""Reusable presentation containers; scientific state remains in the workbench.

Network context is intentionally a presentation-only model. A future FORGE /
ATLAS adapter may map validated context and provenance into it; this module
must not infer identity from group definitions, labels, or geometry counts.
"""
from dataclasses import dataclass
from collections.abc import Mapping
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView, QFormLayout, QFrame, QGroupBox, QHBoxLayout, QLabel,
    QHeaderView, QScrollArea, QTableWidget, QTableWidgetItem, QToolButton,
    QVBoxLayout, QWidget, QSizePolicy,
)
from vorpy.src.results.model import ResultModel
from vorpy.src.results.shared_network import is_provisional_source_network_id


@dataclass(frozen=True)
class NetworkContextPresentation:
    """Read-only values that the Workbench may safely present.

    ``identity_state`` is explicit so provisional or unresolved values cannot
    be styled as validated provenance. ``context_scope`` remains unknown until
    a producer proves global versus isolated ownership.
    """

    source_network_id: str | None = None
    solve_universe_id: str | None = None
    geometry_scheme: str | None = None
    frame: str | None = None
    context_scope: str = "unknown"
    solvent_mode: str | None = None
    solvent_mode_state: str = "unresolved"
    retained_generators: str | None = None
    completeness: str = "UNKNOWN"
    status: str = "NOT_CALCULATED"
    identity_state: str = "unresolved"
    note: str = ""


@dataclass(frozen=True)
class AnalysisCatalogEntryPresentation:
    """One producer-attested group or interface analysis for a network."""

    selection_label: str
    selection_id: str | None = None
    analysis_id: str | None = None
    context_scope: str = "unknown"
    completeness: str = "UNKNOWN"
    status: str = "UNRESOLVED"
    area: str | None = None
    volume: str | None = None
    selection_kind: str = "unknown"
    selection_target: str | None = None


@dataclass(frozen=True)
class AnalysisCatalogPresentation:
    """Read-only catalog data; it never infers analysis from a definition."""

    source_network_id: str | None = None
    identity_state: str = "unresolved"
    entries: tuple[AnalysisCatalogEntryPresentation, ...] = ()
    note: str = ""


class NetworkContextPanel(QGroupBox):
    """Compact, read-only network identity and completeness presentation."""

    _SCOPE_LABELS = {
        "global": "Global network",
        "global_network": "Global network",
        "isolated": "Isolated geometry",
        "isolated_group": "Isolated group",
        "isolated_interface": "Isolated interface",
        "unknown": "Unknown",
    }
    _IDENTITY_LABELS = {
        "validated": "Validated identity",
        "provisional": "Provisional identity",
        "unresolved": "Identity unresolved",
    }
    _IDENTITY_COLORS = {
        "validated": "#5cc98a",
        "provisional": "#e5b94c",
        "unresolved": "#ef8f8f",
    }
    _SOLVENT_MODE_LABELS = {
        "molecule_only": "Molecule-only",
        "solvent_context": "Solvent context",
        "full_solvent": "Full solvent",
    }

    def __init__(self, parent=None):
        super().__init__("Network context", parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 8)
        layout.setSpacing(4)

        self.identity_badge = QLabel()
        self.identity_badge.setObjectName("networkContextBadge")
        self.identity_badge.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.identity_badge)

        self.fields: dict[str, QLabel] = {}
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        for key, label in (
            ("source_network_id", "Source network ID"),
            ("solve_universe_id", "Solve universe ID"),
            ("geometry_scheme", "Geometry scheme"),
            ("frame", "Frame"),
            ("context_scope", "Analysis context"),
            ("solvent_mode", "Solvent mode"),
            ("retained_generators", "Retained generators"),
            ("completeness", "Completeness"),
            ("status", "Status"),
        ):
            value = QLabel("—")
            value.setTextInteractionFlags(Qt.TextSelectableByMouse)
            value.setWordWrap(True)
            self.fields[key] = value
            form.addRow(label, value)
        layout.addLayout(form)

        self.note = QLabel()
        self.note.setObjectName("sectionLabel")
        self.note.setWordWrap(True)
        self.note.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.note)
        self.set_context(None)

    @property
    def context(self) -> NetworkContextPresentation | None:
        return self._context

    def set_context(self, context: NetworkContextPresentation | None) -> None:
        """Render a presentation snapshot; never computes or validates it."""
        self._context = context
        if context is None:
            self.identity_badge.setText("No network loaded")
            self.identity_badge.setStyleSheet("color: #a9b9c9;")
            for value in self.fields.values():
                value.setText("—")
            self.note.setText(
                "Group and interface definitions are not analyzed network contexts."
            )
            return

        identity_state = context.identity_state.lower()
        self.identity_badge.setText(
            self._IDENTITY_LABELS.get(identity_state, "Identity unresolved")
        )
        self.identity_badge.setStyleSheet(
            f"color: {self._IDENTITY_COLORS.get(identity_state, self._IDENTITY_COLORS['unresolved'])};"
        )
        values = {
            "source_network_id": context.source_network_id,
            "solve_universe_id": context.solve_universe_id,
            "geometry_scheme": context.geometry_scheme,
            "frame": context.frame,
            "context_scope": self._SCOPE_LABELS.get(
                context.context_scope.lower(), "Unknown"
            ),
            "solvent_mode": self._solvent_mode_label(
                context.solvent_mode, context.solvent_mode_state
            ),
            "retained_generators": context.retained_generators,
            "completeness": context.completeness,
            "status": context.status,
        }
        for key, value in values.items():
            self.fields[key].setText(str(value) if value not in (None, "") else "—")
        self.note.setText(context.note)

    @classmethod
    def _solvent_mode_label(cls, mode: str | None, state: str) -> str:
        """Only validated metadata may use the canonical solvent labels."""
        if not mode:
            return "Unknown / unresolved"
        label = cls._SOLVENT_MODE_LABELS.get(mode.lower(), str(mode))
        state = state.lower()
        if state == "validated":
            return label
        if state == "legacy":
            return f"Legacy / unresolved ({label})"
        if state == "provisional":
            return f"Provisional ({label})"
        return f"Unresolved ({label})"


_SOLVENT_MODES = {"MOLECULE_ONLY", "SOLVENT_CONTEXT", "FULL_SOLVENT"}


def _valid_generator_membership(value) -> bool:
    """Accept only the five-field F10 local/system membership evidence."""
    if not isinstance(value, (tuple, list)) or not value:
        return False
    for record in value:
        if not isinstance(record, (tuple, list)) or len(record) != 5:
            return False
        local_index, system_index, stable_id, is_solvent, is_boundary = record
        if not isinstance(local_index, int) or local_index < 0:
            return False
        if system_index is not None and not isinstance(system_index, int):
            return False
        if not isinstance(stable_id, str) or not stable_id:
            return False
        if type(is_solvent) is not bool or type(is_boundary) is not bool:
            return False
    return True


def _solvent_presentation(provenance: Mapping | None, legacy_mode=None):
    """Return mode plus evidence state without guessing from labels or counts."""
    if not isinstance(provenance, Mapping):
        if legacy_mode is not None:
            legacy_mode = getattr(legacy_mode, "value", legacy_mode)
            return (
                str(legacy_mode).strip().upper(),
                "legacy",
                "Legacy solvent mode has no producer membership/completeness evidence.",
            )
        return None, "unresolved", "Producer solvent evidence is unavailable."

    raw_mode = provenance.get("solvent_mode")
    mode = str(raw_mode).strip().upper() if raw_mode is not None else ""
    membership = provenance.get("generator_membership")
    if mode not in _SOLVENT_MODES or not _valid_generator_membership(membership):
        return (
            mode or None,
            "unresolved",
            "Solvent mode is present, but producer membership evidence is incomplete.",
        )

    generator_ids = provenance.get("solvent_generator_ids")
    digest = provenance.get("solvent_generator_digest")
    cells_complete = provenance.get("solvent_cells_complete")
    if mode == "MOLECULE_ONLY":
        valid = (
            not generator_ids
            and not digest
            and cells_complete is False
            and not any(record[3] for record in membership)
        )
    elif mode == "SOLVENT_CONTEXT":
        valid = bool(generator_ids) and bool(digest) and cells_complete is False
    else:
        valid = bool(generator_ids) and bool(digest) and cells_complete is True
    if not valid:
        return (
            mode,
            "unresolved",
            "Solvent mode is present, but producer completeness evidence is inconsistent.",
        )
    return mode, "validated", "Solvent mode is backed by producer membership and completeness evidence."


def _retained_generator_presentation(provenance: Mapping | None) -> str | None:
    """Render F12 retention evidence without inferring it from atom counts."""
    if not isinstance(provenance, Mapping):
        return None
    indices = provenance.get("retained_system_indices")
    generator_ids = provenance.get("retained_generator_ids")
    if not isinstance(indices, (tuple, list)) or not isinstance(
        generator_ids, (tuple, list)
    ):
        return "Unresolved / incomplete"
    if len(indices) != len(generator_ids):
        return "Unresolved / inconsistent"
    if any(not isinstance(value, int) or value < 0 for value in indices):
        return "Unresolved / inconsistent"
    if any(not isinstance(value, str) or not value for value in generator_ids):
        return "Unresolved / inconsistent"
    return f"{len(indices)} retained (explicit policy)"


def context_presentation_from_owner(owner, *, frame=None):
    """Adapt a producer-owned Group/Network to the read-only UI model.

    This is intentionally duck-typed: the Workbench does not import or mutate
    FORGE/NEXUS objects, and missing evidence remains unresolved.
    """
    network = getattr(owner, "net", None)
    if network is None:
        return None

    key = None
    get_context = getattr(network, "get_solve_context_key", None)
    if callable(get_context):
        try:
            key = get_context(frame_id=frame)
        except (TypeError, ValueError, RuntimeError):
            key = None

    # Archive object identity is not producer provenance.  Keep the context
    # unresolved until the producer supplies an explicit source namespace.
    source_network_id = getattr(network, "source_network_id", None)
    source_network_id = str(source_network_id) if source_network_id else None
    identity_state = "unresolved"
    if source_network_id and key is not None:
        identity_state = (
            "provisional"
            if source_network_id == f"network:{key.digest}"
            else "validated"
        )
    elif source_network_id:
        identity_state = "provisional"

    scope_mode = getattr(key, "scope_mode", None)
    if getattr(network, "iface_grps", None) is not None:
        context_scope = "isolated_interface"
    elif scope_mode == "full":
        context_scope = "global_network"
    elif scope_mode == "group":
        context_scope = "isolated_group"
    else:
        context_scope = "unknown"

    provenance = getattr(network, "solvent_provenance", None)
    legacy_mode = getattr(network, "solvent_mode", None)
    if legacy_mode is None and isinstance(getattr(network, "settings", None), Mapping):
        legacy_mode = network.settings.get("solvent_mode")
    solvent_mode, solvent_state, solvent_note = _solvent_presentation(
        provenance, legacy_mode
    )
    retained_note = _retained_generator_presentation(
        getattr(network, "retained_generator_provenance", None)
    )
    notes = [solvent_note]
    if retained_note is None:
        notes.append("Retained-generator evidence is unavailable.")
    if not source_network_id:
        notes.append("Source network namespace is not producer-attested.")
        if getattr(network, "archive_id", None):
            notes.append("Archive object identity is not used as source provenance.")
    if key is None:
        notes.append("Solve-universe identity is unavailable.")
    return NetworkContextPresentation(
        source_network_id=source_network_id,
        solve_universe_id=getattr(key, "solve_universe_id", None),
        geometry_scheme=getattr(key, "scheme", None),
        frame=getattr(key, "frame_id", None),
        context_scope=context_scope,
        solvent_mode=solvent_mode,
        solvent_mode_state=solvent_state,
        retained_generators=retained_note,
        completeness="UNKNOWN",
        status="UNRESOLVED",
        identity_state=identity_state,
        note=" ".join(notes),
    )


def _typed_result_list(results):
    if isinstance(results, ResultModel):
        results = (results,)
    try:
        values = tuple(results or ())
    except TypeError:
        return ()
    return tuple(
        result for result in values
        if isinstance(result, ResultModel)
        and result.identity.result_kind in {"group", "interface"}
    )


def _typed_identity_context(result):
    identity = result.identity
    required = (
        identity.source_network_id,
        identity.solve_universe_id,
        identity.analysis_context_kind,
        identity.selection_id,
        identity.solvent_mode,
    )
    if any(value in (None, "") for value in required):
        return None
    return (
        identity.source_network_id,
        identity.solve_universe_id,
        str(identity.frame) if identity.frame is not None else None,
        identity.analysis_context_kind,
        identity.solvent_mode,
    )


def _typed_geometry_scheme(result):
    details = getattr(result.provenance, "details", {})
    shared = details.get("shared_network_identity")
    if isinstance(shared, Mapping) and shared.get("geometry_scheme"):
        return str(shared["geometry_scheme"])
    return None


def _typed_completeness(result):
    details = getattr(result.provenance, "details", {})
    completeness = details.get("completeness")
    if not isinstance(completeness, Mapping):
        return "UNKNOWN"
    proven = completeness.get(
        "proven", completeness.get("completeness_proven")
    )
    if proven is True:
        return "PROVEN"
    if proven is False:
        return "NOT PROVEN"
    return "UNKNOWN"


def _typed_quantity_text(result, name):
    quantity = getattr(result, "metrics", {}).get(name)
    value = getattr(quantity, "value", None)
    if quantity is None or value is None:
        return None
    try:
        rendered = f"{float(value):g}"
    except (TypeError, ValueError):
        rendered = str(value)
    units = getattr(quantity, "units", "")
    return f"{rendered} {units}".strip()


def catalog_presentation_from_results(results) -> AnalysisCatalogPresentation | None:
    """Adapt one coherent set of typed group/interface Results for display."""
    typed = _typed_result_list(results)
    if not typed:
        return None
    contexts = tuple(_typed_identity_context(result) for result in typed)
    if any(context is None for context in contexts):
        return AnalysisCatalogPresentation(
            note="Typed result identity is incomplete; no analysis is presented."
        )
    first = contexts[0]
    if any(context != first for context in contexts[1:]):
        return AnalysisCatalogPresentation(
            note="Results belong to different solve contexts; choose one network."
        )
    source_network_id, solve_universe_id, _, _, _ = first
    if is_provisional_source_network_id(source_network_id, solve_universe_id):
        return AnalysisCatalogPresentation(
            source_network_id=source_network_id,
            identity_state="provisional",
            note="Source network identity is provisional; analyses remain unresolved.",
        )

    entries = []
    for result in typed:
        identity = result.identity
        if identity.result_kind == "group":
            label = f"Group: {identity.group_id}"
            target = identity.group_id
        else:
            label = f"Interface: {identity.interface_id}"
            target = identity.interface_id
        entries.append(AnalysisCatalogEntryPresentation(
            selection_label=label,
            selection_id=identity.selection_id,
            analysis_id=result.result_id,
            context_scope={
                "GLOBAL_NETWORK": "global_network",
                "ISOLATED_GROUP": "isolated_group",
            }.get(identity.analysis_context_kind, "unknown"),
            completeness=_typed_completeness(result),
            status=getattr(result.status, "value", str(result.status)),
            area=_typed_quantity_text(result, "area"),
            volume=_typed_quantity_text(result, "volume"),
            selection_kind=identity.result_kind,
            selection_target=target,
        ))
    return AnalysisCatalogPresentation(
        source_network_id=source_network_id,
        identity_state="validated",
        entries=tuple(entries),
    )


def validated_group_system_ids(result, entry) -> tuple[int, ...] | None:
    """Resolve a catalog group through producer-owned System IDs only.

    The Workbench may use this read-only mapping for an overlay. A missing or
    incomplete producer Group/NetworkView fails closed instead of treating a
    catalog label as atom membership.
    """
    if getattr(entry, "selection_kind", None) != "group":
        return None
    owner = getattr(result, "export_group", None)
    network = getattr(owner, "net", None)
    if network is None:
        return None
    target = getattr(entry, "selection_target", None)
    transient = getattr(result, "typed_selection_system_ids", None)
    if isinstance(transient, Mapping) and target in transient:
        system_ids = tuple(transient[target])
        atom_ids = {
            getattr(atom, "index", None) for atom in getattr(result, "atoms", ())
        }
        if system_ids and set(system_ids).issubset(atom_ids):
            return system_ids
    groups = [owner]
    groups.extend(tuple(getattr(getattr(network, "sys", None), "groups", ()) or ()))
    group = next(
        (
            candidate
            for candidate in groups
            if (getattr(candidate, "group_id", None) or getattr(candidate, "name", None))
            == target
        ),
        None,
    )
    if group is None:
        return None
    try:
        from vorpy.src.group.shared_view import make_group_view

        view = make_group_view(group, network)
    except (AttributeError, KeyError, TypeError, ValueError, RuntimeError):
        return None
    if not getattr(view, "completeness", None) or not view.completeness.proven:
        return None
    system_ids = tuple(view.system_ids)
    atom_ids = {getattr(atom, "index", None) for atom in getattr(result, "atoms", ())}
    if not system_ids or not set(system_ids).issubset(atom_ids):
        return None
    return system_ids


def context_presentation_from_results(results) -> NetworkContextPresentation | None:
    """Project one coherent typed-Results context without certifying metrics."""
    typed = _typed_result_list(results)
    if not typed:
        return None
    contexts = tuple(_typed_identity_context(result) for result in typed)
    if any(context is None for context in contexts) or any(
        context != contexts[0] for context in contexts[1:]
    ):
        return None
    source_network_id, solve_universe_id, frame, context_kind, solvent_mode = contexts[0]
    provisional = is_provisional_source_network_id(
        source_network_id, solve_universe_id
    )
    solvent_evidence = []
    for result in typed:
        details = getattr(result.provenance, "details", {})
        evidence = details.get("source_solvent_provenance")
        shared = details.get("shared_network_identity")
        if not isinstance(evidence, Mapping) or not isinstance(shared, Mapping):
            solvent_evidence = []
            break
        if any(
            evidence.get(name) != (
                result.identity.solvent_mode
                if name == "solvent_mode"
                else shared.get(name)
            )
            for name in (
                "solvent_mode",
                "solvent_generator_digest",
                "solvent_cells_complete",
            )
        ):
            solvent_evidence = []
            break
        solvent_evidence.append(evidence)
    solvent_state = "validated" if solvent_evidence else "unresolved"
    statuses = {
        getattr(result.status, "value", str(result.status)) for result in typed
    }
    completeness = {_typed_completeness(result) for result in typed}
    status = next(iter(statuses)) if len(statuses) == 1 else "MIXED"
    return NetworkContextPresentation(
        source_network_id=source_network_id,
        solve_universe_id=solve_universe_id,
        geometry_scheme=_typed_geometry_scheme(typed[0]),
        frame=frame,
        context_scope={
            "GLOBAL_NETWORK": "global_network",
            "ISOLATED_GROUP": "isolated_group",
        }.get(context_kind, "unknown"),
        solvent_mode=solvent_mode,
        solvent_mode_state=solvent_state,
        completeness=(next(iter(completeness)) if len(completeness) == 1 else "MIXED"),
        status=status,
        identity_state="provisional" if provisional else "validated",
        note=(
            "Typed Results share one solve context; metric status remains per analysis."
            if solvent_state == "validated"
            else "Typed Results omit matching source solvent evidence; solvent mode is unresolved."
        ),
    )


class AnalysisCatalogPanel(QGroupBox):
    """Read-only table for multiple analyses over one producer network."""

    active_selection_changed = Signal(object)

    _HEADERS = (
        "Selection", "Selection ID", "Analysis ID", "Context",
        "Completeness", "Status", "Exposed area (boundary-dependent)",
        "View volume (boundary-dependent)",
    )

    def __init__(self, parent=None):
        super().__init__("Analysis catalog", parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 8)
        layout.setSpacing(4)

        self.network_label = QLabel()
        self.network_label.setObjectName("sectionLabel")
        layout.addWidget(self.network_label)

        self.empty = QLabel()
        self.empty.setObjectName("emptyState")
        self.empty.setWordWrap(True)
        layout.addWidget(self.empty)

        self.table = QTableWidget(0, len(self._HEADERS))
        self.table.setHorizontalHeaderLabels(self._HEADERS)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setSortingEnabled(False)
        self.table.cellClicked.connect(self._activate_row)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        for column in range(1, len(self._HEADERS)):
            self.table.horizontalHeader().setSectionResizeMode(
                column, QHeaderView.ResizeToContents
            )
        layout.addWidget(self.table)
        self.set_catalog(None)

    @property
    def catalog(self) -> AnalysisCatalogPresentation | None:
        return self._catalog

    @property
    def active_entry(self) -> AnalysisCatalogEntryPresentation | None:
        return self._active_entry

    def set_catalog(self, catalog: AnalysisCatalogPresentation | None) -> None:
        """Render only supplied entries; group/interface definitions are not entries."""
        self._catalog = catalog
        self._active_entry = None
        self.table.clearSelection()
        if catalog is None:
            self.network_label.setText("No network loaded")
            self.empty.setText(
                "No analyses associated with this network.\n"
                "Defined groups and interfaces are not analyzed results."
            )
        else:
            state = catalog.identity_state.lower()
            prefix = {
                "validated": "Network",
                "provisional": "Provisional network",
            }.get(state, "Unresolved network")
            self.network_label.setText(
                f"{prefix}: {catalog.source_network_id or '—'}"
            )
            self.empty.setText(catalog.note or (
                "No analyses associated with this network.\n"
                "Select a producer-attested group or interface result to populate this catalog."
            ))
        entries = catalog.entries if catalog is not None else ()
        self.empty.setVisible(not entries)
        self.table.setVisible(bool(entries))
        self.table.setRowCount(len(entries))
        for row, entry in enumerate(entries):
            values = (
                entry.selection_label,
                entry.selection_id,
                entry.analysis_id,
                NetworkContextPanel._SCOPE_LABELS.get(
                    entry.context_scope.lower(), "Unknown"
                ),
                entry.completeness,
                entry.status,
                entry.area,
                entry.volume,
            )
            for column, value in enumerate(values):
                self.table.setItem(
                    row,
                    column,
                    QTableWidgetItem(str(value) if value not in (None, "") else "—"),
                )

    def _activate_row(self, row, _column):
        entries = self._catalog.entries if self._catalog is not None else ()
        if not 0 <= row < len(entries):
            return
        self._active_entry = entries[row]
        self.active_selection_changed.emit(self._active_entry)


def action_button(action, text=None):
    button = QToolButton()
    button.setDefaultAction(action)
    if text:
        button.setText(text)
        action.changed.connect(lambda: button.setText(text))
    button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
    return button


def scroll_panel(widget):
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.NoFrame)
    scroll.setWidget(widget)
    return scroll


class WorkflowSidebar(QWidget):
    def __init__(self, tabs):
        super().__init__()
        self.setObjectName("workflowPanel")
        self.setMinimumWidth(240)
        self.setMaximumWidth(320)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(tabs)


class ViewerPanel(QWidget):
    def __init__(self, viewer, actions):
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 0)
        layout.setSpacing(4)
        bar = QFrame()
        bar.setObjectName("viewerBar")
        tools = QHBoxLayout(bar)
        tools.setContentsMargins(6, 3, 6, 3)
        tools.setSpacing(3)
        for action, label in actions:
            tools.addWidget(action_button(action, label))
        tools.addStretch()
        bar.setToolTip("Drag to rotate · middle-drag to pan · wheel for depth clipping · Shift-click to toggle selection")
        layout.addWidget(bar)
        layout.addWidget(viewer, 1)
        self.metadata = QLabel("No structure loaded")
        self.metadata.setObjectName("sectionLabel")
        self.metadata.setWordWrap(True)
        layout.addWidget(self.metadata)
        self.selection_label = QLabel("Residue selection active")
        self.selection_label.setObjectName("sectionLabel")
        layout.addWidget(self.selection_label)


class ViewInspector(QWidget):
    def __init__(self, tabs):
        super().__init__()
        self.setMinimumWidth(360)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 6)
        title = QLabel("View Settings")
        title.setObjectName("viewSettingsTitle")
        layout.addWidget(title)
        layout.addWidget(tabs, 1)


class ResultsInspector(QWidget):
    def __init__(self, content):
        super().__init__()
        self.setObjectName("resultsInspector")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 0, 8, 4)
        layout.setSpacing(2)
        content.setMinimumHeight(0)
        content.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Ignored)
        self.content = content
        layout.addWidget(content, 1)
