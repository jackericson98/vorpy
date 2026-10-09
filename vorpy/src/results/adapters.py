"""Adapters over cached Python records. No network/geometry imports or file IO.

Aliases below are serialization names for InterfaceGeometryAnalysis fields,
not new scientific computations. Missing ratios, counts, and terms stay missing.
"""

from collections.abc import Mapping

from .model import (
    AtomIdentity,
    AlphaSelectionResult,
    ContactRecord,
    DualContactComplexResult,
    DualGeneratorComplexResult,
    DualSimplexRecord,
    Environment,
    GaussianCurvature,
    GroupResult,
    InterfaceResult,
    GeneratorComplex,
    MolecularContactSurfaceResult,
    MeanCurvature,
    NetworkResult,
    Provenance,
    Quantity,
    RepresentationExtension,
    Status,
    SystemResult,
    Topology,
    normalize_environment,
    normalize_network_scope,
)

MEAN_CONVENTION = "sum(smooth face H) + 1/2 sum(internal seam beta ds); free boundary has no automatic scalar H term"
GAUSSIAN_CONVENTION = "intrinsic surface integral; boundary geodesic/corner terms complete Gauss-Bonnet separately; extrinsic fold beta is not Gaussian curvature"


def _read(source, name, default=None):
    return (
        source.get(name, default)
        if isinstance(source, Mapping)
        else getattr(source, name, default)
    )


def _first(source, names):
    for name in names:
        if name in source:
            return name, source[name]
    return None, None


def _metadata_value(metadata, names):
    for name in names:
        value = metadata.get(name)
        if value is not None:
            return name, value
    return None, None


def _same_value(actual, expected):
    if (
        isinstance(actual, (int, float))
        and not isinstance(actual, bool)
        and isinstance(expected, (int, float))
        and not isinstance(expected, bool)
    ):
        return float(actual) == float(expected)
    return actual == expected


def _check_cached_identity(metadata, identity):
    """Validate known cache identity; absent legacy fields remain unvalidated."""

    def check(label, names, expected, normalizer=None, *, require_identity=False):
        key, actual = _metadata_value(metadata, names)
        if key is None:
            return
        try:
            actual = normalizer(actual) if normalizer else actual
        except ValueError as error:
            raise ValueError(f"Cached {label} is invalid: {actual!r}") from error
        if expected is None:
            if require_identity:
                raise ValueError(f"Cached {label} is present but ResultIdentity is missing it")
            return
        if not _same_value(actual, expected):
            raise ValueError(
                f"Cached identity mismatch for {label}: {actual!r} != {expected!r}"
            )

    check("system ID", ("system_id", "stable_system_id"), identity.system_id or identity.system)
    check("system", ("system", "system_name"), identity.system)
    check("network ID", ("network_id", "network_archive_id"), identity.network_id)
    check(
        "network scope",
        ("network_scope", "scope"),
        identity.network_scope,
        normalize_network_scope,
    )
    # Existing InterfaceGeometryAnalysis.interface_id identifies that cache,
    # not necessarily the logical ResultIdentity interface. New persisted
    # state must use an explicit result/logical interface field.
    check(
        "interface ID",
        ("result_interface_id", "logical_interface_id", "interface_identity_id"),
        identity.interface_id,
    )
    check("group A", ("group_A", "group_a"), identity.group_A)
    check("group B", ("group_B", "group_b"), identity.group_B)
    check(
        "frame",
        ("frame", "frame_id", "frame_index"),
        identity.frame,
        require_identity=True,
    )

    def normalize_partition(value):
        aliases = {"aw": "aw", "power": "power", "pow": "power"}
        try:
            return aliases[value]
        except (KeyError, TypeError):
            raise ValueError(f"unknown partition {value!r}") from None

    check("partition", ("partition", "scheme"), identity.partition, normalize_partition)
    check(
        "radius configuration",
        ("radius_configuration_id", "radii_configuration_id", "radius_config_id"),
        identity.radius_configuration_id,
    )
    check(
        "alpha convention",
        ("alpha_method", "alpha_convention"),
        identity.alpha_method,
        require_identity=True,
    )
    check("alpha value", ("alpha_value",), identity.alpha_value, require_identity=True)
    check("alpha units", ("alpha_units",), identity.alpha_units, require_identity=True)

    probe = metadata.get("probe")
    if isinstance(probe, Mapping):
        metadata = dict(metadata)
        metadata.setdefault("probe_radius", probe.get("radius", probe.get("probe_radius")))
        metadata.setdefault("probe_units", probe.get("units", probe.get("probe_units")))
    elif probe is not None and "probe_radius" not in metadata:
        metadata = dict(metadata, probe_radius=probe)
    check("probe radius", ("probe_radius",), identity.probe_radius, require_identity=True)
    check("probe units", ("probe_units",), identity.probe_units, require_identity=True)

    observed_environment = []
    key, value = _metadata_value(metadata, ("environment", "solvent_environment"))
    if key is not None:
        try:
            observed_environment.append(normalize_environment(value))
        except ValueError as error:
            raise ValueError(f"Cached environment is invalid: {value!r}") from error
    key, value = _metadata_value(metadata, ("solvent_context",))
    if key is not None:
        if not isinstance(value, bool):
            try:
                value = normalize_environment(value) == "solvent_competing"
            except ValueError as error:
                raise ValueError(f"Cached solvent context is invalid: {value!r}") from error
        observed_environment.append("solvent_competing" if value else "dry")
    for name in ("solvent_competing_generator_count", "solvent_interface_competitors"):
        key, value = _metadata_value(metadata, (name,))
        if key is not None:
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"Cached {name} is invalid: {value!r}")
            observed_environment.append("solvent_competing" if value else "dry")
    if observed_environment and any(value != observed_environment[0] for value in observed_environment):
        raise ValueError("Cached solvent metadata is contradictory")
    if observed_environment and observed_environment[0] != identity.environment:
        raise ValueError(
            "Cached identity mismatch for environment: "
            f"{observed_environment[0]!r} != {identity.environment!r}"
        )


def _quantity(
    source,
    name,
    aliases,
    units,
    scope,
    provenance,
    *,
    certified=False,
    support=None,
    total=None,
    missing_status=None,
):
    key, value = _first(source, aliases)
    status = Status.CERTIFIED if certified else Status.PARTIAL
    if value is None:
        status = (
            missing_status
            if missing_status is not None
            else Status.UNRESOLVED
            if key is not None
            else Status.NOT_CALCULATED
        )
    return Quantity(
        name,
        value,
        units,
        status,
        "cached " + (key or aliases[0]),
        scope,
        support,
        total,
        provenance,
    )


def unavailable_interface(identity, *, status=Status.NOT_CALCULATED, provenance=None):
    """Represent an absent/unsupported future representation without geometry."""
    status = Status(status)
    if status not in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED}:
        raise ValueError("Absent representations require unavailable status")
    provenance = provenance or Provenance("cached scientific state unavailable")

    def q(name, units="1"):
        return Quantity(
            name,
            None,
            units,
            status,
            "unavailable cached state",
            "requested_representation",
            provenance=provenance,
        )

    return InterfaceResult(
        identity,
        provenance,
        status,
        Environment(identity.environment),
        metrics={"area": q("area", "Å²")},
        topology=Topology({name: q(name) for name in TOPOLOGY_ALIASES}),
        mean_curvature=MeanCurvature(
            q("integrated_mean_curvature", "Å"),
            q("complete_mean_curvature", "Å"),
            q("smooth_face_contribution", "Å"),
            q("edge_contribution", "Å"),
            q("boundary_contribution", "Å"),
            MEAN_CONVENTION,
        ),
        gaussian_curvature=GaussianCurvature(
            *[q(name) for name in GAUSSIAN_NAMES], GAUSSIAN_CONVENTION
        ),
        extensions=(
            RepresentationExtension(
                "representation.availability",
                "1",
                status,
                {"representation": identity.representation},
                provenance,
            ),
        ),
    )


TOPOLOGY_ALIASES = {
    "vertices": ("topology_vertices",),
    "edges": ("topology_edges",),
    "faces": ("topology_faces",),
    "euler_characteristic": ("euler_characteristic",),
    "connected_components": ("connected_components", "components"),
    "boundary_edges": ("boundary_edges",),
    "boundary_loops": ("boundary_loops",),
    "nonmanifold_edges": ("nonmanifold_edges",),
    "orientable": ("orientable",),
}
GAUSSIAN_NAMES = (
    "integrated_gaussian_curvature",
    "smooth_face_contribution",
    "intrinsic_discrete_contribution",
    "intrinsic_seam_contribution",
    "interior_vertex_defects",
    "boundary_geodesic_contribution",
    "boundary_corner_contribution",
    "gauss_bonnet_completion",
    "gauss_bonnet_expected",
    "gauss_bonnet_residual",
)


def adapt_interface_analysis(
    analysis,
    identity,
    *,
    provenance=None,
    contacts=(),
    extensions=(),
    cached_metrics=None,
):
    """Snapshot existing InterfaceGeometryAnalysis, with an explicit identity.

    Only physical Voronoi caches are understood here. Dual-side result blocks
    must be supplied by their future producer; this adapter never copies primal
    measurements into a dual representation. It never visits Network properties.
    """
    if identity.representation != "voronoi":
        raise ValueError("Physical analysis cannot populate dual-side results")
    if analysis is None:
        return unavailable_interface(identity, provenance=provenance)
    metadata = _read(analysis, "metadata", {})
    if not isinstance(metadata, Mapping):
        raise ValueError("Cached metadata must be a mapping")
    _check_cached_identity(metadata, identity)
    source = (
        "voronoi_side_2"
        if identity.orientation == "group2 -> group1"
        else "voronoi_side_1"
    )
    if identity.orientation not in {"group1 -> group2", "group2 -> group1"}:
        raise ValueError("Use the cached physical orientation explicitly")
    side = _read(analysis, source, {})
    selection = _read(analysis, "alpha_selection", {})
    provenance = provenance or Provenance(
        "InterfaceGeometryAnalysis",
        details={"cache_metadata": metadata, "cache_side": source},
        diagnostics=tuple(_read(analysis, "unresolved", ())),
    )
    curvature_not_requested = (
        Status.NOT_CALCULATED
        if metadata.get("curvature_requested", True) is False
        else None
    )
    full = "full_selected_interface"
    mc = bool(
        side.get(
            "mean_curvature_certified",
            side.get("certified_mean_total", side.get("mean_certified", False)),
        )
    )
    gc = bool(
        side.get(
            "gaussian_curvature_certified",
            side.get(
                "certified_gaussian_total", side.get("gaussian_total_certified", False)
            ),
        )
    )
    topology_certified = bool(side.get("manifold_incidence_certified", False))
    complete_h = _quantity(
        side,
        "complete_mean_curvature",
        ("total_integrated_mean_curvature", "total_mean_curvature", "total_H"),
        "Å",
        full,
        provenance,
        certified=mc,
        missing_status=curvature_not_requested,
    )
    if mc:
        headline_h = Quantity(
            "integrated_mean_curvature",
            complete_h.value,
            "Å",
            complete_h.status,
            complete_h.method,
            full,
            provenance=provenance,
        )
    else:
        headline_h = _quantity(
            side,
            "integrated_mean_curvature",
            ("partial_integrated_mean_curvature", "partial_H"),
            "Å",
            "supported_contributions",
            provenance,
            missing_status=curvature_not_requested,
        )
    face_count = side.get("topology_faces")
    face_support = side.get("mean_supported_faces", side.get("supported_surfaces"))
    smooth = _quantity(
        side,
        "smooth_face_contribution",
        (
            "smooth_integrated_mean_curvature",
            "smooth_mean_curvature",
            "smooth_H_partial",
            "smooth_H",
        ),
        "Å",
        "supported_face_contributions",
        provenance,
        certified=mc,
        missing_status=curvature_not_requested,
        support=face_support,
        total=face_count,
    )
    edge = _quantity(
        side,
        "edge_contribution",
        (
            "internal_seam_integrated_mean_curvature",
            "internal_crease_mean_curvature",
            "edge_H",
        ),
        "Å",
        "supported_internal_seams",
        provenance,
        certified=mc,
        support=side.get("mean_supported_internal_seams"),
        total=side.get("mean_total_internal_seams"),
        missing_status=curvature_not_requested,
    )
    # The cache's boundary_H=0 encodes omission, not a measured boundary integral.
    boundary = Quantity(
        "boundary_contribution",
        None,
        "Å",
        Status.NOT_SUPPORTED,
        "validated free-boundary convention: no scalar mean-curvature term",
        "free_boundary",
        provenance=provenance,
    )
    mean = MeanCurvature(
        headline_h, complete_h, smooth, edge, boundary, MEAN_CONVENTION
    )
    gaussian_scope = side.get("gaussian_accounting_scope", full)
    aliases = (
        (
            "intrinsic_integrated_gaussian_curvature",
            "total_gaussian_curvature",
            "total_K",
        ),
        (
            "smooth_integrated_gaussian_curvature",
            "smooth_gaussian_curvature",
            "smooth_K_partial",
        ),
        ("intrinsic_discrete_gaussian_curvature",),
        ("intrinsic_seam_gaussian_curvature",),
        ("interior_vertex_gaussian_curvature", "intrinsic_vertex_gaussian"),
        ("boundary_geodesic_curvature", "boundary_geodesic_gaussian"),
        ("boundary_corner_turning", "boundary_corner_gaussian"),
        ("gauss_bonnet_lhs",),
        ("gauss_bonnet_expected",),
        ("gauss_bonnet_residual",),
    )
    gaussian = GaussianCurvature(
        *[
            _quantity(
                side,
                name,
                keys,
                "1",
                full
                if name in {"integrated_gaussian_curvature", "gauss_bonnet_expected"}
                else gaussian_scope,
                provenance,
                certified=gc,
                missing_status=curvature_not_requested,
            )
            for name, keys in zip(GAUSSIAN_NAMES, aliases)
        ],
        GAUSSIAN_CONVENTION,
    )
    genus = side.get("genus_by_component")
    topology = Topology(
        {
            name: _quantity(
                side,
                name,
                keys,
                "1",
                full,
                provenance,
                certified=topology_certified,
                missing_status=curvature_not_requested,
            )
            for name, keys in TOPOLOGY_ALIASES.items()
        },
        Status.CERTIFIED if topology_certified else Status.NOT_CALCULATED,
        tuple(genus) if genus is not None and topology_certified else None,
        "cached validated manifold component classification"
        if genus is not None and topology_certified
        else None,
    )
    metrics = {
        "area": _quantity(
            side,
            "area",
            ("area",),
            "Å²",
            side.get("area_scope", full),
            provenance,
            certified=bool(side.get("area_selection_complete", False)),
        ),
        "mean_curvature_per_area": _quantity(
            side,
            "mean_curvature_per_area",
            ("mean_curvature_per_area",),
            "Å^-1",
            headline_h.scope,
            provenance,
            certified=mc,
            missing_status=curvature_not_requested,
        ),
        "gaussian_curvature_per_area": _quantity(
            side,
            "gaussian_curvature_per_area",
            ("gaussian_curvature_per_area",),
            "Å^-2",
            full,
            provenance,
            certified=gc,
            missing_status=curvature_not_requested,
        ),
    }
    for name in (
        "n_atoms_A",
        "n_atoms_B",
        "n_interface_atoms_A",
        "n_interface_atoms_B",
        "n_interface_residues_A",
        "n_interface_residues_B",
        "n_contacts",
        "n_residue_contacts",
    ):
        metrics[name] = _quantity(side, name, (name,), "1", full, provenance)
    metrics["n_interface_atoms"] = _quantity(
        side, "n_interface_atoms", ("interface_atoms",), "1", full, provenance
    )
    for name in (
        "mean_unsupported_faces",
        "mean_unresolved_internal_seams",
        "gaussian_supported_faces",
        "gaussian_unsupported_faces",
    ):
        metrics[name] = _quantity(
            side,
            name,
            (name,),
            "1",
            "cached_support_audit",
            provenance,
            missing_status=curvature_not_requested,
        )
    # Metadata can provide precomputed metrics. No counting/dividing is performed.
    if cached_metrics:
        metrics.update(cached_metrics)
    selection_aliases = {
        "eligible_pairs": ("eligible_bicolor_pairs",),
        "candidate_pairs": ("candidate_bicolor_pairs",),
        "selected_features": ("selected_pairs",),
        "rejected_features": ("certified_rejected_pairs", "rejected_pairs"),
        "unresolved_features": ("unresolved_pairs",),
        "mapped_surfaces": ("selected_physical_pairs",),
        "excluded_pairs": ("excluded_pairs",),
    }
    selections = {
        name: _quantity(
            selection,
            name,
            keys,
            "1",
            "cached_pair_selection",
            provenance,
            certified=bool(selection.get("certified", False)),
        )
        for name, keys in selection_aliases.items()
    }
    # excluded_pairs may include unresolved pairs; never relabel as rejected.
    cache_extension = RepresentationExtension(
        "vorpy.interface_geometry_analysis",
        "1",
        Status(metadata.get("result_status", "PARTIAL")),
        {
            "selection": selection,
            "coverage": _read(analysis, "coverage", {}),
            "unresolved": _read(analysis, "unresolved", ()),
            "physical_accounting": side,
            "dry_solvent_comparison": _read(analysis, "dry_solvent_comparison", {}),
        },
        provenance,
    )
    environment = Environment(
        identity.environment,
        metadata.get("solvent_detected"),
        metadata.get("solvent_competing_generator_count"),
        metadata.get("solvent_interface_competitors"),
        provenance,
    )
    return InterfaceResult(
        identity,
        provenance,
        Status(metadata.get("result_status", "PARTIAL")),
        environment,
        metrics,
        topology,
        mean,
        gaussian,
        tuple(contacts),
        selections,
        (cache_extension, *extensions),
    )


def adapt_cached_result(identity, *, provenance, status=Status.PARTIAL, **blocks):
    """Common composition entry point for cached system/network/group producers.

    Also accepts future dual-side InterfaceResult blocks supplied by Agent #3.
    Producers supply Quantity objects with their own certification; the adapter
    does not infer geometry, scientific totals, or aggregate certification.
    """
    cls = {
        "system": SystemResult,
        "network": NetworkResult,
        "group": GroupResult,
        "interface": InterfaceResult,
        "molecular_contact_surface": MolecularContactSurfaceResult,
    }[identity.result_kind]
    return cls(identity=identity, provenance=provenance, status=status, **blocks)


def adapt_atom_metadata(metadata, *, system_id, network_id=None, generator_index=None):
    """Adapt already resolved molecular metadata; never use a generator as identity."""

    def first(*keys):
        for key in keys:
            value = metadata.get(key)
            if value is not None and str(value).strip() not in {"", "nan", "None"}:
                return str(value)
        return ""

    return AtomIdentity(
        system_id,
        first("chain", "chain_name"),
        first("residue_number", "auth_seq_id", "res_seq"),
        first("insertion_code", "pdb_ins_code", "icode"),
        first("residue_name", "auth_comp_id", "res_name"),
        first("atom_name", "auth_atom_id", "name"),
        first("alternate_location", "altloc"),
        first("system_atom_id", "system_num") or None,
        first("stable_atom_id", "stable_id") or None,
        generator_index,
        network_id,
    )


def adapt_contact_mapping(
    mapping, identity, atom_A, atom_B, *, quantities=None, provenance=None
):
    """Consume an existing selected/status mapping; no filtration query or attribution.

    Curvature attribution and area support must be supplied as Quantity records
    by the producer. In particular, no shared seam integral is assigned to a
    contact automatically. Chemistry side order is explicit, not generator order.
    """
    raw = mapping.get(
        "selection_status", mapping.get("restriction_status", "NOT_CALCULATED")
    )
    status = {
        "selected": "SELECTED",
        "certified_rejected": "REJECTED",
        "unresolved": "UNRESOLVED",
    }.get(raw, raw)
    return ContactRecord(
        identity.result_id,
        identity,
        atom_A,
        atom_B,
        status,
        mapping.get(
            "dual_feature_id",
            mapping.get("dual_simplex_id", mapping.get("dual_edge_id")),
        ),
        mapping.get(
            "physical_feature_id",
            mapping.get("physical_surface_id", mapping.get("surface_id")),
        ),
        quantities or {},
        provenance,
    )


def _require_representation(identity, allowed, label):
    if identity.representation not in allowed:
        raise ValueError(f"{label} identity has an incompatible representation")


def _unavailable_status(status):
    status = Status(status)
    if status not in {Status.NOT_CALCULATED, Status.NOT_SUPPORTED}:
        raise ValueError("Unavailable adapters require NOT_CALCULATED or NOT_SUPPORTED")
    return status


def _dual_provenance(provenance, source):
    if provenance is not None:
        return provenance
    return Provenance(source)


def _source_feature_id(value):
    if isinstance(value, str) and value:
        return value
    if isinstance(value, tuple) and len(value) == 2:
        return f"{value[0]}:{value[1]}"
    key = getattr(value, "key", None)
    if key is not None:
        return _source_feature_id(key)
    kind = getattr(value, "kind", None)
    feature_id = getattr(value, "feature_id", None)
    if kind is not None and feature_id is not None:
        return f"{kind}:{feature_id}"
    raise TypeError("Dual primal features require a stable feature key")


def _dual_record_status(simplex):
    raw = str(_read(simplex, "status", "")).upper()
    if raw in {status.value for status in Status}:
        return Status(raw)
    if "UNRESOLVED" in raw or "UNKNOWN" in raw:
        return Status.UNRESOLVED
    if _read(simplex, "supported", True) is False:
        return Status.NOT_SUPPORTED
    if _read(simplex, "complete", True) is False:
        return Status.PARTIAL
    return Status.CERTIFIED


def _dual_record(simplex, provenance):
    generator_ids = tuple(_read(simplex, "generator_ids", _read(simplex, "generator_tuple")))
    dimension = _read(simplex, "dimension")
    primal_feature_ids = _read(simplex, "primal_feature_ids", None)
    if primal_feature_ids is None:
        primal_feature_ids = tuple(
            _source_feature_id(feature)
            for feature in _read(simplex, "primal_features", ())
        )
    return DualSimplexRecord(
        _read(simplex, "simplex_id"),
        dimension,
        generator_ids,
        tuple(_source_feature_id(feature) for feature in primal_feature_ids),
        _read(simplex, "bounded", None),
        _read(simplex, "complete", None),
        _read(simplex, "supported", None),
        _dual_record_status(simplex),
        provenance,
    )


def _dual_records(dual, provenance, *, dimensions=None, simplex_ids=None):
    simplices = _read(dual, "simplices")
    if not isinstance(simplices, Mapping):
        raise TypeError("Dual complex must expose a dimension-indexed simplices mapping")
    wanted_dimensions = None if dimensions is None else frozenset(dimensions)
    wanted_ids = None if simplex_ids is None else frozenset(simplex_ids)
    if wanted_ids is not None and any(not isinstance(value, str) or not value for value in wanted_ids):
        raise ValueError("simplex_ids must contain non-empty canonical simplex IDs")
    records = []
    for dimension, by_generator_ids in simplices.items():
        if wanted_dimensions is not None and dimension not in wanted_dimensions:
            continue
        source_simplices = (
            by_generator_ids.values()
            if isinstance(by_generator_ids, Mapping)
            else by_generator_ids
        )
        for simplex in source_simplices:
            record = _dual_record(simplex, provenance)
            if wanted_ids is None or record.simplex_id in wanted_ids:
                records.append(record)
    records.sort(key=lambda record: (record.dimension, record.generator_ids))
    if wanted_ids is not None:
        found_ids = {record.simplex_id for record in records}
        missing_ids = wanted_ids - found_ids
        if missing_ids:
            raise ValueError(
                "Requested dual simplex IDs are absent from the supplied dual complex: "
                + ", ".join(sorted(missing_ids))
            )
    return tuple(records)


def _dual_audit(dual):
    audit = getattr(dual, "audit", None)
    return audit() if callable(audit) else None


def _dual_result_status(records, audit):
    if audit is None or not bool(_read(audit, "valid", False)):
        return Status.PARTIAL
    if any(record.status is not Status.CERTIFIED for record in records):
        return Status.PARTIAL
    return Status.CERTIFIED


def _dual_count_quantities(records, status, provenance, audit):
    counts = {}
    for record in records:
        name = f"simplex_count_dim_{record.dimension}"
        counts[name] = counts.get(name, 0) + 1
    quantities = {
        name: Quantity(
            name,
            count,
            "1",
            status,
            "dual-complex traversal",
            "selected dual simplices",
            provenance=provenance,
        )
        for name, count in counts.items()
    }
    if audit is not None:
        generator_count = _read(audit, "generator_count", None)
        if generator_count is not None:
            quantities["dual_generator_count"] = Quantity(
                "dual_generator_count",
                generator_count,
                "1",
                status,
                "dual-incidence audit",
                "dual generator namespace",
                provenance=provenance,
            )
        unmapped = _read(audit, "unmapped_feature_count", None)
        if unmapped is not None:
            quantities["dual_unmapped_primal_feature_count"] = Quantity(
                "dual_unmapped_primal_feature_count",
                unmapped,
                "1",
                status,
                "dual-incidence audit",
                "primal-to-dual incidence audit",
                provenance=provenance,
            )
        feature_counts = _read(audit, "feature_counts", None)
        if isinstance(feature_counts, Mapping) and unmapped is not None:
            total = sum(feature_counts.values())
            if total:
                support = total - unmapped
                if not isinstance(support, int) or support < 0:
                    raise ValueError("Dual incidence audit has invalid feature coverage")
                quantities["dual_incidence_coverage"] = Quantity(
                    "dual_incidence_coverage",
                    support / total,
                    "1",
                    status,
                    "dual-incidence audit",
                    "mapped primal-feature coverage",
                    support_count=support,
                    total_count=total,
                    provenance=provenance,
                )
    return quantities


def unavailable_dual_contact_complex(identity, *, status=Status.NOT_CALCULATED, provenance=None):
    """Return an explicitly unavailable dual-contact result without fake counts."""
    _require_representation(identity, {"dual"}, "Dual-contact")
    return DualContactComplexResult(
        identity,
        _dual_provenance(provenance, "dual contact complex unavailable"),
        _unavailable_status(status),
    )


def adapt_dual_contact_complex(dual, identity, *, provenance=None, simplex_ids=None):
    """Adapt dimension-one dual incidences without querying geometry or curvature.

    ``simplex_ids`` is an optional explicit contact subset.  It is deliberately
    an ID list rather than a group predicate: FORGE owns the side-selection
    policy, while this adapter only snapshots the selected incidence records.
    """
    _require_representation(identity, {"dual"}, "Dual-contact")
    provenance = _dual_provenance(provenance, "cached dual contact complex")
    if dual is None:
        return unavailable_dual_contact_complex(identity, provenance=provenance)
    records = _dual_records(
        dual, provenance, dimensions={1}, simplex_ids=simplex_ids
    )
    audit = _dual_audit(dual)
    status = _dual_result_status(records, audit)
    return DualContactComplexResult(
        identity,
        provenance,
        status,
        metrics=_dual_count_quantities(records, status, provenance, audit),
        contact_features=records,
    )


def unavailable_dual_generator_complex(identity, *, status=Status.NOT_CALCULATED, provenance=None):
    """Return an explicitly unavailable dual-generator result without fake counts."""
    _require_representation(
        identity, {"dual_a", "dual_b", "dual_generator_complex"}, "Dual-generator"
    )
    return DualGeneratorComplexResult(
        identity,
        _dual_provenance(provenance, "dual generator complex unavailable"),
        _unavailable_status(status),
    )


def adapt_dual_generator_complex(
    dual,
    identity,
    *,
    provenance=None,
    simplex_ids=None,
    topology=None,
):
    """Traverse a cached dual complex into typed generator-complex records.

    ``topology`` must already contain typed, producer-supplied quantities.  No
    topology, area, or curvature is inferred from simplex counts here.
    """
    _require_representation(
        identity, {"dual_a", "dual_b", "dual_generator_complex"}, "Dual-generator"
    )
    provenance = _dual_provenance(provenance, "cached dual generator complex")
    if dual is None:
        return unavailable_dual_generator_complex(identity, provenance=provenance)
    records = _dual_records(dual, provenance, simplex_ids=simplex_ids)
    audit = _dual_audit(dual)
    status = _dual_result_status(records, audit)
    simplex_counts = _dual_count_quantities(records, status, provenance, None)
    simplex_counts = {
        name: quantity
        for name, quantity in simplex_counts.items()
        if name.startswith("simplex_count_dim_")
    }
    generator_complex = GeneratorComplex(
        identity.side,
        {
            int(name.rsplit("_", 1)[1]): quantity
            for name, quantity in simplex_counts.items()
        },
        topology or {},
        status,
        provenance,
    )
    return DualGeneratorComplexResult(
        identity,
        provenance,
        status,
        metrics=_dual_count_quantities(records, status, provenance, audit),
        generator_complex=generator_complex,
        simplex_records=records,
    )


def _alpha_payload(source):
    if source is None:
        return None
    payload = _read(source, "alpha_selection", source)
    if not isinstance(payload, Mapping):
        raise TypeError("Alpha selection must be a mapping or an analysis with alpha_selection")
    return payload


def _alpha_status(payload):
    raw = str(payload.get("status", "")).lower()
    if raw in {"unresolved", "unknown"}:
        return Status.UNRESOLVED
    if raw in {"not_calculated", "not calculated"}:
        return Status.NOT_CALCULATED
    if raw in {"not_supported", "not supported"}:
        return Status.NOT_SUPPORTED
    return Status.CERTIFIED if payload.get("certified") is True else Status.PARTIAL


def _alpha_quantity(name, value, status, provenance, *, coverage=False):
    if status in {Status.UNRESOLVED, Status.NOT_CALCULATED, Status.NOT_SUPPORTED}:
        # Legacy unresolved caches use ``coverage: 0.0`` as a sentinel.  A
        # scientific Result must not reinterpret that sentinel as a measured
        # zero, so unavailable alpha quantities remain explicitly null.
        return Quantity(
            name,
            None,
            "1",
            status,
            "cached alpha selection",
            "alpha-selected dual subset",
            provenance=provenance,
        )
    if value is None:
        return Quantity(
            name,
            None,
            "1",
            Status.NOT_CALCULATED,
            "cached alpha selection",
            "alpha-selected dual subset",
            provenance=provenance,
        )
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
        raise ValueError(f"Alpha selection {name} must be a nonnegative numeric value")
    if not coverage and (not isinstance(value, int) or isinstance(value, bool)):
        raise ValueError(f"Alpha selection {name} must be an integer count")
    if coverage and value > 1:
        raise ValueError("Alpha selection coverage must be between zero and one")
    return Quantity(
        name,
        value,
        "1",
        status,
        "cached alpha selection",
        "alpha-selected dual subset",
        provenance=provenance,
    )


def _stable_pair_strings(values, name):
    if values is None:
        return ()
    pairs = []
    for pair in values:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2 or any(value is None for value in pair):
            raise ValueError(f"{name} must contain two-item stable-identifier pairs")
        pairs.append(tuple(str(value) for value in pair))
    return tuple(pairs)


def unavailable_alpha_selection(identity, *, status=Status.NOT_CALCULATED, provenance=None):
    """Return an explicitly unavailable alpha-selection result without zero counts."""
    _require_representation(identity, {"alpha_selection"}, "Alpha-selection")
    return AlphaSelectionResult(
        identity,
        _dual_provenance(provenance, "alpha selection unavailable"),
        _unavailable_status(status),
    )


def adapt_alpha_selection(source, identity, *, dual=None, provenance=None):
    """Adapt a cached alpha selection without treating it as primal geometry."""
    _require_representation(identity, {"alpha_selection"}, "Alpha-selection")
    provenance = _dual_provenance(provenance, "cached alpha selection")
    payload = _alpha_payload(source)
    if payload is None:
        return unavailable_alpha_selection(identity, provenance=provenance)
    status = _alpha_status(payload)
    aliases = {
        "eligible_pairs": "eligible_bicolor_pairs",
        "candidate_pairs": "candidate_bicolor_pairs",
        "selected_features": "selected_pairs",
        "rejected_features": "certified_rejected_pairs",
        "unresolved_features": "unresolved_pairs",
        "mapped_surfaces": "selected_physical_pairs",
        "excluded_pairs": "excluded_pairs",
        "coverage": "coverage",
    }
    selection = {
        name: _alpha_quantity(
            name, payload.get(source_name), status, provenance, coverage=name == "coverage"
        )
        for name, source_name in aliases.items()
    }
    generator_pairs = tuple(
        tuple(pair) for pair in (payload.get("selected_generator_pairs") or ())
    )
    selected_feature_ids = tuple(
        f"1:{first},{second}" for first, second in generator_pairs
    )
    selected_features = ()
    if dual is not None and selected_feature_ids:
        selected_features = _dual_records(
            dual,
            provenance,
            dimensions={1},
            simplex_ids=selected_feature_ids,
        )
    return AlphaSelectionResult(
        identity,
        provenance,
        status,
        selection=selection,
        selected_dual_feature_ids=selected_feature_ids,
        selected_generator_pairs=generator_pairs,
        selected_system_pairs=_stable_pair_strings(
            payload.get("selected_system_pairs"), "selected_system_pairs"
        ),
        selected_physical_system_pairs=_stable_pair_strings(
            payload.get("selected_physical_system_pairs"),
            "selected_physical_system_pairs",
        ),
        selected_dual_features=selected_features,
    )
