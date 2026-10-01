"""Combinatorial straight-line Apollonius dual of a solved VorPy network.

This module is an analysis layer: it reads the network's generator incidence
and never calls or changes the AW solver. Generator IDs are the network-local
``Network.balls`` index values used by the primal incidence tables.
"""

from dataclasses import dataclass, field
from itertools import combinations
from math import isfinite

import numpy as np


@dataclass(frozen=True)
class PrimalFeatureRef:
    """Reference to one primal network record supporting a dual simplex."""

    kind: str
    feature_id: int
    generator_ids: tuple
    bounded: bool
    complete: bool
    status: str = "supported"
    component_index: int | None = None
    generator_cells_complete: bool = True

    @property
    def key(self):
        suffix = "" if self.component_index is None else f"#component-{self.component_index}"
        return f"{self.kind}:{self.feature_id}{suffix}"


@dataclass
class ApolloniusSimplex:
    """One canonical abstract simplex and its straight geometric realization."""

    dimension: int
    generator_ids: tuple
    generator_coordinates: tuple
    generator_radii: tuple
    primal_features: list = field(default_factory=list)
    complete: bool = False
    bounded: bool = True
    status: str = "supported"
    alpha_birth: float | None = None
    alpha_status: str = "not_calculated"
    diagnostics: dict = field(default_factory=dict)

    @property
    def simplex_id(self):
        return f"{self.dimension}:" + ",".join(map(str, self.generator_ids))

    @property
    def primal_feature_ids(self):
        return tuple(feature.key for feature in self.primal_features)

    @property
    def num_primal_components(self):
        return len(self.primal_features)

    @property
    def geometry(self):
        """Return straight dual vertices in R3, one per generator."""
        return self.generator_coordinates


@dataclass(frozen=True)
class UnsupportedFeature:
    """A primal record whose generator cardinality is not supported as a simplex."""

    kind: str
    feature_id: int
    generator_ids: tuple
    reason: str


class ApolloniusComplex:
    """Extract canonical dual simplices from an already-built VorPy network.

    Surfaces with two generators map to edges, network edges with three
    generators map to triangles, and vertices with four generators map to
    tetrahedra. Repeated rows are retained as separate ``primal_features`` on
    the same canonical simplex. Higher/lower-cardinality records are recorded
    as unsupported rather than expanded into guessed subsets.
    """

    _FEATURE_TABLES = ((1, "surfs", 2), (2, "edges", 3), (3, "verts", 4))

    def __init__(self, network, name=None):
        self.network = network
        self.name = str(name or getattr(network, "group_name", None) or "Network")
        self.simplices = {0: {}, 1: {}, 2: {}, 3: {}}
        self.unsupported_features = []
        self._extract()

    @staticmethod
    def _ids(value):
        if value is None:
            return ()
        try:
            return tuple(int(item) for item in value)
        except (TypeError, ValueError):
            return ()

    @staticmethod
    def _row_value(row, key, default=None):
        try:
            value = row.get(key, default)
        except AttributeError:
            value = default
        return default if value is None else value

    def _generator_metadata(self):
        balls = getattr(self.network, "balls", None)
        if balls is None:
            raise ValueError("ApolloniusComplex requires Network.balls.")
        metadata = {}
        for generator_id, row in balls.iterrows():
            generator_id = int(generator_id)
            try:
                coordinates = tuple(float(value) for value in row["loc"])
                radius = float(row["rad"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Invalid generator geometry for ID {generator_id}.") from error
            if len(coordinates) != 3 or not all(map(isfinite, coordinates)) or not isfinite(radius):
                raise ValueError(f"Non-finite generator geometry for ID {generator_id}.")
            # Synthetic and third-party network objects may omit the cell
            # completion flag; in that case the analysis layer has no evidence
            # of incompleteness. Solved VorPy networks always provide it.
            complete = bool(self._row_value(row, "complete", True))
            metadata[generator_id] = (coordinates, radius, complete)
        return metadata

    def _add(self, dimension, generator_ids, feature):
        key = tuple(sorted(generator_ids))
        if len(set(key)) != dimension + 1:
            self.unsupported_features.append(UnsupportedFeature(
                feature.kind, feature.feature_id, tuple(generator_ids),
                "Repeated generator IDs do not define a simplex.",
            ))
            return
        coordinates, radii, complete = zip(*(self._generators[gid] for gid in key))
        simplex = self.simplices[dimension].get(key)
        if simplex is None:
            simplex = ApolloniusSimplex(
                dimension=dimension,
                generator_ids=key,
                generator_coordinates=tuple(coordinates),
                generator_radii=tuple(radii),
                complete=all(complete) and feature.complete and feature.bounded,
                bounded=feature.bounded,
                diagnostics={"atom_numbers": self._atom_numbers(key)},
            )
            if "is_boundary_generator" in self.network.balls.columns:
                simplex.diagnostics["boundary_generator_flags"] = tuple(
                    bool(self.network.balls.at[gid, "is_boundary_generator"])
                    for gid in key
                )
            self.simplices[dimension][key] = simplex
        else:
            simplex.complete = simplex.complete or (
                all(complete) and feature.complete and feature.bounded
            )
            simplex.bounded = simplex.bounded or feature.bounded
        simplex.primal_features.append(feature)

    def _atom_numbers(self, key):
        return tuple(int(self.network.balls.at[gid, "num"])
                     if "num" in self.network.balls.columns else gid for gid in key)

    def _extract(self):
        self._generators = self._generator_metadata()
        for generator_id, (coordinates, radius, complete) in self._generators.items():
            feature = PrimalFeatureRef("cell", generator_id, (generator_id,),
                                       bounded=complete, complete=complete,
                                       status="complete" if complete else "incomplete")
            self.simplices[0][(generator_id,)] = ApolloniusSimplex(
                0, (generator_id,), (coordinates,), (radius,), [feature], complete,
                bounded=complete, status="supported",
            )

        for dimension, table_name, expected_count in self._FEATURE_TABLES:
            table = getattr(self.network, table_name, None)
            if table is None:
                continue
            for feature_id, row in table.iterrows():
                generator_ids = self._ids(self._row_value(row, "balls", ()))
                if len(generator_ids) != expected_count or len(set(generator_ids)) != expected_count:
                    self.unsupported_features.append(UnsupportedFeature(
                        table_name[:-1], int(feature_id), generator_ids,
                        f"Expected {expected_count} distinct generators; found {len(generator_ids)}.",
                    ))
                    continue
                missing = sorted(set(generator_ids) - self._generators.keys())
                if missing:
                    self.unsupported_features.append(UnsupportedFeature(
                        table_name[:-1], int(feature_id), generator_ids,
                        f"References unknown generator IDs {missing}.",
                    ))
                    continue
                if table_name == "surfs":
                    boundary_edges = self._ids(self._row_value(row, "edges", ()))
                    boundary_verts = self._ids(self._row_value(row, "verts", ()))
                    bounded = bool(boundary_edges and boundary_verts)
                elif table_name == "edges":
                    endpoint_ids = self._ids(self._row_value(row, "verts", ()))
                    bounded = len(endpoint_ids) >= 2 and all(
                        index in self.network.verts.index
                        and self._valid_point(self.network.verts.loc[index].get("loc"))
                        for index in endpoint_ids
                    )
                else:
                    bounded = self._valid_point(self._row_value(row, "loc"))
                generator_cells_complete = all(self._generators[gid][2] for gid in generator_ids)
                # Whole-cell completeness and primal-feature completeness
                # are different: an edge/vertex may be finite even when its
                # participating 3-D AW cells are unbounded.
                complete = bounded
                feature = PrimalFeatureRef(
                    table_name[:-1], int(feature_id), tuple(generator_ids),
                    bounded=bounded, complete=complete,
                    status=("unbounded_primal_feature" if not bounded else
                            "supported" if generator_cells_complete else
                            "incomplete_network_geometry"),
                    generator_cells_complete=generator_cells_complete,
                )
                component_count = (
                    self._surface_component_count(row)
                    if table_name == "surfs" else 1
                )
                for component_index in range(component_count):
                    component_feature = feature
                    if component_count > 1:
                        component_feature = PrimalFeatureRef(
                            kind=feature.kind,
                            feature_id=feature.feature_id,
                            generator_ids=feature.generator_ids,
                            bounded=feature.bounded,
                            complete=feature.complete,
                            status=feature.status,
                            component_index=component_index,
                            generator_cells_complete=feature.generator_cells_complete,
                        )
                    self._add(dimension, generator_ids, component_feature)
        del self._generators

    @staticmethod
    def _surface_component_count(row):
        """Count disconnected triangle components in one stored surface patch."""
        try:
            triangles = np.asarray(row.get("tris", ()), dtype=np.int64)
        except (TypeError, ValueError):
            return 1
        if triangles.size == 0:
            return 1
        if triangles.ndim != 2 or triangles.shape[1] != 3:
            return 1
        parents = list(range(len(triangles)))

        def find(index):
            while parents[index] != index:
                parents[index] = parents[parents[index]]
                index = parents[index]
            return index

        def union(first, second):
            first_root, second_root = find(first), find(second)
            if first_root != second_root:
                parents[second_root] = first_root

        owner = {}
        for triangle_id, triangle in enumerate(triangles):
            for vertex_id in set(map(int, triangle)):
                if vertex_id in owner:
                    union(triangle_id, owner[vertex_id])
                else:
                    owner[vertex_id] = triangle_id
        return len({find(index) for index in range(len(triangles))})

    @staticmethod
    def _valid_point(point):
        try:
            values = np.asarray(point, dtype=float)
            return values.shape == (3,) and bool(np.isfinite(values).all())
        except (TypeError, ValueError):
            return False

    def _supported(self, simplex):
        if simplex.dimension == 0:
            # The generator is the dual vertex even when its whole primal
            # cell is unbounded; only positive-dimensional incidence needs
            # bounded feature geometry.
            return True
        for feature in simplex.primal_features:
            if not (feature.bounded and feature.complete):
                continue
            if simplex.dimension == 1:
                return True
            if all(self._face_is_supported(simplex, tuple(face_ids))
                   for face_ids in combinations(simplex.generator_ids, simplex.dimension)):
                return True
        return False

    def _face_is_supported(self, simplex, face_ids):
        face = self.simplices[len(face_ids) - 1].get(face_ids)
        return (face is not None and self._supported(face)
                and self._primal_feature_claims_face(simplex, face_ids))

    def is_supported(self, dimension, generator_ids):
        """Return whether this canonical simplex is in the supported subcomplex."""
        key = tuple(sorted(int(value) for value in generator_ids))
        simplex = self.simplices[int(dimension)].get(key)
        return bool(simplex is not None and self._supported(simplex))

    def _primal_feature_claims_face(self, simplex, face_ids):
        """Check the stored boundary references on the actual primal record."""
        for feature in simplex.primal_features:
            if feature.kind == "edge" and len(face_ids) == 2:
                row = self.network.edges.loc[feature.feature_id]
                if "surfs" not in row:
                    return face_ids in self.simplices[1]
                for surf_id in self._ids(self._row_value(row, "surfs", ())):
                    if surf_id in self.network.surfs.index:
                        if tuple(sorted(self._ids(self.network.surfs.loc[surf_id, "balls"]))) == face_ids:
                            return True
            elif feature.kind == "vert":
                row = self.network.verts.loc[feature.feature_id]
                link_key = "edges" if len(face_ids) == 3 else "surfs"
                if link_key not in row:
                    return face_ids in self.simplices[len(face_ids) - 1]
                table = self.network.edges if link_key == "edges" else self.network.surfs
                for child_id in self._ids(self._row_value(row, link_key, ())):
                    if child_id in table.index:
                        child_ids = self._ids(table.loc[child_id, "balls"])
                        if tuple(sorted(child_ids)) == face_ids:
                            return True
        return False

    def incidence_diagnostics(self):
        """Classify closure gaps against the actual primal records.

        Missing faces are never synthesized. Only closure relations whose
        coface has bounded geometry and complete participating cells count as
        supported-incidence failures.
        """
        causes = {
            "truly_unbounded_primal_feature": [],
            "incomplete_network_geometry": [],
            "unsupported_higher_order_incidence": list(self.unsupported_features),
            "missing_primal_feature_despite_valid_incidence": [],
            "canonical_id_mismatch": [],
            "duplicate_multiple_primal_components": [],
            "implementation_error": [],
        }
        for table in self.simplices.values():
            for simplex in table.values():
                if simplex.num_primal_components > 1:
                    causes["duplicate_multiple_primal_components"].append(simplex.simplex_id)
                for feature in simplex.primal_features:
                    if not feature.bounded:
                        causes["truly_unbounded_primal_feature"].append(feature.key)
                    elif not feature.generator_cells_complete:
                        causes["incomplete_network_geometry"].append(feature.key)
                for gid in simplex.generator_ids:
                    if gid not in self.network.balls.index:
                        causes["canonical_id_mismatch"].append((simplex.simplex_id, gid))

        missing = []
        missing_face_causes = {
            "truly_unbounded_primal_feature": [],
            "incomplete_network_geometry": [],
            "unsupported_higher_order_incidence": [],
            "missing_primal_feature_despite_valid_incidence": [],
            "implementation_error": [],
        }
        for dim in (2, 3):
            for simplex in self.simplices[dim].values():
                face_dimensions = (dim - 1, 1) if dim == 3 else (dim - 1,)
                for face_dimension in face_dimensions:
                    for face_ids in combinations(simplex.generator_ids, face_dimension + 1):
                        face_ids = tuple(face_ids)
                        if face_ids in self.simplices[face_dimension]:
                            continue
                        relation = {
                            "simplex": simplex.simplex_id,
                            "missing_face": face_ids,
                            "face_dimension": face_dimension,
                        }
                        missing.append(relation)
                        if any(not feature.bounded for feature in simplex.primal_features):
                            missing_face_causes["truly_unbounded_primal_feature"].append(relation)
                        elif any(not feature.generator_cells_complete for feature in simplex.primal_features):
                            missing_face_causes["incomplete_network_geometry"].append(relation)
                        elif self._primal_feature_claims_face(simplex, face_ids):
                            missing_face_causes["missing_primal_feature_despite_valid_incidence"].append(relation)
                            missing_face_causes["implementation_error"].append(relation)
                        else:
                            missing_face_causes["unsupported_higher_order_incidence"].append(relation)
        return {
            "causes": causes,
            "missing_face_relations": missing,
            "missing_face_causes": missing_face_causes,
            "supported_counts": {
                dim: sum(self._supported(simplex) for simplex in table.values())
                for dim, table in self.simplices.items()
            },
        }

    def validate(self):
        """Validate simplicial closure of the supported bounded subcomplex."""
        errors = []
        valid_ids = set(int(index) for index in self.network.balls.index)
        for dimension, table in self.simplices.items():
            for key, simplex in table.items():
                if key != simplex.generator_ids or len(key) != dimension + 1:
                    errors.append(f"Noncanonical dimension-{dimension} simplex {key}.")
                if not set(key).issubset(valid_ids):
                    errors.append(f"Simplex {simplex.simplex_id} has invalid generator IDs.")
                if dimension and not simplex.primal_features:
                    errors.append(f"Simplex {simplex.simplex_id} has no primal support.")
        for triangle in self.simplices[2].values():
            if not self._supported(triangle):
                continue
            for edge in combinations(triangle.generator_ids, 2):
                edge = tuple(edge)
                if not self._face_is_supported(triangle, edge):
                    errors.append(f"Triangle {triangle.simplex_id} lacks edge {edge}.")
        for tetrahedron in self.simplices[3].values():
            if not self._supported(tetrahedron):
                continue
            ids = tetrahedron.generator_ids
            for face in combinations(ids, 3):
                if not self._face_is_supported(tetrahedron, face):
                    errors.append(f"Tetrahedron {tetrahedron.simplex_id} lacks triangle {face}.")
            for edge in combinations(ids, 2):
                if not self._face_is_supported(tetrahedron, edge):
                    errors.append(f"Tetrahedron {tetrahedron.simplex_id} lacks edge {edge}.")
            for vertex in combinations(ids, 1):
                if vertex not in self.simplices[0]:
                    errors.append(f"Tetrahedron {tetrahedron.simplex_id} lacks vertex {vertex}.")
        return errors

    def assert_valid(self):
        errors = self.validate()
        if errors:
            raise ValueError("Invalid Apollonius complex:\n" + "\n".join(errors))
        return True

    @property
    def euler_characteristic(self):
        """Euler characteristic of the supported, incidence-closed subcomplex."""
        return sum((-1) ** dimension * sum(self._supported(simplex) for simplex in table.values())
                   for dimension, table in self.simplices.items())

    @property
    def raw_euler_characteristic(self):
        """Alternating count of all retained primal generator tuples."""
        return sum((-1) ** dimension * len(table)
                   for dimension, table in self.simplices.items())

    @property
    def supported_counts(self):
        return {dim: sum(self._supported(simplex) for simplex in table.values())
                for dim, table in self.simplices.items()}

    @property
    def multiple_primal_simplices(self):
        return [simplex for table in self.simplices.values()
                for simplex in table.values() if simplex.num_primal_components > 1]

    @property
    def incomplete_features(self):
        return [feature for table in self.simplices.values()
                for simplex in table.values() for feature in simplex.primal_features
                if (not feature.complete or not feature.bounded
                    or not feature.generator_cells_complete)]

    @property
    def degenerate_feature_count(self):
        return len(self.unsupported_features)

    def summary(self):
        errors = self.validate()
        diagnostics = self.incidence_diagnostics()
        cause_counts = {key: len(value) for key, value in diagnostics["causes"].items()}
        summary = (
            f"Apollonius Complex: {self.name}\n\n"
            f"Generators:                 {len(self.simplices[0]):,}\n"
            f"Dual edges:                 {len(self.simplices[1]):,}\n"
            f"Dual triangles:             {len(self.simplices[2]):,}\n"
            f"Dual tetrahedra:            {len(self.simplices[3]):,}\n\n"
            f"Degenerate/unsupported:     {self.degenerate_feature_count:,}\n"
            f"Unbounded/incomplete:       {len(self.incomplete_features):,}\n"
            f"Multiple-primal simplices:  {len(self.multiple_primal_simplices):,}\n"
            f"Euler characteristic:       {self.euler_characteristic}\n"
            f"Incidence validation:        {'passed' if not errors else f'{len(errors)} issue(s)'}\n"
            f"Supported simplices V/E/Tr/Tet: {tuple(diagnostics['supported_counts'][d] for d in range(4))}\n"
            f"Missing-face relations:       {len(diagnostics['missing_face_relations']):,}\n"
            f"Missing-face causes:          "
            f"unbounded={len(diagnostics['missing_face_causes']['truly_unbounded_primal_feature']):,}, "
            f"incomplete={len(diagnostics['missing_face_causes']['incomplete_network_geometry']):,}, "
            f"unsupported={len(diagnostics['missing_face_causes']['unsupported_higher_order_incidence']):,}, "
            f"missing-primal={len(diagnostics['missing_face_causes']['missing_primal_feature_despite_valid_incidence']):,}, "
            f"implementation={len(diagnostics['missing_face_causes']['implementation_error']):,}\n"
            f"Feature flags: unbounded={cause_counts['truly_unbounded_primal_feature']:,}, "
            f"incomplete={cause_counts['incomplete_network_geometry']:,}, "
            f"higher-order unsupported={cause_counts['unsupported_higher_order_incidence']:,}, "
            f"ID-mismatch={cause_counts['canonical_id_mismatch']:,}, "
            f"multiple-primal={cause_counts['duplicate_multiple_primal_components']:,}\n"
        )
        if hasattr(self, "alpha_unresolved"):
            resolved = sum(
                simplex.alpha_birth is not None
                for table in self.simplices.values() for simplex in table.values()
            )
            filtration_issues = self.validate_filtration()
            unresolved_comparisons = sum(issue["kind"] == "unresolved_face"
                                         for issue in filtration_issues)
            monotonicity_issues = sum(issue["kind"] == "monotonicity_violation"
                                      for issue in filtration_issues)
            summary += (
                f"\nAlpha births resolved:       {resolved:,}\n"
                f"Alpha births unresolved:    {len(self.alpha_unresolved):,}\n"
                f"Filtration violations:       {monotonicity_issues:,}\n"
                f"Unresolved face comparisons:{unresolved_comparisons:>8,}\n"
            )
        return summary

    def calculate_alpha_births(self, tolerance=1e-6):
        """Calculate exact supported AW filtration births in length units.

        A generator's cell birth is ``-r_i`` only when its center is in its
        closed AW cell; otherwise its constrained cell minimum is unresolved.
        A dual triangle uses the analytic AW trisector parameter (common
        clearance ``rho``) over each bounded primal edge branch. A dual
        tetrahedron uses the common clearance at each corresponding Voronoi
        vertex, checked independently against every generator. Minima over
        trimmed pairwise surfaces are left unresolved because the stored
        triangulated patches do not define their exact analytic domain.
        """
        tolerance = float(tolerance)
        if not isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("Alpha tolerance must be finite and positive.")
        settings = getattr(self.network, "settings", None) or {}
        if settings.get("net_type", "aw") != "aw":
            raise ValueError("AW alpha filtration requires an additively weighted network.")
        self.alpha_tolerance = tolerance
        self.alpha_unresolved = []

        balls = self.network.balls
        locs = {int(index): np.asarray(row["loc"], dtype=float)
                for index, row in balls.iterrows()}
        radii = {int(index): float(row["rad"])
                 for index, row in balls.iterrows()}
        all_locs = np.asarray([locs[index] for index in locs], dtype=float)
        all_radii = np.asarray([radii[index] for index in locs], dtype=float)

        # A cell contains its generator center iff its weighted distance there
        # is no greater than that of every other generator. The global lower
        # bound -r_i is then attained and is the exact cell minimum.
        for generator_key, simplex in self.simplices[0].items():
            generator_id = generator_key[0]
            center = locs[generator_id]
            clearances = np.linalg.norm(all_locs - center, axis=1) - all_radii
            owner_clearance = -radii[generator_id]
            if bool(np.all(owner_clearance <= clearances + tolerance)):
                self._set_birth(simplex, owner_clearance, "exact_cell_center",
                                {"center_in_cell": True})
            else:
                self._unresolve(simplex, "cell_minimum_unresolved",
                                "Generator center is outside its AW cell; exact constrained minimum not implemented.")

        # Pairwise surface minima require minimization over the analytically
        # trimmed 2D patch, including its boundary. Do not substitute samples.
        for simplex in self.simplices[1].values():
            self._unresolve(simplex, "surface_minimum_unresolved",
                            "Exact minimum over the trimmed AW surface patch is not currently available.")

        for simplex in self.simplices[2].values():
            values = []
            diagnostics = []
            for feature in simplex.primal_features:
                try:
                    value, diagnostic = _edge_alpha_minimum(
                        self.network, feature.feature_id, tolerance
                    )
                except (TypeError, ValueError, KeyError, IndexError, np.linalg.LinAlgError) as error:
                    values = []
                    diagnostics.append({"feature": feature.key, "reason": str(error)})
                    break
                values.append(value)
                diagnostics.append({"feature": feature.key, **diagnostic})
            if values:
                self._set_birth(simplex, min(values), "exact_aw_edge_minimum",
                                {"components": diagnostics, "minimized_over": "all primal edge components"})
            else:
                reason = diagnostics[0].get("reason") if diagnostics else "No primal edge components."
                self._unresolve(simplex, "edge_minimum_unresolved", reason, {"components": diagnostics})

        for simplex in self.simplices[3].values():
            if not self._supported(simplex):
                self._unresolve(simplex, "unsupported_primal_incidence",
                                "The bounded network incidence is incomplete or lacks a required lower-dimensional primal face.")
                continue
            values = []
            diagnostics = []
            for feature in simplex.primal_features:
                row = self.network.verts.loc[feature.feature_id]
                point = np.asarray(row["loc"], dtype=float)
                clearances = [float(np.linalg.norm(point - locs[generator_id]) - radii[generator_id])
                              for generator_id in simplex.generator_ids]
                spread = max(clearances) - min(clearances)
                if not self._valid_point(point) or spread > tolerance:
                    values = []
                    diagnostics.append({
                        "feature": feature.key,
                        "status": "unequal_generator_clearances",
                        "clearances": clearances,
                        "spread": spread,
                    })
                    break
                value = float(np.mean(clearances))
                stored_radius = row.get("rad")
                stored_error = (None if stored_radius is None
                                else abs(value - float(stored_radius)))
                diagnostics.append({
                    "feature": feature.key,
                    "status": "exact_common_vertex_clearance",
                    "clearances": clearances,
                    "spread": spread,
                    "stored_vertex_radius_error": stored_error,
                })
                values.append(value)
            if values:
                self._set_birth(simplex, min(values), "exact_aw_vertex_clearance",
                                {"components": diagnostics, "minimized_over": "all primal vertex components"})
            else:
                self._unresolve(simplex, "vertex_clearance_unresolved",
                                "Primal vertex does not have a verified common AW clearance.",
                                {"components": diagnostics})

        violations = self.validate_filtration(tolerance=tolerance)
        substantive = [issue for issue in violations if issue["kind"] == "monotonicity_violation"]
        if substantive:
            details = "; ".join(issue["message"] for issue in substantive[:5])
            raise ValueError(f"AW alpha birth values violate filtration monotonicity: {details}")
        return violations

    def _set_birth(self, simplex, value, status, diagnostics):
        simplex.alpha_birth = float(value)
        simplex.alpha_status = status
        simplex.diagnostics["alpha"] = diagnostics

    def _unresolve(self, simplex, status, reason, diagnostics=None):
        simplex.alpha_birth = None
        simplex.alpha_status = status
        simplex.diagnostics["alpha"] = {"reason": str(reason), **(diagnostics or {})}
        self.alpha_unresolved.append(simplex.simplex_id)

    def validate_filtration(self, tolerance=None):
        """Report unresolved face comparisons and substantive monotonicity errors."""
        tolerance = float(tolerance if tolerance is not None
                          else getattr(self, "alpha_tolerance", 1e-6))
        issues = []
        for dimension in range(1, 4):
            for simplex in self.simplices[dimension].values():
                if simplex.alpha_birth is None or not self._supported(simplex):
                    continue
                for face_ids in combinations(simplex.generator_ids, dimension):
                    face = self.simplices[dimension - 1].get(tuple(face_ids))
                    if face is None:
                        issues.append({
                            "kind": "missing_face",
                            "simplex": simplex.simplex_id,
                            "face": face_ids,
                            "message": f"{simplex.simplex_id} lacks face {face_ids}.",
                        })
                    elif face.alpha_birth is None:
                        issues.append({
                            "kind": "unresolved_face",
                            "simplex": simplex.simplex_id,
                            "face": face.simplex_id,
                            "message": f"Cannot verify filtration: {face.simplex_id} birth is unresolved.",
                        })
                    elif face.alpha_birth > simplex.alpha_birth + tolerance:
                        issues.append({
                            "kind": "monotonicity_violation",
                            "simplex": simplex.simplex_id,
                            "face": face.simplex_id,
                            "message": (f"Face {face.simplex_id} birth {face.alpha_birth:.9g} exceeds "
                                        f"coface {simplex.simplex_id} birth {simplex.alpha_birth:.9g}."),
                        })
        return issues

    def alpha_complex(self, alpha, tolerance=None):
        """Extract the birth-filtered subcomplex, omitting unresolved cofaces."""
        if not hasattr(self, "alpha_tolerance"):
            self.calculate_alpha_births()
        alpha = float(alpha)
        if not isfinite(alpha):
            raise ValueError("Alpha must be finite.")
        tolerance = float(tolerance if tolerance is not None else self.alpha_tolerance)
        included = {0: {}, 1: {}, 2: {}, 3: {}}
        blocked = []
        for dimension in range(4):
            for key, simplex in self.simplices[dimension].items():
                if simplex.alpha_birth is None or simplex.alpha_birth > alpha + tolerance:
                    continue
                if dimension and not self._supported(simplex):
                    blocked.append({"simplex": simplex.simplex_id,
                                    "missing_faces": (),
                                    "reason": "primal incidence is unsupported or incomplete"})
                    continue
                faces = [tuple(face) for face in combinations(key, dimension)] if dimension else []
                missing = [face for face in faces if face not in included[dimension - 1]]
                if missing:
                    blocked.append({"simplex": simplex.simplex_id,
                                    "missing_faces": tuple(missing),
                                    "reason": "faces unresolved, above alpha, or absent"})
                    continue
                included[dimension][key] = simplex
        return AlphaComplex(self, alpha, included, blocked, tolerance)


def _edge_alpha_minimum(network, edge_id, tolerance):
    """Minimize AW clearance over one bounded, analytically resolved edge."""
    from vorpy.src.calculations.edge_geometry import (
        AdditivelyWeightedTrisectorBranch,
        AdditivelyWeightedTrisectorConic,
        AffineEdgeGeometry,
        LineEdgeGeometry,
    )
    from vorpy.src.network.edge_geometry_diagnostics import resolve_aw_network_edge

    resolved = resolve_aw_network_edge(network, edge_id, tolerance=tolerance)
    geometry = resolved.geometry
    parameter_map = None
    if isinstance(geometry, AffineEdgeGeometry):
        wrapper = geometry
        geometry = wrapper.source
        source_lower = float(wrapper.source_start)
        source_upper = float(wrapper.source_end)
        parameter_map = lambda value: (float(value) - source_lower) / (source_upper - source_lower)
        lower = min(source_lower, source_upper)
        upper = max(source_lower, source_upper)
    else:
        lower, upper = float(geometry.t_min), float(geometry.t_max)
    if not isfinite(lower) or not isfinite(upper) or upper < lower:
        raise ValueError("Analytic AW edge has invalid or unbounded parameter bounds.")

    if isinstance(geometry, AdditivelyWeightedTrisectorBranch):
        candidates = [lower, upper]
        parameterization = "rho_parameter"
        clearance = lambda value: float(value)
    elif isinstance(geometry, AdditivelyWeightedTrisectorConic):
        candidates = [lower, upper]
        if geometry.conic_type == "ellipse":
            first = int(np.ceil(lower / np.pi))
            last = int(np.floor(upper / np.pi))
            candidates.extend(float(index * np.pi) for index in range(first, last + 1))
        elif geometry.conic_type == "hyperbola" and geometry.hyperbola_orientation == "rho":
            if lower <= 0.0 <= upper:
                candidates.append(0.0)
        elif geometry.conic_type == "parabola" and lower <= 0.0 <= upper:
            candidates.append(0.0)
        parameterization = f"{geometry.conic_type}_rho_extrema"
        clearance = geometry.rho
    elif isinstance(geometry, LineEdgeGeometry):
        candidates = [lower, upper, 0.5 * (lower + upper)]
        parameterization = "equal_radius_straight_edge"
        clearance = None
    else:
        raise ValueError(f"Unsupported analytic edge geometry {type(geometry).__name__}.")

    values = []
    for parameter in candidates:
        point = np.asarray(geometry.point(parameter), dtype=float)
        generator_clearances = np.linalg.norm(
            point - np.asarray(resolved.locations, dtype=float), axis=1
        ) - np.asarray(resolved.radii, dtype=float)
        spread = float(np.max(generator_clearances) - np.min(generator_clearances))
        if spread > tolerance:
            raise ValueError(f"AW edge generators disagree in clearance by {spread:.6g} Å.")
        if clearance is None:
            values.append(float(np.mean(generator_clearances)))
        else:
            expected = clearance(parameter)
            if float(np.max(np.abs(generator_clearances - expected))) > tolerance:
                raise ValueError("Analytic edge parameter disagrees with generator clearance.")
            values.append(float(expected))
    return min(values), {
        "status": resolved.status,
        "parameterization": parameterization,
        "parameter_bounds": (lower, upper),
        "candidate_clearances": tuple(values),
        "affine_wrapper": parameter_map is not None,
    }


class AlphaComplex:
    """Subcomplex induced by supported AW births no greater than ``alpha``."""

    def __init__(self, source, alpha, simplices, blocked_simplices, tolerance):
        self.source = source
        self.alpha = float(alpha)
        self.simplices = simplices
        self.blocked_simplices = list(blocked_simplices)
        self.tolerance = float(tolerance)

    @property
    def counts(self):
        return {dimension: len(table) for dimension, table in self.simplices.items()}

    def validate(self):
        errors = []
        for dimension in range(1, 4):
            for simplex in self.simplices[dimension].values():
                for face in combinations(simplex.generator_ids, dimension):
                    if tuple(face) not in self.simplices[dimension - 1]:
                        errors.append(f"{simplex.simplex_id} is missing face {face}.")
        return errors

    def summary(self):
        counts = self.counts
        return (
            f"AW Alpha Complex: {self.source.name}\n"
            f"Alpha (Å): {self.alpha:.6g}\n"
            f"Effective radius: r_i + alpha\n"
            f"Vertices: {counts[0]:,}\nEdges: {counts[1]:,}\n"
            f"Triangles: {counts[2]:,}\nTetrahedra: {counts[3]:,}\n"
            f"Blocked by unavailable faces: {len(self.blocked_simplices):,}\n"
            f"Closure validation: {'passed' if not self.validate() else 'failed'}\n"
        )
