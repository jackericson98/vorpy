"""Closed, versioned schemas for cached scientific records.

Load allocates these explicitly allowed classes without calling constructors.
No module/class names from an archive are imported or executed.
"""
from dataclasses import fields
from uuid import uuid4


def stable_id(obj):
    value = getattr(obj, 'archive_id', None)
    if value is None:
        value = str(uuid4())
        object.__setattr__(obj, 'archive_id', value)
    return value


_ENVIRONMENT_ALIASES = {
    'dry': 'dry',
    'wet': 'solvent_competing',
    'solvent': 'solvent_competing',
    'solvent-competing': 'solvent_competing',
    'solvent_competing': 'solvent_competing',
    'solvent competing': 'solvent_competing',
}
_NETWORK_SCOPE_ALIASES = {
    'system': 'system',
    'system network': 'system',
    'whole-system': 'system',
    'whole_system': 'system',
    'whole system': 'system',
    'interface': 'interface',
    'interface/dedicated': 'interface',
    'dedicated': 'interface',
    'dedicated interface': 'interface',
    'interface network': 'interface',
}


def canonical_environment(value):
    """Map old solvent labels to Agent #1's session vocabulary."""
    if value is None:
        return 'unknown'
    text = str(value).strip().lower()
    return _ENVIRONMENT_ALIASES.get(text, text or 'unknown')


def canonical_network_scope(value):
    """Map old network-scope labels without rejecting future labels."""
    if value is None:
        return 'unknown'
    text = str(value).strip().lower()
    return _NETWORK_SCOPE_ALIASES.get(text, text or 'unknown')


def normalize_vocabulary(obj):
    """Normalize persisted labels; leave scientific values untouched."""
    if hasattr(obj, 'network_scope'):
        object.__setattr__(obj, 'network_scope', canonical_network_scope(obj.network_scope))
    for name in ('metadata', 'archive_provenance'):
        value = getattr(obj, name, None)
        if not isinstance(value, dict):
            continue
        if 'environment' in value:
            value['environment'] = canonical_environment(value['environment'])
        if 'network_scope' in value:
            value['network_scope'] = canonical_network_scope(value['network_scope'])


def registry():
    from vorpy.src.interface.geometry_analysis import InterfaceGeometryAnalysis
    from vorpy.src.geometry.duals.model import DualComplex, DualSimplex, DualIncidenceAudit
    from vorpy.src.geometry.filtrations.alpha import AlphaFiltration, AlphaSimplex, AlphaSubcomplex
    from vorpy.src.geometry.interfaces.model import (
        InterfaceRepresentation, InterfaceSurface, InterfaceEdge, InterfaceVertex,
    )
    from vorpy.src.geometry.interfaces.curvature import CurvatureSummary
    from vorpy.src.analyze.apollonius import (
        PrimalFeatureRef, UnsupportedFeature, ApolloniusSimplex, ApolloniusComplex,
    )
    from vorpy.src.geometry.aw_alpha.complex import AWAlphaRecord, AWAlphaFiltration, AWAlphaSubcomplex
    from vorpy.src.geometry.molecular_contact import (
        MolecularContactSurface, SphericalPatch, CircularArc, JunctionVertex,
        MolecularContactSurfaceCurvature, RefinedSeamPiece, TopologicalCut,
    )
    from vorpy.src.boundary import BoundaryConfig, BoundaryGenerator
    from vorpy.workbench.domain import AnalysisResult, Atom, Bond, GeometryLayer
    from vorpy.workbench.project import Project, StructureSource, GroupDefinition, InterfaceDefinition, AtomKey
    from vorpy.workbench.services.info_parser import Measurement, ChainRecord, ResidueRecord, NetworkSummary
    dataclasses = {
        'interface_analysis': InterfaceGeometryAnalysis,
        'dual_simplex': DualSimplex, 'dual_incidence_audit': DualIncidenceAudit,
        'alpha_simplex': AlphaSimplex,
        'physical_interface': InterfaceRepresentation,
        'physical_surface': InterfaceSurface, 'physical_edge': InterfaceEdge,
        'physical_vertex': InterfaceVertex, 'curvature_summary': CurvatureSummary,
        'primal_reference': PrimalFeatureRef, 'unsupported_feature': UnsupportedFeature,
        'apollonius_simplex': ApolloniusSimplex, 'aw_alpha_record': AWAlphaRecord,
        'boundary_config': BoundaryConfig, 'boundary_generator': BoundaryGenerator,
        'workbench_result': AnalysisResult, 'display_atom': Atom, 'display_bond': Bond,
        'display_layer': GeometryLayer, 'workbench_project': Project,
        'structure_source': StructureSource, 'group_definition': GroupDefinition,
        'interface_definition': InterfaceDefinition, 'atom_key': AtomKey,
        'measurement': Measurement, 'chain_summary': ChainRecord,
        'residue_summary': ResidueRecord, 'network_summary': NetworkSummary,
        'molecular_contact_surface': MolecularContactSurface,
        'molecular_contact_surface_curvature': MolecularContactSurfaceCurvature,
        'molecular_spherical_patch': SphericalPatch,
        'molecular_circular_arc': CircularArc,
        'molecular_junction_vertex': JunctionVertex,
        'molecular_refined_seam_piece': RefinedSeamPiece,
        'molecular_topological_cut': TopologicalCut,
    }
    result = {kind: (cls, tuple(field.name for field in fields(cls)))
              for kind, cls in dataclasses.items()}
    result.update({
        'dual_complex': (DualComplex, ('network', 'scheme', 'generators', 'simplices',
            'feature_to_simplices', 'unresolved_features', 'metadata')),
        'alpha_filtration': (AlphaFiltration, ('network', 'scheme', 'records', 'metadata',
            'tolerance', 'blocked', '_aw_source')),
        'alpha_subcomplex': (AlphaSubcomplex, ('scheme', 'alpha', 'records', 'blocked', 'metadata')),
        'aw_alpha_filtration': (AWAlphaFiltration, ('network', 'name', 'tolerance',
            'incidence', 'records', 'issues', '_built', '_max_dimension_built')),
        'aw_alpha_subcomplex': (AWAlphaSubcomplex, ('source', 'alpha', 'records', 'blocked',
            'tolerance', 'truncate_negative')),
        'apollonius_complex': (ApolloniusComplex, ('network', 'name', 'simplices', 'unsupported_features')),
    })
    return result


def scope(network):
    return ('interface' if getattr(network, 'iface_grps', None) is not None
            or getattr(network, 'network_mode', '') == 'interface' else 'system')


def restore_runtime_keys(iface):
    """Rebind disposable cache keys to restored objects, never recompute results."""
    analysis = getattr(iface, 'geometry_analysis', None)
    if analysis is None or getattr(iface, 'archive_cache_state', None) == 'STALE':
        iface._geometry_analysis_key = iface._geometry_source_key = None
        return
    from vorpy.src.interface.geometry_analysis import geometry_cache_key
    curvature_requested = bool(analysis.metadata.get('curvature_requested', True))
    key = (*geometry_cache_key(iface, analysis.metadata.get('alpha_value', 0.)), curvature_requested)
    iface._geometry_analysis_key = key
    iface._geometry_source_key = (key[0:4], key[5], curvature_requested)
