import json
from dataclasses import asdict

from vorpy.src.geometry.molecular_contact import (
    CircularArc,
    JunctionVertex,
    MolecularContactSurface,
    RefinedSeamPiece,
    SphericalPatch,
    TopologicalCut,
)
from vorpy.src.geometry.visualization import export_molecular_contact_bundle
from vorpy.src.geometry.visualization.molecular_contact import export_molecular_contact_surface


def _surface(side):
    carrier = f'carrier-{side}'
    points = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-1.0, 0.0, 0.0))
    junctions = tuple(
        JunctionVertex(f'junction-{index}', carrier, point, point,
                       (f'circle-{index}',), (f'arc-{index}',), 'CERTIFIED')
        for index, point in enumerate(points)
    )
    spans = (1.5707963267948966, 1.5707963267948966, 3.141592653589793)
    arcs = tuple(
        CircularArc(
            f'arc-{index}', carrier, f'circle-{index}',
            f'junction-{index}', f'junction-{(index + 1) % 3}',
            points[index], points[(index + 1) % 3], (0.0, 0.0, 1.0), 0.0,
            spans[index], (carrier, f'partner-{side}'), (f'contact-{side}',), 1,
        )
        for index in range(3)
    )
    patch = SphericalPatch(
        f'patch-{side}', carrier, (0.0, 0.0, 0.0), 1.0, side,
        (f'partner-{side}',), (), (f'contact-{side}',),
        'source-network', 'source-interface', 'mcs-v1',
        tuple(arc.arc_id for arc in arcs), tuple(junction.vertex_id for junction in junctions),
        f'component-{side}', 1.0, 0.0, 'arrangement_component', 'PARTIAL',
        f'selection-{side}',
    )
    seam = RefinedSeamPiece(
        f'seam-{side}', (carrier, f'same-group-{side}'), 'circle-seam',
        (0.0, 0.0, 0.0), 1.0, (0.0, 0.0, 1.0), (0.0, 1.0),
        ('junction-0', 'junction-1'), points[0], points[1],
        (patch.patch_id,), ('arc-0',), (f'contact-{side}',),
        ((patch.patch_id, 1),), 'internal', 'CERTIFIED',
    )
    cut = TopologicalCut(
        f'cut-{side}', patch.patch_id, carrier, (0, 1),
        ('junction-0', 'junction-1'), points[0], points[1],
    )
    return MolecularContactSurface(
        side, 'expanded_union_of_balls', 'mcs-v1', (patch,), arcs, junctions,
        (carrier,), (f'contact-{side}',), 1.0, 1.0, 1, 5, None, 0, (),
        'cached analytic area', 0.0, 'PARTIAL', ('topology pending',),
        f'selection-{side}', (seam,), 3, 4, 1, 2, 1,
        (('junction-0', 'boundary_path'),), None, None, 'CERTIFIED', (cut,), 1,
        (('patch-{}'.format(side), 1),), 4, (1, 0, 0), 'UNRESOLVED',
    )


def _files(root):
    return {
        path.relative_to(root): path.read_bytes()
        for path in root.rglob('*') if path.is_file()
    }


def test_cache_only_surface_export_variants_are_deterministic(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('visualization attempted scientific reconstruction')

    import vorpy.src.geometry.molecular_contact as geometry
    for name in ('build_molecular_contact_surface', 'build_molecular_contact_surface_pair',
                 '_build_carrier', '_face_area', '_arc_integral', '_pair_cap_area'):
        monkeypatch.setattr(geometry, name, forbidden)

    surfaces = {'A': _surface('A'), 'B': _surface('B')}
    before = {side: json.dumps(asdict(surface), sort_keys=True) for side, surface in surfaces.items()}

    only_a = export_molecular_contact_bundle({'molecular_contact_surface:A': surfaces['A']}, tmp_path / 'a')
    only_b = export_molecular_contact_bundle({'molecular_contact_surface:B': surfaces['B']}, tmp_path / 'b')
    both_1 = export_molecular_contact_bundle(surfaces, tmp_path / 'both-1', interaction='AB')
    both_2 = export_molecular_contact_bundle(surfaces, tmp_path / 'both-2', interaction='AB')

    assert only_a['capabilities']['molecular_contact_surface:B']['state'] == 'NOT_CALCULATED'
    assert only_b['capabilities']['molecular_contact_surface:A']['state'] == 'NOT_CALCULATED'
    assert both_1['capabilities']['molecular_contact_surface:A']['state'] == 'PARTIAL'
    assert both_1['capabilities']['molecular_contact_surface:B']['state'] == 'PARTIAL'
    assert both_1['layers']['molecule_A']['state'] == 'NOT_CALCULATED'
    assert both_1['layers']['physical_partition_interface']['state'] == 'NOT_CALCULATED'
    assert both_1['layers']['dual_contact_complex']['state'] == 'NOT_CALCULATED'
    assert both_1['layers']['molecular_contact_surface:A']['metadata']['topology_status'] == 'UNRESOLVED'
    assert both_1['layers']['molecular_contact_surface:A']['metadata']['counts']['refined_seams'] == 1
    assert both_1['layers']['molecular_contact_surface:A']['metadata']['counts']['topological_cuts'] == 1
    assert both_1['layers']['molecular_contact_surface:A']['metadata']['topology']['topological_cut_count'] == 1
    mapping = (tmp_path / 'both-1/interfaces/AB/molecular_contact_surface_A/molecular_contact_surface_A_mapping.csv').read_text()
    assert 'topological_cut' not in mapping
    assert _files(tmp_path / 'both-1') == _files(tmp_path / 'both-2')
    assert (tmp_path / 'a/interfaces/interaction/molecular_contact_surface_A/molecular_contact_surface_A.obj').exists()
    assert (tmp_path / 'b/interfaces/interaction/molecular_contact_surface_B/molecular_contact_surface_B_boundary.obj').exists()
    assert all(json.dumps(asdict(surfaces[side]), sort_keys=True) == before[side] for side in surfaces)


def test_partial_topology_surface_still_renders_known_geometry(tmp_path):
    surface = _surface('A')
    metadata = export_molecular_contact_surface(surface, tmp_path)
    assert metadata['geometry_status'] == 'CERTIFIED'
    assert metadata['topology_status'] == 'UNRESOLVED'
    assert metadata['counts']['mesh_triangles'] > 0
    assert metadata['topology']['euler_characteristic'] is None
