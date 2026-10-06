import csv
import json
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np

from vorpy.src.system import System
from vorpy.src.command.vpy_cmnd import Command
from vorpy.src.network import Network
from vorpy.src.interface.selected_export import export_selected_interface
from vorpy.src.boundary import SOLVENT_RESIDUES

PDB = Path(__file__).resolve().parents[2] / 'data' / '2KAI.pdb'


def test_archive_omits_only_runtime_water_lookup():
    from vorpy.src.io.network_archive import _archive_field
    context = {'topology_to_parent': lambda i: (i, 'identity'), 'topology_to_system': {0: 0}}
    topology = {'context': context, 'waters': {}}
    iface = SimpleNamespace(water_topology=topology)
    snapshot = _archive_field(iface, 'interface', 'water_topology')
    assert snapshot == {'context': {'topology_to_system': {0: 0}}, 'waters': {}}
    assert 'topology_to_parent' in topology['context']
    iface.water_geometries = [{'topology_analysis': topology}]
    assert _archive_field(iface, 'interface', 'water_geometries') == [{'topology_analysis': snapshot}]


def test_cached_selected_export_without_selection_or_solve(tmp_path, monkeypatch):
    from vorpy.tests.interface.test_geometry_analysis import iface_for, install_filtration
    from vorpy.tests.geometry.test_alpha_interfaces import _network
    net = _network()
    net.surfs['points'] = [np.array([[1., 0., 0.], [1., 1., 0.], [1., 0., 1.]])] * 3
    net.surfs['tris'] = [np.array([[0, 1, 2]])] * 3
    install_filtration(monkeypatch, net)
    iface = iface_for(net)
    iface.dir = str(tmp_path)
    result = iface.analyze_geometry()
    def prohibited(*args, **kwargs):
        pytest.fail('Exporter recomputed selection or solved a network')
    monkeypatch.setattr(Network, 'build', prohibited)
    monkeypatch.setattr(type(iface._geometry_filtration), 'interface_at', prohibited)
    metadata = export_selected_interface(iface)
    assert metadata['selected_generator_pairs'] == [(10, 13)]
    assert metadata['physical_generator_pairs'] == [(10, 13)]
    assert result.alpha_selection['selected_generator_pairs'] == [[10, 13]]
    assert (tmp_path / 'alpha_interface' / 'selected_surfaces.off').exists()


def test_2kai_loads_with_mixed_author_chains():
    system = System('2kai', make_dir=False)
    assert len(system.balls) == 2246
    assert [(c.name, len(c.atoms)) for c in system.chains] == [('A', 645), ('B', 1154), ('I', 437)]
    assert len(system.residues) == 287
    assert len(system.sol.atoms) == len(system.sol.residues) == 10
    assert system.has_explicit_solvent
    assert all(system.balls.loc[i, 'chn'] is system.sol for i in system.sol.atoms)


@pytest.mark.parametrize('records', [('HOH',), ('HOH', 'NA', 'CL'), ('NA', 'SO4'), ('LIG',)])
def test_solvent_ownership_on_existing_chain(tmp_path, records):
    path = tmp_path / 'mixed.pdb'
    lines = [f'ATOM  {1:5d}  CA  ALA A   1    {0.:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C\n']
    for serial, res in enumerate(records, 2):
        lines.append(f'HETATM{serial:5d}  O   {res:3s} A{serial:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           O\n')
    path.write_text(''.join(lines))
    system = System(str(path), make_dir=False)
    expected = sum(res in SOLVENT_RESIDUES for res in records)
    assert len(system.sol.atoms) == len(system.sol.residues) == expected
    assert system.has_explicit_solvent == bool(expected)
    assert len(system.chains) == 1
    assert len(system.chains[0].atoms) == len(lines) - expected
    assert len(system.residues) == len(lines) - expected
    from vorpy.src.interface.water import _build_interface_index_context
    context = _build_interface_index_context(SimpleNamespace(
        sys=system, net=SimpleNamespace(balls=system.balls),
        group1_indices={0}, group2_indices=set()))
    water_ids = {i for i, row in system.balls.iterrows() if row['res_name'] == 'HOH'}
    assert set(context['water_by_system']) == water_ids


@pytest.mark.slow
@pytest.mark.timeout(1800) if importlib.util.find_spec('pytest_timeout') else pytest.mark.slow
def test_real_dry_cli_selection_and_exports(tmp_path, monkeypatch):
    path = tmp_path / 'dry_2kai.pdb'
    path.write_text(''.join(line for line in PDB.read_text().splitlines(True)
                            if line[:6].strip() not in {'ATOM', 'HETATM'}
                            or line[17:20].strip().upper() not in SOLVENT_RESIDUES))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('sys.argv', ['vorpy', str(path), '-g', 'c', '0', 'and', 'c', '1', '-g', 'c', '2', '-i'])
    calls = []
    original = Network.build
    def observed(net, *args, **kwargs):
        calls.append(net)
        return original(net, *args, **kwargs)
    monkeypatch.setattr(Network, 'build', observed)
    command = Command()
    command.run()
    system = command.sys
    assert len(system.balls) == 2236
    assert not system.has_explicit_solvent
    assert system.sol.atoms == system.sol.residues == []
    iface = system.ifaces[0]
    assert len(iface.group1_indices) == 1799 and len(iface.group2_indices) == 437
    assert calls == [iface.net]
    analysis = iface.geometry_analysis
    assert analysis.metadata['solvent_context'] is False
    assert analysis.metadata['additional_network_solves'] == 0
    assert analysis.alpha_selection['selected_pairs'] > 0
    assert analysis.voronoi_side_1['topology_faces'] > 0
    def prohibited(*args, **kwargs):
        pytest.fail('Export attempted selection or an additional solve')
    monkeypatch.setattr(Network, 'build', prohibited)
    monkeypatch.setattr(type(iface._geometry_filtration), 'interface_at', prohibited)
    metadata = export_selected_interface(iface)
    assert metadata['selected_generator_pairs'] == sorted(map(tuple, analysis.alpha_selection['selected_generator_pairs']))
    assert metadata['physical_generator_pairs'] == sorted(tuple(sorted(s.generator_ids)) for s in iface._geometry_representation.selected_surfaces)
    root = Path(iface.dir) / 'alpha_interface'
    with (root / 'mapping.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert {(int(r['generator_i']), int(r['generator_j'])) for r in rows} == set(metadata['selected_generator_pairs'])
    assert {int(r['surface_id']) for r in rows if r['selected'] == 'True'} == set(metadata['selected_surface_ids'])
    assert (root / 'selected_surfaces.off').is_file()
    mesh_lines = [line for line in (root / 'selected_surfaces.off').read_text().splitlines() if line.strip()]
    point_count, face_count, _ = map(int, mesh_lines[1].split())
    selected_rows = [iface.net.surfs.loc[s.feature_id] for s in iface._geometry_representation.selected_surfaces]
    assert point_count == sum(len(row['points']) for row in selected_rows)
    assert face_count == sum(len(row['tris']) for row in selected_rows)
    np.testing.assert_allclose(
        np.array([list(map(float, line.split()[:3])) for line in mesh_lines[2:2+point_count]]),
        np.concatenate([row['points'] for row in selected_rows]), rtol=0, atol=5.1e-5)
    assert (root / 'selected_interface.pml').is_file()
    (tmp_path / 'dry-2kai-verification.json').write_text(json.dumps(analysis.log_record(), indent=2))
