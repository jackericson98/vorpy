"""Integration-only regressions against the frozen static physical geometry."""
import json
from io import StringIO

import pytest

from scripts.export_open_interface_geometry_examples import ROOT, MATCHED, fixture_analysis
from vorpy.src.interface.geometry_analysis import write_geometry_info, write_geometry_log
from vorpy.src.inputs.logs import read_logs
from vorpy.src.analyze.tools.compare.read_logs2 import read_logs2
from vorpy.src.log_columns import LogWriter
from vorpy.tests.inputs.test_log_column_order import write_log


@pytest.fixture(scope='module')
def static_analyses():
    if not (MATCHED/'power_geometry_cache.pkl').exists():
        pytest.skip('Frozen local matched-probe physical geometry is unavailable')
    from vorpy.src.network.network import Network
    from vorpy.src.analyze import open_interface_curvature as kernel
    calls = []
    original = kernel.integrated_mean_curvature
    def observed(*args, **kwargs):
        calls.append(kwargs)
        return original(*args, **kwargs)
    def prohibited(*args, **kwargs):
        raise AssertionError('Integration attempted an additional network solve')
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(Network, 'build', prohibited)
        patch.setattr(kernel, 'integrated_mean_curvature', observed)
        results = {scheme:fixture_analysis(scheme) for scheme in ('pow','aw')}
    assert len(calls) == 2
    return results


def test_power_preserves_precise_intrinsic_and_completed_values(static_analyses):
    v = static_analyses['pow'].voronoi_side_1
    assert (v['topology_vertices'],v['topology_edges'],v['topology_faces']) == (605,955,350)
    assert v['euler_characteristic'] == 0
    assert v['connected_components'] == 2
    assert v['boundary_loops'] == 4
    assert v['genus_by_component'] == [0,0]
    assert v['manifold_incidence_certified']
    assert v['area'] == pytest.approx(898.1185793538328,abs=1e-10)
    assert v['total_integrated_mean_curvature'] == pytest.approx(-64.45049624200274,abs=1e-12)
    assert v['mean_curvature_certified']
    assert v['mean_supported_internal_seams'] == 835
    assert v['mean_unresolved_internal_seams'] == 0
    assert v['boundary_integrated_mean_curvature'] == 0
    assert v['intrinsic_integrated_gaussian_curvature'] == pytest.approx(5.188337486733097,abs=1e-12)
    assert v['boundary_corner_turning'] == pytest.approx(-5.1883374867332295,abs=1e-12)
    assert v['gauss_bonnet_lhs'] == pytest.approx(0,abs=2e-12)
    assert v['gauss_bonnet_residual'] == pytest.approx(-1.3233858453531866e-13,abs=2e-13)
    assert v['gaussian_curvature_certified']
    assert v['intrinsic_integrated_gaussian_curvature'] != v['gauss_bonnet_expected']


def test_aw_retains_partial_and_subcomplex_scope(static_analyses):
    v = static_analyses['aw'].voronoi_side_1
    assert (v['topology_vertices'],v['topology_edges'],v['topology_faces']) == (584,924,339)
    assert v['euler_characteristic'] == -1
    assert v['connected_components'] == 2 and v['boundary_loops'] == 5
    assert v['genus_by_component'] == [0,0]
    assert v['area'] == pytest.approx(787.0925777403216,abs=1e-10)
    assert v['partial_integrated_mean_curvature'] == pytest.approx(-54.150160911300546,abs=1e-10)
    assert v['total_integrated_mean_curvature'] is None
    assert not v['mean_curvature_certified']
    assert v['mean_unsupported_faces'] == 5
    assert v['mean_unresolved_internal_seams'] == 14
    assert v['intrinsic_integrated_gaussian_curvature'] is None
    assert not v['gaussian_curvature_certified']
    assert v['gaussian_accounting_scope'] == 'supported_subcomplex'
    assert v['gauss_bonnet_lhs'] is None and v['gauss_bonnet_residual'] is None
    subset = v['supported_subcomplex_gauss_bonnet']
    assert subset['residual'] == pytest.approx(-0.0006654952459843599,abs=1e-10)


@pytest.mark.parametrize('scheme',['pow','aw'])
def test_reversal_and_zero_solves(static_analyses, scheme):
    result = static_analyses[scheme]
    forward, reverse = result.voronoi_side_1,result.voronoi_side_2
    for key in ('smooth_integrated_mean_curvature','internal_seam_integrated_mean_curvature',
                'partial_integrated_mean_curvature','total_integrated_mean_curvature'):
        assert reverse[key] == (None if forward[key] is None else -forward[key])
    for key in ('intrinsic_integrated_gaussian_curvature','gauss_bonnet_lhs',
                'gauss_bonnet_residual','area','genus_by_component','boundary_loops'):
        assert reverse[key] == forward[key]
    assert result.metadata['additional_network_solves'] == 0


@pytest.mark.parametrize('reader',[read_logs,read_logs2])
@pytest.mark.parametrize('scheme',['pow','aw'])
def test_extended_and_old_logs(reader,scheme,static_analyses,tmp_path):
    path = tmp_path/'logs.csv'
    write_log(path,LogWriter)
    old = reader(str(path),return_dict=True)
    with path.open('a',newline='') as stream:
        write_geometry_log(LogWriter(stream),static_analyses[scheme])
    parsed = reader(str(path),return_dict=True)
    extension = parsed.pop('interface geometry')
    assert parsed == old
    assert extension == json.loads(json.dumps(static_analyses[scheme].log_record()))
    assert 'vor_intrinsic_integrated_gaussian_curvature' in extension
    assert 'vor_gauss_bonnet_lhs' in extension


@pytest.mark.parametrize('scheme',['pow','aw'])
def test_text_distinguishes_intrinsic_and_completion(static_analyses,scheme):
    stream = StringIO(); write_geometry_info(stream,static_analyses[scheme])
    text = stream.getvalue()
    assert 'Intrinsic Gaussian curvature' in text and 'Gauss-Bonnet' in text
    if scheme == 'pow':
        assert 'Intrinsic total: 5.188337' in text
    else:
        assert 'Intrinsic total: UNRESOLVED' in text
        assert 'Supported-subcomplex diagnostic' in text


def test_real_edta_cached_pipeline_regression():
    old_path = ROOT/'output/interface_geometry_phase1/EDTA_interface/geometry_analysis.json'
    new_path = ROOT/'output/open_interface_integration/EDTA/geometry_analysis.json'
    if not old_path.exists() or not new_path.exists():
        pytest.skip('Run the documented real EDTA cached-export reproduction first')
    old,new = (json.loads(path.read_text()) for path in (old_path,new_path))
    assert new['alpha_selection']['eligible_bicolor_pairs'] == 6
    assert new['alpha_selection']['selected_pairs'] == 4
    for key in ('interface_atoms','components','mean_certified','gaussian_total_certified'):
        assert new['voronoi_side_1'][key] == old['voronoi_side_1'][key]
    for key in ('area','smooth_K_partial'):
        assert new['voronoi_side_1'][key] == pytest.approx(old['voronoi_side_1'][key],abs=1e-12,rel=0)
    assert new['metadata']['additional_network_solves'] == 0
