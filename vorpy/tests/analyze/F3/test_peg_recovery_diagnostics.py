from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'pyproject.toml').is_file())
script = root / 'scripts/figure_reorganization/F3/investigate_peg_first_layer.py'
if not script.is_file():
    script = root / 'scripts/investigate_peg_first_layer.py'
spec = spec_from_file_location('peg_diagnostics', script)
module = module_from_spec(spec)
spec.loader.exec_module(module)


def test_triclinic_cryst1_is_not_replaced_by_an_orthogonal_box():
    vectors = module.lattice('CRYST1   50.946   50.946   50.946  60.00  60.00  90.00 P 1           1')
    lengths = np.linalg.norm(vectors, axis=1)
    np.testing.assert_allclose(lengths, [50.946] * 3)
    cosines = vectors @ vectors.T / np.outer(lengths, lengths)
    np.testing.assert_allclose(cosines, [[1, 0, .5], [0, 1, .5], [.5, .5, 1]], atol=1e-14)


def test_layer_identity_is_not_inferred_from_a_serial_offset():
    source = ['ATOM    202  OW  SOL    2       1.000   2.000   3.000  1.00  0.00           O\n']
    assert module.layer_membership(source, source) == [0]
    shifted = source[0][:6] + '  201' + source[0][11:]
    with pytest.raises(ValueError, match='serial disagrees'):
        module.layer_membership(source, [shifted])


def test_checkpoint_reader_uses_declared_root_not_last_embedded_network(monkeypatch):
    class Reader:
        def __init__(self, archive):
            pass

        def json(self, member):
            if member == 'metadata.json':
                return {'root_network': 0, 'scheme': 'aw'}
            return {'objects': [{'fields': {k: {'table': 'selected/' + k} for k in ['balls', 'verts', 'edges', 'surfs']}},
                                {'fields': {k: {'table': 'wrong/' + k} for k in ['balls', 'verts', 'edges', 'surfs']}}]}

    monkeypatch.setattr(module, 'Reader', Reader)
    monkeypatch.setattr(module, 'table', lambda reader, path, columns: pd.DataFrame({'selected': [path]}))
    net, _ = module.read_geometry(None)
    assert net.balls['selected'].iloc[0] == 'selected/balls'
