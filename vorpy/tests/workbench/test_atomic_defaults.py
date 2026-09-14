import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from types import SimpleNamespace

import pandas as pd
import pytest
from PySide6.QtWidgets import QApplication

from vorpy.workbench.atomic_defaults import apply_atom_defaults, apply_system_defaults, validate_defaults
from vorpy.workbench.domain import Atom, AnalysisResult
from vorpy.workbench.project import result_from_json, result_to_json
from vorpy.workbench.ui.atomic_defaults_dialog import AtomicDefaultsDialog


def atom(**kwargs):
    return Atom(0, 1, 'CA', 'C', (0., 0., 0.), residue_name='GLY', radius=1.8,
                mass=12.011, **kwargs)


def test_independent_property_precedence_and_reset_round_trip():
    original = atom()
    settings = {'C': {'radius': 2., 'mass': 13., 'charge': 0.2},
                'RES:GLY:CA': {'radius': 1.9, 'charge': 0.}}
    edited = apply_atom_defaults(original, settings)
    assert (edited.radius, edited.mass, edited.charge) == (1.9, 13., 0.)
    result = result_from_json(result_to_json(AnalysisResult(None, 'test', atoms=[edited])))
    reset = apply_atom_defaults(result.atoms[0], {})
    assert (reset.radius, reset.mass, reset.charge) == (1.8, 12.011, None)


def test_ion_overrides_do_not_change_residue_alpha_carbon():
    calcium = Atom(1, 2, 'CA', 'CA', (0., 0., 0.), residue_name='CA', radius=1.)
    values = {'CA': {'mass': 40.}, 'ION:CA': {'radius': 1.2, 'charge': 2.}}
    assert apply_atom_defaults(calcium, values).charge == 2.
    assert apply_atom_defaults(atom(), values).radius == 1.8


def test_system_properties_applied_before_build_and_spatial_cache_invalidated():
    invalidations = []
    system = SimpleNamespace(
        balls=pd.DataFrame([{'element': 'C', 'res_name': 'GLY', 'name': 'CA', 'rad': 1.8, 'mass': 12.011, 'charge': ''}]),
        invalidate_spatial_index=lambda: invalidations.append(True),
    )
    apply_system_defaults(system, {'RES:GLY:CA': {'radius': 2.1, 'mass': 13., 'charge': -0.5}})
    assert system.balls.loc[0, ['rad', 'mass', 'charge']].tolist() == [2.1, 13., -0.5]
    assert system.max_atom_rad == 2.1
    assert invalidations == [True]


@pytest.mark.parametrize('name,value', [('radius', 0), ('mass', -1), ('charge', float('nan'))])
def test_invalid_properties_are_rejected(name, value):
    with pytest.raises(ValueError):
        validate_defaults({'C': {name: value}})


def test_dialog_draft_survives_scope_changes_and_cancel_does_not_mutate_input():
    QApplication.instance() or QApplication([])
    source = {'O': {'radius': 1.6}}
    dialog = AtomicDefaultsDialog(source, [atom()])
    dialog.element_editors['mass'].enabled.setChecked(True)
    dialog.element_editors['mass'].value.setValue(13.)
    dialog.tabs.setCurrentIndex(1)
    index = next(i for i in range(dialog.residue_list.count()) if dialog.residue_list.item(i).data(256) == 'GLY')
    dialog.residue_list.setCurrentRow(index)
    row = next(i for i in range(dialog.residue_table.rowCount()) if dialog.residue_table.item(i, 0).text() == 'CA')
    assert dialog.residue_table.item(row, 1).text() == 'C'
    radius = dialog.residue_table.cellWidget(row, 2)
    assert radius.value.value() == 1.8
    radius.enabled.setChecked(True)
    radius.value.setValue(1.95)
    dialog.residue_list.setCurrentRow(0)
    dialog.residue_list.setCurrentRow(index)
    assert dialog.values()['RES:GLY:CA']['radius'] == 1.95
    assert dialog.values()['C']['mass'] == 13.
    assert dialog.values()['O']['radius'] == 1.6
    dialog.reject()
    assert source == {'O': {'radius': 1.6}}


def test_periodic_cells_do_not_overlap_at_small_window_size():
    QApplication.instance() or QApplication([])
    dialog = AtomicDefaultsDialog({})
    dialog.resize(850, 540)
    dialog.show()
    QApplication.processEvents()
    buttons = list(dialog.element_buttons.values())
    for i, first in enumerate(buttons):
        for second in buttons[i + 1:]:
            assert not first.geometry().intersects(second.geometry())
    dialog.close()
