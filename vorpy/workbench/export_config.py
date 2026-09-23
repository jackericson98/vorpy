"""One export configuration shared by controls, summary, and export dispatch."""
from dataclasses import dataclass, field
from pathlib import Path

from vorpy.workbench.services.export import PRESETS

PRESET_LABELS = {'Small': 'Small', 'Medium': 'Standard', 'Large': 'Full', 'Custom': 'Custom'}
GENERAL_OPTIONS = (('info', 'Network info'), ('logs', 'Logs'))
ATOMIC_OPTIONS = (('atoms', 'Group atoms'), ('surr_atoms', 'Surface neighbors'),
                  ('surr_resids', 'Surrounding residues'), ('full_molecule', 'Full molecule'))
GEOMETRY_ROWS = (('Surfaces', 'surfs'), ('Edges', 'edges'), ('Vertices', 'verts'))
GEOMETRY_MODES = (('All', ''), ('Shell', 'shell_'), ('Separate', 'sep_'), ('Cell', 'atom_'))
OPTION_LABELS = dict(GENERAL_OPTIONS + ATOMIC_OPTIONS)
OPTION_LABELS.update({prefix + key: f'{mode} {label.lower()}'
                      for label, key in GEOMETRY_ROWS for mode, prefix in GEOMETRY_MODES})


@dataclass
class ExportConfig:
    options: dict[str, bool] = field(default_factory=lambda: dict.fromkeys(PRESETS['Medium'], True))
    directory: str = ''
    precision: int = 3
    file_type: str = 'off'
    atom_format: str = 'pdb'
    color_overrides: dict[str, str] = field(default_factory=dict)

    @property
    def preset(self):
        selected = {key for key, enabled in self.options.items() if enabled}
        return next((name for name, keys in PRESETS.items() if selected == set(keys)), 'Custom')

    def select_preset(self, preset):
        self.options = dict.fromkeys(PRESETS[preset], True)

    def output_directory(self, result):
        if self.directory.strip():
            return str(Path(self.directory.strip()).expanduser().resolve())
        if result is not None and result.source:
            return str((result.source.parent / 'VorPy_Output').resolve())
        return ''


def config_from_json(state):
    """Restore known export fields from a saved Workbench session."""
    config = ExportConfig()
    if not isinstance(state, dict):
        return config
    if isinstance(state.get('options'), dict):
        config.options = {key: value for key, value in state['options'].items()
                          if key in OPTION_LABELS and isinstance(value, bool)}
    if isinstance(state.get('directory'), str):
        config.directory = state['directory']
    if isinstance(state.get('precision'), int):
        config.precision = max(0, min(12, state['precision']))
    if state.get('file_type') in {'off', 'ply', 'vtp'}:
        config.file_type = state['file_type']
    if state.get('atom_format') in {'pdb', 'xyz'}:
        config.atom_format = state['atom_format']
    if isinstance(state.get('color_overrides'), dict):
        from matplotlib.colors import is_color_like
        config.color_overrides = {key: value for key, value in state['color_overrides'].items()
                                  if key in {'surfaces', 'edges', 'vertices'}
                                  and isinstance(value, str) and is_color_like(value)}
    return config
