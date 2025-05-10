"""
GUI module for VorPy's graphical user interface.
"""

from vorpy.src.GUI.vorpy_gui import VorPyGUI
from vorpy.src.GUI.system.system_frame import SystemFrame
from vorpy.src.GUI.system.system_exports import SystemExportsWindow
from vorpy.src.GUI.group.groups_frame import GroupsFrame
from vorpy.src.GUI.group.build.build_frame import BuildFrame
from vorpy.src.GUI.group.build.color_settings_window import ColorSettingsWindow
from vorpy.src.GUI.group.export.export_frame import ExportFrame
from vorpy.src.GUI.help.help_window import HelpWindow
from vorpy.src.GUI.progress_window import ProgressWindow

__all__ = [
    # Main GUI Classes
    'VorPyGUI',
    'ProgressWindow',
    
    # System GUI Components
    'SystemFrame',
    'SystemExportsWindow',
    
    # Group GUI Components
    'GroupsFrame',
    'BuildFrame',
    'ColorSettingsWindow',
    'ExportFrame',
    
    # Help and Documentation
    'HelpWindow'
]


