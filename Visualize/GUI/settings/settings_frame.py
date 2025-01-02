from Visualize.GUI.settings.export.export_frame import export_frame
from Visualize.GUI.settings.build.build_frame import build_frame


def create_settings_section(root, parent):
    build_frame(root, parent)
    export_frame(root, parent)
