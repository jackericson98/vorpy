from Visualize.GUI.settings.export.export_frame import ExportFrame
from Visualize.GUI.settings.build.build_frame import BuildFrame


def create_settings_section(root, parent):
    build_frame = BuildFrame(root, parent)
    export_frame = ExportFrame(root, parent)
