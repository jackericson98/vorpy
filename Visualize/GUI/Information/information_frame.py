from Visualize.GUI.Information.group.GrpInfo import group_frame
from Visualize.GUI.Information.System.SysInfo import system_frame


def create_information_section(gui, parent):
    system_frame(gui, parent)
    group_frame(gui, parent)



