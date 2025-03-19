from Visualize.GUI.Information.group.GrpInfo import group_frame
from Visualize.GUI.Information.System.SysInfo import SystemFrame


def create_information_section(gui, parent):
    SystemFrame(gui, parent)
    group_frame(gui, parent)



