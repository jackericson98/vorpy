import os
from System.sys_funcs.calcs.sorting import ndx_search
from System.sys_funcs.output.verts import write_off_verts
from System.sys_funcs.output.edges import write_edges
from System.sys_funcs.output.surfs import write_surfs


class Interface:
    def __init__(self, group1, group2, name=None, balls=None, verts=None, edges=None, surfs=None, sa=0, curv=None):
        """ Initialize the input parameters for an interface"""

        # Interfacial attributes
        self.name = name
        self.group1 = group1
        self.group2 = group2

        # Network Attributes
        self.balls = balls
        self.verts = verts
        self.edges = edges
        self.surfs = surfs

        # Interface Data
        self.sa = sa
        self.curv = curv

        # Run the set up so everything is in place
        self.setup()

    def setup(self):
        """ Make sure everything is as is should be before starting the investigation"""
        # Make sure one of the groups has a network solved
        if self.group1.net is None or self.group1.net.verts is None or len(self.group1.net.verts) == 0:
            if self.group2.net is None or self.group2.net.verts is None or len(self.group2.net.verts) == 0:
                self.group1.build()
            else:
                self.group1, self.group2 = self.group2, self.group1
        # If there is no name assigned the default based on the groups
        if self.name is None:
            self.name = self.group1.name + '_' + self.group2.name + '_iface'

    def get_iface(self, group1, group2):
        # Set up the lists and
        i_balls, i_verts, i_edges, i_surfs, i_sa = [], [], [], [], 0
        group1.iface_sa = 0
        iface_curvs = []
        # Go through the balls in the group
        self.balls = [_ for _ in self.group1.ball_ndxs if _ in self.group2.ball_ndxs]

        # Get the interface elements

        # Function to check if any item in the list is in another list
        def any_in_list(check_list, target_list):
            return any(item in target_list for item in check_list)

        # Find the vertices
        self.verts = self.group1.net.verts[self.group1.net.verts['balls'].apply(lambda x: any_in_list(x, self.balls))]
        # Find the edges
        self.edges = self.group1.net.edges[self.group1.net.edges['balls'].apply(lambda x: any_in_list(x, self.balls))]
        # Find the surfaces
        self.surfs = self.group1.net.surfs[self.group1.net.surfs['balls'].apply(lambda x: any_in_list(x, self.balls))]

    def export_iface_verts(self, directory=None):
        """
        Exports the interfacial vertices between the group and its bff
        :param grp: Group object for exporting
        :param directory: Directory to export to
        """
        # Move to the directory
        if directory is not None and os.path.exists(directory):
            os.chdir(directory)
        # write the vertices
        write_off_verts(self.group1.net, self.verts, directory=directory, file_name=self.name + "_verts")

    def export_iface_edges(self, directory=None):
        """
        Exports the edges of the interface
        :param grp: Group to pull the interface from
        :param directory: Output directory for the interface
        """
        # Move to the directory
        if directory is not None and os.path.exists(directory):
            os.chdir(directory)
        # write the vertices
        write_edges(self.group1.net, self.edges, directory=directory, file_name=self.name + "_edges")

    def export_iface_info(grp, directory=None):
        """
        Exports the information for an interface
        :param grp: The group that holds the interface information
        :param directory: Output directory for the group interface info
        """
        # Move to the directory
        if directory is not None and os.path.exists(directory):
            os.chdir(directory)
        # Create the file
        with open("info.txt", 'w', encoding='utf-8') as info:
            # Write the main header
            info.write(grp.name + " - " + grp.bff.name + " interface \n\n")
            # Information sub header
            info.write("Interface:\n\n")
            # Write the information
            info.write(
                "  {} Surfaces, {} {} atoms, {} {} atoms\n".format(len(grp.iface_surfs), len(grp.atoms), grp.name,
                                                                   len(grp.bff.atoms), grp.bff.name))
            # Network counts
            info.write("  {} Vertices, {} Edges\n\n".format(len(grp.iface_verts), len(grp.iface_edges),
                                                            len(grp.iface_surfs)))
            # Write the analysis header
            info.write("\nAnalysis:\n\n")
            # Write the analysis
            info.write(u"  Surface Area = {:.5f} \u212B\u00B2, Average Curvature = {:.5}\n\n"
                       .format(grp.iface_sa, grp.iface_curv))
            # Surfaces header
            info.write("\nSurfaces:\n\n")
            # Go through each of the surfaces in the group
            for surf in grp.iface_surfs:
                info.write("  Surface {} - \n".format(grp.sys.net.surfs['satoms'][surf]))
                info.write("    Surface Area = {:.5f} \u212B\u00B2\n".format(grp.sys.net.surfs['sa'][surf]))
                info.write("    Volume contributions = {:.5f}, {:.5f} \u212B\u00B3\n"
                           .format(grp.sys.net.surfs['vols'][surf][0], grp.sys.net.surfs['vols'][surf][0]))
                info.write("    Gaussian Curvature = {:.5f}\n".format(grp.sys.net.surfs['curv'][surf]))

    def export(self, balls=True, verts=True, edges=True, surfs=True, info=True, cells=True):
        pass