import os
from System.sys_funcs.output.atoms import write_pdb
from System.sys_funcs.output.surfs import write_surfs1
from System.sys_funcs.output.edges import write_edges1
from System.sys_funcs.output.verts import write_off_verts1


class Interface:
    def __init__(self, group1, group2, balls=None, surfs=None, edges=None, verts=None, directory=None, name=None,
                 build=True):

        """
        Interface Object.
        """
        # Interface Groups
        self.g1 = group1
        self.g2 = group2

        # Network components
        self.balls = balls
        self.surfs = surfs
        self.edges = edges
        self.verts = verts

        # Interface attributes
        self.dir = directory
        self.name = name

        # Check if the interface wants to be built
        if build:
            self.get_components()

    def get_components(self):
        """
        Gathers the components of the network to make the interface components
        """
        # First name the interface
        if self.name is None:
            # Set the name to be a mix of the two group names
            self.name = self.g1.name + '__' + self.g2.name + '_Interface'

        # Gather the surfaces shared by the two groups
        if self.surfs is None or len(self.surfs) == 0:
            # Create a set out of the group2 ball ndxs
            g2_bndxs = set(self.g2.ball_ndxs)
            # Get the overlapping surfaces
            self.surfs = self.g1.net.surfs[self.g1.net.surfs['balls'].apply(lambda balls: any(ball in g2_bndxs for ball in balls))]

        # Gather the edges associated with the surfaces of the interface
        if self.balls is None or len(self.balls) == 0:
            # Flatten the list of lists and extract unique indices using set for uniqueness
            self.balls = [idx for sublist in self.surfs['balls'] for idx in sublist]

        # Gather the edges associated with the surfaces of the interface
        if self.edges is None or len(self.edges) == 0:
            # Flatten the list of lists and extract unique indices using set for uniqueness
            unique_edges = [idx for sublist in self.surfs['edges'] for idx in sublist]
            # Gather the edges
            self.edges = self.g1.net.edges.loc[unique_edges]

        # Gather the verts associated with the surfaces of the interface
        if self.verts is None or len(self.verts) == 0:
            # Flatten the list of lists and extract unique indices using set for uniqueness
            unique_verts = [idx for sublist in self.surfs['verts'] for idx in sublist]
            # Gather the edges
            self.verts = self.g1.net.verts.loc[unique_verts]

    def export_info(self):
        """
        Exports a txt file with the information for the interface. Surface area, contributing balls, outer edge length
        """

    def export(self, all_=False, balls=False, surfs=False, sep_surfs=False, edges=False, sep_edges=False, verts=False,
               sep_verts=False, info=False):
        """
        Interface exports method. Controls what goes out of the interface
        """
        # Set the interface export directory
        if self.dir is None:
            i = 1
            my_dir = self.g1.sys.files['dir'] + "/" + self.name
            first = True
            while os.path.exists(my_dir):
                if first:
                    my_dir += "__"
                    first = False
                my_dir = my_dir[:-(1 + len(str(i)))] + '_' + str(i)
                i += 1
            self.dir = my_dir
            os.mkdir(self.dir)
        os.chdir(self.dir)

        # Export the balls of the interface
        if balls or all_:
            # Write the pdb
            write_pdb(atoms=self.balls, file_name=self.name + '_atoms', sys=self.g1.sys, directory=self.dir)

        # Export the surfaces
        if sep_surfs or all_:
            # Make the surface exporting directory
            os.mkdir(self.dir + '/surfs')
            # loop through the surfaces
            for i, surf in self.surfs.iterrows():
                # Export the surface
                write_surfs1([surf], '_'.join(surf['balls']), self.g1.settings, directory=self.dir + '/surfs')

        # Export the surfaces as one big group
        if surfs or all_:
            # Export the surfaces file
            write_surfs1(self.surfs, self.name + '_surface', self.g1.settings, directory=self.dir)

        # Export the separate edges
        if sep_edges or all_:
            # Make the edges directory
            os.mkdir(self.dir + '/edges')
            # Loop through the edges
            for i, edge in self.edges.iterrows():
                # Export the edge
                write_edges1([edge], '_'.join(edge['balls']), directory=self.dir + '/edges')

        # Export the full edge file
        if edges or all_:
            # Export the edges in one single file
            write_edges1(self.edges, self.name + '_edges', directory=self.dir)

        # Export the separate vertices
        if sep_verts or all_:
            # Make the edges directory
            os.mkdir(self.dir + '/verts')
            # Loop through the edges
            for i, vert in self.verts.iterrows():
                # Export the edge
                write_off_verts1([vert], '_'.join(vert['balls']), directory=self.dir + '/edges')

        # Export the full edge file
        if verts or all_:
            # Export the edges in one single file
            write_off_verts1(self.verts, self.name + '_verts', directory=self.dir)

        # Export the info
        if info or all_:
            pass
