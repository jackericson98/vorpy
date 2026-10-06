"""Presentation order for sectioned logs; column names remain the schema."""

import csv


ATOM_ORDER = [
    'Index', 'Name', 'Chain', 'Residue', 'Residue Sequence',
    'X', 'Y', 'Z', 'Radius', 'Mass',
    'Volume', 'Van Der Waals Volume', 'Surface Area', 'Sphericity', 'Isometric Quotient',
    'Complete Cell?', 'Inner Ball?',
    'Number of Neighbors', 'Closest Neighbor', 'Closest Neighbor Distance', 'neighbors',
    'Number of Overlaps', 'Contact Area', 'Non-Overlap Volume', 'Overlap Volume',
    'Layer Distance Average', 'Layer Distance RMSD', 'Minimum Point Distance', 'Maximum Point Distance',
]
GROUP_ORDER = ['Name', 'Volume', 'VDW Volume', 'Surface Area', 'Mass', 'Density']
SURFACE_ORDER = ['Index', 'Ball 1', 'Ball 2', 'Surface Area',
                 'Ball 1 Volume Contribution', 'Ball 2 Volume Contribution', 'Contact Area', 'Overlap']
BUILD_ORDER = ['Name', 'Location', 'Completion Date', 'vorPy version', 'Network Type',
               'Surface Resolution', 'Box Size', 'Maximum Allowable Vertex', 'Maximum Found Vertex',
               'Total Time', 'Vertex Time', 'Connect Time', 'Surface Building Time', 'Analysis time']
SPATIAL_ORDER = ['Center of Mass', 'VDW Center of Mass', 'Bounding Box',
                 'Moment of Inertia', 'Moment of Inertia Tensor', 'Spatial Moment of Inertia']


def ordered_columns(section, headers):
    """Group related measurements without dropping aliases or unknown columns."""
    first = {'Atoms': ATOM_ORDER, 'Surfaces': SURFACE_ORDER,
             'group information': GROUP_ORDER, 'build information': BUILD_ORDER}.get(section)
    if first is None:
        return list(headers)  # Edge and vertex geometry already has a compact order.
    mean = [h for h in headers if 'Mean' in h and 'Curvature' in h]
    gaussian = [h for h in headers if ('Gaussian' in h or 'Gauss' in h) and 'Curvature' in h]
    preferred = first + mean + gaussian + ['Representative Surface Energy']
    preferred += [h for h in headers if h not in preferred and h not in SPATIAL_ORDER]
    preferred += SPATIAL_ORDER
    return list(dict.fromkeys(h for h in preferred if h in headers))


class LogWriter:
    """Reorder each header and its data with the same cached permutation."""

    def __init__(self, stream, **kwargs):
        self.writer = csv.writer(stream, **kwargs)
        self.section = None
        self.expect_header = False
        self.indices = None

    def writerow(self, row):
        if len(row) == 1 and row[0] in {
            'build information', 'build informaiton', 'group information',
            'Atoms', 'Surfaces', 'Edges', 'Vertices', 'Interface Geometry',
        }:
            self.section = 'build information' if row[0] == 'build informaiton' else row[0]
            self.expect_header = True
            self.indices = None
            return self.writer.writerow([self.section])
        if self.expect_header:
            headers = ordered_columns(self.section, row)
            self.indices = [row.index(h) for h in headers]
            self.expect_header = False
        if self.indices is not None:
            if len(row) != len(self.indices):
                raise ValueError(f'{self.section}: header and row widths differ')
            row = [row[i] for i in self.indices]
        return self.writer.writerow(row)
