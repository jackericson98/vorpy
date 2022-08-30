# Calculate critical angle function. Finds the critical angle for the surface.
def calc_crit_ang(self):
    # Set the root difference to infinity and the critical angle to pi (exclusive case)
    ang_inc = np.pi / 4
    ang = np.pi / 2
    min_diff = 1e+6
    counter = 0
    # Keep looping until the difference in magnitudes is less than half the minimum distance
    while counter < 50:
        # Get the plus/minus angs
        ang0, ang1 = ang + ang_inc, ang - ang_inc
        proj_points, roots, root_diffs = [], [], []
        # Go to each point projection from the plus/minus angles and find the roots of the surface points
        for ang in ang0, ang1:
            # Get the point using the angle and then find the magnitudes of the projections (i.e. the roots)
            proj_points.append([self.atoms[0].rad * np.cos(ang), self.atoms[0].rad * np.sin(ang), 0])
            roots.append(self.calc_surf_point(proj_points[-1], return_roots=True))
            root_diffs.append(abs(abs(roots[-1][0]) - abs(roots[-1][1])))
        better_ang = [ang0, ang1][root_diffs.index(min(root_diffs))]
        if min(root_diffs) < self.min_dist / 2:
            self.crit_ang = better_ang
            return
        # Narrow down the best estimate of the critical angle
        x = min(root_diffs) / min_diff
        if min(root_diffs) < min_diff:
            min_diff = min(root_diffs)
            ang_inc = ang_inc * x / 2
            ang = better_ang
        else:
            ang_inc = ang_inc * x * 2
        counter += 1
    self.crit_ang = ang

# Go through each triangle on the surface
for i in range(len(self.tris)):
# Make a copy of the triangle for
tri = self.tris[i].copy()
tri.sort()
points = [[nps2d[0][tri[i]], nps2d[1][tri[i]]] for i in range(len(tri))]
com = Point(calc_com(points=points))
print(com, polygon)
if not polygon.contains(com):
    remove_ndxs.append(i)


# Set the counter to 0
counter = 0
# Go through each point on the triangle checking to see if it is an edge point
for j in range(3):
    # If the triangles jth point index is less than the number of vertex & edge points increment the counter
    if tri[j] < len(self.perimeter):
        counter += 1
# If all three of the points are on an edge we need to check it
if counter == 3:
    # Calculate the side distances for the triangle
    side_dists = []
    for k in range(3):
        side_dists.append(calc_dist(self.points[tri[k]], self.points[tri[(k+1) % 3]]))
    # Check the length of one of the longest legs
    if min(side_dists) * 10 < max(side_dists) and max(side_dists) > 3 * self.min_dist:
        remove_ndxs.append(i)
    else:
        # Set up the pass triangle boolean
        keep_tri = False
        # If the triangle has points on either side of a vertex, we exclude it
        for vert_ndx in self.vert_ndxs:
            # Check to see if there is a vertex index between any of the points
            if tri[0] <= vert_ndx <= tri[1] or tri[1] <= vert_ndx <= tri[0] or \
                    tri[1] <= vert_ndx <= tri[2]:
                keep_tri = True
        if tri[0] in self.vert_ndxs or tri[1] in self.vert_ndxs or tri[2] in self.vert_ndxs:
            if tri[2] - tri[0] == 2 or (tri[1] == 1 and tri[2] == len(self.perimeter) - 1):
                keep_tri = True
            else:
                keep_tri = False
        if not keep_tri:
            remove_ndxs.append(i)