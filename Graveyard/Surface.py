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