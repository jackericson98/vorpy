# if b0[0] <= loc[0] <= b1[0] and b0[1] <= loc[1] <= b1[1] and b0[2] <= loc[2] <= b1[2]:
#     self.loc, self.rad = loc, rad
# else:
#     return
# # Otherwise, choose the first
# else:
# loc, rad = verts[0][0], abs(verts[0][1])
# # Check to see if the vertex is in the box or not
# if b0[0] <= loc[0] <= b1[0] and b0[1] <= loc[1] <= b1[1] and b0[2] <= loc[2] <= b1[2]:
#     self.loc, self.rad = loc, rad
# else:
#     Worst
# case
# scenario
# try the Hu Method
# pass
# loc = self.fv2()
# if len(loc) > 0:
#     self.loc = loc
#     self.rad = np.linalg.norm(self.loc - self.atoms[0].loc) - self.atoms[0].rad