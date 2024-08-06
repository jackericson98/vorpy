import pandas as pd


def mark_doublets(verts):
    # Check to see if the input for vertices is in dataframe format or not
    if isinstance(verts, pd.DataFrame):
        # Instantiate the vdubs list and averts dictionary
        vdubs, b_verts = [], {}
        # Go through each of the vertices to see of the vertex has been found yet
        for i, vert in verts.iterrows():
            # Set the standard doublet value to 0
            vdub = 0
            # We only need to add to and check the first atom in each vertex bc the indices are sorted
            ball = vert['vatoms'][0]
            # If the atom has a list we are going to check that for the vertex and append to it
            if ball in b_verts:
                # Go through the vertices in the list for the atom
                for vert_check in b_verts[ball]:
                    # If found it is a doublet
                    if vert['vatoms'] == verts['vatoms'][vert_check]:
                        vdub = 1
                # Add the vertex to the atom verts list for that atom
                b_verts[ball].append(i)
            # Create a new list for the atom with the vertex's index
            else:
                b_verts[ball] = [i]
            # Add the designation for the vertex
            vdubs.append(vdub)

    else:
        vdubs = None
        print("Gimme that sweet sweet dataframe please")
    return vdubs
