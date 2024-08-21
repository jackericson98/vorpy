

# def CDT(points, edge_ndxs, alpha=0.3):
#     # Find minimum and maximum x and y values
#     min_x = min(points, key=lambda x: x[0])[0]
#     max_x = max(points, key=lambda x: x[0])[0]
#     min_y = min(points, key=lambda x: x[1])[1]
#     max_y = max(points, key=lambda x: x[1])[1]
#     # Sort the points in terms of the y values
#     sorted_points = zip(sorted(enumerate(points), key=lambda x: (-x[1][1], x[1][0])))
#     point_dict = {_[0]: _[1] for _ in sorted_points}
#     new_edge_ndxs = [(sorted_points.index(_[0], key=lambda x: x[0]))]
#     for edge in edge_ndxs:
#         edge_mapping[]
#     print(sorted_points)
#     print(index_map)
#     # Adjust the edge_ndxs for the sorted points
#     # Note: Need to handle indices that might be out of range or incorrect as given in the example
#     updated_edge_ndxs = [(index_map.get(idx1 - 1), index_map.get(idx2 - 1)) for idx1, idx2 in edge_ndxs if
#                          idx1 - 1 in index_map and idx2 - 1 in index_map]
#
#     # Calculate the delta values for x and y points
#     del_x, del_y = alpha * (max_x - min_x), alpha * (max_y - min_y)
#     # Get the far out projection points
#     p_1, p_2 = [min_x - del_x, min_y - del_y], [max_x - del_x, min_y - del_y]
#     # Figure out which points coincide with the edge events
#     print(sorted_points)
#     print(updated_edge_ndxs)


if __name__ == '__main__':
    # points = [[0.5, 0], [0.1, -1], [0, 0], [0.5, -1], [0.5, 1], [1, 0], [0.25, 0.1], [0.75, -0.1]]
    # edge_ndxs = [[2, 3], [3, 5], [5, 4], [4, 2]]
    # CDT(points, edge_ndxs)

    # Example data: a list of points (x, y) and a list of tuples referring to indices in this list
    points = [(3, 2), (5, 5), (1, 1), (7, 2), (4, 5), (9, 9)]
    edge_ndxs = [(4, 5), (3, 1), (5, 2)]  # Assuming indices start at 0 and are within the correct range

    # Step 1: Sort the points by y-values descending, and x-values ascending on tie
    sorted_points = sorted(enumerate(points), key=lambda x: (-x[1][1], x[1][0]))
    index_map = {original_index: new_index for new_index, (original_index, _) in enumerate(sorted_points)}

    # Adjust the edge_ndxs for the sorted points
    updated_edge_ndxs = [(index_map[idx1], index_map[idx2]) for idx1, idx2 in edge_ndxs]

    # Display results
    print("Original Points:", points)
    print("Sorted Points:", [x for _, x in sorted_points])
    print("Original Edge Indices:", edge_ndxs)
    print("Updated Edge Indices:", updated_edge_ndxs)

