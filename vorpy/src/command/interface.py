from itertools import combinations


def build_interfaces(sys, num_requested_groups):
    interface_pairs = get_interface_pairs(
        sys=sys,
        num_requested_groups=num_requested_groups,
    )

    for group_1, group_2 in interface_pairs:
        if group_2 is None:
            print(
                f'\nInterface mode: "{group_1.name}" '
                f"against surrounding atoms"
            )
        else:
            print(
                f'\nInterface mode: "{group_1.name}" '
                f'against "{group_2.name}"'
            )


    return interface_pairs


def get_interface_pairs(sys, num_requested_groups):
    """
    Determine which interfaces should be built.

    Returns
    -------
    list[tuple]
        Each tuple contains two interface sides.

        A value of None for the second side means:
        build the group against its surrounding molecules.
    """

    groups = sys.groups

    if groups is None or len(groups) == 0:
        raise ValueError("Interface mode requires at least one system group.")

    # No explicit -g command:
    # default/main group against surrounding molecules.
    if num_requested_groups == 0:
        return [(groups[0], None)]

    # One explicit group:
    # that group against surrounding molecules.
    if num_requested_groups == 1:
        return [(groups[0], None)]

    # Two or more explicit groups:
    # build every unique pairwise interface.
    return list(combinations(groups, 2))