from Visualize.commands.interpret import *


def print_list(names, list_name=None, width=150, height=30, cutoff=15):
    """
    Prints a long list in columns with numbers and allows the user to scroll through the list
    :param names:
    :param list_name:
    :param width:
    :param height:
    :param cutoff:
    :return:
    """
    # Check to see if a list name was provided
    if list_name is None:
        list_name = "My Objects"
    # First find the longest input in the list
    max_len = 0
    for name in names:
        if len(name) > max_len:
            max_len = len(name)
    if max_len > cutoff:
        max_len = cutoff

    # Figure out the columns. num cols = width / # of digits in index of last element + 2 ('. ') + max_len + 2 spaces
    num_cols = int(width / (2 + len(str(len(names) - 1)) + max_len + 2))
    # Print the first set of numbers
    i, row = 0, 0
    # Go through the names row by row, also print the header
    print(list_name)
    while row < height:
        row_str = ""
        for col in range(num_cols):
            if i >= len(names):
                row = 100000000000
            else:
                row_str += str(i) + ". " + " " * (len(str(len(names) - 1)) - len(str(i))) + names[i] + " " * (
                            max_len - len(names[i])) + "  "
                i += 1
        print(row_str)
    # If that is all the data we are done and able to quit
    if len(names) < num_cols * height:
        return
    # In the case where the user wants to see a really long list, allow them to scroll
    scrolling = True
    while scrolling:
        my_response = input("Enter an index or a range or type 'q' to quit. (\'356\' or \'400-600\')\nindex >>>   ")
        if my_response.lower() in quits:
            return
        elif my_response.lower() in helps:
            help_()
        nums = None
        for i in range(len(my_response)):
            if my_response[i] in splitters:
                try:
                    nums = [int(my_response[:i], int(my_response[i + 1:]))]
                except ValueError:
                    nums = None
        # Check to see if a single number has been entered
        if nums is None:
            try:
                nums = [int(my_response)]
            except ValueError:
                nums = None
        # Print the lists
        if nums is not None and len(nums) == 1 and nums[0] < len(names):
            print_list(names[nums[0]:nums[0] + num_cols * height],
                       list_name=list_name + ": Elements " + str(nums[0]) + "-" + str(nums[0] + num_cols * height),
                       height=height, width=width, cutoff=max_len)
        elif nums is not None and nums[0] < nums[1] < len(names):
            # Check to see if the height needs to be changed
            new_height = height
            if nums[1] - nums[0] > num_cols * height:
                new_height = (nums[1] - nums[0]) // num_cols + 1
            print_list(names[nums[0]:nums[1]], list_name=list_name + ": Elements " + str(nums[0]) + "-" + str(nums[1]),
                       height=new_height, width=width, cutoff=max_len)
        else:
            invalid_input(my_response)



def show(sys, usr_npt=None):
    """
    Shows the input group type
    :return:
    """
    # If the user types 'Show' have a catch for it
    if usr_npt is None or len(usr_npt) == 1:
        show_var = get_obj(sys=sys, return_ndx=False).lower()
    # Get the list that the user wants to be shown if none was provided
    elif len(usr_npt) == 2 and usr_npt[1] in my_objects:
        show_var = usr_npt[1].lower()
    else:
        invalid_input(usr_npt)
        return

    # Get the correct list to show the user
    if show_var in mol_objs:
        show_name = "{} Molecules".format(sys.name)
        show_list = sys.mol_names
    elif show_var in res_objs:
        show_name = "{}".format(sys.name)
        show_list = sys.res_names
    elif show_var in atom_objs:
        show_name = "{}".format(sys.name)
        show_list = sys.atom_names
    elif show_var in ndx_objs:
        show_name = "{}".format(sys.name)
        show_list = sys.ndx_names
    else:
        show_name = ""
        show_list = []

    # Show the list
    if len(show_list) == 0:
        print("no objects to show")
        return
    else:
        print_list(show_list, show_name)
