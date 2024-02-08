
################################################  Load Grants Vertices #################################################

my_check_file = open("C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/vdos results/EDTA_Mg_generic_atoms/EDTA_Mg_generic.txt", 'r')

my_check_file = my_check_file.readlines()
check_verts = []
for i in range(1, len(my_check_file)):
    line = my_check_file[i].split(',')
    # Create the array
    my_ndxs = line[0:4]

    my_tested_ndxs = []
    # Go through the indices testing if they are integers or not
    for ndx in my_ndxs:
        try:
            my_int = int(ndx) - 1
            my_tested_ndxs.append(my_int)
        except IndexError:
            my_tested_ndxs = None
            break
    if my_tested_ndxs is not None:
        # Add it to the list
        check_verts.append(my_tested_ndxs + [round(float(_), 3) for _ in line[4:7]])


######################################## Load My Vertices #############################################################

my_file = open("C:/Users/jacke/OneDrive - Georgia State University/GSU NSC/Jack/Vorpy/test_data/Outputs/EDTA_Mg/EDTA_Mg_verts.txt", 'r')
my_file = my_file.readlines()
my_verts, my_check_verts, missing_verts = [], [], []
for i in range(len(my_file)):
    line = my_file[i].split()
    # Create the array
    if line[0].lower() == 'vert':
        my_ndxs = line[1:5]
        my_tested_ndxs = []
        # Go through the indices testing if they are integers or not
        for ndx in my_ndxs:
            try:
                my_tested_ndxs.append(int(ndx))
            except IndexError:
                my_tested_ndxs = None
                break
        if my_tested_ndxs is not None:
            for ndx in check_verts:
                if ndx[:4] == my_tested_ndxs:
                    my_check_verts.append(ndx)
                    my_verts.append(my_tested_ndxs + [round(float(_), 3) for _ in line[5:9]])
                else:
                    missing_verts.append(ndx)
for i in range(len(my_check_verts)):

    print(my_check_verts[i], my_verts[i])
missing_verts = []
for vert in check_verts:
    if vert not in my_check_verts:
        missing_verts.append(vert)
print(len(my_check_verts))
print(missing_verts)


