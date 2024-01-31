import csv
import pandas as pd


def check_verts(complete_logs, test_logs):
    with open(complete_logs, 'r') as good:
        good_logs = csv.reader(good)
        checking_verts = False
        my_good_verts, good_vert_atoms = [], []
        for line in good_logs:
            # Check to see if we are in the vertices part of the logs
            if line[0] == 'Vertices':
                checking_verts = True
                continue
            # If we arent cheking vertices keep going
            if not checking_verts:
                continue
            my_good_verts.append({'num': line[0], 'atoms': line[1:5], 'loc': line[5:8], 'rad': line[8]})
            good_vert_atoms.append(line[1:5])

    with open(test_logs, 'r') as test:
        check_logs = csv.reader(test)
        checking_verts = False
        my_check_verts, test_vert_atoms = [], []
        for line in check_logs:
            # Check to see if we are in the vertices part of the logs
            if line[0] == 'Vertices':
                checking_verts = True
                continue
            # If we arent cheking vertices keep going
            if not checking_verts:
                continue
            my_check_verts.append({'num': line[0], 'atoms': line[1:5], 'loc': line[5:8], 'rad': line[8]})
            test_vert_atoms.append(line[1:5])
    # Go through the vertices and see if there are missing or extra vertices
    extra_verts = [my_check_verts[i] for i, _ in enumerate(test_vert_atoms) if _ not in good_vert_atoms]
    missing_verts = [my_good_verts[i] for i, _ in enumerate(good_vert_atoms) if _ not in test_vert_atoms]

    print("Missing Vertics: \n\n")
    for _ in missing_verts:
        print(_)
    print('\n Extra Vertices: \n\n')
    for _ in extra_verts:
        print(_)


if __name__ == '__main__':
    check_verts(complete_logs='/Users/jackericson/PycharmProjects/vorpy/Data/user_data/cambrin_Complete/sys/cambrin_logs.csv',
                test_logs='/Users/jackericson/PycharmProjects/vorpy/Data/user_data/cambrin_19/sys/cambrin_logs.csv')

