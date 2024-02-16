import os

rootdir ='C:/Users/i7-8700/PycharmProjects/foam_gen/Data/user_data'


my_file_types = {}
for subdir, dirs, files in os.walk(rootdir):
    for dir_ in dirs:
        dir_info = dir_.split('_')
        try:
            my_num = int(dir_info[-1])
            dir_string = '_'.join(dir_info[:-1])
        except ValueError:
            dir_string = '_'.join(dir_info)

        if dir_string in my_file_types:
            my_file_types[dir_string] += 1
        else:
            my_file_types[dir_string] = 1
count = 0
reverse_list = []
for _ in my_file_types:
    for i in range(20 - my_file_types[_]):
        count += 1
        poopy = _.split('_')
        print('python3 foam_gen.py', *poopy)
        reverse_list.insert(0, ['python3 foam_gen.py', *poopy])
# for _ in reverse_list:
#     print(*_)
print("{} Runs Left".format(count))


cvs = {}
for _ in my_file_types:
    name = _.split('_')
    cv = name[1]
    dens = name[3]
    if cv in cvs:
        if dens in cvs[cv]:
            cvs[cv][dens] += my_file_types[_]
        else:
            cvs[cv][dens] = my_file_types[_]
    else:
        cvs[cv] = {dens: my_file_types[_]}

print()
