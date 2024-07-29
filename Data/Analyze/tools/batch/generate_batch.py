import os


# for j in range(20):
#     if j == 0:
#         extra = ''
#     else:
#         extra = '_' + str(j)
#     for i in range(1, 21):
#         print('python3 vorpy.py C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data/1.0_0.01_300_{}_False_physical1{}/1.0_0.01_300_{}_False_physical1.pdb -s nt compare -s mv 1000'.format(round(i*0.025, 5), extra, round(i*0.025, 5)))


# for i in range(18):
#     i_val = round((i+1)*0.025, 3)
#     for j in range(21):
#         j_val = round((j+6)*0.25, 3)
#         print('python3 vorpy.py C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data/1.0_{}_300_{}_False_gamma/1.0_{}_300_{}_False_gamma.pdb -s nt compare -s mv 1000'.format(j_val, i_val, j_val, i_val))
#
#
# for method in {'pow', 'del'}:
#     for i in range(1, 12):
#         with open('../../../../cambrin_frames.bat', 'a') as write_file:
#             write_file.write('py vorpy.py c{} -s nt {} -e dir C:/Users/jacke/GSU Dropbox Dropbox/John Ericson/GSU lab/Jack/Vorpy/Data/IV_Molecular/Cambrin/Frames/c_{}/atomistic\n'.format(i, method, i))



# main_dir = 'C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data'
# file_names = {}
# my_sub_dirs = [_[0] for _ in os.walk(main_dir)]
# for file_dir in my_sub_dirs:
#     my_info = file_dir.split("/")[-1]
#     dir_name = my_info[10:]
#     dir_info = dir_name.split('_')
#     if len(dir_info) >= 2 and dir_info[-1] == 'physical1':
#         if dir_name in file_names:
#             file_names[dir_name].append((dir_name, dir_info))
#         else:
#             file_names[dir_name] = [(dir_name, dir_info)]
#     if len(dir_info) >= 2 and dir_info[-2] == 'physical1':
#         my_dir_name = '_'.join(dir_info[:-1])
#         if my_dir_name in file_names:
#             file_names[my_dir_name].append((dir_name, dir_info))
#         else:
#             file_names[my_dir_name] = [(dir_name, dir_info)]
#
# # for file in file_names:
# #     if len(file_names[file]) < 10:
# #         for _ in range(10 - len(file_names[file])):
# #             print('python3 foam_gen.py ' + ' '.join(file_names[file][0][1]))
# #
# # for file in file_names:
# #     for ind_file in file_names[file]:
# #
# #         if ind_file[1][-1] != 'physical1':
# #             pdb_file = '_'.join(ind_file[1][:-1]) + ' physical1.pdb'
# #         else:
# #             pdb_file = ind_file[0] + ' physical1.pdb'
# #         print('python3 vorpy.py ' + main_dir + '/' + ind_file[0] + '/' + pdb_file + ' -s nt compare -s mv 1000')
#
# # extra = ''
# # for j in range(17):
# #     num1 = round((j+4)*0.025, 3)
# #     for k in range(20):
#         num2 = round((k+1)*0.025, 3)
#         print('python3 vorpy.py C:/Users/jacke/PycharmProjects/foam_gen/Data/user_data/1.0_{}_300_{}_False_lognormal'.format(num1, num2) + extra + '/1.0_{}_300_{}_False_lognormal.pdb -s nt compare -s mv 1000'.format(num1, num2))
# print("\n")

for i in range(1, 12):
    for net_type in ['vor', 'pow', 'del']:
        print('py vorpy.py hairpin_0{} -s nt {}'.format(i, net_type))

