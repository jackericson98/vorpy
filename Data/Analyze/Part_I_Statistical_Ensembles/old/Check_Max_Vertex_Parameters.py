import csv
from Data.Analyze.tools.batch.compile_logs import get_logs_and_pdbs

data = {}

loggan_paul = get_logs_and_pdbs(False)
for directory in loggan_paul:
    if 'aw' not in loggan_paul[directory] or 'pow' not in loggan_paul[directory]:
        continue
    vals = directory.split('_')
    cv, den = float(vals[1]), float(vals[3])
    # Open the aw logs and get the maximum vertex
    with open(loggan_paul[directory]['aw'], 'r') as aw_logs:
        aw_read = csv.reader(aw_logs)
        reading = False
        vert_rads = []
        for line in aw_read:
            if reading:
                vert_rads.append(float(line[8]))
            if len(line) == 9 and line[8] == 'r':
                reading = True
            else:
                reading = False
    # Record the maximum vertex for aw
    max_aw = max(vert_rads)
    len_aw_rads = len(vert_rads)
    # Open the aw logs and get the maximum vertex
    with open(loggan_paul[directory]['pow'], 'r') as pow_logs:
        pow_read = csv.reader(pow_logs)
        reading = False
        vert_rads = []
        for line in pow_read:
            if reading:
                vert_rads.append(float(line[8]))
            if len(line) == 9 and line[8] == 'r':
                reading = True
            else:
                reading = False
    # Record the maximum vertex for aw
    max_pow = max(vert_rads)
    len_pow_rads = len(vert_rads)

    if (cv, den) in data:
        data[(cv, den)]['pow'].append(max_pow)
        data[(cv, den)]['aw'].append(max_aw)
    else:
        data[(cv, den)] = {'pow': [max_pow], 'aw': [max_aw], 'pow_verts': len_pow_rads}

for _ in data:
    print(_, 'pow', data[_]['pow'])
    print(_, 'aw', data[_]['aw'])
