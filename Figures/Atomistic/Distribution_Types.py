import matplotlib.pyplot as plt
import numpy as np


names =        ['Na5',  'EDTA_Mg',  'DB1976',     'Hairpin', 'Cambrin', '1BNA',    'Hammerhead',       'p53tet',  'NCP']
types =        ['Ion',  'Molecule', 'Molecule',   'Nucleic', 'Protein', 'Nucleic', 'Nucleic',          'Protein', 'Complex']
poly_cnt =     [1,       32,        41,       570,       674,       773,       2032,      2088,               25086]
ion_cnt =      [0,       1,         10,       0,         11,        0,         30,        8,                  62]
tot_cnt =      [61,      612,       4440,     3518,      4897,      2965,      52879,     26069,              92704]
poly_res_cnt = [1,       2,         1,        46,        18,        24,        20,        355,                1262]
sol_res_cnt =  [20,      192,       1463,     948,       1439,      730,       16939,     7991,               22539]


# Sort indices based on poly_cnt
sort_indices = np.argsort(poly_cnt)

# Sort data based on poly_cnt
sorted_names = [names[i] for i in sort_indices]
sorted_types = [types[i] for i in sort_indices]
sorted_poly_cnt = [poly_cnt[i] for i in sort_indices]

# Color mapping
colors = {'Protein': 'blue', 'Nucleic': 'green', 'Complex': 'orange', 'Ion': 'black', 'Molecule': 'grey'}

# Plotting
fig, ax = plt.subplots()

# Scatter plot with colored points based on atom type
for i, (name, t, count) in enumerate(zip(sorted_names, sorted_types, sorted_poly_cnt)):
    ax.scatter(i, count, label=name, color=colors[t], s=50)

# Log scale on the y-axis
ax.set_yscale('log')

# Labeling
ax.set_xticks(np.arange(len(sorted_names)))
ax.set_xticklabels(sorted_names)
ax.set_xlabel('Names')
ax.set_ylabel('Atom Counts (log scale)')

# Legend with unique labels for each type
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors[t], markersize=10, label=t) for t in set(types)]
ax.legend(handles=legend_elements, title='Atom Type', loc='upper left')

plt.title('Atom Counts in Different Atom Files')
plt.show()


