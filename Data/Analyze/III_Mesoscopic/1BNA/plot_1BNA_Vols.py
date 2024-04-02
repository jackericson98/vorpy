import matplotlib.pyplot as plt
from Data.Analyze.tools.compare.compare_files import compare_files


if __name__ == '__main__':
    prefix = 'C:/Users/jacke/Documents/data1/'
    my_info = compare_files(pdb_files=[prefix + '1BNA.pdb',
                                       prefix + '1BNA.pdb',
                                       prefix + '1BNA_coarse_ad.pdb',
                                       prefix + '1BNA_coarse_ad.pdb',
                                       prefix + '1BNA_coarse_ncap.pdb',
                                       prefix + '1BNA_coarse_ncap.pdb',
                                       prefix + '1BNA_coarse_scbb_ad.pdb',
                                       prefix + '1BNA_coarse_scbb_ad.pdb',
                                       prefix + '1BNA_coarse_scbb_ncap.pdb',
                                       prefix + '1BNA_coarse_scbb_ncap.pdb'],
                            log_files=[prefix + '1BNA_atom_vor_logs.csv',
                                       prefix + '1BNA_atom_pow_logs.csv',
                                       prefix + '1BNA_coarse_ad_vor_logs.csv',
                                       prefix + '1BNA_coarse_ad_pow_logs.csv',
                                       prefix + '1BNA_coarse_ncap_vor_logs.csv',
                                       prefix + '1BNA_coarse_ncap_pow_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ad_vor_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ad_pow_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ncap_vor_logs.csv',
                                       prefix + '1BNA_coarse_scbb_ncap_pow_logs.csv'], totals=True)

    # Sample data
    labels = ['Atoms', 'Avg Dist', 'Encapsulate', 'SC/BB AD', 'SC/BB Encap.']
    data = [round(my_info['totals'][_]['vol'], 2) for _ in my_info['totals']]  # Sample data for the first set
    data1 = data[::2]
    data2 = data[1::2]
    ymax = max(data)
    # Bar width
    bar_width = 0.35

    # Index for the x-axis
    x = range(len(labels))

    # Create the bar graph
    plt.bar(x, data1, width=bar_width, label='Additively Weighted', color='skyblue', edgecolor='black')
    plt.bar([i + bar_width for i in x], data2, width=bar_width, label='Power', color='orange', edgecolor='black')

    # Add labels and title
    plt.ylabel('Volume', fontdict=dict(size=15))
    plt.title('1BNA Volume by Scheme', fontdict=dict(size=20))

    # Angle the labels and add values at the top of the bars
    plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right')
    for i, v in enumerate(data1):
        plt.text(i, v / 2, str(v) + ' \u212B\u00B3', ha='center', va='center', rotation=90)
    for i, v in enumerate(data2):
        plt.text(i + bar_width, v / 2, str(v) + ' \u212B\u00B3', ha='center', va='center', rotation=90)
    plt.ylim(0, 1.25 * ymax)
    # Add legend with appropriate layout
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.97), shadow=True, ncol=2)

    # Show the plot
    plt.tight_layout()
    plt.show()
