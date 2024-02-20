import matplotlib.pyplot as plt
from Data.Analyze.compare.compare_files import compare_files


if __name__ == '__main__':
    prefix = 'C:/Users/jacke/Documents/data/'
    my_info = compare_files(pdb_files=[prefix + '181L.pdb',
                                       prefix + '181L.pdb',
                                       prefix + '181L_coarse_ad.pdb',
                                       prefix + '181L_coarse_ad.pdb',
                                       prefix + '181L_coarse_ncap.pdb',
                                       prefix + '181L_coarse_ncap.pdb',
                                       prefix + '181L_coarse_scbb_ad.pdb',
                                       prefix + '181L_coarse_scbb_ad.pdb',
                                       prefix + '181L_coarse_scbb_ncap.pdb',
                                       prefix + '181L_coarse_scbb_ncap.pdb',
                                       prefix + '181L_martini.pdb',
                                       prefix + '181L_martini.pdb'],
                            log_files=[prefix + '181L_atom_vor_logs.csv',
                                       prefix + '181L_atom_pow_logs.csv',
                                       prefix + '181L_coarse_ad_logs.csv',
                                       prefix + '181L_coarse_ad_pow_logs.csv',
                                       prefix + '181L_coarse_ncap_vor_logs.csv',
                                       prefix + '181L_coarse_ncap_pow_logs.csv',
                                       prefix + '181L_coarse_scbb_ad_vor_logs.csv',
                                       prefix + '181L_coarse_scbb_ad_pow_logs.csv',
                                       prefix + '181L_coarse_scbb_ncap_vor_logs.csv',
                                       prefix + '181L_coarse_scbb_ncap_pow_logs.csv',
                                       prefix + '181L_martini_vor_logs.csv',
                                       prefix + '181L_martini_pow_logs.csv'], totals=True)

    # Sample data
    labels = ['Atomistic', 'Average Distance', 'Encapsulate', 'Side-Chain/Backbone AD', 'Side-Chain/Backbone Encap', 'CG Martini']
    data = [round(my_info['totals'][_]['vol']) for _ in my_info['totals']]  # Sample data for the first set
    data1 = data[::2]
    data2 = data[1::2]
    ymax = max(data)
    # Bar width
    bar_width = 0.35

    # Index for the x-axis
    x = range(len(labels))

    # Create the bar graph
    plt.bar(x, data1, width=bar_width, label='Additively Weighted')
    plt.bar([i + bar_width for i in x], data2, width=bar_width, label='Power')

    # Add labels and title
    plt.ylabel('Volume')
    plt.title('181L Volume by Scheme')

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
