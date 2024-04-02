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
                                       prefix + '1BNA_coarse_scbb_ncap_pow_logs.csv'], totals=True, curv=True)

    # Sample data
    labels = ['Atoms', 'Avg Dist', 'Encapsulate', 'SC/BB AD', 'SC/BB Encap']
    data = [round(my_info['totals'][_]['avg curv'], 3) for _ in my_info['totals']]  # Sample data for the first set
    data1 = data[::2]
    data2 = data[1::2]
    max_height = max(data1)
    # Bar width
    bar_width = 0.35

    # Index for the x-axis
    x = range(len(labels))

    # Create the bar graph
    plt.bar(x, data1, width=bar_width, label='Additively Weighted')
    # plt.bar([i + bar_width for i in x], data2, width=bar_width, label='Power')

    # Add labels and title
    plt.ylabel('Average Curvature (Gaussian)')
    plt.title('1BNA Average Curvature by Scheme')

    # Angle the labels and add values at the top of the bars
    plt.xticks([i + bar_width / 2 for i in x], labels, rotation=45, ha='right')
    for i, v in enumerate(data1):
        plt.text(i, max_height / 2, str(v), ha='center', va='center', rotation=90)
    # for i, v in enumerate(data2):
    #     plt.text(i + bar_width, v / 2, str(v), ha='center', va='center', rotation=90)

    # Add legend
    # plt.legend()

    # Show the plot
    plt.tight_layout()
    plt.show()
