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
                                       prefix + '181L_martini_pow_logs.csv'])

