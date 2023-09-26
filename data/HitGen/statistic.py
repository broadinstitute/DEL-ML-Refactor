import seaborn as sns
import copy
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import yaml
import os
palette = {"Allosteric":"tab:blue",
           "Orthosteric":"tab:green",
           "Cryptic":"tab:purple"}
def get_quantiles(df, col):
    """Get the quantiles of the dataframe for the columns of interest
    Args:
        df (pd.DataFrame): Dataframe contains effect size and counts
        col (str): Column of interest
    Returns:
        pd.DataFrame: Dataframe of quantiles for the column of interest
    """
    return df[col].quantile([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])

def get_all_quantiles(df, condition, count_thresh, blank_thresh):
    """Get the quantiles for each experiment condition in different filtering approach
    Args:
        df (pd.DataFrame): Dataframe contains effect size and counts
        condition (str): Experiment condition (CK1a, CK1a-inh, CK1d, CK1d-inh)
        count_thresh (int): Count threshold
        blank_thresh (int): Blank effect size threshold
    Returns:
        dict: Dictionary of quantiles for each experiment condition in different filtering approach
    """
    results = {'HC quantile':None , 'ES prefilter quantile':None, 'ES postfilter quantile':None}
    quantile_count_col = {'CK1a': 'A counts', 'CK1a-inh': 'A-inh counts', 'CK1d': 'D counts', 'CK1d-inh': 'D-inh counts'}
    quantile_effect_col = {'CK1a': 'A Effect Size', 'CK1a-inh': 'A-inh Effect Size', 'CK1d': 'D Effect Size', 'CK1d-inh': 'D-inh Effect Size'}
    df_tmp = df[df[quantile_count_col[condition]] > count_thresh]
    df_tmp_hc_quantile = get_quantiles(df_tmp, [quantile_count_col[condition]])
    results['HC quantile'] = df_tmp_hc_quantile

    df_tmp_es_prefilter_quantile = get_quantiles(df_tmp, [quantile_effect_col[condition]])
    results['ES prefilter quantile'] = df_tmp_es_prefilter_quantile

    df_tmp = df_tmp[df_tmp['blank Effect Size'] <= blank_thresh]
    df_tmp_es_postfilter_quantile = get_quantiles(df_tmp, [quantile_effect_col[condition]])
    results['ES postfilter quantile'] = df_tmp_es_postfilter_quantile

    return results

def get_sublibary_stat(sublib_of_interest, sublib_mol_count, experiment_conditoin):
    """Get the statistics of sublibrary of interest
    Args:
        sublib_of_interest (pd.DataFrame): Dataframe of sublibrary of interest
        sublib_mol_count (dict): Dictionary of sublibrary and its corresponding molecule count
        experiment_conditoin (str): Experimental condition
    Returns:
        pd.DataFrame: Dataframe of sublibrary statistics
    """
    sublib_mol_count_interest = sublib_of_interest.groupby('library.name').count()['SMILES'].to_dict()
    output_dict = {'library.name':[], 'fraction':[] ,'pass_filter_rate':[], 'experiment_condition':[]}
    for sublib,  mol_count in sublib_mol_count.items():
        output_dict['library.name'].append(sublib)
        output_dict['experiment_condition'].append(experiment_conditoin)
        if sublib not in sublib_mol_count_interest:
            output_dict['fraction'].append(0)
            output_dict['pass_filter_rate'].append(0)
        else:
            output_dict['pass_filter_rate'].append(sublib_mol_count_interest[sublib]/mol_count)
            output_dict['fraction'].append(sublib_mol_count_interest[sublib]/len(sublib_of_interest))
    return pd.DataFrame(output_dict)

def plot_effect_size_whole_lib(df, effect_size_columns, output_path, effect_size_cutoff=100):
    """Plot the effect size for the whole library
    Args:
        df (pd.DataFrame): Dataframe contains effect size and counts
        effect_size_columns (list): List of effect size columns
        output_path (str): Output path
        effect_size_cutoff (int, optional): Effect size cutoff. Defaults to 100.
    """
    if not os.path.exists(os.path.join(output_path, 'lib_stat', 'all')):
        os.makedirs(os.path.join(output_path, 'lib_stat', 'all'), exist_ok=True)
    save_path = os.path.join(output_path, 'lib_stat', 'all')
    print("Plotting effect size boxenplot...")
    sns.boxenplot(data=df[effect_size_columns])
    plt.ylabel('log Effect Size')
    plt.yscale('symlog')
    plt.xticks(rotation=90)
    plt.savefig(os.path.join(save_path, 'effect_size_boxen.png'), bbox_inches='tight')
    plt.clf()
    
    print(f"Plotting effect size distribution with {effect_size_cutoff} cutoff...")
    for es in effect_size_columns:
        sns.kdeplot(df[df[es] <= effect_size_cutoff][es], label=es)
    plt.legend()
    plt.savefig(os.path.join(save_path, f'effect_size_dist_cutoff_{effect_size_cutoff}.png'), bbox_inches='tight')
    plt.clf()

    print("Plotting effect size diff boxenplot...")
    effect_size_diff_columns = [f"{col} Diff" for col in effect_size_columns] 
    for col in effect_size_columns:
        df[f'{col} Diff'] = df[col] - df['blank Effect Size']

    sns.boxenplot(df[effect_size_diff_columns])
    plt.yscale('symlog')
    plt.ylabel('log Effect Size Diff from the blank')
    plt.xticks(rotation=90)
    plt.savefig(os.path.join(save_path, 'effect_size_diff_boxen.png'), bbox_inches='tight')
    plt.clf()

    print(f"Plotting effect size diff distribution with {effect_size_cutoff} cutoff...")
    for es in effect_size_diff_columns:
        if es == 'blank Effect Size Diff':
            continue
        sns.kdeplot(df[(df[es] <= effect_size_cutoff) & (df[es] >= - effect_size_cutoff)][es], label=es)
    plt.legend()
    plt.savefig(os.path.join(save_path, f'effect_size_diff_dist_cutoff_{effect_size_cutoff}.png'), bbox_inches='tight')
    plt.clf()


if __name__ == "__main__":
    # Load data
    print("Loading data...\n")
    with open('./config.yml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    molecular_enrichments_df = pd.read_csv(os.path.join(config['data_path'], 'merged_molecule_enrichments.csv'))
    df_CK1a_all_labels = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1a_all_labels.csv'),
                                     usecols=['A Effect Size', 'A-inh Effect Size','library.name','customlabel','SMILES'])
    df_CK1d_all_labels = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1d_all_labels.csv'),
                                     usecols=['D Effect Size', 'D-inh Effect Size','library.name','customlabel','SMILES'])
    df_CK1a = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1a_filtered.csv'), usecols=['SMILES','library.name'])
    df_CK1a_inh = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1a_inh_filtered.csv'), usecols=['SMILES','library.name'])
    df_CK1d = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1d_filtered.csv'), usecols=['SMILES','library.name'])
    df_CK1d_inh = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1d_inh_filtered.csv'), usecols=['SMILES','library.name'])

    # Rename and define the columns used for analysis
    molecular_enrichments_df.rename(columns={'product_smiles': 'SMILES', 'mol.id':'CompoundIndex'}, inplace=True)
    effect_size_columns = [f'{x} Effect Size' for x in ['A', 'A-inh', 'D', 'D-inh', 'blank']]
    counts_columns = [f'{l} counts' for l in ['A', 'A-inh', 'D', 'D-inh', 'blank']]
    relevant_columns = ['CompoundIndex', 'SMILES', 'library.name'] + effect_size_columns + counts_columns
    count_thresh = config['count_threshold']
    blank_effect_size_thresh = config['blank_effect_size_threshold']
    df_CK1a_orthosteric = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Orthosteric'])
    df_CK1d_orthosteric = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Orthosteric'])

    # Sublibrary count analysis for each experiment condition
    # The goal of this part is to get the fraction rate of each sublibrary 
    # in each experiment condition after filtering
    print("Calculating sublibrary count statistics...\n")
    sublib_mol_count = molecular_enrichments_df.groupby('library.name').count()['SMILES'].to_dict()
    experiments = [df_CK1a_orthosteric, df_CK1d_orthosteric, df_CK1a, df_CK1d, df_CK1a_inh, df_CK1d_inh, molecular_enrichments_df]
    experiment_conditions = ['CK1a_orthosteric','CK1d_orthosteric', 'CK1a', 'CK1d', 'CK1a_inh','CK1d_inh', 'all']
    for df_exp, exp_cond in zip(experiments, experiment_conditions):
        print(f"Calculating sublibrary statistics for {exp_cond}")

        if 'orthosteric' in exp_cond:
            os.makedirs(os.path.join(config['output_path'], 'lib_stat', exp_cond.split("_")[0]), exist_ok=True)
        else:
            os.makedirs(os.path.join(config['output_path'], 'lib_stat', exp_cond), exist_ok=True)
        sublibrary_stat = get_sublibary_stat(df_exp, sublib_mol_count, exp_cond)
        sublibrary_frac_bar = sns.catplot(data=sublibrary_stat,
                                        x='library.name',
                                        y='fraction',
                                        kind='bar',
                                        col='experiment_condition')
        plt.xticks(rotation=90)
        sublibrary_pass_bar = sns.catplot(data=sublibrary_stat,
                                        x='library.name',
                                        y='pass_filter_rate',
                                        kind='bar',
                                        col='experiment_condition')
        plt.xticks(rotation=90)
        if 'orthosteric' in exp_cond:
            sublibrary_frac_bar.savefig(os.path.join(config['output_path'], 'lib_stat', exp_cond.split("_")[0], 'orthosteric_sublibrary_fraction.png'))
            sublibrary_pass_bar.savefig(os.path.join(config['output_path'], 'lib_stat', exp_cond.split("_")[0], 'orthosteric_sublibrary_pass_filter_rate.png'))
        else:
            sublibrary_frac_bar.savefig(os.path.join(config['output_path'], 'lib_stat', exp_cond, 'sublibrary_fraction.png'))
            sublibrary_pass_bar.savefig(os.path.join(config['output_path'], 'lib_stat', exp_cond, 'sublibrary_pass_filter_rate.png'))
        plt.close()

    # Quantile analysis for each experiment condition in different filtering approach
    # We analyze the quantile by count and effect size threshold for each experiment condition
    # then generate a merged quantile table
    print("")
    print("Generating quantile tables...\n")
    # Step 1: Get quantiles for each experiment condition in different filtering approach
    df_CK1a_quantiles = get_all_quantiles(molecular_enrichments_df, 'CK1a', count_thresh, blank_effect_size_thresh)
    df_CK1a_inh_quantiles = get_all_quantiles(molecular_enrichments_df, 'CK1a-inh', count_thresh, blank_effect_size_thresh)
    df_CK1d_quantiles = get_all_quantiles(molecular_enrichments_df, 'CK1d', count_thresh, blank_effect_size_thresh)
    df_CK1d_inh_quantiles = get_all_quantiles(molecular_enrichments_df, 'CK1d-inh', count_thresh, blank_effect_size_thresh)

    df_blank = molecular_enrichments_df[molecular_enrichments_df['blank counts'] > count_thresh]
    df_blank_hc_quantiles = get_quantiles(df_blank, ['blank counts'])
    df_blank_es_prefilter_quantiles = get_quantiles(df_blank, ['blank Effect Size'])

    # Step 2: Merge the quantile tables of each experiment condition then output
    es_prefilter_quantiles_df = pd.concat([df_CK1a_quantiles['ES prefilter quantile'],\
                                           df_CK1a_inh_quantiles['ES prefilter quantile'],\
                                           df_CK1d_quantiles['ES prefilter quantile'],\
                                           df_CK1d_inh_quantiles['ES prefilter quantile'],\
                                           df_blank_es_prefilter_quantiles['blank Effect Size']], axis=1)

    es_quantiles_df = pd.concat([df_CK1a_quantiles['ES postfilter quantile'],\
                                 df_CK1a_inh_quantiles['ES postfilter quantile'],\
                                 df_CK1d_quantiles['ES postfilter quantile'],\
                                 df_CK1d_inh_quantiles['ES postfilter quantile']], axis=1)

    hc_quantiles_df = pd.concat([df_CK1a_quantiles['HC quantile'],\
                                 df_CK1a_inh_quantiles['HC quantile'],\
                                 df_CK1d_quantiles['HC quantile'],\
                                 df_CK1d_inh_quantiles['HC quantile'],\
                                 df_blank_hc_quantiles['blank counts']], axis=1)
    es_prefilter_quantiles_df.to_csv(os.path.join(config['output_path'], 'lib_stat', 'all', 'es_prefilter_quantiles.csv'), index=False)
    es_quantiles_df.to_csv(os.path.join(config['output_path'], 'lib_stat', 'all', 'es_quantiles.csv'), index=False)
    hc_quantiles_df.to_csv(os.path.join(config['output_path'], 'lib_stat', 'all', 'hc_quantiles.csv'), index=False)

    # Plot effect size for the whole library and each experiment condition
    # Whole library effectt size
    print("Plotting effect size for the whole library...\n")
    plot_effect_size_whole_lib(molecular_enrichments_df, effect_size_columns, config['output_path'])
    print("Plotting effect size for each experiment condition...\n")
    # Get effect size for each experiment condition in CK1a and transform to log scale
    df_CK1a_orthosteric = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Orthosteric'])
    df_CK1a_allosteric = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Allosteric'])
    df_CK1a_cryptic = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Cryptic'])
    df_CK1a_orthosteric['log 10 Effect Size'] = np.log10(df_CK1a_orthosteric['A Effect Size'])
    df_CK1a_allosteric['log 10 Effect Size'] = np.log10(np.mean(df_CK1a_allosteric[['A Effect Size', 'A-inh Effect Size']], axis=1))
    df_CK1a_cryptic['log 10 Effect Size'] = np.log10(df_CK1a_cryptic['A-inh Effect Size'])
    df_effectsize_boxenplot = pd.concat([df_CK1a_orthosteric, df_CK1a_allosteric, df_CK1a_cryptic])

    # Plot and save
    sns.boxenplot(data=df_effectsize_boxenplot.sort_values('customlabel'), x='customlabel', y='log 10 Effect Size', palette=palette)
    plt.savefig(os.path.join(config['output_path'], 'lib_stat', 'CK1a', 'effective_size.png'), bbox_inches='tight')
    plt.clf()
    sns.boxenplot(data=df_CK1a_orthosteric.sort_values('library.name'), x="library.name", y="log 10 Effect Size")
    plt.xticks(rotation=90)
    plt.savefig(os.path.join(config['output_path'], 'lib_stat', 'CK1a', 'orthosteric_effective_size_perlib.png'), bbox_inches='tight')
    plt.clf()
    # Get effect size for each experiment condition in CK1d and transform to log scale
    df_CK1d_orthosteric = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Orthosteric'])
    df_CK1d_allosteric = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Allosteric'])
    df_CK1d_cryptic = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Cryptic'])
    df_CK1d_orthosteric['log 10 Effect Size'] = np.log10(df_CK1d_orthosteric['D Effect Size'])
    df_CK1d_allosteric['log 10 Effect Size'] = np.log10(np.mean(df_CK1d_allosteric[['D Effect Size', 'D-inh Effect Size']], axis=1))
    df_CK1d_cryptic['log 10 Effect Size'] = np.log10(df_CK1d_cryptic['D-inh Effect Size'])
    df_effectsize_boxenplot = pd.concat([df_CK1d_orthosteric, df_CK1d_allosteric, df_CK1d_cryptic])

    # Plot and save
    sns.boxenplot(data=df_CK1d_orthosteric.sort_values('library.name'), x="library.name", y="log 10 Effect Size")
    plt.xticks(rotation=90)
    plt.savefig(os.path.join(config['output_path'], 'lib_stat', 'CK1d', 'orthosteric_effective_size_perlib.png'), bbox_inches='tight')
    plt.clf()
    sns.boxenplot(data=df_effectsize_boxenplot.sort_values('customlabel'), x='customlabel', y='log 10 Effect Size', palette=palette)
    plt.savefig(os.path.join(config['output_path'], 'lib_stat', 'CK1d', 'effective_size.png'), bbox_inches='tight')
    plt.clf()