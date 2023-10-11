import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import os
import yaml
import copy

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

if __name__ == "__main__":
    # Load data
    with open('./config.yml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    df_all = pd.read_csv(os.path.join(config['data_path'], "Broad_10M_Screen_Results_2022-07-12.csv"))
    df_CK1a_all_labels = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1a_all_labels.csv'))
    df_CK1d_all_labels = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1d_all_labels.csv'))
    df_CK1a = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1a_filtered.csv'))
    df_CK1a_inh = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1a_inh_filtered.csv'))
    df_CK1d = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1d_filtered.csv'))
    df_CK1d_inh = pd.read_csv(os.path.join(config['output_path'], 'stratified', 'CK1d_inh_filtered.csv'))

    # Prepare output path
    os.makedirs(os.path.join(config['output_path'],'lib_stat', 'all'), exist_ok=True)
    os.makedirs(os.path.join(config['output_path'],'lib_stat', 'CK1a'), exist_ok=True)
    os.makedirs(os.path.join(config['output_path'],'lib_stat', 'CK1d'), exist_ok=True)
    
    # Calculate quantile table for the whole library
    df_all_quantiles= get_quantiles(df_all, 'ZScore')
    df_all_quantiles.to_csv(os.path.join(config['output_path'], 'lib_stat', 'all', 'zscore_quantiles.csv'), index=False)
    
    # Plot Zscore of CK1a
    df_CK1a_orthosteric = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Orthosteric'])
    df_CK1a_allosteric = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Allosteric'])
    df_CK1a_cryptic = copy.deepcopy(df_CK1a_all_labels[df_CK1a_all_labels['customlabel'] == 'Cryptic'])
    df_CK1a_orthosteric['log 10 ZScore'] = np.log10(df_CK1a_orthosteric['ZScore_A'])
    df_CK1a_allosteric['log 10 ZScore'] = np.log10(np.mean(df_CK1a_allosteric[['ZScore_A', 'ZScore_A_inh']], axis=1))
    df_CK1a_cryptic['log 10 ZScore'] = np.log10(df_CK1a_cryptic['ZScore_A_inh'])
    df_effectsize_boxenplot = pd.concat([df_CK1a_orthosteric, df_CK1a_allosteric, df_CK1a_cryptic])
    sns.boxenplot(data=df_effectsize_boxenplot.sort_values('customlabel'), x='customlabel', y='log 10 ZScore', palette=palette)    
    plt.savefig(os.path.join(config['output_path'], 'lib_stat', 'CK1a', 'ZScore.png'), bbox_inches='tight')
    plt.clf()

    # Plot Zscore of CK1d
    df_CK1d_orthosteric = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Orthosteric'])
    df_CK1d_allosteric = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Allosteric'])
    df_CK1d_cryptic = copy.deepcopy(df_CK1d_all_labels[df_CK1d_all_labels['customlabel'] == 'Cryptic'])
    df_CK1d_orthosteric['log 10 ZScore'] = np.log10(df_CK1d_orthosteric['ZScore_D'])
    df_CK1d_allosteric['log 10 ZScore'] = np.log10(np.mean(df_CK1d_allosteric[['ZScore_D', 'ZScore_D_inh']], axis=1))
    df_CK1d_cryptic['log 10 ZScore'] = np.log10(df_CK1d_cryptic['ZScore_D_inh'])
    df_effectsize_boxenplot = pd.concat([df_CK1d_orthosteric, df_CK1d_allosteric, df_CK1d_cryptic])
    sns.boxenplot(data=df_effectsize_boxenplot.sort_values('customlabel'), x='customlabel', y='log 10 ZScore', palette=palette)    
    plt.savefig(os.path.join(config['output_path'], 'lib_stat', 'CK1d', 'ZScore.png'), bbox_inches='tight')
    plt.clf()