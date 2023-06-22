import yaml
import os
import numpy as np
import pandas as pd

def get_filtered_hc_df(df, hit_count_col_0, hit_count_col_1=None, thresh=0.0, blank_thresh=0.0,):
    """Filter data based on hit count threshold and blank hit count threshold
    Args:
        df (pd.DataFrame): Dataframe to be filtered
        hit_count_col_0 (str): Column name of hit count 0
        hit_count_col_1 (str): Column name of hit count 1
        thresh (float): Hit count threshold
        blank_thresh (float): Blank hit count threshold
    Returns:
        filt_df (pd.DataFrame): Filtered dataframe
        frac (float): Fraction of data that passed the filter
    """
    if hit_count_col_1:
        filt_df = df[(df[hit_count_col_0] > thresh) & (df[hit_count_col_1] > thresh) & 
                     (df['blank_hit_counts_0'] <= blank_thresh) & 
                     (df['blank_hit_counts_1'] <= blank_thresh)
                    ]
    else:
        filt_df = df[(df[hit_count_col_0] > thresh) & (df['blank_hit_counts_0'] <= blank_thresh)]
        
    print(f"{len(filt_df)} data ({len(filt_df) / len(df) * 100:.2f}%) passed the hit count filter with threshold ({thresh})")
    return filt_df, len(filt_df) / len(df)

def label_generation(df_exp, df_exp_inh, relevant_columns, exp_code):
    """Generate labels for each compound based on the experimental condition
    Args:
        df_exp (pd.DataFrame): Dataframe of experimental condition
        df_exp_inh (pd.DataFrame): Dataframe of experimental condition (inhibitor)
        relevant_columns (list): List of relevant columns
        exp_code (str): Experimental condition code
    Returns:
        df_common_exp_and_inh (pd.DataFrame): Dataframe of compounds that are common in both experimental conditions
        df_exclusive_exp_inh (pd.DataFrame): Dataframe of compounds that are exclusive in experimental condition (inhibitor)
        df_exclusive_exp (pd.DataFrame): Dataframe of compounds that are exclusive in experimental condition
    """
    assert exp_code in ['CK1a', 'CK1d']
    # Select relevant columns for df_exp and df_exp_inh
    df_exp = df_exp[relevant_columns]
    df_exp_inh = df_exp_inh[relevant_columns]

    # Merge df_exp and df_exp_inh based on relevant columns
    df_common_exp_and_inh = df_exp.merge(df_exp_inh, how='inner', on=relevant_columns)
    df_common_exp_and_inh['customlabel'] = 'Allosteric'
    df_common_exp_and_inh[f'{exp_code}_er_merge'] = np.mean(df_common_exp_and_inh[[f'{exp_code}_er', f'{exp_code}_inh_er']], axis=1)
    df_common_exp_and_inh[f'{exp_code}_er_lb_merge'] = np.mean(df_common_exp_and_inh[[f'{exp_code}_er_lb', f'{exp_code}_inh_er_lb']], axis=1)

    # Find exclusive rows in df_exp_inh
    df_exclusive_exp_inh = df_exp_inh[~df_exp_inh['SMILES'].isin(df_exp['SMILES'])].copy()
    df_exclusive_exp_inh['customlabel'] = 'Cryptic'
    #df_exclusive_exp_inh.rename(columns={f'{exp_code}_inh_er': 'Enrichment Ratio', f'{exp_code}_inh_er_lb': 'Enrichment Ratio LB'}, inplace=True)


    # Find exclusive rows in df_exp
    df_exclusive_exp = df_exp[~df_exp['SMILES'].isin(df_exp_inh['SMILES'])].copy()
    df_exclusive_exp['customlabel'] = 'Orthosteric'
    #df_exclusive_exp.rename(columns={f'{exp_code}_er': 'Enrichment Ratio', f'{exp_code}_er_lb': 'Enrichment Ratio LB'}, inplace=True)

    # Return the processed dataframes
    return df_common_exp_and_inh, df_exclusive_exp_inh, df_exclusive_exp


if __name__ == "__main__":
    # Load config and defined columns used in stratification
    with open('./config.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    thresh = config['positive_threshold']
    blank_thresh = config['blank_threshold']

    suffix = ['hit_counts_0', 'hit_counts_1', 'er', 'er_lb']
    blank_column = ['blank_hit_counts_0', 'blank_hit_counts_1']
    meta_column = ['CompoundIndex', 'SMILES', 'library.name']
    experiment_columns = [f'{exp_cond}_{sfx}' for exp_cond in config['experimental_condition'] for sfx in suffix]
    relevant_columns = experiment_columns + meta_column + blank_column

    # Load data
    print("Loading data...")
    df = pd.read_csv(os.path.join(config['output_path'],'preprocessed', 'preprocessed.csv'))

    # Filtering data with low hit counts
    print("Filtering data with low hit counts...")
    negative_condition = (df['blank_hit_counts_0'] > blank_thresh) & (df['blank_hit_counts_1'] > blank_thresh) &\
                         (df['CK1a_hit_counts_0'] <= thresh) & (df['CK1a_hit_counts_1'] <= thresh) &\
                         (df['CK1a_inh_hit_counts_0'] <= thresh) & (df['CK1a_inh_hit_counts_1'] <= thresh) &\
                         (df['CK1d_hit_counts_0'] <= thresh) & (df['CK1d_hit_counts_1'] <= thresh) &\
                         (df['CK1d_inh_hit_counts_0'] <= thresh) & (df['CK1d_inh_hit_counts_1'] <= thresh)
    df_negative = df[negative_condition]
    df_CK1a, frac_CK1a = get_filtered_hc_df(df, 'CK1a_hit_counts_0', "CK1a_hit_counts_1", thresh, blank_thresh)
    df_CK1a_inh, frac_CK1a_inh = get_filtered_hc_df(df, 'CK1a_inh_hit_counts_0', "CK1a_inh_hit_counts_1", thresh, blank_thresh)
    df_CK1d, frac_CK1d = get_filtered_hc_df(df, 'CK1d_hit_counts_0', "CK1d_hit_counts_1", thresh, blank_thresh)
    df_CK1d_inh, frac_CK1d_inh = get_filtered_hc_df(df, 'CK1d_inh_hit_counts_0', "CK1d_inh_hit_counts_1", thresh, blank_thresh)

    # Stratify data into three groups: Allosteric, Cryptic, Orthosteric
    print("Stratifying data...")
    df_common_CK1a_and_inh, df_exclusive_CK1a_inh, df_exclusive_CK1a = label_generation(df_CK1a, df_CK1a_inh, relevant_columns, exp_code='CK1a')
    df_common_CK1d_and_inh, df_exclusive_CK1d_inh, df_exclusive_CK1d = label_generation(df_CK1d, df_CK1d_inh, relevant_columns, exp_code='CK1d')
    df_CK1a_all_labels = pd.concat([df_exclusive_CK1a, df_exclusive_CK1a_inh,df_common_CK1a_and_inh])
    df_CK1d_all_labels = pd.concat([df_exclusive_CK1d, df_exclusive_CK1d_inh,df_common_CK1d_and_inh])

    # Write data to csv
    print("Writing data to csv...")
    output_path_stratified = os.path.join(config['output_path'], 'stratified')
    if not os.path.exists(output_path_stratified):
        os.makedirs(output_path_stratified)

    df_exclusive_CK1a.to_csv(os.path.join(output_path_stratified,'CK1a_orthosteric_153k.csv'), index=False)
    df_exclusive_CK1d.to_csv(os.path.join(output_path_stratified, 'CK1d_orthosteric_58k.csv'), index=False)
    df_negative.to_csv(os.path.join(output_path_stratified, 'negative_99k.csv'), index=False)
    df_CK1a_all_labels.to_csv(os.path.join(output_path_stratified, 'CK1a_all_labels.csv'), index=False)
    df_CK1d_all_labels.to_csv(os.path.join(output_path_stratified, 'CK1d_all_labels.csv'), index=False)
