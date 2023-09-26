import numpy as np
import pandas as pd
import copy
import os
import yaml
def filtering(df, condition, count_thresh, blank_thresh, relevant_columns):
    filtering_col = {'CK1a': 'A counts', 'CK1a-inh': 'A-inh counts', 'CK1d': 'D counts', 'CK1d-inh': 'D-inh counts'}
    assert condition in filtering_col.keys(), 'condition must be one of {}'.format(filtering_col.keys())

    df_tmp = df[df[filtering_col[condition]] > count_thresh]
    df_tmp = df_tmp[df_tmp['blank Effect Size'] <= blank_thresh]
    return df_tmp[relevant_columns]

def label_generation(df_exp, df_exp_inh, relevant_columns, exp_code, config):
    """
    df_exp (pd.DataFrame): Dataframe of experimental condition
    df_exp_inh (pd.DataFrame): Dataframe of experimental condition (inhibitor)
    relevant_columns (list): List of relevant columns
    exp_code (str): Experimental condition code (CK1a or CK1d)
    config (dict): Configuration dictionary
    """
    assert exp_code in ['CK1a', 'CK1d']
    result = {'Allosteric': None, 'Orthosteric': None, 'Cryptic': None, 
              'Allosteric_competitive': None, 'Allosteric_non_competitive': None}
    competive_measure = config[f'competitive_measure_{exp_code}']
    
    # Stratify to allosteric
    df_common_exp_and_inh = df_exp.merge(df_exp_inh, how='inner', on=relevant_columns)
    df_common_exp_and_inh['customlabel'] = 'Allosteric'
    result['Allosteric'] = df_common_exp_and_inh

    # Stratify to orthosteric
    df_exclusive_exp = copy.deepcopy(df_exp[~df_exp.SMILES.isin(df_exp_inh.SMILES)])
    df_exclusive_exp['customlabel'] = 'Orthosteric'
    result['Orthosteric'] = df_exclusive_exp

    # Stratify to cryptic
    df_exclusive_exp_inh = copy.deepcopy(df_exp_inh[~df_exp_inh.SMILES.isin(df_exp.SMILES)])
    df_exclusive_exp_inh['customlabel'] = 'Cryptic'
    result['Cryptic'] = df_exclusive_exp_inh

    # Stratify to allosteric competitive and non-competitive
    if exp_code == 'CK1a':
        df_common_exp_and_inh['competitive'] = np.where((df_common_exp_and_inh['A Effect Size'] 
                                                        > df_common_exp_and_inh['A-inh Effect Size'] * competive_measure), 
                                                        f'{exp_code}_common_competitive_to_inh', f'{exp_code}_common_non_competitive_to_inh')
    elif exp_code == 'CK1d':
        df_common_exp_and_inh['competitive'] = np.where((df_common_exp_and_inh['D Effect Size'] 
                                                        > df_common_exp_and_inh['D-inh Effect Size'] * competive_measure), 
                                                        f'{exp_code}_common_competitive_to_inh', f'{exp_code}_common_non_competitive_to_inh')
    df_common_competitive_exp_to_inh = df_common_exp_and_inh.loc[df_common_exp_and_inh['competitive'] == f'{exp_code}_common_competitive_to_inh']
    df_common_non_competitive_exp_to_inh = df_common_exp_and_inh.loc[df_common_exp_and_inh['competitive'] == f'{exp_code}_common_non_competitive_to_inh', ]
    df_common_competitive_exp_to_inh.reset_index(drop=True, inplace=True)
    df_common_non_competitive_exp_to_inh.reset_index(drop=True, inplace=True)
    result['Allosteric_competitive'] = df_common_competitive_exp_to_inh
    result['Allosteric_non_competitive'] = df_common_non_competitive_exp_to_inh

    print(exp_code)
    for key, value in result.items():
        print(key, len(value))
    return result

if __name__ == '__main__':
    with open('./config.yml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    # Loaind data
    print('Loading data...')
    molecular_enrichments_df = pd.read_csv(os.path.join(config['data_path'], 'merged_molecule_enrichments.csv'))
    molecular_enrichments_df.rename(columns={'product_smiles': 'SMILES', 'mol.id':'CompoundIndex'}, inplace=True)
    # Define columns for analysis and threshold for filtering
    effect_size_columns = [f'{x} Effect Size' for x in ['A', 'A-inh', 'D', 'D-inh', 'blank']]
    counts_columns = [f'{l} counts' for l in ['A', 'A-inh', 'D', 'D-inh', 'blank']]
    relevant_columns = ['CompoundIndex', 'SMILES', 'library.name'] + effect_size_columns + counts_columns
    count_thresh = config['count_threshold']
    blank_effect_size_thresh = config['blank_effect_size_threshold']
    print("")
    # Filter data
    print('Filtering data...')
    df_CK1a = filtering(molecular_enrichments_df, 'CK1a', count_thresh, blank_effect_size_thresh, relevant_columns)
    df_CK1a_inh = filtering(molecular_enrichments_df, 'CK1a-inh', count_thresh, blank_effect_size_thresh, relevant_columns)
    df_CK1d = filtering(molecular_enrichments_df, 'CK1d', count_thresh, blank_effect_size_thresh, relevant_columns)
    df_CK1d_inh = filtering(molecular_enrichments_df, 'CK1d-inh', count_thresh, blank_effect_size_thresh, relevant_columns)
    print("")
    # Stratify data into three major categories: Allosteric, Orthosteric, and Cryptic
    # The purpose of CK1x exclusive competitive is unclear. feel free to leave them out
    print('Stratifying data...')
    results_CK1a = label_generation(df_CK1a, df_CK1a_inh, relevant_columns, 'CK1a', config)
    df_CK1a_all_labels = pd.concat([results_CK1a['Allosteric'], results_CK1a['Orthosteric'], results_CK1a['Cryptic']])
    df_CK1a_orthosteric = results_CK1a['Orthosteric']
    df_CK1a_exclusive_competitive = pd.concat([results_CK1a['Allosteric'], results_CK1a['Allosteric_competitive']])
    df_CK1a_exclusive_competitive.reset_index(drop=True, inplace=True)
    print("")
    results_CK1d = label_generation(df_CK1d, df_CK1d_inh, relevant_columns, 'CK1d', config)
    df_CK1d_all_labels = pd.concat([results_CK1d['Allosteric'], results_CK1d['Orthosteric'], results_CK1d['Cryptic']])
    df_CK1d_orthosteric = results_CK1d['Orthosteric']
    df_CK1d_exclusive_competitive = pd.concat([results_CK1d['Allosteric'], results_CK1d['Allosteric_competitive']])
    df_CK1d_exclusive_competitive.reset_index(drop=True, inplace=True)
    print("")
    print("Writing data to csv...")
    output_path_stratified = os.path.join(config['output_path'], 'stratified')
    if not os.path.exists(output_path_stratified):
        os.makedirs(output_path_stratified)

    df_CK1a.to_csv(os.path.join(output_path_stratified, 'CK1a_filtered.csv'), index=False)
    df_CK1a_inh.to_csv(os.path.join(output_path_stratified, 'CK1a_inh_filtered.csv'), index=False)
    df_CK1a_all_labels.to_csv(os.path.join(output_path_stratified, 'CK1a_all_labels.csv'), index=False)
    df_CK1a_orthosteric.to_csv(os.path.join(output_path_stratified, 'CK1a_orthosteric.csv'), index=False)
    df_CK1a_exclusive_competitive.to_csv(os.path.join(output_path_stratified, 'CK1a_exclusive_competitive.csv'), index=False)

    df_CK1d.to_csv(os.path.join(output_path_stratified, 'CK1d_filtered.csv'), index=False)
    df_CK1d_inh.to_csv(os.path.join(output_path_stratified, 'CK1d_inh_filtered.csv'), index=False)
    df_CK1d_all_labels.to_csv(os.path.join(output_path_stratified, 'CK1d_all_labels.csv'), index=False)
    df_CK1d_orthosteric.to_csv(os.path.join(output_path_stratified, 'CK1d_orthosteric.csv'), index=False)
    df_CK1d_exclusive_competitive.to_csv(os.path.join(output_path_stratified, 'CK1d_exclusive_competitive.csv'), index=False)