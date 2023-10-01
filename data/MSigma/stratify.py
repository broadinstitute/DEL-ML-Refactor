import pandas as pd
import copy
import os
import yaml

def label_generation(df_exp, df_exp_inh, relevant_columns, exp_code):
    """
    df_exp (pd.DataFrame): Dataframe of experimental condition
    df_exp_inh (pd.DataFrame): Dataframe of experimental condition (inhibitor)
    relevant_columns (list): List of relevant columns
    exp_code (str): Experimental condition code (CK1a or CK1d)
    config (dict): Configuration dictionary
    """
    assert exp_code in ['CK1a', 'CK1d']
    result = {'Allosteric': None, 'Orthosteric': None, 'Cryptic': None}
    
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
    print(exp_code)
    for key, value in result.items():
        print(key, len(value))
    return result

if __name__ == "__main__":
    # Load data
    with open('./config.yml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    output_path_stratified = os.path.join(config['output_path'], 'stratified')
    if not os.path.exists(output_path_stratified):
        os.makedirs(output_path_stratified, exist_ok=True)
    relevant_columns = ['CompoundIndex','SMILES', 'Sample', 'ZScore']
    df_all = pd.read_csv(os.path.join(config['data_path'], "Broad_10M_Screen_Results_2022-07-12.csv"), usecols=relevant_columns)

    df_CK1a = copy.deepcopy(df_all.loc[df_all["Sample"] == "Sigma-A"])
    df_CK1a_inh = copy.deepcopy(df_all.loc[df_all["Sample"] == "Sigma-A-inh"])
    df_CK1a.rename(columns={'ZScore': 'ZScore_A'}, inplace=True)
    df_CK1a_inh.rename(columns={'ZScore': 'ZScore_A_inh'}, inplace=True)
    results_CK1a = label_generation(df_CK1a, df_CK1a_inh, ['SMILES'], 'CK1a')
    df_CK1a_all_labels = pd.concat([results_CK1a['Allosteric'], results_CK1a['Orthosteric'], results_CK1a['Cryptic']])
    df_CK1a_orthosteric = results_CK1a['Orthosteric']
    print("")

    df_CK1d = copy.deepcopy(df_all.loc[df_all["Sample"] == "Sigma-D"])
    df_CK1d_inh = copy.deepcopy(df_all.loc[df_all["Sample"] == "Sigma-D-inh"])
    df_CK1d.rename(columns={'ZScore': 'ZScore_D'}, inplace=True)
    df_CK1d_inh.rename(columns={'ZScore': 'ZScore_D_inh'}, inplace=True)
    results_CK1d = label_generation(df_CK1d, df_CK1d_inh, ['SMILES'], 'CK1d')
    df_CK1d_all_labels = pd.concat([results_CK1d['Allosteric'], results_CK1d['Orthosteric'], results_CK1d['Cryptic']])
    df_CK1d_orthosteric = results_CK1d['Orthosteric']

    # Following the naming convention of other two libraries (HitGen and DOS-DEL),
    # We append a filtered suffix to the file name but there is no filtering happens in MSigma
    df_CK1a.to_csv(os.path.join(output_path_stratified, 'CK1a_filtered.csv'), index=False)
    df_CK1a_inh.to_csv(os.path.join(output_path_stratified, 'CK1a_inh_filtered.csv'), index=False)
    df_CK1a_all_labels.to_csv(os.path.join(output_path_stratified, 'CK1a_all_labels.csv'), index=False)
    df_CK1a_orthosteric.to_csv(os.path.join(output_path_stratified, 'CK1a_orthosteric.csv'), index=False)

    df_CK1d.to_csv(os.path.join(output_path_stratified, 'CK1d_filtered.csv'), index=False)
    df_CK1d_inh.to_csv(os.path.join(output_path_stratified, 'CK1d_inh_filtered.csv'), index=False)
    df_CK1d_all_labels.to_csv(os.path.join(output_path_stratified, 'CK1d_all_labels.csv'), index=False)
    df_CK1d_orthosteric.to_csv(os.path.join(output_path_stratified, 'CK1d_orthosteric.csv'), index=False)