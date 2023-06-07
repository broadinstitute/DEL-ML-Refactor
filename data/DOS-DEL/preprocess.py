import os
import pandas as pd
import yaml
from tqdm import tqdm

def pad(st):
    """Pad a string with zeros to make it 3 digits long"""
    return '0'*(3-len(str(st))) + str(st)

def load_library(config):
    """Load the library of structures. Each (sub)library contain SMIELS string of the molecules."""
    lib_dfs = []
    for lib in tqdm(config['sublibrary']):
        lib_path = os.path.join(config['data_path'], "libraries", f"lib{lib}.csv")
        lib_df = pd.read_csv(lib_path)
        lib_df['sublibrary'] = lib
        lib_dfs.append(lib_df)

    lib_dfs = pd.concat(lib_dfs)
    lib_dfs.rename(columns={'structure': 'SMILES'}, inplace=True)
    #lib_dfs.to_csv(os.path.join(output_path_preprocessed, "lib_dfs.csv"), index=False)
    return lib_dfs


def load_sample(config):
    """Load the sample of structures. Each sample file contains the screening results (counts) of each sublibrary 
    under a specific experimental condition. See the config file for the mapping between sample number and
    experimental condition.
    """
    sample_dfs = []
    sample_to_exp_condition = config['sample_to_exp_condition']
    sample_list = list(sample_to_exp_condition.keys())
    for sample in tqdm(range(min(sample_list), max(sample_list) + 1)):
        tmp_dfs = []
        for lib in config['sublibrary']:
            sample_path_per_lib = os.path.join(config['data_path'], 'samples', f'run038_samp000{sample}_lib{lib}.csv')
            sample_df_per_lib = pd.read_csv(sample_path_per_lib)
            tmp_dfs.append(sample_df_per_lib)
        
        sample_df = pd.concat(tmp_dfs)
        sample_df['experimental_cond'] = sample_to_exp_condition[sample]
        sample_df['sample'] = sample
        sample_df['run'] = sample % 2
        
        sample_dfs.append(sample_df)
        
    all_sample_df = pd.concat(sample_dfs)
    return all_sample_df

def merge_sample_and_experiment(config, all_sample_df, output_path):
    """Merge the sample and experiment results"""
    tmp = []
    for exp_cond in tqdm(config['experimental_condition']):
        # Retrieve the sample results for the current experiment condition
        sample_df = all_sample_df[all_sample_df['experimental_cond'] == exp_cond]
        sample_df_0 = sample_df[sample_df['run'] == 0].reset_index(drop=True)
        sample_df_1 = sample_df[sample_df['run'] == 1].reset_index(drop=True)

        # Retrieve the screening results for the current experiment condition across all sublibraries
        dfs = []
        for lib in config['sublibrary']:
            lib_exp_cond_info = pd.read_csv(os.path.join(config['data_path'], f"{exp_cond}_lib{lib}.csv"))
            dfs.append(lib_exp_cond_info)
        exp_cond_df = pd.concat(dfs)

        # Naming the compound by concatenating the library ID and the DEL synthesis cycle number
        exp_cond_df['CompoundIndex'] = exp_cond_df['lib_id'].astype(str).str.cat([exp_cond_df['cycle1'].astype(str), \
                                                                                  exp_cond_df['cycle2'].astype(str), \
                                                                                  exp_cond_df['cycle3'].astype(str)], sep='.')
        # Concatenate the corresponding molecule with screening results
        exp_cond_df['experimental_cond'] = exp_cond
        exp_cond_df_full = pd.concat([lib_dfs, exp_cond_df], axis=1)
        exp_cond_df_full.reset_index(drop=True, inplace=True)
        exp_cond_df_full['hit_counts_0'] = sample_df_0['value']
        exp_cond_df_full['hit_counts_1'] = sample_df_1['value']
        tmp.append(exp_cond_df_full)
        exp_cond_df_full.to_csv(os.path.join(output_path, f'{exp_cond}.csv'), index=False)

    # Concate the results for all experiment conditions
    all_exp_cond_df = pd.concat(tmp)
    all_exp_cond_df['CompoundIndex'] = all_exp_cond_df['cycle1'].apply(pad).str.cat([all_exp_cond_df['cycle2'].apply(pad), \
                                                                                    all_exp_cond_df['cycle3'].apply(pad)], sep='.')
    all_exp_cond_df.reset_index(inplace=True)
    all_exp_cond_df.drop(columns=['index'], inplace=True)
    return all_exp_cond_df

def merge_library_for_all_experiment(experimental_condition, lib_dfs, all_exp_cond_df):
    for exp_cond in tqdm(experimental_condition):
        exp_cond_df = all_exp_cond_df[all_exp_cond_df['experimental_cond'] == exp_cond].copy().reset_index()
        exp_cond_df = exp_cond_df.rename(columns={
            'er': f'{exp_cond}_er',
            'er_lb': f'{exp_cond}_er_lb',
            'er_ub': f'{exp_cond}_er_ub',
            'hit_counts_0': f'{exp_cond}_hit_counts_0',
            'hit_counts_1': f'{exp_cond}_hit_counts_1'
        })

        lib_dfs = pd.concat([lib_dfs, exp_cond_df[[f'{exp_cond}_er', f'{exp_cond}_er_ub', f'{exp_cond}_er_lb', f'{exp_cond}_hit_counts_0', 
                                            f'{exp_cond}_hit_counts_1']]], axis=1)
    return lib_dfs

if __name__ == "__main__":
    # Load config file
    with open("config.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    output_path_preprocessed = os.path.join(config['output_path'], "preprocessed")
    os.makedirs(output_path_preprocessed, exist_ok=True)

    # Step 1: Load library
    print("Loading library...")
    lib_dfs = load_library(config)
    lib_dfs.to_csv(os.path.join(output_path_preprocessed, "lib_dfs.csv"), index=False)

    # Step2: Load sample
    print("Loading sample...")
    all_sample_df = load_sample(config)
    all_sample_df.to_csv(os.path.join(output_path_preprocessed, "all_sample.csv"), index=False)

    # Step 3: Merge sample with screeing results from different experiment conditions and save the merged results individually
    print("Merging sample and experiment...")
    all_exp_cond_df = merge_sample_and_experiment(config, all_sample_df, output_path_preprocessed)
    all_exp_cond_df.to_csv(os.path.join(output_path_preprocessed, "all_exp_cond.csv"), index=False)
    smile2compoundIndex = all_exp_cond_df[['SMILES', 'CompoundIndex']].set_index('SMILES')['CompoundIndex'].to_dict()

    # Step 4: Merge the library and all readouts of all experiment conditions (except blank)
    # Potential enhancement in the future: Step 3 and Step 4 could be merged together
    print("Merging library and experiment...")
    lib_dfs_w_exp_cond = lib_dfs.reset_index().copy()
    lib_dfs_w_exp_cond['library.name'] = lib_dfs_w_exp_cond['sublibrary']
    lib_dfs_w_exp_cond = merge_library_for_all_experiment(config['experimental_condition'], lib_dfs_w_exp_cond, all_exp_cond_df)
    
    # Step 5: Merge the blank results
    print("Merging blank...")
    blank_sample_df = all_sample_df[all_sample_df['experimental_cond'] == 'blank']
    blank_sample_0_df = blank_sample_df[blank_sample_df['run'] == 0].reset_index(drop=True).copy()
    blank_sample_0_df['blank_hit_counts_0'] = blank_sample_0_df['value']
    blank_sample_1_df = blank_sample_df[blank_sample_df['run'] == 1].reset_index(drop=True).copy()
    blank_sample_1_df['blank_hit_counts_1'] = blank_sample_1_df['value']
    lib_dfs_w_exp_cond = pd.concat([lib_dfs_w_exp_cond, blank_sample_0_df[['blank_hit_counts_0']], 
                                    blank_sample_1_df['blank_hit_counts_1']], axis=1)
    
    # Clean and save the results
    print("Saving final results...")
    lib_dfs_w_exp_cond['CompoundIndex'] = lib_dfs_w_exp_cond['SMILES'].map(smile2compoundIndex)
    lib_dfs_w_exp_cond.drop(columns=['index'], inplace=True)
    lib_dfs_w_exp_cond.to_csv(os.path.join(output_path_preprocessed, "preprocessed.csv"), index=False)
    
    




    
    