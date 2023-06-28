import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import os

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

def get_enrichment_ratio_plot(df, experiment_condition, output_path):
    """Get the enrichment ratio plot
    Args:
        df (pd.DataFrame): Dataframe of preprocessed data for the experiment condition
        experiment_condition (str): Experimental condition
        output_path (str): Output path
    """
    palette = {"Allosteric":"tab:blue",
            "Orthosteric":"tab:green",
            "Cryptic":"tab:purple"}

    # Get the columns of interest
    df_ortho_er = df[df['customlabel'] == 'Orthosteric'].loc[:, ['library.name', 'customlabel',f'{experiment_condition}_er', f'{experiment_condition}_er_lb']]
    df_allo_er = df[df['customlabel'] == 'Allosteric'].loc[:, ['library.name', 'customlabel',f'{experiment_condition}_er_merge', f'{experiment_condition}_er_lb_merge']]
    df_cryptic_er = df[df['customlabel'] == 'Cryptic'].loc[:, ['library.name', 'customlabel',f'{experiment_condition}_inh_er', f'{experiment_condition}_inh_er_lb']]

    # Rename the columns to facilitate seaborn plotting
    df_ortho_er.rename(columns={f'{experiment_condition}_er':'Enrichment Ratio', f'{experiment_condition}_er_lb':'Enrichment Ratio Lower Bound'}, inplace=True)
    df_allo_er.rename(columns={f'{experiment_condition}_er_merge':'Enrichment Ratio', f'{experiment_condition}_er_lb_merge':'Enrichment Ratio Lower Bound'}, inplace=True)
    df_cryptic_er.rename(columns={f'{experiment_condition}_inh_er':'Enrichment Ratio', f'{experiment_condition}_inh_er_lb':'Enrichment Ratio Lower Bound'}, inplace=True)

    # Concatenate the dataframes for seaborn plotting and take log10 of the enrichment ratio
    df_er_boxenplot = pd.concat([df_ortho_er, df_allo_er, df_cryptic_er])
    df_er_boxenplot['log 10 Enrichment Ratio'] = df_er_boxenplot['Enrichment Ratio'].apply(lambda x: np.log10(x))
    df_er_boxenplot['log 10 Enrichment Ratio Lower Bound'] = df_er_boxenplot['Enrichment Ratio Lower Bound'].apply(lambda x: np.log10(x))

    # Plot the enrichment ratio for each label (Allosteric, Orthosteric, Cryptic)
    fig_er, ax_er = plt.subplots()
    sns.boxenplot(data=df_er_boxenplot.sort_values('customlabel'), x='customlabel', y='log 10 Enrichment Ratio', palette=palette, ax=ax_er)
    fig_er.savefig(os.path.join(output_path, experiment_condition, 'enrichment_ratio.png'))
    
    # Plot the enrichment ratio lower bound for each label (Allosteric, Orthosteric, Cryptic)
    fig_er_lb, ax_er_lb = plt.subplots()
    sns.boxenplot(data=df_er_boxenplot.sort_values('customlabel'), x='customlabel', y='log 10 Enrichment Ratio Lower Bound', palette=palette, ax=ax_er_lb)
    fig_er_lb.savefig(os.path.join(output_path, experiment_condition, 'enrichment_ratio_lower_bound.png'))

    # Plot the enrichment ratio of orthosteric for each library
    fig_ortho_er_per_lib, ax_ortho_er_per_lib = plt.subplots()
    sns.boxenplot(data=df_er_boxenplot[df_er_boxenplot['customlabel'] == 'Orthosteric'].sort_values('library.name'), x='library.name', y='log 10 Enrichment Ratio', ax=ax_ortho_er_per_lib)
    ax_ortho_er_per_lib.set_title(f'{experiment_condition} Orthosteric Enrichment Ratio')
    fig_ortho_er_per_lib.savefig(os.path.join(output_path, experiment_condition, 'orthosteric_enrichment_ratio_perlib.png'))

    # Plot the enrichment ratio lower bound of orthosteric for each library
    fig_ortho_erlb_per_lib, ax_ortho_erlb_per_lib = plt.subplots()
    sns.boxenplot(data=df_er_boxenplot[df_er_boxenplot['customlabel'] == 'Orthosteric'].sort_values('library.name'), x='library.name', y='log 10 Enrichment Ratio Lower Bound', ax=ax_ortho_erlb_per_lib)
    ax_ortho_erlb_per_lib.set_title(f'{experiment_condition} Orthosteric Enrichment Ratio Lower Bound')
    fig_ortho_erlb_per_lib.savefig(os.path.join(output_path, experiment_condition, 'orthosteric_enrichment_ratio_lower_bound_perlib.png'))
    plt.close()

if __name__ == "__main__":
    print("Loading data...")
    print("")
    df = pd.read_csv("./output/preprocessed/preprocessed.csv")
    df_negative = pd.read_csv("./output/stratified/negative_99k.csv")
    df_CK1a = pd.read_csv("./output/stratified/CK1a_filtered.csv")
    df_CK1a_inh = pd.read_csv("./output/stratified/CK1a_inh_filtered.csv")
    df_CK1a_all_label = pd.read_csv("./output/stratified/CK1a_all_labels.csv")
    df_CK1d = pd.read_csv("./output/stratified/CK1d_filtered.csv")
    df_CK1d_inh = pd.read_csv("./output/stratified/CK1d_inh_filtered.csv")
    df_CK1d_all_label = pd.read_csv("./output/stratified/CK1d_all_labels.csv")
    df_exclusive_CK1a = pd.read_csv("./output/stratified/CK1a_orthosteric_153k.csv")
    df_exclusive_CK1d = pd.read_csv("./output/stratified/CK1d_orthosteric_58k.csv")

    output_path = "./output/lib_stat/"
    os.makedirs(output_path)
    
    experiments = [df_exclusive_CK1a, df_exclusive_CK1d, df_CK1a, df_CK1d, df_CK1a_inh, df_CK1d_inh, df_negative, df]
    experiment_conditions = ['CK1a_orthosteric','CK1d_orthosteric', 'CK1a', 'CK1d', 'CK1a_inh','CK1d_inh', 'negative', 'all']
    quantiles = [0.0 ,0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    sublib_mol_count = df.groupby('library.name').count()['SMILES'].to_dict()

    for df_exp, exp_cond in zip(experiments, experiment_conditions):
        print(f"Calculating sublibrary statistics for {exp_cond}")

        if 'orthosteric' in exp_cond:
            os.makedirs(os.path.join(output_path, exp_cond.split("_")[0]), exist_ok=True)
        else:
            os.makedirs(os.path.join(output_path, exp_cond), exist_ok=True)
        sublibrary_stat = get_sublibary_stat(df_exp, sublib_mol_count, exp_cond)
        sublibrary_frac_bar = sns.catplot(data=sublibrary_stat,
                                        x='library.name',
                                        y='fraction',
                                        kind='bar',
                                        col='experiment_condition')
        sublibrary_pass_bar = sns.catplot(data=sublibrary_stat,
                                        x='library.name',
                                        y='pass_filter_rate',
                                        kind='bar',
                                        col='experiment_condition')
        if 'orthosteric' in exp_cond:
            sublibrary_frac_bar.savefig(os.path.join(output_path, exp_cond.split("_")[0], 'orthosteric_sublibrary_fraction.png'))
            sublibrary_pass_bar.savefig(os.path.join(output_path, exp_cond.split("_")[0], 'orthosteric_sublibrary_pass_filter_rate.png'))
        else:
            sublibrary_frac_bar.savefig(os.path.join(output_path, exp_cond, 'sublibrary_fraction.png'))
            sublibrary_pass_bar.savefig(os.path.join(output_path, exp_cond, 'sublibrary_pass_filter_rate.png'))
        plt.close()
    print("")
    print("Calculating enrichment ratio quantile statistics...")
    enrichment_ratio = []
    enrichment_ratio_lb = []
    for pfx, df in [('CK1a', df_CK1a), ('CK1a_inh', df_CK1a_inh), ('CK1d', df_CK1d), ('CK1d_inh', df_CK1d_inh)]:
        enrichment_ratio.append(df[f'{pfx}_er'].quantile(quantiles))
        enrichment_ratio_lb.append(df[f'{pfx}_er_lb'].quantile(quantiles))
    enrichment_ratio_quantile = pd.concat(enrichment_ratio, axis=1)
    enrichment_ratio_lb_quantile = pd.concat(enrichment_ratio_lb, axis=1)
    enrichment_ratio_quantile.to_csv(os.path.join(output_path, 'enrichment_ratio_quantile.csv'))
    enrichment_ratio_lb_quantile.to_csv(os.path.join(output_path, 'enrichment_ratio_lb_quantile.csv'))

    print("Plotting enrichment ratio")
    get_enrichment_ratio_plot(df_CK1a_all_label, 'CK1a', output_path)
    get_enrichment_ratio_plot(df_CK1d_all_label, 'CK1d', output_path)



