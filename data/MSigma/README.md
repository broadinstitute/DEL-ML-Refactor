# Millipore sigma (Msigma)
This folder contains all data and results related to Millipore Sigma library from raw preprocessing, analysing and machine learning results.

## Dataset overview
The data in raw folder contains csv file (Broad_10M_Screen_Results_2022-07-12) of provided by Millipore Sigma. Following attributes are used in our analysis:
   - `SMILES`: The smiles string of the final DEL compound
   - `Zscore`: The binding affinity of the compound to under different experiment condition
   - `Sample`: Experiment condition (Sigma-A: CK1a, Sigma-A-inh: CK1a-inh, Sigma-D: CK1d, Sigma-D-inh: CK1d-inh)

## Run analyses
### Step 0: Download the data
You can download the full version of `Broad_10M_Screen_Results_2022-07-12.csv` from [here] (link TBD)

### Step 1: Stratifying
Simply run:

```
python ./stratify.py
```
The script stratifies the compound by their effective size into three categories, which are allosteric, orthosteric and cryptic binders. It produces `{CK1a,CK1d}_orthosteric.csv`, `{CK1a,CK1d}_all_labels.csv` and `{CK1a, CK1a-inh, CK1d, CK1d-inh}_filtered.csv`. In this project, we are only interested in the orthosteric binders for each protein (i.e., `{CK1a,CK1d}_orthosteric.csv`). There is no filtering and sublibrary here but we still follow the same format in other two DEL libraries. We use `{CK1a,CK1d}_all_labels.csv` directly to generate Zscore plot.

### Step 2: Generating statistics for the libraries
Simply run:
```
python ./statistics.py
```
Since there is no filtering and sublibrary in MSigma, it only generates the quantile tables of Zscore for the whole library and Zscore plot for CK1a and CK1d.

All the output files are stored under `output/lib_stat`.
