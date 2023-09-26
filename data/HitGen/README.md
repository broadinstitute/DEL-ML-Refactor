# HitGen
This folder contains all data and results related to HitGen library from raw preprocessing, analysing and machine learning results.

## Dataset overview

The data in `raw` folder contains csv file (`merged_molecule_enrichments_sample.csv`) of HitGen 3-synthon DNA-ecnode library provided by HitGen. Following attributes are used in our analysis:
  - `library.name`: The sub-library where the molecules come from
  - `product_smiles`: SMILES string of the final DEL molecule (renamed as `SMILES` in our script)
  - `mol.id (depreciated)`: Molecular ID. We do not use this anymore as we found this is not unique
  - ` {A, A-inh, D, D-inh, blank} Effectt Size`: A number to tell the enrichment score of DEL molecule (i.e., how strong the molecule can bind to the specific target). This is the main attribute we used to stratify the molecules and do exploration analysis
  - ` {A, A-inh, D, D-inh, blank} counts`: Similar idea like effect size but measured in raw molecular count (i.e., the higher mean stronger binding). We only use this for filtering low affinity molecules.


## Run analyses

### Step 0: Download the data
You can download the full version of `merged_molecule_enrichments.csv` from [here] (link TBD)

### Step 1: Stratifying
Simply run:

```
python ./stratify.py
```
The script stratifies the compound by their effective size into three categories, which are allosteric, orthosteric and cryptic binders. It produces `{CK1a,CK1d}_orthosteric.csv`, `{CK1a,CK1d}_all_labels.csv`, `{CK1a, CK1a-inh, CK1d, CK1d-inh}_filtered.csv` and `{CK1a,CK1d}_exclusive_competitive.csv`. In this projects, we are only interested in the orthosteric binders for each protein (i.e., `{CK1a,CK1d}_orthosteric.csv`). `{CK1a,CK1d}_all_labels.csv` and `{CK1a, CK1a-inh, CK1d, CK1d-inh}_filtered.csv` are only used for understanding how the different binder looks like in the sub-libriaries. 

All the output files are stored under `output/stratified/`. The meaning of columns are illustrated in dataset overview 

### Step 3: Generating statistics for the libraries
Simply run:
```
python ./statistics.py
```

The script will generate the plots for each experiment condition (except the last one):
  - rawe molecular ratio of each sublibrary for each experiment conditoin 
  - molecular ratio of each sublibrary after filtering low hit count for each experiment condition
  - The effective size boxenplot of CK1a and CK1d
  - The effective size boxenplot, density plot with default threshold of the whole library
  - The quantile table of effective size in the whole library under different filtering condition (before filtering, after filtering by effective size and after filtering by counts)
  - The effective size of binder of interested (orthosteric) for CK1d and CK1a (i.e., experiment of interest)

All the output files are stored under `output/lib_stat`.
