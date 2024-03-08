# DOS-DEL

This folder contains all data and results related to DOS-DEL (Diversity Oriented Synthesis-DNA Encoded Library) library from raw preprocessing, analysing and machine learning results.

## Dataset overview
The data in `raw` folder contain raw readouts of DOS-DEL library. It contains two sub-folders (`libraries` and `sample`) and several files. In summary, there are 6 sublibraies (002,004,005,006,007,100), each one is screened under 4 experimental conditions (CK1a, CK1a_inh,CK1d and CK1d_inh) and each condition are run two times (results are under `sample`)
We use lib002 as an example to illustrate the detail of each table (All examples provided here are sampled version. To run the experiment, please refer to **Run experiment** section)

- `libraries/lib002_sample.csv`: The file contains all SMILES strings of small molecules generated from DEL in sublibaray 002 (column: structure)
- `samples/run038_samp0000{experimental condition code}_lib002_sample.csv`: These files contain the sequencing readout (columns: value) under different experimental condition. The code of all experimental conditoins are described in `config.yaml` 
- `{experimental condition code}_lib002_sample.csv`: These files mainly contain the enrichment score and corresponding lower and upper bound.  The meaning of each column are expplained below
  - `lib_id`: sublibary id (2 in our example)
  - `cycle{x} x=1,2,3`: DOS-DEL is a three-cycle DEL pipeline. The value in the columns represents scaffold used in each synthetic cycle.
  - `ctrl_ct`:
  - `targ_ct`:
  - `er`: Enrichment score. Roughly speaking, it measures how strong the molecule can bind to the target
  - `er_lb`: The lower bound of enrichment score
  - `er_up`: The upper bound of enrichment score


## Run experiments
### Step 0: Download the data
The data is downloaded from [here](http://chembio-dev-02:3838/del-app/). (Detail instruction: TBD)

### Step 1: Preprocessing
Simply run:
```
python ./preprocessing.py
```
It produces `preprocessed.csv` under the folder `outputs/preprocessed`. The meaning of each columns are exaplained below
- `SMILES`: The SMILES string of the small molecule
- `sublibrary/library.name`: The id of sublibrary (TODO: Possibily duplicate columns)
- `{experimental_condition}_er`: Enrichment score of the molecule in the given experimental condition
- `{experimental_condition}_er_ub`: upper bound of `{experimental_condition}\_er`
- `{experimental_condition}_lb`: lower bound of `{experimental_condition}_er`
- `{experimental_condition}_hit_count_{number of run}`: The readout of each run in the given experimental condition. The number of run are either 0 or 1
- `blank_hit_count_{number of run}`: The readout of each run in the blank condition (i.e., no target). The number of run are either 0 or 1

### Step 2: Stratifying
Simply run:
```
python ./stratify.py
```
The script stratifies the compound by their enrichment scores into three categories, which are allosteric, orthosteric and cryptic binders. It produces `CK1a_orthosteric_153k.csv`, `CK1d_orthosteric_58k.csv`, `CK1a_all_labels.csv` and `CK1d_all_labels.csv`. In this projects, we are only interested in the orthosteric binders for each protein (`*_orthosteric_*.csv`). `*_all_lables.csv` are only used for understanding how the different binder looks like in the libriary. All the output files are stored under `output/stratified/`. The meaning of columns are the same in the **preprocessing** step.

### Step 3: Generate statistics for the library
Simply run:
```
python ./statistic.py
```
The script will generate the plots for each experiment condition (except the last one):
  - rawe molecular ratio of each sublibrary for each experiment conditoin 
  - molecular ratio of each sublibrary after filtering low hit count for each experiment condition
  - The enrichment ratio and its lower bound of each experiment condition
  - The quantile table of enrichment ratio and its lower bound for each
  - The enrichment ratio and its lower bound of binder of interested (orthosteric) for CK1d and CK1a (i.e., experiment we are interested in)

All the output files are stored under `output/lib_stat`. If you just want to reproduce the figure in this step, please just run this step directly from `/home/unix/chengkua/project/DEL-ML-Refactor_dev/data/DOS-DEL` in workstation (gpbd5-975). Due to the large file size, the files used to generate figures is not uploaded here

## TODO
- [ ] As reflected in the issue log, The [merge](https://github.com/broadinstitute/DEL-ML-Refactor/blob/6d81c047161e7b2f679352e3ae29a40d7db69b6c/data/DOS-DEL/stratify.py#L29) in label generation function use relevant columns to ensure there is no duplidcate, which is pretty dirty. Refactor it by de-duplicate first then use primary key to join (i.e., SMILES)
- [ ] Consider to move some hardcoded variable in the code to `config.yaml`
- [x] Sorting out the missing code piece that produce effective size fileds in DOS-DEL (Sumaiya and Kuan)
- [x] Split exploratory analysis from DOS-DEL-Analysis.ipynb
- [x] Split stratify from DOS-DEL-Analysis.ipynb

## Refactoring issue log
- 06/21/2023:  Found the supposed primary key (CompoundIndex) is a fake primary key. This issue is highly likely to affect the subsequent analysis since the binder type are stratified by CompoundIndex. The current repo already uses SMILE as PK instead. 
