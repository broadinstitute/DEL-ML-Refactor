# DOS-DEL

This folder contains all data and results related to DOS-DEL (Diversity Oriented Synthesis-DNA Encoded Library) library from raw preprocessing, analysing and machine learning results.

## Dataset overview
The data in `raw` folder contain raw readouts of DOS-DEL library. It contains two sub-folders (`libraries` and `sample`) and several files. In summary, there are 6 sublibraies (002,004,005,006,007,100), each one is screened under 4 experimental conditions (CK1a, CK1a_inh,CK1d and CK1d_inh) and each condition are run two times (results are under `sample`)
We use lib002 as an example to illustrate the detail of each table (All examples provided here just a sample version. To run the experiment, please refer to **Run experiment** section)

- `libraries/lib002_sample.csv`: The file contains all SMILES strings of small molecules generated from DEL in sublibaray 002 (column: structure)
- `samples/run038_samp0000{experimental condition code}_lib002_sample.csv`: These files contain the sequencing readout (columns: value) under different experimental condition. The code of all experimental conditoins are described in `config.yaml` 
- `{experimental condition code}_lib002_sample.csv`: These files mainly contain the enrichment score and corresponding lowe and upper bound.  The meaning of each column are expplained below
  - `lib_id`: sublibary id (2 in our example)
  - `cycle{x} x=1,2,3`: DOS-DEL is a three-cycle DEL pipeline. The value in the columns represents scaffold used in each synthetic cycle.
  - `ctrl_ct`:
  - `targ_ct`:
  - `er`: Enrichment score. Roughly speaking, it measures how strong the molecule can bind to the target
  - `er_lb`: The lower bound of enrichment score
  - `er_up`: The upper bound of enrichment score


## Run experiments
### Step 0: Download the data
TBD

### Step 1: Preprocessing
Simply run:
```
python ./preprocessing.py
```
It produce `preprocessed.csv` under the folder `outputs/preprocessed`. The meaning of each columns are exaplained below
- `SMILES`: The SMILES string of the small molecule
- `sublibrary/library.name`: The id of sublibrary (TODO: Possibily duplicate columns)
- `{experimental_condition_code}_er`: Enrichment score of the molecule in the given experimental condition
- `{experimental_condition_code}_er_ub: upper bound of `{experimental_condition_code}\_er`
- `{experimental_condition_code}_lb`: lower bound of `{experimental_condition_code}_er`
- `{experimental_condition_code}_hit_count_{number of run}`: The readout of each run in the given experimental condition. The number of run are either 0 or 1
- blank_hit_count_{number of run}: The readout of each run in the blank condition (i.e., no target). The number of run are either 0 or 1
