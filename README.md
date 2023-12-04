# Paper Title
Abstract:
"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

This repository contains pretrained models and scripts used for prediction mentioned in the paper (link)

## Pre-requisites:
- Linux (Tested on Ubuntu 22.04)
- NVIDIA GPU (optional, tested on NVIDIA RTX A6000)
- Python (3.10)
Please refer to installation guide for more details'

## Installation guide
- Installation guide (GPU)
- Installation guide (CPU-only)

## Step 0: Data preparation
Prepare your data in the format like `example/compound.csv`. In Summary, you can combine any metada of compounds but there must be a column named **SMILES**. We will use `compound.csv` as example to demonstrate the usage of other scripts


## Step 1: Feature extraction
In our paper, we use Morgan Fingerprint from RDkit with `nbits=2048, radius=2, useChirality=True`. Simply run
```
python feature_extractor.py --input_file ./example/compound.csv --save_path ./example/ --experiment compound
```
Feel free to customized the `--save_path` and `--experiment` flag to suit your need. For simplicity, we save the extracted feature to the same folder and call it compound

## Step 2: Binder/Non-binder Prediction
We released the best two type of models (MLP and GNN) in each DEL librabry. Simply run
### MLP
```
python prediction.py --input_file ./example/compound.h5 --save_path ./example/ --experiment compound --checkpoint ./data/HitGen/models/CK1a/MLP.keras
```

We do not use gpu by default as we observe it does not provide a clear speed up. If you want to run on GPU, add `--use_gpu` flag:
```
python prediction.py --input_file ./example/compound.h5 --save_path ./example/ --experiment compound --checkpoint ./data/HitGen/models/CK1a/MLP.keras --use_gpu
```
Note that we do not see a clear speed up of using GPU (it is somtimes slower). We believe the reason is the overhead of moving data from CPU to GPU dominates the speed up of very small model (our case).

### GNN
We use the graph neural network (GNN) implemented in chemprop. Following the [chemprop instruction](https://github.com/chemprop/chemprop#predicting),  run
```
chemprop_predict --smiles_columns SMILES --test_path ./example/compound.csv --checkpoint_path data/HitGen/models/CK1a/chemprop.pt --preds_path ./example/compound_prediction.csv
```
The above command uses the GNN models pretrained on HitGen CK1a molecules to predict how likely the molecules in the `example/compound.csv` are binders.

## t-SNE visualization
After **Step 1**, you can run the following script to visualize the high dimension data in 2d space
```
python tsne.py --input_folder ./example --experiment compound --save_path ./example
```











