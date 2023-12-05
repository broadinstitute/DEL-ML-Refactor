Paper Title
====
Abstract:
"Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."

This repository contains pretrained models and scripts used for prediction mentioned in the paper (link)

## Pre-requisites:
- Linux (Tested on Ubuntu 22.04)
- NVIDIA GPU (Tested on NVIDIA RTX A6000 with cuda version 12.1)
- Python (3.10)
- Tensorflow (2.14)
- chemprop (1.6.1)
- RDkit (2023.9.2)
- Pytorch (2.1.1)
  
Please refer to [installation](docs/INSTALLATION.md) guide for how to set up the working environment



## Step 0: Data preparation
Prepare your data in the format like `example/compound.csv`. In summary, you can combine any metada of compounds but there must be a column named **SMILES**. We will use `compound.csv` as example to demonstrate the usage of other scripts


## Step 1: Feature extraction
In our paper, we use Morgan Fingerprint from RDkit with `nbits=2048, radius=2, useChirality=True`. Simply run
```
python feature_extractor.py --input_file ./example/compound.csv --save_path ./example/ --experiment compound_feature
```
Feel free to customized the `--save_path` and `--experiment` flag to suit your need. For simplicity, we save the extracted feature to the same folder and call it compound feature

## Step 2: Binder/Non-binder Prediction
We released the best two type of models (MLP and GNN) in each DEL librabry. Simply run
### Multi-layer perceptron (MLP)
```
python prediction.py --input_file ./example/compound_feature.h5 --save_path ./example/ --experiment compound_pred_mlp --checkpoint ./data/HitGen/models/CK1a/MLP.keras
```
The above command uses the MLP models pretrained on HitGen CK1a molecules to predict how likely the molecules in the `example/compound.csv` are binders.

We do not use gpu by default as we observe it does not provide a clear speed up. We believe the reason is the overhead of moving data from CPU to GPU dominates the speed up of very small model (our case). If you want to run on GPU, add `--use_gpu` flag:
```
python prediction.py --input_file ./example/compound_feature.h5 --save_path ./example/ --experiment compound_pred_mlp --checkpoint ./data/HitGen/models/CK1a/MLP.keras --use_gpu
```


### Graph neural network (GNN)
We use the graph neural network (GNN) implemented in chemprop. Following the [chemprop instruction](https://github.com/chemprop/chemprop#predicting),  run
```
chemprop_predict --smiles_columns SMILES --test_path ./example/compound.csv --checkpoint_path data/HitGen/models/CK1a/chemprop.pt --preds_path ./example/compound_prediction_chemprop.csv
```
The above command uses the GNN models pretrained on HitGen CK1a molecules to predict how likely the molecules in the `example/compound.csv` are binders.

## t-SNE visualization
After **Step 1**, you can run the following script to visualize the high dimension data in 2d space:
```
python tsne.py --input_file ./example/compound.h5 --save_path ./example/ --experiment compound --perplexity 2
```
By default, we set `n_jobs=-1` (i.e., Using all CPUs in the computer). Note `perplexity=2` in the above script is just for this example dataset. The default value used in the script is `perplexity=30`. You may need to tune this parameter to fit your dataset by running:
```
python tsne.py --input_file ./example/compound.h5 --save_path ./example/ --experiment compound --perplexity YOUR_VALUE
```

Regarding the best practice to use t-SNE and more dicussions about the method, we recommend users to read this [blog post](https://distill.pub/2016/misread-tsne/) and this [video](https://www.youtube.com/watch?v=CsUqmug7ZMc)

## Reference
If you find our work useful in your research or if you use parts of this code please consider citing our paper:

Author list, Paper Title, Journal, Year. Paper link










