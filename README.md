DEL+ML paradigm for finding actionable discovery – a cross DEL and cross ML model assessment
====
Preprint: [DEL+ML paradigm for actionable hit discovery – a cross DEL and cross ML model assessment](https://chemrxiv.org/engage/chemrxiv/article-details/66a05468c9c6a5c07aae574d).

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


### Output
Here is the output format of the prediction
```
SMILES,prediction:
CCN1CCCC1Cn1cnc2c3ccc(OC)cc3nc-2c1O,0.678
CCC1(c2ccccc2)CC(=O)C(C2CC(c3ccc(OCc4ccc(C(F)(F)F)cc4)cc3)Cc3ccccc32)=C(O)O1,0.111
CCOC(=O)C(C)(OC(C)=O)c1cc(C)c(/N=C/N(C)C)c(C)c1,0.281
Cc1cc(C)n2s/c(=N\C(=O)C(c3ccc(Cl)cc3)C(C)C)nc2n1,0.112
OC1=C(Cl)/C(=N\Cc2ccccc2)C(O)O1,0.342
```
It is a .csv file that contain the prediction score from model output of each molecule. 

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










