# Installaiton guide

## With CUDA GPU 
Please make sure your NVIDIA GPU and cuda toolkit are properly intalled and configured. Should you have any question for setting and configuring GPU, please reach out to your IT for supports, or see CPU-only section

First, create a conda environment with RDkit pre-installed
```
conda create -c conda-forge -n del-ml-gpu rdkit
conda activate del-ml-gpu
```
Then simply using the `environment_gpu.yml` in to set up the environment
```
conda env create --name del-ml-gpu -f ./environment_gpu.yml
```

## CPU-only



