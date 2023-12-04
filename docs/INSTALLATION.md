# Installaiton guide

Please make sure your NVIDIA GPU and cuda toolkit are properly intalled and configured. Should you have any question for setting and configuring GPU, please reach out to your IT for supports. We use Anaconda to manage our packages and environemnts, please follow the instructions [here](https://docs.anaconda.com/free/anaconda/install/linux/) to download and install anaconda.

Once GPU, cuda and conda are properly configured and installed, simply run the following to set up and activate the environemtn:
```
conda env create --name del-ml-gpu -f ./environment_gpu.yml
conda activate del-ml-gpu
```



