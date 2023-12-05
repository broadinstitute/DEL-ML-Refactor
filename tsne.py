from sklearn.manifold import TSNE
import numpy as np
import time
import h5py
import glob
import os
import pickle
import argparse
np.random.seed(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="folder to store features to be reduced")
    parser.add_argument("--save_path", type=str, required=True, help="folder to store features to be reduced")
    parser.add_argument("--experiment", type=str, required=True, help="experiment name")
    parser.add_argument("--perplexity", type=int, default=30, help="perplexity of TSNE")
    parser.add_argument("--n_jobs", type=int, default=-1, help="number of CPU to run TSNE")
    args = parser.parse_args()

    data = []
    print("Loading data...")
    with h5py.File(args.input_file, "r") as f:
        feature = f["feature"][:]

    print("Running TSNE...")
    t_start = time.time()
    data_embedded = TSNE(n_components=2, n_jobs=args.n_jobs, 
                         perplexity=args.perplexity, random_state=0).fit_transform(feature)
    print('Time elapsed: {} seconds'.format(time.time()-t_start))

    print("Saving data...")
    with open(os.path.join(args.save_path, f"data_embedded_{args.experiment}.pkl"), "wb") as f:
        pickle.dump(data_embedded, f)

