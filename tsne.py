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
    parser.add_argument("--input_folder", type=str, required=True, help="folder to store features to be reduced")
    parser.add_argument("--experiment", type=str, required=True, help="experiment name")
    args = parser.parse_args()

    data_len = []
    data = []
    print("Loading data...")
    for file in glob.glob(os.path.join(args.input_folder, "*.h5")):
        if args.experiment in file or 'broad_cp_140k' in file or 'negative' in file:
            print(file)
            file_name = os.path.basename(file).replace(".h5", "")
            lib_name = os.path.basename(os.path.dirname(file))
            key = "_".join([lib_name, file_name])
            with h5py.File(file, "r") as f:
                feature = f["feature"][:]
            data.append(feature)
            data_len.append((key, feature.shape[0]))
    data = np.concatenate(data, axis=0)
    print("Running TSNE...")
    t_start = time.time()
    data_embedded = TSNE(n_components=2, n_jobs=-1, random_state=0).fit_transform(data)
    print('Time elapsed: {} seconds'.format(time.time()-t_start))

    print("Saving data...")
    with open(f"data_embedded_{args.experiment}.pkl", "wb") as f:
        pickle.dump(data_embedded, f)
    with open(f"data_embedded_meta_{args.experiment}.pkl", "wb") as f:
        pickle.dump(data_len, f)

