import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import tensorflow as tf
import h5py
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="feature file")
    parser.add_argument("--save_path", type=str, required=True, help="path to store results")
    parser.add_argument("--checkpoint", type=str, required=True, help="model chekcpoint")
    parser.add_argument("--experiment", type=str, required=True, help="experiment name")
    parser.add_argument("--use_gpu", action='store_true', default=False, help="use gpu")
    args = parser.parse_args()
    if not args.use_gpu:
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
        print("Not using GPU\n")

    print("Loading data...")
    with h5py.File(args.input_file, 'r') as hf:
        smiles = hf["SMILES"][:]
        feature = hf["feature"][:]
    print("Predicting...")
    model = tf.keras.models.load_model(args.checkpoint)
    prediction = model.predict(feature)

    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
    with open(os.path.join(args.save_path, f"{args.experiment}.csv"), 'w') as fw:
        fw.write("smiles,prediction\n")
        for i in range(len(smiles)):
            fw.write(f"{smiles[i]},{prediction[i][0]}\n")



