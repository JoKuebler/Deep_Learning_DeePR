from src.A_file_reader import FileReader
from src.B_encoding import Encoder
from src.C_conv_net import ConvolutionalNetwork
from src.D_refine_net import RefinementNetwork
import argparse
import numpy as np
import os

if __name__ == '__main__':

    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--retrain", help="Retrain and store model")
    parser.add_argument("-l", "--load", help="Load model")
    parser.add_argument("-i", "--input", help="File to predict")

    args = parser.parse_args()

    # Init FileReader object for Training Data
    file_read = FileReader()
    # Init convolutional network
    conv_net = ConvolutionalNetwork()
    # Init Refinement network
    ref_net = RefinementNetwork()
    # Init Encoder object
    encoder = Encoder()

    # When retrain parameter is set to True
    if args.retrain:

        # Get Training Data as list
        pos_data, neg_data = file_read.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/positive_data.txt',
                                                          '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/negative_data.txt')
        pos_test, neg_test = file_read.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/test_set_pos.txt',
                                                          '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/test_set_neg.txt')
        # Encode Training Data
        enc_data, target = encoder.encode(pos_data, neg_data)
        enc_test, target_test = encoder.encode(pos_test, neg_test)

        # Train network
        conv_net.train_network(enc_data, target, enc_test, target_test)

        # Store model in given directory
        conv_net.save_model(args.retrain)

    else:

        # Load model from given directory
        conv_net.load_model(args.load)

        refine_training = []
        refine_training_labels = []
        seq_id_to_map = []

        # Predict sequences to get prob profiles for refine network training
        for file in os.listdir(args.input):
            # Read in protein and cut into windows
            pred_data, seq_id, chain_id = file_read.read_pred_data(args.input + file, 34, 1)

            # Encode input
            enc_pred, target = encoder.encode(pred_data)

            predictions, refine_data, refine_target = conv_net.predict(pred_data, enc_pred, seq_id, chain_id)

            # Store prob profile for next network to train
            refine_training.append(refine_data)
            refine_training_labels.append(refine_target)
            seq_id_to_map.append(seq_id + ' ' + chain_id)

        # This takes the data of the CNN predictions and trains next network
        ref_net.train_predict(np.asarray(refine_training), np.asarray(refine_training_labels), seq_id_to_map)

