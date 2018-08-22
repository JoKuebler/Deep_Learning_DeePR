from src.A_file_reader import FileReader
from src.B_encoding import Encoder
from src.C_conv_net import ConvolutionalNetwork
from src.D_refine_net import RefinementNetwork
import argparse
import numpy as np

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
    # ref_net = RefinementNetwork()
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

        # Read in protein and cut into windows
        pred_data = file_read.read_pred_data(args.input, 34, 1)

        # Encode input
        enc_pred, target = encoder.encode(pred_data)

        conv_net.predict(pred_data, enc_pred)



