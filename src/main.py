from src.A_file_reader import FileReader
from src.B_encoding import Encoder
from src.C_conv_net import ConvolutionalNetwork
from src.D_refine_net import RefinementNetwork
from src.DataPreprocessing.A_query_aligner import Aligner
from src.DataPreprocessing.B_file_parser import MatchParser
import argparse
import numpy as np
import os


def preprocess(align_object, matcher_object):
    """
    Runs several preprocessing steps
    and makes it easy to comment out if not needed
    :param align_object:
    :return:
    """

    # align_object.pairwise_align()

    matcher_object.parse_info()

    # matcher_object.tprpred_pvalues('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/match_files/third_set/match_dict_beauty.json')

    # matcher_object.tprpred_plot('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/tprpred_res.fa')


def network_training(reader_object, encoder_object, conv_object, ref_object):
    """
    Trains the network and evaluates parameters
    :param conv_object:
    :return:
    """

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

        # conv_net.load_model(args.retrain)

        # Read in protein and cut into windows
        pred_data, seq_id = file_read.read_pred_data(args.input, 34, 1)

        # Encode input
        enc_pred, target = encoder.encode(pred_data)

        conv_net.predict(pred_data, enc_pred)

    else:

        # Load model from given directory
        conv_net.load_model(args.load)

        refine_training = []
        refine_training_labels = []
        seq_id_to_map = []

        # Predict sequences to get prob profiles for refine network training
        for file in os.listdir(args.input_dir):
            # Read in protein and cut into windows
            pred_data, seq_id, chain_id = file_read.read_pred_data(args.input_dir + file, 34, 1)

            # Encode input
            enc_pred, target = encoder.encode(pred_data)

            predictions, refine_data, refine_target = conv_net.predict(pred_data, enc_pred, seq_id, chain_id)

            # Store prob profile for next network to train
            refine_training.append(refine_data)
            refine_training_labels.append(refine_target)
            seq_id_to_map.append(seq_id + ' ' + chain_id)

        print(np.asarray(refine_training).shape)
        print(np.asarray(refine_training_labels).shape)

        # This takes the data of the CNN predictions and trains next network
        ref_net.train_predict(np.asarray(refine_training), np.asarray(refine_training_labels), seq_id_to_map)


if __name__ == '__main__':

    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--retrain", help="Retrain and store model")
    parser.add_argument("-l", "--load", help="Load model")
    parser.add_argument("-i", "--input", help="File to predict")
    parser.add_argument("-d", "--input_dir", help="Input directory to predict multiple files")

    args = parser.parse_args()

    # Preprocessing
    # Path to queries
    aligner = Aligner('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/queries/third_set/')
    # Path to MASTER match files
    matcher = MatchParser('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/match_files/third_set/',
                          '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/all_matches/')

    # Network
    # Init FileReader object for Training Data
    file_read = FileReader()
    # Init Encoder object
    encoder = Encoder()
    # Init convolutional network
    conv_net = ConvolutionalNetwork()
    # Init Refinement network
    ref_net = RefinementNetwork()

    # Running
    # Preprocess Data
    preprocess(aligner, matcher)

    # Train Network
    # network_training(file_read, encoder, conv_net, ref_net)




