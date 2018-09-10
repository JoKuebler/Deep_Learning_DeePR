from src.A_file_reader import FileReader
from src.B_encoding import Encoder
from src.C_conv_net import ConvolutionalNetwork
from src.D_refine_net import RefinementNetwork
from src.DataPreprocessing.A_query_aligner import Aligner
from src.DataPreprocessing.B_data_getter import DataGetter
import argparse
import numpy as np
import os


def preprocess(align_object, data_getter_object):
    """
    Runs several preprocessing steps
    and makes it easy to comment out if not needed
    :param align_object:
    :param data_getter_object
    :return:
    """

    full_length_dir = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/full_length_fasta/'
    single_chain_dir = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/single_chain_fasta/'

    # Containing 6462 Hits from 4057 unique PDB structures
    match_data, pos_data = data_getter_object.read_match_json('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/match_files/third_set/match_dict_query_beauty.json')

    data_getter_object.write_pos_data(pos_data)

    # align_object.pairwise_align()

    # data_getter_object.parse_info()

    # data_getter_object.download_fasta(match_data, full_length_dir)

    # data_getter_object.single_chains_fasta(match_data, full_length_dir, single_chain_dir)


def network_training(reader_object, encoder_object, conv_object, ref_object):
    """
    Trains the network and evaluates parameters
    :param conv_object:
    :param encoder_object
    :param reader_object
    :param ref_object
    :return:
    """

    # When retrain parameter is set to True
    if args.retrain:

        # Get Training Data as list
        pos_data, neg_data = reader_object.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/positive_data.txt',
                                                              '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/negative_data.txt')
        pos_test, neg_test = reader_object.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/test_set_pos.txt',
                                                              '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/test_set_neg.txt')
        # Encode Training Data
        enc_data, target = encoder_object.encode(pos_data, neg_data)
        enc_test, target_test = encoder_object.encode(pos_test, neg_test)

        # Train network
        conv_object.train_network(enc_data, target, enc_test, target_test)

        # Store model in given directory
        conv_object.save_model(args.retrain)

        # conv_object.load_model(args.retrain)

        # Read in protein and cut into windows
        pred_data, seq_id = reader_object.read_pred_data(args.input, 34, 1)

        # Encode input
        enc_pred, target = encoder_object.encode(pred_data)

        conv_object.predict(pred_data, enc_pred)

    else:

        # Load model from given directory
        conv_object.load_model(args.load)

        refine_training = []
        refine_training_labels = []
        seq_id_to_map = []

        # Predict sequences to get prob profiles for refine network training
        for file in os.listdir(args.input_dir):
            # Read in protein and cut into windows
            pred_data, seq_id, chain_id = reader_object.read_pred_data(args.input_dir + file, 34, 1)

            # Encode input
            enc_pred, target = encoder_object.encode(pred_data)

            predictions, refine_data, refine_target = conv_object.predict(pred_data, enc_pred, seq_id, chain_id)

            # Store prob profile for next network to train
            refine_training.append(refine_data)
            refine_training_labels.append(refine_target)
            seq_id_to_map.append(seq_id + ' ' + chain_id)

        print(np.asarray(refine_training).shape)
        print(np.asarray(refine_training_labels).shape)

        # This takes the data of the CNN predictions and trains next network
        ref_object.train_predict(np.asarray(refine_training), np.asarray(refine_training_labels), seq_id_to_map)


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
    data_getter = DataGetter('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/match_files/third_set/',
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
    preprocess(aligner, data_getter)

    # Train Network
    # network_training(file_read, encoder, conv_net, ref_net)




