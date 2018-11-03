from src.A_file_reader import FileReader
from src.B_encoding import Encoder
from src.C_conv_net import ConvolutionalNetwork
from src.D_refine_net import RefinementNetwork
from src.E_SVM import SVM
from src.DataPreprocessing.A_query_aligner import Aligner
from src.DataPreprocessing.B_data_getter import DataGetter
import argparse
import os
import time
import numpy as np
from keras.utils import plot_model


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
    hhpred_result_dir = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/HHpred/results/'
    hhpred_querydb_results = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/HHpred/results_querydb/'
    psiblast_result_dir = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/PSIBlast/results/'

    # Containing 1561 Hits from 800 unique PDB structures
    match_data, pos_data = data_getter_object.read_match_json('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/match_dict_post_templ_filter.json')

    # data_getter_object.write_pos_data(pos_data)

    # align_object.pairwise_align()

    # data_getter_object.parse_info()

    # data_getter_object.download_fasta(match_data, full_length_dir)

    # data_getter_object.single_chains_fasta(match_data, full_length_dir, single_chain_dir)

    # data_getter_object.hhpred_init_filter(hhpred_result_dir, match_data)

    # confirmed, unconfirmed, final_seqs, final_entries = data_getter_object.check_range(hhpred_querydb_results, match_data)

    # Get enriched BLAST data
    # data_getter_object.get_blast_seqs(psiblast_result_dir, match_data)

    # Get negative data from scope
    data_getter_object.get_neg_data(600000)


def network_training(reader_object, encoder_object, conv_object, ref_object, svm):
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
        pos_data, neg_data = reader_object.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/training_sets/JanekData/train_pos.txt',
                                                              '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/training_sets/JanekData/train_neg.txt')
        pos_val, neg_val = reader_object.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/training_sets/JanekData/valid_pos.txt',
                                                            '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/training_sets/JanekData/valid_neg.txt')
        pos_test, neg_test = reader_object.read_training_data('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/training_sets/JanekData/test_pos.txt',
                                                              '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/training_sets/JanekData/test_neg.txt')
        # Encode Training Data
        enc_data, target = encoder_object.encode(pos_data, neg_data)
        enc_val, target_val = encoder_object.encode(pos_val, neg_val)
        enc_test, target_test = encoder_object.encode(pos_test, neg_test)

        # Shuffle training but keep labels in order
        indices = np.arange(enc_data.shape[0])
        np.random.shuffle(indices)
        shuf_enc = enc_data[indices]
        shuf_target = target[indices]

        # Train network or SVM
        conv_object.train_network(shuf_enc, shuf_target, enc_val, target_val, enc_test, target_test)
        # svm.train_svm(enc_data, target)

        # conv_object.cross_validate(enc_data, target)

        # Store model in given directory
        conv_object.save_model(args.retrain)

        # conv_object.load_model(args.retrain)

        # Read in protein and cut into windows
        pred_data, seq_ids = reader_object.read_pred_data(args.input, 34, 1)

        for idx, seq_fragments in enumerate(pred_data):

            # Encode input
            enc_pred, target = encoder_object.encode(seq_fragments)

            conv_object.predict(seq_fragments, enc_pred, seq_ids[idx], target)
            # svm.predict(seq_fragments, enc_pred, seq_ids[idx], target)

    else:

        start_time = time.time()
        network_data = [x[:-5] for x in os.listdir(args.load) if x.endswith('.json')]

        for file_name in network_data:

            if 'sgd-75' in file_name:

                hits = 0

                print(file_name)

                # Load model from given directory
                conv_object.load_model(args.load, file_name)
                print("--- %s seconds ---" % (time.time() - start_time))

                # Read in protein and cut into windows
                pred_data, seq_ids = reader_object.read_pred_data(args.input, 34, 1)
                print(" READING IN --- %s seconds ---" % (time.time() - start_time))

                # to save memory input proteins are predicted in chunks of 100
                chunks = [pred_data[x:x+100] for x in range(0, len(pred_data), 100)]

                for chunk in chunks:

                    # Encode input, first argument is pos data and will return pos label vector with same length and second
                    # argument is negative data and will do the same, is None if no negative data is given
                    enc_pred, pro_len = encoder_object.encode_predictions(chunk)
                    print("ENCODING --- %s seconds ---" % (time.time() - start_time))

                    # To predict F1 score give target as parameter
                    no_found = conv_object.predict(chunk, enc_pred, pro_len, seq_ids)
                    print("PREDICTION --- %s seconds ---" % (time.time() - start_time))
                    hits += no_found
                print('TOTAL FOUND: ', hits)


def visualize_model(model):

    plot_model(model, to_file='model.png', show_shapes=True, show_layer_names=True)


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
    # Init SVM
    svm = SVM()

    # Running
    # Preprocess Data
    # preprocess(aligner, data_getter)

    # Train Network
    network_training(file_read, encoder, conv_net, ref_net, svm)

    # Plot model
    # plot_model(conv_net.model)




