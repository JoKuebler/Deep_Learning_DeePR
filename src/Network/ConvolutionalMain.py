from src.Network.ConvolutionalNetwork import ConvolutionalNetwork
from src.Preprocessing.PreprocessorConv import PreprocessorConv
from src.Helpers.HHR_Parser import HhrParser
import numpy as np


class Convolutional:

    def __init__(self):

        # Match search result .match file
        self.match_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/results/query1qqe228.match'
        self.second_match_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/results/query1fch366.match'
        self.total_matches = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/total_matches.match'
        self.rmsd_treshold = 2
        self.padded_length = 750

        # Predictions
        self.test_predict = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/1kt1.fasta'

    def init_preprocessor(self):

        # Preprocessing Protein to predict
        preprocessor_object = PreprocessorConv(self.padded_length)

        return preprocessor_object

    # This takes query MASTER hits and filters out duplicates as well as identical chains
    # More filtering according to length and amino acids
    # Writes proteins into single chains for hhpred
    def get_training_sequences(self, preprocessor_object):

        # Filter out duplicates in match file
        matches_dict = preprocessor_object.filter_duplicates(self.total_matches, self.rmsd_treshold)

        print('Unique Sequences with matches: ' + str(len(matches_dict)))
        #
        # Download each Fasta from PDB ID in the match file
        # preprocessor_object.download_fasta(matches_dict, '/ebio/abt1_share/update_tprpred/data/PDB_Approach/'
        #                                                  'All_at_once/new_download_fasta/')

        # Filter out hits with identical chains
        chain_filtered = preprocessor_object.filter_chains('/ebio/abt1_share/update_tprpred/data/'
                                                           'PDB_Approach/All_at_once/new_download_fasta/', matches_dict)

        print('Unique Sequences with matches after filtering identical chains: ' + str(len(matches_dict)))

        # Filter out sequences which are too long (returned in BioPython format)
        length_filtered = preprocessor_object.length_filter(chain_filtered, self.padded_length, matches_dict)

        print('Unique Sequences with matches after filtering too long chains: ' + str(len(matches_dict)))

        # Filter out unknown amino acids
        aa_filtered = preprocessor_object.aa_filter(length_filtered, matches_dict)

        print('Unique Sequences with matches after filtering unknown amino acids: ' + str(len(matches_dict)))

        overlap_filtered_dict = preprocessor_object.filter_overlaps(matches_dict)

        print(len(overlap_filtered_dict))

        # Write all files to single chains
        preprocessor_object.single_chains_fasta(aa_filtered, '/ebio/abt1_share/update_tprpred'
                                                             '/data/PDB_Approach/All_at_once/single_chains_all/')
        return overlap_filtered_dict

    # Takes matches dict and looks at hhr files produced by hhpred
    @staticmethod
    def eval_hhpred_results(matches_dict):

        # Filter out sequences which were not a hit in a TPR related scope class when run with hhpred
        # Filter out sequences where hit is not in template range of hhpred hit
        hhr_parser = HhrParser('/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/all_hhr/')
        # returns directory where remaining sequences are stored
        hhr_parser.filter_files(matches_dict)

    # Once final sequences are validated with HHpred they get encoded here
    def encode_data(self, preprocessor_object, final_sequences_dir, ):

        # this returns sequences as list of raw strings and as list of biopython records
        sequences = preprocessor_object.read_final_sequences(final_sequences_dir)

        # One hot encode each sequence and create numpy array use string list
        encoded_sequences = preprocessor_object.one_hot_encode(sequences[0], self.padded_length)

        hhr_parser = HhrParser()
        matches_dict = hhr_parser.read_matches_json('/ebio/abt1_share/update_tprpred/data/'
                                                    'PDB_Approach/All_at_once/tpr_info.json')

        # Create target vectors (labels) use records list
        target_vectors = preprocessor_object.create_target_vector(matches_dict, sequences[1])

        encoded_target_vector = np.asarray(target_vectors)

        print(encoded_sequences.shape)
        print(encoded_target_vector.shape)

        return [encoded_sequences, encoded_target_vector]

    @staticmethod
    def init_network():

        # Train network with Training and Test Set
        convolutional_network = ConvolutionalNetwork()

        return convolutional_network

    @staticmethod
    def train_network(network_object, training_data):

        network_object.compile_network()

        network_object.train_network(training_data[0], training_data[1])

    def predict(self, network, preprocessor):

        records = preprocessor.parse_prediction_input(self.test_predict)

        encoded = preprocessor.one_hot_encode(records, self.padded_length)

        results = network.predict(encoded)

        new_file = open('/ebio/abt1_share/update_tprpred/data/PDB_Approach/convnet_predictions' + '.txt', 'w')

        for i in range(0, len(results[0])):
            new_file.write(str(results[0][i]) + '\n')

        new_file.close()

        return results