from src.Network.ConvolutionalNetwork import ConvolutionalNetwork
from src.Preprocessing.PreprocessorConv import PreprocessorConv
import numpy as np


class Convolutional:

    def __init__(self):

        # Match search result .match file
        self.match_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/Results/query.match'
        self.rmsd_treshold = 2.0

    @staticmethod
    def init_preprocessor():

        # Preprocessing Protein to predict
        preprocessor_object = PreprocessorConv()

        return preprocessor_object

    def init_training_data(self, preprocessor_object):

        padded_length = 500

        # Filter out duplicates in match file
        matches_dict = preprocessor_object.filter_duplicates(self.match_file, self.rmsd_treshold)

        # Download each Fasta from PDB ID in the match file
        # preprocessor_object.download_fasta(matches_dict, '/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta/')

        # Filter out sequences which are too long (returned in BioPython format)
        sequence_records = preprocessor_object.length_filter('/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta/',
                                                             matches_dict, padded_length)

        # One hot encode each sequence
        encoded_sequences = preprocessor_object.one_hot_encode(sequence_records, padded_length)

        encoded_array = np.asarray(encoded_sequences)

        print(encoded_array.shape)

    @staticmethod
    def init_network():

        # Train network with Training and Test Set
        convolutional_network = ConvolutionalNetwork()

        return convolutional_network

    @staticmethod
    def train_network(network_object, training_data):

        network_object.compile_network()

        # network_object.train_network(training_data[0], training_data[1], training_data[2], training_data[3])