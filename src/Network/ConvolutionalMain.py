from src.Network.ConvolutionalNetwork import ConvolutionalNetwork
from src.Preprocessing.PreprocessorConv import PreprocessorConv
import numpy as np


class Convolutional:

    def __init__(self):

        # Match search result .match file
        self.match_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/Results/query.match'
        self.rmsd_treshold = 2.5
        self.padded_length = 500

        # Predictions
        self.test_predict = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/3elr.fasta'

    def init_preprocessor(self):

        # Preprocessing Protein to predict
        preprocessor_object = PreprocessorConv(self.padded_length)

        return preprocessor_object

    def init_training_data(self, preprocessor_object):

        # Filter out duplicates in match file
        matches_dict = preprocessor_object.filter_duplicates(self.match_file, self.rmsd_treshold)

        # Download each Fasta from PDB ID in the match file
        # preprocessor_object.download_fasta(matches_dict, '/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta/')

        # Filter out sequences which are too long (returned in BioPython format)
        sequence_records = preprocessor_object.length_filter('/ebio/abt1_share/update_tprpred/data/PDB_Approach/Fasta/',
                                                             matches_dict, self.padded_length)
        # One hot encode each sequence and create numpy array
        encoded_sequences = preprocessor_object.one_hot_encode(sequence_records, self.padded_length)

        # Create target vectors (labels)
        target_vectors = preprocessor_object.create_target_vector(matches_dict, sequence_records)

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