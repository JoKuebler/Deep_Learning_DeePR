from src.Network.ConvolutionalNetwork import ConvolutionalNetwork
from src.Preprocessing.PreprocessorConv import PreprocessorConv
from src.Helpers.HHR_Parser import HhrParser
from src.Preprocessing.PreprocessorFF import PreprocessorFeedForward
import numpy as np


class Convolutional:

    def __init__(self):

        # Match search result .match file
        self.match_file = '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/results/updated/total.match'
        self.total_matches = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/total_matches.match'
        self.rmsd_treshold = 2.5
        self.padded_length = 750

        # Predictions
        self.test_predict = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/test_predict.txt'

    def init_preprocessor(self):

        # Preprocessing Protein to predict
        preprocessor_object = PreprocessorConv(self.padded_length)

        return preprocessor_object

    # This takes query MASTER hits and filters out duplicates as well as identical chains
    # More filtering according to length and amino acids
    # Writes proteins into single chains for hhpred
    def get_training_sequences(self, preprocessor_object):

        # PASS MATCH FILE OBTAINED FROM MASTER TO BUILD DICT AND DOWNLOAD FASTAS
        # -------------------------------------------- #
        # Filter out duplicates in match file
        # matches_dict = preprocessor_object.filter_duplicates(self.match_file, self.rmsd_treshold)
        # print('Unique Sequences with matches: ' + str(len(matches_dict)))
        # Download each Fasta from PDB ID in the match file
        # preprocessor_object.download_fasta(matches_dict, '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/updated_fasta_download/')
        # ---------------------------------------------#

        # UNCOMMENT ONE OF BOTH BLOCKS SEPARATELY

        # ONCE FASTAS ARE DOWNLOADED SPLIT AND FILTER THEM TO GET FINAL FASTAS
        # -------------------------------------------- #
        # Write all files to single chains
        # preprocessor_object.single_chains_fasta('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/updated_fasta_download/',
        #                                        '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/update_single_chains/')
        # preprocessor_object.filter_chains_new('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/update_single_chains/')
        # preprocessor_object.length_filter_new('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/update_single_chains/')
        # preprocessor_object.aa_filter_new('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/update_single_chains/')
        # matches_dict_new = preprocessor_object.build_matches_dict(self.match_file, '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/update_single_chains/')
        matches_dict_new = HhrParser.read_matches_json('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/match_dict.json')
        # filter_overlap = preprocessor_object.filter_overlaps(matches_dict_new)
        # HhrParser.write_matches_json('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/', filter_overlap)
        # -------------------------------------------- #

        return matches_dict_new

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

    # cuts input sequence in 34 windows and stores positive data
    def encode_as_window(self, preprocessor_object, final_sequences_dir, matches_dict):

        # this returns sequences as list of raw strings and as list of biopython records
        sequences = preprocessor_object.read_final_sequences(final_sequences_dir)

        new_file = open('/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/positive_data' + '.txt', 'w')

        hhr_parser = HhrParser()
        matches_json = hhr_parser.read_matches_json(matches_dict)

        # sequences[1] stores biopython records
        for record in sequences[1]:

            pdb_id = str(record.id[0:4])
            chain_id = str(record.id[5])

            for entry in matches_json[pdb_id]:

                if entry['chain'] == chain_id:
                    fragment = record.seq[int(entry['tpr_start']):int(entry['tpr_end'])+1]
                    if len(fragment) < 33:
                        print(pdb_id)
                        print(entry)
                    else:
                        new_file.write(str(record.seq[int(entry['tpr_start']):int(entry['tpr_end'])+1]) + '\n')

        new_file.close()

    @staticmethod
    def init_network():

        # Train network with Training and Test Set
        convolutional_network = ConvolutionalNetwork()

        return convolutional_network

    @staticmethod
    def train_network(network_object, training_data):

        network_object.compile_network()

        network_object.train_network(training_data[0], training_data[1], training_data[2], training_data[3])

    def predict(self, network, preprocessor):

        records = preprocessor.parse_prediction_input(self.test_predict)

        encoded = preprocessor.one_hot_encode(records, self.padded_length)

        results = network.predict(encoded)

        new_file = open('/ebio/abt1_share/update_tprpred/data/PDB_Approach/convnet_predictions' + '.txt', 'w')

        for i in range(0, len(results[0])):
            new_file.write(str(results[0][i]) + '\n')

        new_file.close()

        return results

    def predict_fragments(self, network, preprocessor):

        with open(self.test_predict) as input_file:
            content = input_file.readlines()

        stripped = [x.rstrip() for x in content]

        prediction_set = preprocessor.one_hot_encode(stripped, 34)

        results = network.predict(prediction_set)

        print(results)

