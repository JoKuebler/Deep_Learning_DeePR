import numpy as np
from Bio import SeqIO

from src.Helpers.Result_Evaluator import Result_Evaluator
from src.Network.FeedForwardNetwork import FeedForwardNetwork
from src.Preprocessing.PreprocessorFF import PreprocessorFeedForward


class FeedForward:

    def __init__(self):
        # Positive training set
        self.training_set_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/positive_data.txt'

        # Negative training set
        self.training_negative_set_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/All_at_once/negative_data.txt'

        # Positive test set
        self.test_set_file = '/ebio/abt1/jkuebler/Desktop/windows.txt'

        # Negative test set
        self.test_negative_set_file = '/ebio/abt1_share/update_tprpred/data/Training_Data/test_set_neg.txt'

        # Different Test Files
        self.test_one = '/ebio/abt1_share/update_tprpred/data/Test_Proteins/ecoli_test.fa'
        self.test_two = '/ebio/abt1_share/update_tprpred/data/Proteomes/HumanExampleSeqs.txt'
        self.test_three = '/ebio/abt1_share/update_tprpred/data/Training_Data/test_single.txt'
        self.test_four = '/ebio/abt1_share/update_tprpred/data/Test_Proteins/RPAP3.fa'

        # Scope70 Fasta
        self.scope_predictions = '/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/scope70.fas'

        # EColi proteom
        self.e_coli = '/ebio/abt1_share/update_tprpred/data/Proteomes/Escherichia_coli_K12.fas'

        self.evaluator = Result_Evaluator()

    @staticmethod
    def init_preprocessor(input):

        # Preprocessing Protein to predict
        preprocessor_object = PreprocessorFeedForward(input)

        return preprocessor_object

    @staticmethod
    def init_network():

        # Train network with Training and Test Set
        neural_network_object = FeedForwardNetwork()

        return neural_network_object

    def init_training_data(self, preprocessor_object):

        # Read in training Data
        print('Positive Training Data:', self.training_set_file)
        pos_set = preprocessor_object.read_line_files(self.training_set_file)
        print(len(pos_set))

        print('Negative Training Data:', self.training_negative_set_file)
        neg_set = preprocessor_object.read_line_files(self.training_negative_set_file)

        print('Positive Test Data:', self.test_set_file)
        pos_test_set = preprocessor_object.read_line_files(self.test_set_file)

        print('Negative Test Data:', self.test_negative_set_file)
        neg_test_set = preprocessor_object.read_line_files(self.test_negative_set_file)

        # Multidimensional Approach encodes Data as (x, 34, 20) Array
        encoded_pos_set = preprocessor_object.encoding_training_helper(pos_set)
        encoded_neg_set = preprocessor_object.encoding_training_helper(neg_set)
        encoded_pos_test_set = preprocessor_object.encoding_training_helper(pos_test_set)
        encoded_neg_test_set = preprocessor_object.encoding_training_helper(neg_test_set)

        # Add Pos and Neg Training set to one Array
        encoded_total_set = encoded_pos_set + encoded_neg_set
        encoded_test_total_set = encoded_pos_test_set + encoded_neg_test_set

        # Convert to np Array
        encoded_array = np.asarray(encoded_total_set)
        encoded_array_test = np.asarray(encoded_test_total_set)

        # Create Labels for Training
        labels = preprocessor_object.create_labels(len(encoded_pos_set), len(encoded_neg_set))

        # Create Labels for Test Set
        labels_test = preprocessor_object.create_labels(len(encoded_pos_test_set), len(encoded_neg_test_set))

        return [encoded_array, labels, encoded_array_test, labels_test]

    @staticmethod
    def train_network(network_object, training_data):

        network_object.compile_network()

        network_object.train_network(training_data[0], training_data[1], training_data[2], training_data[3])

    @staticmethod
    def load_network(network_object):

        network_object.load_model('/ebio/abt1_share/update_tprpred/code/PycharmProjects/TrainingData/src/NetworkData/model.json',
                                  '/ebio/abt1_share/update_tprpred/code/PycharmProjects/TrainingData/src/NetworkData/model.h5')

        network_object.compile_network()

    @staticmethod
    def save_network(network_object):

        network_object.save_model()

    @staticmethod
    def cross_validate(network_object, training_data):

        network_object.cross_validate(training_data[0], training_data[1])

    @staticmethod
    def predict(network_object, preprocessor_object):

        # Make prediction
        # Cut Input Protein into windows
        fragments_per_protein = preprocessor_object.cutInWindows(34)

        prediction_set = preprocessor_object.encoding_helper(fragments_per_protein)

        new_file = open('/ebio/abt1_share/update_tprpred/data/network_predictions' + '.txt', 'w')
        counter = 0
        foundTPR = False

        for protein in prediction_set:

            prediction_set_array = np.asarray(protein)
            results = network_object.predict(prediction_set_array)

            initial_findings_index = preprocessor_object.filter_print_results(results, 0.7)

            for elem in initial_findings_index:
                foundTPR = True
                new_file.write('Protein: ' + str(preprocessor_object.ids[counter]) + '\n')
                new_file.write("Probability of being a TPR starting at position " + str(elem) + ' '
                               ':' + str(preprocessor_object.probabilites[elem]) + '\n')

                new_file.write('Sequence: ' + str(fragments_per_protein[counter][elem - 1].upper()) + '\n')

            counter += 1

        if not foundTPR:
            new_file.write('No TPR repeats found')

        new_file.close()

    # predict multiple sequences in one fasta file
    def single_predict(self, network_object, preprocessor_object, threshold):

        print('File to predict', preprocessor_object.inputFile)

        # Run Time improvements
        new_file = open('/ebio/abt1_share/update_tprpred/data/network_predictions' + '.txt', 'w')
        counter = 0
        found_tpr = True
        final_props = []

        for record in SeqIO.parse(preprocessor_object.inputFile, 'fasta'):

            if len(record.seq) >= 34:

                fragments = preprocessor_object.cutInWindows_single(record, 34)

                prediction_set = preprocessor_object.encoding_helper_single(fragments)

                prediction_set_array = np.asarray(prediction_set)

                results = network_object.predict(prediction_set_array)

                initial_findings_index = preprocessor_object.filter_print_results(results, threshold)

                for elem in initial_findings_index:
                    new_file.write('Protein: ' + str(preprocessor_object.ids[counter]) + '\n')
                    new_file.write("Probability of being a TPR starting at position " + str(elem) + ' '
                                   ':' + str(preprocessor_object.probabilites[elem]) + '\n')

                    final_props.append(preprocessor_object.probabilites[elem])

                    new_file.write('Sequence: ' + str(fragments[elem - 1].upper()) + '\n')

                counter += 1

                if not found_tpr:
                    new_file.write('No TPR repeats found')

        print("Average probability: " + str(self.evaluator.average_probability(final_props)))

        new_file.close()

    def predict_training_data(self, network_object, preprocessor_object, threshold):

        new_file = open('/ebio/abt1_share/update_tprpred/data/network_predictions' + '.txt', 'w')
        final_props = []
        found_tpr = True

        with open(preprocessor_object.inputFile) as fileContent:

            for line in fileContent:

                prediction_set = preprocessor_object.encoding_helper_single(line)

                prediction_set_array = np.asarray(prediction_set)

                results = network_object.predict(prediction_set_array)

                initial_findings_index = preprocessor_object.filter_print_results(results, threshold)

                for elem in initial_findings_index:

                    new_file.write("Probability of being a TPR starting at position " + str(elem) + ' '
                                   ':' + str(preprocessor_object.probabilites[elem]) + '\n')

                    final_props.append(preprocessor_object.probabilites[elem])

                    new_file.write('Sequence: ' + str(line.upper()) + '\n')

                if not found_tpr:
                    new_file.write('No TPR repeats found')

        print("Average probability: " + str(self.evaluator.average_probability(final_props)))
        print("Maximal probability: " + str(max(final_props)))
        print("Minimal probability: " + str(min(final_props)))

        new_file.close()
