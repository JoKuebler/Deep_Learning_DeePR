import time

import numpy as np

from src.Network.NeuralNetwork import NeuralNetwork
from src.Preprocessing.Preprocessor import Preprocessor

start_time = time.time()
# TRAINING
# training Set
training_set_file = '/ebio/abt1_share/update_tprpred/data/Training_Data/positive_set.txt'
training_negative_set_file = '/ebio/abt1_share/update_tprpred/Training_Data/negative_set.txt'
test_set_file = '/ebio/abt1_share/update_tprpred/data/Training_Data/test_set_pos.txt'
test_negative_set_file = '/ebio/abt1_share/update_tprpred/data/Training_Data/test_set_neg.txt'

multiple_predictions = '/ebio/abt1_share/update_tprpred/data/Test_Proteins/vikram_more.fa'

#scope_predictions = '/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/scope70.fas'

human_proteins = '/ebio/abt1_share/update_tprpred/data/Proteomes/25.H_sapiens.fasta'

multiple_proteins = '/ebio/abt1_share/update_tprpred/data/old_maybe_gold_later/testRandom.txt'

# Test Protein for Prediction
inputProtein = '/ebio/abt1_share/update_tprpred/data/Test_Proteins/1W3B.fa'

# Preprocessing Protein to predict
PreprocessorObject = Preprocessor(inputProtein)
#
# pos_set = PreprocessorObject.read_line_files(training_set_file)
# neg_set = PreprocessorObject.read_line_files(training_negative_set_file)
# pos_test_set = PreprocessorObject.read_line_files(test_set_file)
# neg_test_set = PreprocessorObject.read_line_files(test_negative_set_file)
#
# # Multidimensional Approach encodes Data as (x, 34, 20) Array
# encoded_pos_set = PreprocessorObject.encoding_training_helper(pos_set)
# encoded_neg_set = PreprocessorObject.encoding_training_helper(neg_set)
# encoded_pos_test_set = PreprocessorObject.encoding_training_helper(pos_test_set)
# encoded_neg_test_set = PreprocessorObject.encoding_training_helper(neg_test_set)
#
# # Add Pos and Neg Training set to one Array
# encoded_total_set = encoded_pos_set + encoded_neg_set
# encoded_test_total_set = encoded_pos_test_set + encoded_neg_test_set
#
# # Convert to np Array
# encoded_array = np.asarray(encoded_total_set)
# encoded_array_test = np.asarray(encoded_test_total_set)
#
# # Create Labels for Training
# labels = PreprocessorObject.create_labels(len(encoded_pos_set), len(encoded_neg_set))
#
# # Create Labels for Test Set
# labels_test = PreprocessorObject.create_labels(len(encoded_pos_test_set), len(encoded_neg_test_set))

# Train network with Training and Test Set
NeuralNetworkObject = NeuralNetwork()
#
NeuralNetworkObject.load_model('/ebio/abt1_share/update_tprpred/code/PycharmProjects/TrainingData/src/NetworkData/model.json',
                               '/ebio/abt1_share/update_tprpred/code/PycharmProjects/TrainingData/src/NetworkData/model.h5')

NeuralNetworkObject.compile_network()

# NeuralNetworkObject.train_network(encoded_array, labels, encoded_array_test, labels_test)

# NeuralNetworkObject.cross_validate(encoded_array, labels, encoded_array_test, labels_test)

# Save and load model to/from json
# NeuralNetworkObject.save_model()


# Make prediction
# Cut Input Protein into windows
fragments_per_protein = PreprocessorObject.cutInWindows(34)

print('Cutted in fragments in: ' + "%s seconds" % (str(time.time() - start_time)))
start_time = time.time()

prediction_set = PreprocessorObject.encoding_helper(fragments_per_protein)

print('Encoded in : ' + "%s seconds" % (str(time.time() - start_time)))
start_time = time.time()

new_file = open('/ebio/abt1_share/update_tprpred/data/network_predictions' + '.txt', 'w')
counter = 0
foundTPR = False

for protein in prediction_set:

    prediction_set_array = np.asarray(protein)
    results = NeuralNetworkObject.predict(prediction_set_array)

    initial_findings_index = PreprocessorObject.filter_print_results(results, 0.7)

    for elem in initial_findings_index:
        foundTPR = True
        new_file.write('Protein: ' + str(PreprocessorObject.ids[counter]) + '\n')
        new_file.write("Probability of being a TPR starting at position " + str(elem) + ' '
                       ':' + str(PreprocessorObject.probabilites[elem]) + '\n')

        new_file.write('Sequence: ' + str(fragments_per_protein[counter][elem - 1].upper()) + '\n')

    counter += 1

if not foundTPR:
    new_file.write('No TPR repeats found')

print('Results written in : ' + "%s seconds" % (str(time.time() - start_time)))

new_file.close()



# Run Time improvements
# for record in SeqIO.parse(PreprocessorObject.inputFile, 'fasta'):
#
#     fragments = PreprocessorObject.cutInWindows_single(record, 34)
#
#     prediction_set = PreprocessorObject.encoding_helper_single(fragments)
#
# print('Single cutting and encoding: ' + "%s seconds" % (str(time.time() - start_time)))
# start_time = time.time()

