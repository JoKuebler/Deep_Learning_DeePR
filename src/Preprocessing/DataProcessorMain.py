# Created by Jonas Kuebler on 05/10/2018
# MAIN Class which is called when running Data preprocessing config
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from src.Preprocessing.FastaParser import FastaParser
from src.Tools.TPRPred import TPRPredRunner


# Take homologue fasta file
def get_pos_training_data():

    # Scope70 fasta file
    input_scope = '/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/scope70.fas'

    # Fasta File to be input for TPRpred
    input_file = 'FASTA TO PREDICT'

    # Create Fasta Parser with input File, output directory and file and optional identifier to filter sequences out
    fasta_parser = FastaParser(input_file,
                               'OUTPUT_DIRECTORY', 'positive_data_sequences.fa', 'a.118.8')

    # Filter sequences out by identifier
    filtered_input = fasta_parser.filter_by_identifier()

    # Directorry to safe TPR predictions
    output_directory = '/ebio/abt1_share/update_tprpred/data/tprpred_predictions/'
    output_file = 'predictions.txt'

    # initialize TPRPredRunner object with input directory and output directory
    tprpred_runner = TPRPredRunner(filtered_input, output_directory, output_file)

    # runs TPR pred for each file in input directory and stores predictions in output directory
    tprpred_runner.predict_tpr()


def get_neg_training_data():

    # Get negative samples for training set
    input_fasta = '/ebio/abt1_share/update_tprpred/data/Training_Data/SCOPE/scope70.fas'
    output_directory = '/ebio/abt1_share/update_tprpred/Data/Training_Data/'
    output_file = 'new_negative_set'
    identifier = ['a.7.14', 'a.7.16', 'a.23.4', 'a.24.24', 'a.118.5', 'a.118.7', 'a.118.8',
                  'a.118.20', 'a.246.2', 'd.157.1', 'e.61.1']

    # create Biopython Fasta Parser
    fasta_parser = FastaParser(input_fasta, output_directory, output_file, identifier)

    # write fragments of size 34 to file
    fasta_parser.filter_backwards_cut(34, output_file)


if __name__ == '__main__':

    # print('Nothing activated')
    # get_pos_training_data()

    get_neg_training_data()

