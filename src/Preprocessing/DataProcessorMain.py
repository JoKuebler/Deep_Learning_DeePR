# Created by Jonas Kuebler on 05/10/2018
# MAIN Class which is called when running Data preprocessing config
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from src.Preprocessing.FastaParser import FastaParser
from src.Tools.TPRPred import TPRPredRunner


def get_pos_training_data_scope():

    # Extract all TPRs from scope
    output_directory_tpr_fasta = '/ebio/abt1_share/update_tprpred/data/tpr_fasta_files_new/'

    fasta_parser = FastaParser('/ebio/abt1_share/toolkit_sync/databases/hh-suite/scope70/scope70.fas',
                               'a.118.8', output_directory_tpr_fasta)

    fasta_parser.write_to_files_by_identifier()


def get_pos_training_data_homologues():

    # Get homologue sequences by Blast Search against nr90 and then searching Hits for TPRs with TPRpred
    output_directory_tpr_fasta = '/ebio/abt1_share/update_tprpred/data/tpr_fasta_files/'

    # After manually getting homologues with PSI Blast write each sequence to one file in order to make TPR prediction
    for file in os.listdir('/ebio/abt1_share/update_tprpred/data/homologues/'):

        # create Biopython Fasta Parser
        fasta_parser = FastaParser('/ebio/abt1_share/update_tprpred/data/homologues/' + file,
                                   '', output_directory_tpr_fasta)

        # split homologue sequences to one file per sequence
        fasta_parser.write_to_files_all()

    # Take homologue sequences and predict TPRs with TPRpred
    # Directory where tpr Fastas are stored
    input_directory = '/ebio/abt1_share/update_tprpred/data/tpr_fasta_files/'
    # Directorry to safe TPR predictions
    output_directory = '/ebio/abt1_share/update_tprpred/data/tprpred_predictions/'

    # initialize TPRPredRunner object with input directory and output directory
    tprpred_runner = TPRPredRunner(input_directory, output_directory)

    # runs TPR pred for each file in input directory and stores predictions in output directory
    tprpred_runner.predict_tpr()


def get_neg_training_data():

    # Get negative samples for training set
    input_fasta = '/ebio/abt1_share/update_tprpred/data/Training_Data/SCOPE/scope70.fas'
    output_directory = '/ebio/abt1_share/update_tprpred/data/Training_Data/'
    output_file = 'new_negative_set'
    identifier = 'a.'

    # create Biopython Fasta Parser
    fasta_parser = FastaParser(input_fasta, identifier, output_directory)

    # write fragments of size 34 to file
    fasta_parser.write_to_files_fragments(34, output_file)


if __name__ == '__main__':

    # get_pos_training_data_scope()

    get_neg_training_data()

