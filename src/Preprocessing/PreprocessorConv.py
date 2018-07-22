# Created by Jonas Kuebler on 07/15/2018
import subprocess
import os
from Bio import SeqIO


class PreprocessorConv:

    def __init__(self):
        self.amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    @staticmethod
    def filter_duplicates(match_file, rmsd_treshold):

        print('Match file used: ' + str(match_file))

        # files with matches
        open_file = open(match_file, 'r')
        
        # list to hold unique matches
        unique_split = []

        unique_dict = {}

        # parse each line
        for line in open_file:

            # split at white space
            split = line.split()

            # get rid of path just store filename
            split[1] = split[1].split('/')[-1]

            # PDB Id in file
            pdb_id = split[1][0:4].upper()

            # Check if match is from a chain by looking at the name (xxxx_A.pds vs xxxx.pds)
            if len(split[1]) > 8:
                chain_id = split[1][5:6]
            else:
                chain_id = ''

            # RMSD in file
            rmsd = float(split[0])

            # replace unnecessary chars
            for ch in ['[', ']', '(', ')']:
                split[2] = split[2].replace(ch, '')

            tpr_pos = split[2]

            # only keep the ones lower as the treshold
            if rmsd < rmsd_treshold:

                # if filename with identical match positions already exists then filter out
                if (pdb_id, tpr_pos) not in unique_split:

                    unique_split.append((pdb_id, tpr_pos))

                    # dictionary fill
                    if pdb_id in unique_dict:
                        unique_dict[pdb_id].append((tpr_pos, rmsd, chain_id))
                    else:
                        unique_dict[pdb_id] = [(tpr_pos, rmsd, chain_id)]
        print("Number of PDB entries containing unique matches: " + str(len(unique_dict)))
        return unique_dict

    # Download Fasta Files for PDB IDs
    @staticmethod
    def download_fasta(unique_dict, download_dir):

        print('FASTA Files to be downloaded: ' + str(len(unique_dict)))
        for key in unique_dict:

            subprocess.run(['curl', '-o', download_dir + str(key) +
                            '.fasta', 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=' +
                            str(key) + '&compressionType=uncompressed'])

    # Filters sequences according to length so the ones too long are removed
    def length_filter(self, download_dir, matches_dict, length):

        print("Desired length for Training: " + str(length))

        training_sequences = []

        for filename in os.listdir(download_dir):

            records = list(SeqIO.parse(download_dir + filename, 'fasta'))

            # just use the chain the match was in by checking the dictionary from the parser
            record = self.filter_chain_from_fasta(filename, matches_dict, records)

            if len(record.seq) < length:
                training_sequences.append(record)

        print("Number of sequences after length filtering: " + str(len(training_sequences)))
        return training_sequences

    @staticmethod
    def unknown_aa_filter(training_sequences):

        final_seqs = []

        not_allowed = ['X', 'U']

        for record in training_sequences:

            if not any(x in record.seq for x in not_allowed):
                final_seqs.append(record)

        print('Sequences disregarded due unknown amino acids: ' + str(len(training_sequences) - len(final_seqs)))
        return final_seqs

    # Checks fasta records and only uses the desired chain from the unique dictionary
    @staticmethod
    def filter_chain_from_fasta(filename, match_dict, records):

        filename_no_extension = filename[0:4].upper()
        # Triple (RMSD, POSMATCH, CHAIN) in the dictionary
        triple = match_dict[filename_no_extension]

        # Chain ID
        chain_id = triple[0][2]

        for record in records:
            # Even though I stored all chains which were matches I decided
            # here to use only one chain match of each PDB entry in the match file
            # That is why after it finds a record it gets immediately returned
            if str(filename_no_extension + ':' + chain_id) in record.name:
                return record

    # One hot encode sequences
    def one_hot_encode(self, sequences, padded_length):

        encoded_sequences = []

        for sequence in sequences:

            # Initial Zero Vector
            zero_vector = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

            # Initial maxtrix where each vector will be stored in
            sequence_matrix = []

            # Copy variable
            current_vector = zero_vector.copy()

            # For each amino acid in the sequence
            for amino_acid in sequence:

                # Set correct index to 1.0 according to amino acid
                current_vector[self.amino_acids.index(amino_acid.upper())] = 1.0

                # Append vector to matrix
                sequence_matrix.append(current_vector)

                # Reset vector after each amino acid
                current_vector = zero_vector.copy()

            # Zero pad sequences to make them all of the same length
            while len(sequence_matrix) < padded_length:
                sequence_matrix.append(zero_vector)

            # return matrix with (length, 20) dim
            encoded_sequences.append(sequence_matrix)

        return encoded_sequences

    # create zero one encoded target vector for each sequence
    def create_target_vector(self, matches_dict, sequence_records):

        # TODO get positions out of matchtes_dict and encode sequence as 1 inside a TPR and 0 as when no TPR

        return 0