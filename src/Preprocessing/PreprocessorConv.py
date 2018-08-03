# Created by Jonas Kuebler on 07/15/2018
import subprocess
import os
import numpy as np
from Bio import SeqIO


class PreprocessorConv:

    def __init__(self, padded_length):
        self.amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.padded_length = padded_length

    @staticmethod
    # creates a dictionary which stores all matches with sequence, position, chain and RMSD
    def filter_duplicates(match_file, rmsd_treshold):

        print('Match file used: ' + str(match_file))

        # files with matches
        open_file = open(match_file, 'r')
        
        # list to hold unique matches
        unique_split = []

        unique_dict = {}
        i = 0

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
                chain_id = 'A'

            # RMSD in file
            rmsd = float(split[0])

            # replace unnecessary chars
            for ch in ['[', ']', '(', ')']:
                split[2] = split[2].replace(ch, '')

            tpr_start = split[2].split(',')[0]
            tpr_end = split[2].split(',')[1]

            match_entry = {"rmsd": rmsd, "chain": chain_id, "tpr_start": tpr_start, "tpr_end": tpr_end}
            match_entry_array = [match_entry]

            # only keep the ones lower as the treshold
            if rmsd < rmsd_treshold:

                # if filename with identical match positions already exists then filter out
                if (pdb_id, tpr_start) not in unique_split:

                    unique_split.append((pdb_id, tpr_start))

                    # dictionary fill
                    if pdb_id in unique_dict:
                        unique_dict[pdb_id].append(match_entry)
                    else:
                        unique_dict[pdb_id] = match_entry_array
        print("Number of different PDB entries containing matches: " + str(len(unique_dict)))

        return unique_dict

    # Takes the matches dictionary and filters chains out if they are identical
    # Sequences provides via download directory
    @staticmethod
    def filter_chains(download_dir, matches_dict):

        unique_chains_dict = {}
        final_records = []

        # For all fasta files in directory
        for filename in os.listdir(download_dir):

            pdb_id = filename[0:4].upper()
            records = list(SeqIO.parse(download_dir + filename, 'fasta'))
            sequence_store = []
            unique_chains = []

            # Look at all chains
            for record in records:

                chain_id = record.id[5]
                # Look at all matches in match dict
                for entry in matches_dict[pdb_id]:

                    # Only take the chains which match
                    if chain_id == entry['chain']:
                        # and only if they are not equal in sequence
                        if record.seq not in sequence_store:
                            sequence_store.append(record.seq)
                            unique_chains.append(chain_id)
                            final_records.append(record)
                else:
                    continue

            unique_chains_dict[pdb_id] = unique_chains

            # Delete identical sequence chains from match dict
            i = 0
            while i < len(matches_dict[pdb_id]):
                if matches_dict[pdb_id][i]['chain'] not in unique_chains_dict[pdb_id]:
                    del matches_dict[pdb_id][i]
                else:
                    i += 1

            # if all matches for one pdb are removed delete the dictionary key
            if len(matches_dict[pdb_id]) == 0:
                del matches_dict[pdb_id]
                continue

        return final_records

    @staticmethod
    def length_filter(chain_filtered, padded_length, matches_dict):

        length_filtered = []

        for record in chain_filtered:

            pdb_id = record.id[0:4]
            chain = record.id[5]

            if len(record.seq) > padded_length:
                index = 0
                for entry in matches_dict[pdb_id]:
                    if entry['chain'] == chain:
                        del matches_dict[pdb_id][index]
                    index = index + 1

                    # if all matches for one pdb are removed delete the dictionary key
                    if len(matches_dict[pdb_id]) == 0:
                        del matches_dict[pdb_id]

                continue
            else:
                length_filtered.append(record)

        return length_filtered

    def aa_filter(self, length_filtered, matches_dict):

        aa_filtered = []

        for record in length_filtered:

            pdb_id = record.id[0:4]
            chain = record.id[5]

            if self.unknown_aa_filter(record.seq):
                aa_filtered.append(record)
            else:
                index = 0
                for entry in matches_dict[pdb_id]:
                    if entry['chain'] == chain:
                        del matches_dict[pdb_id][index]
                    index = index + 1

                    # if all matches for one pdb are removed delete the dictionary key
                    if len(matches_dict[pdb_id]) == 0:
                        del matches_dict[pdb_id]

                continue

        return aa_filtered

    @staticmethod
    def unknown_aa_filter(sequence):

        not_allowed = ['X', 'U']

        if not any(x in str(sequence) for x in not_allowed):
            return True
        else:
            return False

    # Download Fasta Files for PDB IDs
    @staticmethod
    def download_fasta(unique_dict, download_dir):

        print('FASTA Files to be downloaded: ' + str(len(unique_dict)))
        for key in unique_dict:

            subprocess.run(['curl', '-o', download_dir + str(key) +
                            '.fasta', 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=' +
                            str(key) + '&compressionType=uncompressed'])

    @staticmethod
    def single_chains_fasta(final_records, output_directory):

        for record in final_records:

            SeqIO.write(record, output_directory + record.name[0:6].replace(':', '_') + '.fasta', "fasta")

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

        encoded_as_array = np.asarray(encoded_sequences)

        return encoded_as_array

    # Reads in input file and returns sequences as Biopython records
    @staticmethod
    def parse_prediction_input(input_file):

        records = list(SeqIO.parse(input_file, 'fasta'))

        return records

    # create zero one encoded target vector for each sequence
    def create_target_vector(self, match_dict, final_records):

        target_vector_list = []
        # TODO Built target vector from match dict and final records

        # # For all training sequences
        # for record in final_records:
        #     print(record)
        #     # Get raw PDB ID from the record ID
        #     pdb_id = record.id[0:4]
        #     chain_id = record.id[5]
        #     print(chain_id)
        #
        #     # create zero vector with length of sequence
        #     target_vector = [[0.0, 1.0]] * len(record.seq)
        #
        #     # get all matches for particular PDB id
        #     for entry in match_dict[pdb_id]:
        #
        #         # get start and end position out of dictionary
        #         tpr_pos = entry[0].split(',')
        #         # Fix index
        #         start_pos = int(tpr_pos[0])-1
        #         end_pos = int(tpr_pos[1])-1
        #
        #         print('PDB', pdb_id)
        #         print('ENTRY', entry)
        #         print('LEN', len(target_vector))
        #         print(start_pos)
        #         print(end_pos)
        #
        #         # Set correct positions to 1
        #         for i in range(start_pos, end_pos+1):
        #             target_vector[i] = [1.0, 0.0]
        #     # Padding
        #     while len(target_vector) < self.padded_length:
        #         target_vector.append([0.0, 1.0])
        #
        #     target_vector_list.append(target_vector)

        return target_vector_list
