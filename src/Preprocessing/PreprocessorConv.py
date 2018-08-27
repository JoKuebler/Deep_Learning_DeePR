# Created by Jonas Kuebler on 07/15/2018
import subprocess
import os
import numpy as np
from Bio import SeqIO
from src.Helpers.HHR_Parser import HhrParser


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
        open_file = open(match_file, 'r', encoding='utf-8-sig')
        
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

    def build_matches_dict(self, match_file, fasta_dir):

        match_dict = {}
        # list to hold unique matches
        unique_split = []

        open_file = open(match_file, 'r', encoding='utf-8-sig')
        counter = 0

        for file in os.listdir(fasta_dir):

            filename = str(file[:-6])

            print(filename)
            counter += 1
            print(counter)

            # parse each line
            for line in open_file:

                # split at white space
                split = line.split()
                split[1] = split[1].split('/')[-1]
                # PDB Id in file
                pdb_id = split[1][0:4].upper()
                # RMSD in file
                rmsd = float(split[0])

                # Check if match is from a chain by looking at the name (xxxx_A.pds vs xxxx.pds)
                if len(split[1]) > 8:
                    chain_id = split[1][5:6]
                else:
                    chain_id = 'A'

                match_str = pdb_id + '_' + chain_id

                if filename == match_str:

                    # replace unnecessary chars
                    for ch in ['[', ']', '(', ')']:
                        split[2] = split[2].replace(ch, '')

                    tpr_start = split[2].split(',')[0]
                    tpr_end = split[2].split(',')[1]

                    match_entry = {"rmsd": rmsd, "chain": chain_id, "tpr_start": tpr_start, "tpr_end": tpr_end}
                    match_entry_array = [match_entry]

                    # if filename with identical match positions already exists then filter out
                    if (pdb_id, tpr_start) not in unique_split:
                        unique_split.append((pdb_id, tpr_start))

                        # dictionary fill
                        if pdb_id in match_dict:
                            match_dict[pdb_id].append(match_entry)
                        else:
                            match_dict[pdb_id] = match_entry_array

            open_file.seek(0)

        parser = HhrParser
        print(len(match_dict))
        overlap_filtered = self.filter_overlaps(match_dict)
        print(len(overlap_filtered))

        parser.write_matches_json('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/', overlap_filtered)

    # Takes the matches dictionary and filters chains out if they are identical
    # Sequences provides via download directory
    @staticmethod
    def filter_chains(download_dir, matches_dict):

        unique_chains_dict = {}
        final_records = []

        skipped = 0

        # For all fasta files in directory
        for filename in os.listdir(download_dir):

            pdb_id = filename[0:4].upper()
            records = list(SeqIO.parse(download_dir + filename, 'fasta'))
            sequence_store = []
            unique_chains = []

            if pdb_id in matches_dict:

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
                                skipped += 1
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

        print('SKIPPED ', skipped, 'SEQUENCES BECAUSE OF IDENTICAL CHAINS')
        print('FINAL RECORDS ', len(final_records))

        sum_records = 0
        for key in matches_dict:
            sum_records += len(matches_dict[key])

        print('MATCHES DICT ENTRIES ', sum_records)

        return final_records, matches_dict

    @staticmethod
    def filter_chains_new(directory):

        seen_seq = []

        for file in os.listdir(directory):

            record = list(SeqIO.parse(directory + file, 'fasta'))

            if record[0].seq not in seen_seq:
                seen_seq.append(record[0].seq)
            else:
                print(record[0].id)
                subprocess.run(['rm', directory + file])

    @staticmethod
    def length_filter_new(directory):

        for file in os.listdir(directory):

            record = list(SeqIO.parse(directory + file, 'fasta'))

            if len(record[0].seq) > 750:
                subprocess.run(['rm', directory + file])
            elif len(record[0].seq) < 34:
                subprocess.run(['rm', directory + file])

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

    def aa_filter_new(self, directory):

        for file in os.listdir(directory):

            record = list(SeqIO.parse(directory + file, 'fasta'))

            if not self.unknown_aa_filter(record[0].seq):
                subprocess.run(['rm', directory + file])

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
    def filter_overlaps(matches_dict):

        all_starts = []

        for pdb_id in matches_dict:

            for entry in matches_dict[pdb_id]:

                all_starts.append(int(entry['tpr_start']))

            all_starts.sort()
            print(pdb_id)
            print('ALL', all_starts)

            to_delete = []
            for i in range(0, len(all_starts)):
                if i+1 < len(all_starts):
                    if all_starts[i+1] - all_starts[i] <= 34:

                        to_delete.append(all_starts[i+1])
                        all_starts.remove(all_starts[i + 1])

            all_starts = []

            print('DEL', to_delete)

            for index, entry in enumerate(matches_dict[pdb_id]):

                if int(entry['tpr_start']) in to_delete:
                    print('DELETED', matches_dict[pdb_id][index])
                    del matches_dict[pdb_id][index]

        return matches_dict

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
    def single_chains_fasta(download_directory, output_directory):

        final_records = []

        for file in os.listdir(download_directory):

            records = list(SeqIO.parse(download_directory + file, 'fasta'))

            for rec in records:
                final_records.append(rec)

        for record in final_records:

            SeqIO.write(record, output_directory + record.name[0:6].replace(':', '_') + '.fasta', "fasta")

    @staticmethod
    def read_final_sequences(directory):

        seqs = []
        records_seqs = []

        for file in os.listdir(directory):
            records = list(SeqIO.parse(directory + file, 'fasta'))
            for record in records:
                seqs.append(record.seq)
                records_seqs.append(record)

        return [seqs, records_seqs]

    # Cuts single protein in windows and returns fragments
    # This may improves the runtime since each protein is cutted and stored to the memory one by one
    def cutInWindows_single(self, record, window_size):

        sequence = record.seq

        # declare start of window
        window_start = int()
        protein_length = len(record.seq)

        # calculate end of first window by adding size to start position
        window_end = window_start + window_size

        fragments = []

        # slide window over file until end is reached
        while window_end <= protein_length:

            current_fragment = sequence[window_start:window_end]

            # increment window start pos
            window_start += 1

            # increment window end pos
            window_end += 1

            fragments.append(current_fragment)

        return fragments

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

        # For all training sequences
        for record in final_records:

            # Get raw PDB ID from the record ID
            pdb_id = record.id[0:4]
            chain_id = record.id[5]

            # create zero vector with length of sequence
            target_vector = [[0.0, 1.0]] * len(record.seq)

            # get all matches for particular PDB id
            for entry in match_dict[pdb_id]:

                if chain_id == entry['chain']:
                    # get start and end position out of dictionary
                    tpr_start = int(entry['tpr_start'])
                    tpr_end = int(entry['tpr_end'])
                    print(len(target_vector))
                    print(tpr_end)
                    # Set correct positions to 1
                    for i in range(tpr_start, tpr_end+1):

                        target_vector[i] = [1.0, 0.0]

            # Padding
            while len(target_vector) < self.padded_length:
                target_vector.append([0.0, 1.0])

            target_vector_list.append(target_vector)

        return target_vector_list
