from Bio import SeqIO
import os
import json


class FileReader:

    @staticmethod
    def read_training_data(positive_set, negative_set):
        """
        Reads in training data
        :return: positive and negative set as list
        """
        with open(positive_set, 'r') as pos_file, open(negative_set, 'r') as neg_file:
            pos_data = pos_file.readlines()
            neg_data = neg_file.readlines()

        pos_file.close()
        neg_file.close()

        # Remove new line symbol at end of each line
        pos_data = [x.rstrip('\n') for x in pos_data]
        neg_data = [x.rstrip('\n') for x in neg_data]

        return pos_data, neg_data

    @staticmethod
    def read_pred_data(input_file, window_size, step_size):

        records = list(SeqIO.parse(input_file, 'fasta'))

        for record in records:

            sequence = record.seq
            seq_id = record.id[:4]
            chain_id = record.id[5]

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
                window_start += step_size

                # increment window end pos
                window_end += step_size

                fragments.append(current_fragment)

            return fragments, seq_id, chain_id

    @staticmethod
    def cut_match_dict(final_fasta_dir, match_json):

        # Open match dict
        with open(match_json) as file:
            data = json.load(file)

        print(len(data))

        ids = []

        for file in os.listdir(final_fasta_dir):

            pdb_id = file[0:4]
            ids.append(pdb_id)

        for key in list(data):
            if key in ids:
                continue
            else:
                del data[key]

        # store hhr in json file in the directory
        output_file = open(final_fasta_dir + 'merged_cut.json', 'w+')

        output_file.write(str(data).replace('\'', '"'))

        output_file.close()

    @staticmethod
    def pdb_pos_mapper(map_file_dir, match_json):
        """
        Maps position in MASTER match files to real positions in PDB files
        :param map_file_dir: directory containing each master match
        :param match_json: dictionary with old/wrong positions to be replaced
        :return: new dictionary
        """
        # Open match dict
        with open(match_json) as file:
            data = json.load(file)

        # Look at each file in the match directory
        for file in os.listdir(map_file_dir):

            # Get first two lines of each match file
            with open(map_file_dir + file) as content:
                head = [next(content) for x in range(2)]

                # Define Line 1 & 2 and split it into list
                line1, line2 = head[0], head[1]
                one_list, two_list = line1.split(), line2.split()

                # parse relevant information out of file
                fake_pos = int(one_list[-1].split(',')[0].replace('[', '').replace('(', ''))
                pdb_id = one_list[-2].split('/')[-1][0:4].upper()
                chain = two_list[4]
                true_pos = two_list[5]

                # Check if PDB ID is in match dict
                if pdb_id in data:
                    if len(two_list[4]) > 4:
                        continue
                    else:
                        for entry in data[pdb_id.upper()]:
                            if int(entry["tpr_start"]) == fake_pos:
                                if chain == entry["chain"]:
                                    entry["tpr_start"] = str(true_pos)
                                    entry["tpr_end"] = str(int(true_pos) + 34)
                                    index = data[pdb_id].index(entry)

                                    data[pdb_id][index] = entry

        # store hhr in json file in the directory
        output_file = open('updated.json', 'w+')

        output_file.write(str(data).replace('\'', '"'))

        output_file.close()


if __name__ == '__main__':

    reader = FileReader()
    reader.pdb_pos_mapper('/tmp/jonas/run/results/all_matches/', '/tmp/jonas/run/results/merged_cut.json')

    # reader.cut_match_dict('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/second_update_single_chains/final_fasta_sets_combined/',
    #                       '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/second_update_single_chains/match_dict_merged.json')








