from Bio import SeqIO
import os


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
    def pdb_pos_mapper(map_file_dir, match_dict):
        """
        Maps position in MASTER match files to real positions in PDB files
        :param map_file_dir: directory containing each master match
        :param match_dict: dictionary with old/wrong positions to be replaced
        :return: new dictionary
        """

        for file in os.listdir(map_file_dir):

            with open(map_file_dir + file) as content:
                head = [next(content) for x in range(2)]

                line1, line2 = head[0], head[1]

                one_list, two_list = line1.split(), line2.split()
                print(one_list)
                print(two_list)

                fake_pos = int(one_list[-1].split(',')[0].replace('[', '').replace('(', ''))
                true_pos = two_list[5]


if __name__ == '__main__':

    reader = FileReader()
    reader.pdb_pos_mapper('/tmp/jonas/run/results/all_matches/', '/tmp/jonas/run/results/match_dict_merged.json')








