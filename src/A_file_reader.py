from Bio import SeqIO


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

        prediction_seqs = []
        seq_ids = []

        # If no header in fasta/txt file then this is the way data gets read in
        if len(records) == 0:
            with open(input_file) as file:
                content = file.readlines()

                stripped = [x.strip() for x in content]
                seq_ids = [x for x in content]

                prediction_seqs.append(stripped)

                return prediction_seqs, seq_ids
        else:
            for record in records:

                sequence = record.seq
                seq_id = record.id
                # chain_id = record.id[5]

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

                prediction_seqs.append(fragments)
                seq_ids.append(seq_id)

            return prediction_seqs, seq_ids






