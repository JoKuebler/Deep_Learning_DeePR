import numpy as np
import json


class Encoder:

    def __init__(self):

        self.amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
                            'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        self.side_chain_volume = self.normalize_values([105.0, 79.0, 100.0, 40.0, 62.0, 29.3, 51.3, 58.7, 80.7, 44.6,
                                                        0.0, 41.9, 27.5, 71.5, 93.5, 93.5, 94.1, 115.5, 117.3, 145.5])
        self.hydrophobicity = self.normalize_values([0.60, 0.61, 1.15, 0.46, 0.47, 0.05, 0.05, 0.06, 0.0, 1.07,
                                                     0.07, 1.95, 0.61, 1.32, 2.22, 1.53, 1.18, 2.02, 1.88, 2.65])

        self.unknown_aa = ['X', 'x', 'B', 'b', 'Z', 'z', 'U', 'u']

        # TPR JSON which stores information about all training sequences
        self.tpr_info = self.read_json('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/MASTER_sets/first_set/tpr_info.json')

    def encode(self, pos_data, neg_data=None):
        """
        Call Helper function to encode raw list of sequences into 34x20 vectors
        :param pos_data: positive training sequences as list
        :param neg_data: negative training sequences as list
        :return: numpy array of encoded sequences
        """
        neg_data = [] if not neg_data else neg_data
        # List comprehension to encode fragments after checking for unknown amino acids
        enc_fragments = [self.enc_positions(fragment.rstrip()) if not any(x in fragment for x in self.unknown_aa)
                         else print('ENCODING FAILED') for fragment in pos_data + neg_data]

        target_vector = self.create_labels(len(pos_data), len(neg_data))

        return np.asarray(enc_fragments), target_vector

    def enc_positions(self, fragment):
        """
        Encodes sequences into 34x20 vector
        :param fragment: single fragment
        :return: vector for sequence
        """

        # define zero vector
        zero_vector = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        current_vector = zero_vector.copy()

        # array to store vector for each amino acid
        fragment_vector = []

        # encode each amino acid
        for aminoAcid in fragment:

            # find index in amino acid array and set char at index to 1
            current_vector[self.amino_acids.index(aminoAcid.upper())] = 1.0
            # current_vector[len(zero_vector) - 1] = self.side_chain_volume[self.amino_acids.index(aminoAcid.upper())]
            # current_vector[len(zero_vector) - 2] = self.hydrophobicity[self.amino_acids.index(aminoAcid.upper())]

            # append to vector for whole fragment
            fragment_vector.append(np.array(current_vector))

            # reset vector for next amino acid
            current_vector = zero_vector.copy()

        return fragment_vector

    @staticmethod
    def create_labels(amount_pos, amount_neg):
        """
        Creats target vector for sequence being/not a TPR
        :param amount_pos: Amount of positive labels
        :param amount_neg: Amount of negative labels
        :return: Target vector in 34x2 format ([prob being TPR, prob no TPR])
        """
        # Create positive target vector
        target_pos = np.zeros((amount_pos, 2))
        target_pos[:, 0] = 1

        # Create negative target vector
        target_neg = np.zeros((amount_neg, 2))
        target_neg[:, 1] = 1

        # Final target vector
        target = np.concatenate((target_pos, target_neg), axis=0)

        return target

    def create_refine_data(self, predictions, seq_id, chain_id):

        data1 = np.array([x[0] for x in predictions])
        data2 = np.array([(1 - x[0]) for x in predictions])

        data = []

        for i in range(len(data1)):
            current = [data1[i], data2[i]]
            data.append(current)

        while len(data) < 717:
            data.append([0, 1])

        # padded_data = np.pad(data, (0, 717-len(data)), 'constant', constant_values=(0, 0))

        padded_data = data

        target_vector = (717, 2)
        target_vector = np.zeros(target_vector)

        for entry in self.tpr_info[seq_id]:
            if entry["chain"] == chain_id:
                for i in range(len(target_vector)):
                    if i == int(entry["tpr_start"])-1:
                        target_vector[i] = [1, 0]
                    else:
                        target_vector[i] = [0, 1]

        # if all([v == 0 for v in target_vector]):
        #     print(seq_id, ' LEADS TO A ZERO TARGET VECTOR')

        return padded_data, target_vector

    @ staticmethod
    def read_json(json_file):

        with open(json_file) as file:
            data = json.load(file)

        return data

    @staticmethod
    def normalize_values(array_to_norm):

        max_val, min_val = max(array_to_norm), min(array_to_norm)

        array_normed = [(elem - min_val)/(max_val-min_val) for elem in array_to_norm]

        return array_normed

