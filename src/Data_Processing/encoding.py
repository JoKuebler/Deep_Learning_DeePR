import numpy as np


class Encoder:

    def __init__(self):

        # Positive training set
        self.amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
                            'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']

        self.unknown_aa = ['X', 'x', 'B', 'b', 'Z', 'z', 'U', 'u']

    def encode(self, pos_data, neg_data):
        """
        Call Helper function to encode raw list of sequences into 34x20 vectors
        :param pos_data: positive training sequences as list
        :param neg_data: negative training sequences as list
        :return: numpy array of encoded sequences
        """
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
        target = np.concatenate((target_pos, target_pos), axis=0)

        return target

