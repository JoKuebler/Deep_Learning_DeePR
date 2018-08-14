class Encoder:

    def __init__(self):

        # Positive training set
        self.amino_acids = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
                            'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']

        self.unknown_amino_acids = ['X', 'x', 'B', 'b', 'Z', 'z', 'U', 'u']

    def encode(self, pos_data, neg_data):

        enc_data = []
        enc_lab = []

        return enc_data, enc_lab
