# Created by Jonas Kuebler on 06/10/2018
# Class to add amino acid properties to the data
import numpy as np


class Properties:

    def __init__(self):
        self.positional = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
                           'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        self.polarityValues = [10.5, 10.4, 11.3, 13.0, 12.3, 9.2, 8.6, 11.6, 10.5, 5.5,
                               9.0, 8.0, 8.1, 5.9, 5.2, 4.9, 5.7, 5.2, 6.2, 5.4]
        self.hydrophobicityValues = [0.60, 0.61, 1.15, 0.46, 0.47, 0.05, 0.05, 0.06, 0.0, 1.07,
                                     0.07, 1.95, 0.61, 1.32, 2.22, 1.53, 1.18, 2.02, 1.88, 2.65]
        self.helixAbundance = [1.21, 1.05, 1.23, 0.99, 1.59, 0.57, 0.76, 0.76, 1.27, 0.66,
                               0.43, 0.34, 1.41, 0.90, 1.09, 1.34, 1.30, 1.16, 0.61, 1.02]

    # Encode amino acid sequence to multidimensional array with dimension (x, 34,23)
    def encode_positions(self, fragment):

        # vector without properties
        # empty_vector = [False, False, False, False, False, False, False, False, False, False,
        #                 False, False, False, False, False, False, False, False, False, False]

        # vector with properties
        empty_vector = [False, False, False, False, False, False, False, False, False, False,
                        False, False, False, False, False, False, False, False, False, False, 0.0]

        # only properties
        # empty_vector = [0.0]

        # need copy
        current_vector = empty_vector.copy()

        # array to store vector for each amino acid
        fragment_vector = []

        # encode each amino acid
        for aminoAcid in fragment:
            # find index in amino acid array and set char at index to 1
            current_vector[self.positional.index(aminoAcid.upper())] = True

            # add polarity, hydrophobicity and helix abundance
            # current_vector[len(empty_vector)-3] = self.polarityValues[self.positional.index(aminoAcid.upper())]
            # current_vector[len(empty_vector)-2] = self.hydrophobicityValues[self.positional.index(aminoAcid.upper())]
            current_vector[len(empty_vector)-1] = self.helixAbundance[self.positional.index(aminoAcid.upper())]

            #current_vector[0] = self.polarityValues[self.positional.index(aminoAcid.upper())]
            #current_vector[1] = self.hydrophobicityValues[self.positional.index(aminoAcid.upper())]
            #current_vector[0] = self.helixAbundance[self.positional.index(aminoAcid.upper())]

            # append to vector for whole fragment
            fragment_vector.append(np.array(current_vector))

            # reset current vector for next amino acid
            current_vector = empty_vector.copy()

        return fragment_vector

