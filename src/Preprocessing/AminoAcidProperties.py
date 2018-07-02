# Created by Jonas Kuebler on 06/10/2018
# Class to add amino acid properties to the data
import numpy as np


class Properties:
    """Class to hold properties/features for amino acids
       Also provides method to encode amino acid at each position

     :param self.positional: List, amino acid one letter representation
     :param self.polarityValues: List, amino acid polarity values
     :param self.hydrophobicityValues: List, amino acid hydrophobicity values
     :param self.helixAbundance: List, amino acid values for favoring helices
     """

    def __init__(self):
        self.positional = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C',
                           'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        self.polarityValues = self.normalize_values([10.5, 10.4, 11.3, 13.0, 12.3, 9.2, 8.6, 11.6, 10.5, 5.5,
                                                    9.0, 8.0, 8.1, 5.9, 5.2, 4.9, 5.7, 5.2, 6.2, 5.4])
        self.hydrophobicityValues = self.normalize_values([0.60, 0.61, 1.15, 0.46, 0.47, 0.05, 0.05, 0.06, 0.0, 1.07,
                                                           0.07, 1.95, 0.61, 1.32, 2.22, 1.53, 1.18, 2.02, 1.88, 2.65])
        self.helixAbundance = self.normalize_values([1.21, 1.05, 1.23, 0.99, 1.59, 0.57, 0.76, 0.76, 1.27, 0.66,
                                                    0.43, 0.34, 1.41, 0.90, 1.09, 1.34, 1.30, 1.16, 0.61, 1.02])
        self.betaSheetAbundance = self.normalize_values([0.93, 0.87, 0.74, 0.54, 0.37, 0.75, 1.19, 0.89, 1.10, 1.19,
                                                        0.75, 0.55, 0.83, 1.70, 1.60, 1.30, 1.05, 1.38, 1.47, 1.37])
        self.betaTurnAbundance = self.normalize_values([1.01, 0.95, 1.19, 1.52, 0.95, 1.43, 0.98, 1.46, 0.96, 0.96,
                                                        1.56, 1.56, 0.74, 0.59, 0.47, 0.50, 0.60, 0.66, 1.14, 0.60])
        self.sidechainVolume = self.normalize_values([105.0, 79.0, 100.0, 40.0, 62.0, 29.3, 51.3, 58.7, 80.7, 44.6,
                                                      0.0, 41.9, 27.5, 71.5, 93.5, 93.5, 94.1, 115.5, 117.3, 145.5])

    # Encode amino acid sequence to multidimensional array with dimension (x, 34,23)
    def encode_positions(self, fragment):

        # vector without properties
        # empty_vector = [False, False, False, False, False, False, False, False, False, False,
        #                 False, False, False, False, False, False, False, False, False, False]

        # vector with properties
        empty_vector = [False, False, False, False, False, False, False, False, False, False,
                        False, False, False, False, False, False, False, False, False, False,
                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        # only properties
        # empty_vector = [0.0, 0.0, 0.0, 0.0]

        # need copy
        current_vector = empty_vector.copy()

        # array to store vector for each amino acid
        fragment_vector = []

        # encode each amino acid
        for aminoAcid in fragment:
            # find index in amino acid array and set char at index to 1
            current_vector[self.positional.index(aminoAcid.upper())] = True

            # add polarity, hydrophobicity and helix abundance
            current_vector[len(empty_vector)-1] = self.polarityValues[self.positional.index(aminoAcid.upper())]
            current_vector[len(empty_vector)-2] = self.hydrophobicityValues[self.positional.index(aminoAcid.upper())]
            current_vector[len(empty_vector)-3] = self.helixAbundance[self.positional.index(aminoAcid.upper())]
            current_vector[len(empty_vector)-4] = self.betaSheetAbundance[self.positional.index(aminoAcid.upper())]
            current_vector[len(empty_vector)-5] = self.betaTurnAbundance[self.positional.index(aminoAcid.upper())]
            current_vector[len(empty_vector)-6] = self.sidechainVolume[self.positional.index(aminoAcid.upper())]

            # current_vector[0] = self.polarityValues[self.positional.index(aminoAcid.upper())]
            # current_vector[1] = self.hydrophobicityValues[self.positional.index(aminoAcid.upper())]
            # current_vector[2] = self.helixAbundance[self.positional.index(aminoAcid.upper())]
            # current_vector[3] = self.betaSheetAbundance[self.positional.index(aminoAcid.upper())]

            # append to vector for whole fragment
            fragment_vector.append(np.array(current_vector))

            # reset current vector for next amino acid
            current_vector = empty_vector.copy()

        return fragment_vector

    @staticmethod
    def normalize_values(array_to_norm):

        array_normed = []

        for elem in array_to_norm:
            numerator = elem - min(array_to_norm)
            denominator = max(array_to_norm) - min(array_to_norm)
            array_normed.append(numerator/denominator)

        return array_normed

