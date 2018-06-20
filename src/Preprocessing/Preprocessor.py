# Created by Jonas Kuebler on 05/15/2018
from Bio import SeqIO
from keras.preprocessing.text import one_hot
from src.Preprocessing.AminoAcidProperties import Properties
import numpy as np


# This class provides multiple functions to preprocess data such as
# Cutting input Protein into windows with size of choice
# Encode Amino Acids and Fragments to Zero-One Vector
# Decode Zero-One Vector to Amino Acid sequence
# Helper functions for creating labels, reading lines etc.
class Preprocessor:

    def __init__(self, inputFile):
        self.inputFile = inputFile
        self.ids = []
        self.probabilites = {}
        self.property_adder = Properties()

    # to find most likely fragment to be TPR protein has to be divided first
    def cutInWindows(self, window_size):

        fragments_for_each_protein = []

        for record in SeqIO.parse(self.inputFile, "fasta"):

            sequence = record.seq

            self.ids.append(record.id)

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

            fragments_for_each_protein.append(fragments)

        return fragments_for_each_protein

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

    # Text encoding fragments
    @staticmethod
    def text_encode_fragment(fragment):

        result = one_hot(" ".join(fragment), 21)

        return result

    # Helper to read in File Line by line
    # returns list containing the lines
    @staticmethod
    def read_line_files(file):

        with open(file) as input_file:
            content = input_file.readlines()

        input_file.close()

        return content

    # Helper to encode Fragments to multidimensional Array
    def encoding_helper(self, fragments_per_protein):

        encoded_fragments = []
        encoded_fragments_per_protein = []

        for fragment_list in fragments_per_protein:

            for fragment in fragment_list:

                if 'x' not in fragment:
                    if 'X' not in fragment:
                        encoded = self.property_adder.encode_positions(fragment.rstrip())
                        encoded_fragments.append(encoded)

            encoded_fragments_per_protein.append(encoded_fragments)
            encoded_fragments = []

        return encoded_fragments_per_protein

    # Helper to encode Fragments to multidimensional Array single protein to safe memory
    def encoding_helper_single(self, fragments):

        encoded_fragments = []

        for fragment in fragments:

            if 'x' not in fragment:
                if 'X' not in fragment:
                    if 'u' not in fragment:
                        if 'U' not in fragment:
                            encoded = self.property_adder.encode_positions(fragment.rstrip())
                            encoded_fragments.append(encoded)

        return encoded_fragments

    def encoding_training_helper(self, file_content_list):

        encoded_fragments = []

        for fragment in file_content_list:
            encoded = self.property_adder.encode_positions(fragment.rstrip())
            encoded_fragments.append(encoded)

        return encoded_fragments

    # Helper to create labels
    @staticmethod
    def create_labels(amount_pos, amount_neg):

        pos_labels = [1] * amount_pos
        neg_labels = [0] * amount_neg
        labels = np.array(pos_labels + neg_labels).reshape(amount_pos + amount_neg, 1)

        return labels

    # filter and print results so no overlaps are predicted
    def filter_print_results(self, results, threshold):

        # counter to keep index
        counter = 0
        # index for filtering
        index = 1

        # index of TPRs found before filtering
        initial_findings_index = []

        # figure out indexes and probabilities and store them
        for result in results:
            counter += 1
            if result > threshold:

                initial_findings_index.append(counter)
                self.probabilites[counter] = result

        # print(initial_findings_index)
        # print(prob_dict)

        # in order to modify list in the loop we need a while loop here
        while index < len(initial_findings_index):

            # if element is deleted do not increment index
            deleted = False

            # elements to compare
            current = initial_findings_index[index]
            before = initial_findings_index[index-1]

            # when they are closer then 34 AA together
            if current - before < 34:

                # remove only if propability is worse
                if self.probabilites[current] < self.probabilites[before]:
                    initial_findings_index.remove(current)
                    deleted = True
                elif self.probabilites[current] > self.probabilites[before]:
                    initial_findings_index.remove(before)
                    deleted = True

            if not deleted:
                index += 1

        return initial_findings_index


