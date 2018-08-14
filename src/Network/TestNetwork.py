import numpy as np
from keras.layers import Conv1D, MaxPool1D, GlobalMaxPooling1D, LSTM, MaxPooling1D
from keras.layers.core import Dense, Dropout, Flatten
from keras.models import Sequential
from keras.optimizers import Adam
from keras.regularizers import l2
from Bio import SeqIO

fragments_for_each_protein = []

for record in SeqIO.parse('/ebio/abt1/jkuebler/Desktop/1xnf.fasta', "fasta"):

    sequence = record.seq
    new_file = open('/ebio/abt1/jkuebler/Desktop/windows' + '.txt', 'w')

    # declare start of window
    window_start = int()
    protein_length = len(record.seq)

    # calculate end of first window by adding size to start position
    window_end = window_start + 34

    fragments = []

    # slide window over file until end is reached
    while window_end <= protein_length:
        current_fragment = sequence[window_start:window_end]

        # increment window start pos
        window_start += 1

        # increment window end pos
        window_end += 1

        new_file.write(str(current_fragment + '\n'))


