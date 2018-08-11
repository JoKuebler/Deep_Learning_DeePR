from keras.layers import Conv1D, MaxPooling1D, GlobalMaxPooling1D, LSTM
from keras.layers.core import Dense, Dropout, Flatten
from keras.models import Sequential
from keras.optimizers import Adam
from keras.regularizers import l2
import numpy as np


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(64, 3, padding='same', kernel_regularizer=l2(0.01), input_shape=(34, 20))

        # Define optimizer
        self.optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            Conv1D(64, 5, padding='same', kernel_regularizer=l2(0.01)),
            Conv1D(64, 7, padding='same', kernel_regularizer=l2(0.01)),
            GlobalMaxPooling1D(),
            Dense(2, activation='softmax', name='output_layer')
        ])

    def compile_network(self):

        self.model.compile(loss='binary_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])

    # Trains the network
    def train_network(self, training_samples, training_labels, test_samples, test_labels):

        self.model.summary()

        # Fit network to data with parameters: Batch Size, Epochs
        self.model.fit(training_samples, training_labels, batch_size=100, epochs=50, shuffle=True, verbose=2)

        results = self.predict(test_samples[:242, :])

        new_file = open('/ebio/abt1_share/update_tprpred/data/PDB_Approach/convnet_predictions' + '.txt', 'w')

        for i in range(0, len(results)):
            new_file.write(str(results[i]) + '\n')

        new_file.close()

        acc = np.sum(np.argmax(results, axis=1) == np.argmax(test_labels[:242, :], axis=1))/242

        print(acc)

        # Evaluate model
        # scores = self.model.evaluate(test_samples, test_labels)
        #
        # print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
        #
        # print('[ERROR,ACCURACY]', scores)

    def predict(self, sequences):

        print(str(len(sequences)) + ' Sequences predicted')

        predictions = self.model.predict(sequences, batch_size=20, verbose=2)

        return predictions


