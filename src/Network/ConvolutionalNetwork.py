from keras.layers import Conv1D
from keras.layers.core import Dense, Dropout
from keras.models import Sequential
from keras.optimizers import Adam
from keras.regularizers import l2


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(68, 34, padding='same', kernel_regularizer=l2(0.01), input_shape=(750, 20))

        # Define optimizer
        self.optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            Dropout(0.5),
            Dense(128),
            Dense(2, activation='softmax', name='output_layer')
        ])

    def compile_network(self):

        self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])

    # Trains the network
    def train_network(self, training_samples, training_labels):

        self.model.summary()

        # Fit network to data with parameters: Batch Size, Epochs
        self.model.fit(training_samples, training_labels, batch_size=64, epochs=20, shuffle=True, verbose=2)

    def predict(self, sequences):

        print(str(len(sequences)) + ' Sequences predicted')

        predictions = self.model.predict(sequences, batch_size=64, verbose=2)

        return predictions


