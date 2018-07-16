from keras.layers import Conv1D
from keras.layers.core import Dense, Dropout
from keras.models import Sequential
from keras.optimizers import adam
from keras.regularizers import l2


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(68, 34, padding='same', kernel_regularizer=l2(0.01), input_shape=(500,20))

        # Define optimizer
        self.optimizer = adam(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            Dropout(0.5),
            Dense(128),
            Dense(2, activation='softmax', name='output_layer')
        ])

    def compile_network(self):

        self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer, metrics=['total_accuracy'])

