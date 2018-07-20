from keras.layers import Conv1D
from keras.layers.core import Dense, Dropout
from keras.models import Sequential
from keras.optimizers import Adam
from keras.regularizers import l2


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(68, 34, padding='same', kernel_regularizer=l2(0.01), input_shape=(500, 20))

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

        self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer, metrics=['total_accuracy'])

