from keras.layers import Dropout, Flatten, Conv1D
from keras.layers.core import Dense
from keras.models import Sequential
from keras.optimizers import Adam, SGD
from keras.regularizers import l2
from keras.models import model_from_json
import numpy as np


class RefinementNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(64, 3, padding='same', kernel_regularizer=l2(0.01), input_shape=(717, 2))

        # Define optimizer
        self.optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
        # self.optimizer = SGD(lr=0.0001, decay=1e-6, momentum=0.9, nesterov=True)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            Conv1D(64, 5, padding='same', kernel_regularizer=l2(0.01)),
            Conv1D(64, 7, padding='same', kernel_regularizer=l2(0.01)),
            Conv1D(64, 9, padding='same', kernel_regularizer=l2(0.01)),
            Conv1D(64, 11, padding='same', kernel_regularizer=l2(0.01)),
            Dense(64),
            Dense(28),
            Dense(2, activation='softmax')
        ])

    def train_network(self, data, target, test_data=None, test_target=None):
        """
        Training of the network
        :param data: Positive and negative training data
        :param target: Target vector
        :param test_data: Positive and negative test data
        :param test_target: Target vector
        """

        # Model needs to be compiled before training
        self.model.compile(loss='binary_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])
        # Overview over model
        self.model.summary()

        # Train network to data with parameters: Batch Size, Epochs
        self.model.fit(data, target, validation_split=0.1, batch_size=100, epochs=30, shuffle=True, verbose=2)

        # Evaluate model and print results
        # scores = self.model.evaluate(test_data, test_target)
        # print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
        # print('[ERROR,ACCURACY]', scores)

    def predict(self, cnn_probs, seq_id_to_map):
        """
        Predict encoded sequences
        :param cnn_probs:
        :param seq_id_to_map:
        :return:
        """

        # Make predictions
        predictions = self.model.predict(cnn_probs, verbose=2)

        # Make results more readable
        for i in range(len(predictions)):
            print(seq_id_to_map[i])
            print('SUM ', sum(predictions[i]))
            print(predictions[i])
            for x in range(len(predictions[i])):
                print(predictions[i][x])
            # max_value = max(predictions[i])
            # print('MAX ', max_value)
            # for x in range(len(predictions[i])):
            #     if predictions[i][x] > 0.07:
            #         print(list(predictions[i]).index(predictions[i][x]))

        return predictions

    def save_model(self, directory):
        """
        Stores the model architecture and weights to json and h5 file
        :param directory: Directory to store given via script parameter
        """

        # Write model to json
        with open(directory + 'model.json', 'w') as json_file:
            json_file.write(self.model.to_json())

        # Write weights to h5 file
        self.model.save_weights(directory + 'model.h5')
        print('Model Saved to Disk!')

    def train_predict(self, conv_predictions, conv_labels, seq_id_to_map):

        self.train_network(conv_predictions[:820], conv_labels[:820])

        self.predict(conv_predictions[821:], seq_id_to_map[821:])






