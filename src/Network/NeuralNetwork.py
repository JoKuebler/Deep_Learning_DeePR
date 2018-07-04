# Created by Jonas Kuebler on 05/18/2018
import numpy
from keras import optimizers
from keras.layers.core import Dense, Flatten, Dropout
from keras.models import Sequential
from keras.models import model_from_json
from sklearn.model_selection import StratifiedKFold


# This class defines the neural network and provides training and prediction methods
# Adjust Parameters and Layers for Network used here
class NeuralNetwork:

    def __init__(self):

        # Define input layer
        self.inputLayer = Dense(34, input_shape=(34, 26), activation='relu')

        # Define optimizer
        self.optimizer = optimizers.SGD(lr=0.0001, decay=1e-6, momentum=0.9, nesterov=True)

        # self.optimizer = optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.inputLayer,
            Dropout(0.2),
            Dense(34, activation='relu'),
            Flatten(),
            Dense(1, activation='sigmoid')
        ])

        # Define model including layers and activation functions
        # self.model = Sequential([
        #     LSTM(34,  input_shape=(34, 21), return_sequences=True),
        #     Flatten(),
        #     Dense(1, activation='sigmoid')
        # ])

    # Compiles the network
    def compile_network(self):

        self.model.compile(loss='binary_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])

    # Trains the network
    def train_network(self, training_samples, training_labels, test_samples, test_labels):

        self.model.summary()

        # Fit network to data with parameters: Batch Size, Epochs
        self.model.fit(training_samples, training_labels, validation_split=0.2, batch_size=250, epochs=15, shuffle=True, verbose=2)

        # Evaluate model
        scores = self.model.evaluate(test_samples, test_labels)

        print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))

        print('[ERROR,ACCURACY]', scores)

    def predict(self, fragments):

        print(str(len(fragments)) + ' Sequences predicted')

        predictions = self.model.predict(fragments, batch_size=50, verbose=2)

        return predictions

    def save_model(self):

        model_json = self.model.to_json()

        with open('/ebio/abt1_share/update_tprpred/code/PycharmProjects/TrainingData/src/NetworkData/model.json', 'w') as json_file:
            json_file.write(model_json)

        self.model.save_weights('/ebio/abt1_share/update_tprpred/code/PycharmProjects/TrainingData/src/NetworkData/model.h5')
        print('Model Saved to Disk!')

    # Loads model data from json and weight file
    def load_model(self, model_json, weights_file):

        # load json and create model
        json_file = open(model_json, 'r')

        loaded_model_json = json_file.read()

        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        # load weights into new model
        loaded_model.load_weights(weights_file)

        self.model = loaded_model

    def cross_validate(self, training_samples, training_labels):

        kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
        cvscores = []

        for train, test in kfold.split(training_samples, training_labels):

            # Fit network to data with parameters: Batch Size, Epochs
            self.model.fit(training_samples[train], training_labels[train], validation_split=0.1, batch_size=75, epochs=15,
                           shuffle=True, verbose=2)

            # Evaluate model
            scores = self.model.evaluate(training_samples[test], training_labels[test])

            print("%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
            cvscores.append(scores[1] * 100)

        print('FINISHED')
        print("%.2f%% (+/- %.2f%%)" % (numpy.mean(cvscores), numpy.std(cvscores)))





