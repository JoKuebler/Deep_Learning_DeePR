from keras.layers import Conv1D, GlobalMaxPooling1D
from keras.layers.core import Dense
from keras.models import Sequential
from keras.optimizers import Adam
from keras.regularizers import l2
from keras.models import model_from_json


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
        scores = self.model.evaluate(test_data, test_target)
        print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
        print('[ERROR,ACCURACY]', scores)

    def predict(self, seqs, enc_seq):
        """
        Predict encoded sequences
        :param seqs: Raw sequences for output later
        :param enc_seq: Encoded sequences to predict
        :return:
        """

        print(str(len(seqs)) + ' Sequences predicted')

        # Make predictions
        predictions = self.model.predict(enc_seq, batch_size=20, verbose=2)

        # show the inputs and predicted outputs
        for i in range(len(seqs)):
            print(i+1, seqs[i], predictions[i])

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

    def load_model(self, directory):
        """
        Loads model from json and h5 file in given directory
        :param directory: Directory to load from given via script parameter
        """

        # load json and create model
        json_file = open(directory + 'model.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)

        # load weights into new model
        loaded_model.load_weights(directory + 'model.h5')

        # Set model
        self.model = loaded_model


