from keras.layers import Conv1D, GlobalMaxPooling1D, Input, BatchNormalization, merge, TimeDistributed, LSTM, Bidirectional
from keras.layers.core import Dense
from keras.models import Sequential
from keras.optimizers import Adam
from keras.regularizers import l2
from keras.models import model_from_json, Model
from keras.layers.merge import concatenate
from src.B_encoding import Encoder


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(64, 3, padding='same', kernel_regularizer=l2(0.01), input_shape=(34, 20))
        input_one = Input(shape=(34, 20), name='inp')

        # Define optimizer
        self.optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)

        # Define model including layers and activation functions
        a = Conv1D(64, 3, activation='sigmoid', padding='same', kernel_regularizer=l2(0.0001))(input_one)
        a_b = BatchNormalization()(a)
        b = Conv1D(64, 5, activation='sigmoid', padding='same', kernel_regularizer=l2(0.0001))(input_one)
        b_b = BatchNormalization()(b)
        c = Conv1D(64, 7, activation='sigmoid', padding='same', kernel_regularizer=l2(0.0001))(input_one)
        c_b = BatchNormalization()(c)
        d = Conv1D(64, 9, activation='sigmoid', padding='same', kernel_regularizer=l2(0.0001))(input_one)
        d_b = BatchNormalization()(d)
        e = Conv1D(64, 9, activation='sigmoid', padding='same', kernel_regularizer=l2(0.0001))(input_one)
        e_b = BatchNormalization()(e)
        x = concatenate([a_b, b_b, c_b, d_b, e_b], axis=-1)

        t = Dense(180, activation='relu', kernel_regularizer=l2(0.0001))(x)

        glob = GlobalMaxPooling1D()(t)

        out = Dense(2, activation='softmax', name='output_layer')(glob)

        self.model = Model(inputs=input_one, outputs=out)

        self.encoder = Encoder()

    def train_network(self, data, target, test_data=None, test_target=None):
        """
        Training of the network
        :param data: Positive and negative training data
        :param target: Target vector
        :param test_data: Positive and negative test data
        :param test_target: Target vector
        """

        # Model needs to be compiled before training
        self.model.compile(loss='categorical_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])
        # Overview over model
        self.model.summary()

        # Train network to data with parameters: Batch Size, Epochs
        self.model.fit(data, target, validation_split=0.1, batch_size=100, epochs=50, shuffle=True, verbose=2)

        # Evaluate model and print results
        scores = self.model.evaluate(test_data, test_target)
        print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
        print('[ERROR,ACCURACY]', scores)

    def predict(self, seqs, enc_seq, seq_id=None, chain_id=None):
        """
        Predict encoded sequences
        :param seqs: Raw sequences for output later
        :param enc_seq: Encoded sequences to predict
        :param seq_id: Sequence ID to find the true TPR in the tpr json for target vector
        :param chain_id: Chain ID to find the true TPR in the tpr json for target vector
        :return: Return predictions of CNN as well as padded data vector and target vector for refinement network
        """

        # Print how many sequences are predicted (windows)
        print(str(len(seqs)) + ' Sequences predicted')

        # Make predictions
        predictions = self.model.predict(enc_seq, batch_size=20, verbose=2)

        # show the inputs and predicted outputs
        for i in range(len(seqs)):
            print(i+1, seqs[i], predictions[i])

        if seq_id is not None:
            # for each prediction made by the CNN pass the data to encode it for the refinement network
            padded_refine_data, target_vector = self.encoder.create_refine_data(predictions, seq_id, chain_id)
        else:
            padded_refine_data, target_vector = [], []

        return predictions, padded_refine_data, target_vector

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


