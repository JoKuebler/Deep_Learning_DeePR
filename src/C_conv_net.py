from keras.layers import Conv1D, GlobalMaxPooling1D
from keras.layers.core import Dense
from keras.models import Sequential
from keras.optimizers import SGD
from keras.regularizers import l2
from keras.models import model_from_json
from src.B_encoding import Encoder
from sklearn.model_selection import StratifiedKFold
import numpy as np
import matplotlib.pyplot as plt


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(96, 9, padding='valid', kernel_regularizer=l2(0.01), input_shape=(34, 20))

        # Define optimizer
        self.optimizer = SGD(lr=0.009, momentum=0.9, nesterov=False)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            Conv1D(128, 9, padding='valid', kernel_regularizer=l2(0.01)),
            Conv1D(256, 9, padding='valid', kernel_regularizer=l2(0.01)),
            Conv1D(512, 8, padding='valid', kernel_regularizer=l2(0.01)),
            Conv1D(1024, 3, padding='valid', kernel_regularizer=l2(0.01)),
            Conv1D(144, 1, padding='valid', kernel_regularizer=l2(0.01)),
            GlobalMaxPooling1D(),
            Dense(2, activation='softmax', name='output_layer')
        ])

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
        self.model.compile(loss='binary_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])
        # Overview over model
        self.model.summary()

        # Train network to data with parameters: Batch Size, Epochs
        history = self.model.fit(data, target, batch_size=100, epochs=70, shuffle=True, verbose=2)

        # self.plot_loss(history)
        # self.plot_acc(history)

        outputs = [layer.output for layer in self.model.layers]
        print(outputs)

        # Evaluate model and print results
        scores = self.model.evaluate(test_data, test_target)
        print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
        print('[ERROR,ACCURACY]', scores)

    def predict(self, seqs, enc_seq, seq_id, target=None, chain_id=None):
        """
        Predict encoded sequences
        :param seqs: Raw sequences for output later
        :param enc_seq: Encoded sequences to predict
        :param seq_id: Sequence ID to find the true TPR in the tpr json for target vector
        :param chain_id: Chain ID to find the true TPR in the tpr json for target vector
        :return: Return predictions of CNN as well as padded data vector and target vector for refinement network
        """

        # Print how many sequences are predicted (windows)
        print('>', seq_id)
        print(str(len(seqs)) + ' Sequences predicted')

        # Make predictions
        predictions = self.model.predict(enc_seq, batch_size=20, verbose=2)

        # show the inputs and predicted outputs
        for i in range(len(seqs)):
            print(i+1, seqs[i], predictions[i])

        if target is not None:
            # To calculate F1 score target has to be given
            self.f1_score(0.9, predictions, target)

        # print max 10 probs
        self.max_n(predictions, 10)

        # For Refine LSTM
        # if seq_id is not None:
        #     # for each prediction made by the CNN pass the data to encode it for the refinement network
        #     padded_refine_data, target_vector = self.encoder.create_refine_data(predictions, seq_id, chain_id)
        # else:
        #     padded_refine_data, target_vector = [], []

        return predictions
        # padded_refine_data, target_vector

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

    def cross_validate(self, data, target):
        """
        k-fold cross validation of model with given training set
        :param data: training samples
        :param target: training labels
        :return:
        """

        kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
        cvscores = []

        one_dim_target = [1 if x[0] == 1 else 0 for x in target]
        print(one_dim_target)

        for train, test in kfold.split(data, one_dim_target):

            # Fit network to data with parameters: Batch Size, Epochs
            self.model.fit(data[train], target[train], validation_split=0.1, batch_size=75, epochs=55,
                           shuffle=True, verbose=2)

            # Evaluate model
            scores = self.model.evaluate(data[test], target[test])

            print("%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
            cvscores.append(scores[1] * 100)

        print('FINISHED')
        print("%.2f%% (+/- %.2f%%)" % (np.mean(cvscores), np.std(cvscores)))

    @staticmethod
    def plot_loss(history):

        plt.plot(history.history['loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train'], loc='upper left')
        plt.show()

    @staticmethod
    def plot_acc(history):

        plt.plot(history.history['acc'])
        plt.title('model accuracy')
        plt.ylabel('accuracy')
        plt.xlabel('epoch')
        plt.legend(['train'], loc='upper left')
        plt.show()

    # Implemented so that the validation set is given via input to predict
    def f1_score(self, threshold, predictions, labels):
        """
        Calculates the F1 score for validation set
        :param threshold:
        :param predictions:
        :param labels:
        :return:
        """

        tp = predictions[:, 0] > threshold

        one_dim_target = [1 if x[0] == 1 else 0 for x in labels]

        tp_count = 0
        fp_count = 0
        fn_count = 0

        for index, elem in enumerate(tp):
            if tp[index] == one_dim_target[index] and elem:
                tp_count += 1
            elif tp[index] != one_dim_target[index] and not elem:
                fn_count += 1
            elif tp[index] != one_dim_target[index] and elem:
                fp_count += 1

        f1 = (2 * tp_count)/(2 * tp_count + fn_count)

        print('F1 Score: ', f1)

    @staticmethod
    def max_n(predictions, top):
        """
        prints out the top n probabilities with given position in the input sequence
        :param predictions: per window probs
        :param top: desired amount of elements to output
        :return:
        """

        final_list = []
        pos_probs = predictions[:, [0][0]].tolist()
        index_getter = pos_probs.copy()

        for i in range(0, top):
            max1 = 0

            for j in range(len(pos_probs)):
                if pos_probs[j] > max1:
                    max1 = pos_probs[j]

            final_list.append((max1, index_getter.index(max1) + 1))
            pos_probs.remove(max1)

        print(final_list)

