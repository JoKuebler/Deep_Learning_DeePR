from keras.layers import Conv1D, GlobalMaxPooling1D, BatchNormalization, MaxPooling1D, AveragePooling1D, Dropout
from keras.layers.core import Dense
from keras.models import Sequential
from keras.optimizers import SGD, Adam
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
        # self.optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
        self.optimizer = SGD(lr=0.005, momentum=0.9, nesterov=False)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            MaxPooling1D(pool_size=2),
            Conv1D(128, 9, padding='same', activation='relu', kernel_regularizer=l2(0.01)),
            MaxPooling1D(pool_size=2),
            Conv1D(256, 9, padding='same', activation='relu', kernel_regularizer=l2(0.01)),
            MaxPooling1D(pool_size=2),
            Conv1D(512, 8, padding='same', activation='relu', kernel_regularizer=l2(0.01)),
            MaxPooling1D(pool_size=2),
            Conv1D(1024, 3, padding='same', activation='relu', kernel_regularizer=l2(0.01)),
            MaxPooling1D(pool_size=1),
            Conv1D(144, 1, padding='same', activation='relu', kernel_regularizer=l2(0.01)),
            GlobalMaxPooling1D(),
            Dropout(0.1),
            Dense(2, activation='softmax', name='output_layer')
        ])

        # self.model = Sequential([
        #     self.input_layer,
        #     Conv1D(64, 9, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     MaxPooling1D(pool_size=3, strides=2, padding='same'),
        #     Conv1D(128, 9, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     MaxPooling1D(pool_size=3, strides=2, padding='same'),
        #     Conv1D(256, 8, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     Conv1D(256, 3, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     MaxPooling1D(pool_size=3, strides=2, padding='same'),
        #     Conv1D(512, 3, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     Conv1D(512, 3, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     MaxPooling1D(pool_size=3, strides=2, padding='same'),
        #     Conv1D(1024, 3, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     Conv1D(1024, 1, strides=1, padding='same', activation='relu'),
        #     BatchNormalization(),
        #     MaxPooling1D(pool_size=3, strides=2, padding='same'),
        #     AveragePooling1D(pool_size=1, strides=1, padding='valid'),
        #     GlobalMaxPooling1D(),
        #     Dense(2, activation='softmax', name='output_layer')
        # ])

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
        self.model.fit(data, target, batch_size=16, epochs=75, shuffle=True, verbose=2)

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

        # sort according to max N probs
        # len(predictions) stores all predictions in a nice way (prob, pos, seq)
        # gets passed to final prediction filter
        top_n = self.max_n(predictions, seqs, len(predictions))

        final_predictions, cut_out = self.filter_predictions(top_n, 0.65)

        print('Start\t\t', 'Sequence\t\t\t\t', 'End\t\t', 'Probability')
        # show the inputs and predicted outputs
        for i in range(len(final_predictions)):
            print(final_predictions[i][1], final_predictions[i][2], final_predictions[i][1]+33, final_predictions[i][0])

        print('Filtered out due to overlapping with higher probabilities')
        for i in range(len(cut_out)):
            print(cut_out[i][1], cut_out[i][2], cut_out[i][1] + 33, cut_out[i][0])

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
        plt.plot(history.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'], loc='upper left')
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
    def max_n(predictions, seqs, top):
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

            final_list.append((max1, index_getter.index(max1) + 1, str(seqs[index_getter.index(max1)])))
            pos_probs.remove(max1)

        return final_list

    def filter_predictions(self, predictions, threshold):
        """
        Complex method to filter output predictions so no overlapping TPRs are predicted and always the ones with
        the highest probabilities are shown to the user
        :param predictions: probabilities for each window
        :param threshold: probability threshold
        :return:
        """

        init_filter = []
        cut_out = []
        index = 1
        after = ''

        # Apply threshold
        [init_filter.append(x) if x[0] > 0.65 else '' for x in predictions]

        # Sort by position
        pos_sort = sorted(init_filter, key=lambda triple: triple[1])

        print(pos_sort)

        # in order to modify list in the loop we need a while loop here
        while index < len(pos_sort):

            # if element is deleted do not increment index
            deleted = False
            double_check = True

            # elements to compare
            current, before = pos_sort[index], pos_sort[index-1]

            if index < len(pos_sort)-1:
                after = pos_sort[index+1]

            # when only two predictions are made there cant be a comparison with three elements
            if len(pos_sort) > 2:
                # Special case of three windows are relative close to each other all three have to be considered
                if after[1] - before[1] <= 34 and current[1] - before[1] < 34 and current is not after:
                    # Case 1 the middle one is removed
                    if before[0] < current[0] < after[0]:
                        pos_sort, cut_out, deleted = self.rem_app_element(pos_sort, cut_out, current)
                        double_check = False
                    # Case 2 the right one is removed
                    elif after[0] < before[0] < current[0]:
                        pos_sort, cut_out, deleted = self.rem_app_element(pos_sort, cut_out, after)
                        double_check = False
                    else:
                        double_check = True

            # If no element is removed in the previous check
            if double_check:
                # when they are closer then 34 AA together
                if current[1] - before[1] < 34:
                    # remove only if propability is worse
                    if current[0] < before[0]:
                        pos_sort, cut_out, deleted = self.rem_app_element(pos_sort, cut_out, current)
                    elif current[0] > before[0]:
                        pos_sort, cut_out, deleted = self.rem_app_element(pos_sort, cut_out, before)

            # increment index only if no element is deleted
            if not deleted:
                index += 1

        # could also use initial findings but new variable makes it more clear here
        final = pos_sort
        return final, cut_out

    @staticmethod
    def rem_app_element(filter_list, cut_out, element):

        filter_list.remove(element)
        cut_out.append(element)

        return filter_list, cut_out, True

