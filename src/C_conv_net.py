from keras.layers import Conv1D, GlobalMaxPooling1D, BatchNormalization, MaxPooling1D, AveragePooling1D, Dropout
from keras.layers.core import Dense
from keras.models import Sequential
from keras.optimizers import SGD, Adam
from keras.regularizers import l2
from keras.models import model_from_json
from src.B_encoding import Encoder
from sklearn.model_selection import StratifiedKFold
from src.F_callbacks import Histories
import numpy as np
import matplotlib.pyplot as plt


class ConvolutionalNetwork:

    def __init__(self):

        # Define input layer
        self.input_layer = Conv1D(96, 9, padding='same', kernel_regularizer=l2(0.01), input_shape=(34, 20))

        # Define optimizer
        # self.optimizer = Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
        self.optimizer = SGD(lr=0.0001, momentum=0.9, nesterov=False)

        # Define model including layers and activation functions
        self.model = Sequential([
            self.input_layer,
            Dropout(0.5),
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
            Dense(128),
            Dropout(0.25),
            Dense(2, activation='softmax', name='output_layer')
        ])

        self.encoder = Encoder()

    def train_network(self, data, target, val_data, val_target, test_data=None, test_target=None):
        """
        Training of the network
        :param data: Positive and negative training data
        :param target: Target vector
        :param test_data: Positive and negative test data
        :param test_target: Target vector
        """

        # Model needs to be compiled before training
        self.model.compile(loss='binary_crossentropy', optimizer=self.optimizer, metrics=['accuracy'])

        # prepare callback
        histories = Histories()

        # Overview over model
        self.model.summary()

        # Train network to data with parameters: Batch Size, Epochs
        history = self.model.fit(data, target, batch_size=256, epochs=70, shuffle=True, verbose=2, validation_data=[val_data, val_target], callbacks=[histories])

        self.plot_loss(history)
        self.plot_acc(history)

        outputs = [layer.output for layer in self.model.layers]
        print(outputs)

        # Evaluate model and print results
        scores = self.model.evaluate(x=test_data, y=test_target)
        print("\n%s: %.2f%%" % (self.model.metrics_names[1], scores[1] * 100))
        print('[ERROR,ACCURACY]', scores)

    def predict(self, chunk, enc_seq, pro_length, seq_id):
        """
        Predict encoded sequences
        :param chunk: First part of raw sequences for output later
        :param enc_seq: Encoded sequences to predict
        :param pro_length: lengths of input proteins in order to map back to sequences
        :param seq_id: to output name/id of sequence predicted
        :return: Return predictions of CNN as well as padded data vector and target vector for refinement network
        """

        # Make predictions
        predictions = self.model.predict(enc_seq, batch_size=512, verbose=2)
        pred_per_chunk = 0

        start = 0

        # show the inputs and predicted outputs
        for i in range(len(pro_length)):

            cur_prot = predictions[start:start+pro_length[i]]

            print('>', seq_id[i])
            print('PREDICTED:', len(cur_prot))

            top_n = self.max_n(cur_prot, chunk[i], len(cur_prot))
            final_predictions, cut_out = self.filter_predictions(top_n, 0.9)
            pred_per_chunk += len(final_predictions)

            print('Start\t\t', 'Sequence\t\t\t\t', 'End\t\t', 'Probability')
            for x in range(len(final_predictions)):
                print(final_predictions[x][1], final_predictions[x][2], final_predictions[x][1]+33, final_predictions[x][0])

            print('Filtered out due to overlapping with higher probabilities')
            for c in range(len(cut_out)):
                print(cut_out[c][1], cut_out[c][2], cut_out[c][1] + 33, cut_out[c][0])
            print('\n\n')

            # To print all probabilities
            for idx, elem in enumerate(cur_prot):
                print(idx+1, chunk[i][idx], elem)

            start = start + pro_length[i]

        return pred_per_chunk

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

    def load_model(self, directory, file_name):
        """
        Loads model from json and h5 file in given directory
        :param directory: Directory to load from given via script parameter
        """

        # load json and create model
        json_file = open(directory + file_name + '.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)

        # load weights into new model
        loaded_model.load_weights(directory + file_name + '.h5')

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
            self.model.fit(data[train], target[train], validation_split=0.1, batch_size=128, epochs=75,
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

    @staticmethod
    def max_n(predictions, seqs, top):
        """
        prints out the top n probabilities with given position in the input sequence
        :param predictions: per window probs
        :param top: desired amount of elements to output
        :return:
        """

        final_list = []
        index_seen = []
        pos_probs = predictions[:, [0][0]].tolist()
        index_getter = pos_probs.copy()

        for i in range(0, top):
            max1 = 0

            for j in range(len(pos_probs)):
                if pos_probs[j] > max1:
                    max1 = pos_probs[j]

            indices = [i for i, x in enumerate(pos_probs) if x == max1]

            # This solves issue if there are two exact same probabilities
            # caused by same sequences
            if len(indices) > 1:

                # To get every position not always the first one
                for ind in indices:
                    if ind not in index_seen:
                        final_list.append((max1, ind + 1, str(seqs[index_getter.index(max1)])))
                        index_seen.append(ind)
                    else:
                        continue
                pos_probs.remove(max1)

            # normal case
            else:
                final_list.append((max1, index_getter.index(max1) + 1, str(seqs[index_getter.index(max1)])))
                pos_probs.remove(max1)

        return final_list

    @staticmethod
    def filter_predictions(predictions, threshold):
        """
        Complex method to filter output predictions so no overlapping TPRs are predicted and always the ones with
        the highest probabilities are shown to the user
        :param predictions: probabilities for each window
        :param threshold: probability threshold
        :return:
        """

        init_filter = []
        cut_out = []
        oldfit_counter = 0

        # Apply threshold
        [init_filter.append(x) if x[0] > threshold else '' for x in predictions]

        # Sort by position
        pos_sort = sorted(init_filter, key=lambda triple: triple[1])

        tmp, final = [], []

        # For all hits which survive the threshold
        for idx, elem in enumerate(pos_sort):

            old = tmp
            tmp, oldfit = [], []

            # List of a span of 34 aa starting with each element
            [tmp.append(x) if x[1] - 34 < elem[1] else '' for x in pos_sort[idx:]]

            # hit with highest prob in current list
            max_prob = max(tmp, key=lambda item: item[0])

            # if final list is empty
            if not final:
                final.append(max_prob)
            else:

                # This leads to no overlapping by checking the last element of the final list which is always sorted
                if max_prob not in final:
                    # if max prob is higher then last element of final and in same 34 range then delete old element in final and add new one
                    if final[-1][1] + 34 > max_prob[1] and final[-1][0] < max_prob[0]:
                        del final[-1]
                        final.append(max_prob)
                    # if max element is more then 34 amino acids distant of the last element in final it can be added
                    elif final[-1][1] + 34 <= max_prob[1]:
                        final.append(max_prob)

            # oldfit is a list which checks which of the hits in previous iteration fit in after one is deleted in previous step
            # if final has more then one element the second to last has also to be considered
            if len(final) > 1:
                [oldfit.append(x) if x[1] + 34 <= final[-1][1] and final[-2][1] + 34 < x[1] and x not in final else '' for x in old]
            else:
                [oldfit.append(x) if x[1] + 34 <= final[-1][1] and x not in final else '' for x in old]

            # if multiple elements could now fit take the one with the highest prob and add it to final
            if oldfit:
                max_prob_old_fitting = max(oldfit, key=lambda item: item[0])
                final.append(max_prob_old_fitting)
                oldfit_counter += 1
                #print('ADDED BY OLDFIT', max_prob_old_fitting)

            # Sort by position befor next iteration
            final = sorted(final, key=lambda triple: triple[1])

        # For printing the hits cutted out
        [cut_out.append(cut) if cut not in final else '' for cut in pos_sort]
        #print(oldfit_counter, 'TPRs ADDED BY MAGIC FILTER METHOD APPROACH')

        # Sort by position
        final = sorted(final, key=lambda triple: triple[1])

        # print('HITS:', len(final))
        # print('THRESHOLD:', threshold)

        return final, cut_out

    @staticmethod
    def rem_app_element(filter_list, cut_out, element):

        filter_list.remove(element)
        cut_out.append(element)

        return filter_list, cut_out, True

