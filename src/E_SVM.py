from sklearn.svm import SVC


class SVM:

    def __init__(self):

        self.svm = SVC(gamma=0.001, C=100.)

    def train_svm(self, training_data, target):

        d1_target = [1 if x[0] == 1 else 0 for x in target]

        self.svm.fit(self.d3_to_d2(training_data), d1_target)

    def predict(self, seqs, enc_seq, seq_id, target=None):

        predictions = self.svm.predict(self.d3_to_d2(enc_seq))

        # Print how many sequences are predicted (windows)
        print('>', seq_id)
        print(str(len(seqs)) + ' Sequences predicted')

        # show the inputs and predicted outputs
        for i in range(len(seqs)):
            print(i+1, seqs[i], predictions[i])

    @staticmethod
    def d3_to_d2(d3_array):

        nsamples, nx, ny = d3_array.shape
        d2_array = d3_array.reshape((nsamples, nx * ny))

        return d2_array

