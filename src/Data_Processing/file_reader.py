class FileReader:

    def __init__(self, positive_set, negative_set):

        # Positive training set
        self.positive_set = positive_set

        # Negative training set
        self.negative_set = negative_set

    def read_training_data(self):
        """
        Reads in training data
        :return: positive and negative set as list
        """
        with open(self.positive_set, 'r') as pos_file, open(self.negative_set, 'r') as neg_file:
            pos_data = pos_file.readlines()
            neg_data = neg_file.readlines()

        pos_file.close()
        neg_file.close()

        # Remove new line symbol at end of each line
        pos_data = [x.rstrip('\n') for x in pos_data]
        neg_data = [x.rstrip('\n') for x in neg_data]

        return pos_data, neg_data