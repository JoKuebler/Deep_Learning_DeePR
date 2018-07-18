from src.Network.ConvolutionalNetwork import ConvolutionalNetwork
from src.Preprocessing.PreprocessorConv import PreprocessorConv


class Convolutional:

    def __init__(self):

        # Match search result .match file
        self.match_file = '/ebio/abt1_share/update_tprpred/data/PDB_Approach/Results/query.match'

    @staticmethod
    def init_preprocessor():

        # Preprocessing Protein to predict
        preprocessor_object = PreprocessorConv()

        return preprocessor_object

    def init_training_data(self, preprocessor_object):

        preprocessor_object.filter_duplicates(self.match_file)

    @staticmethod
    def init_network():

        # Train network with Training and Test Set
        convolutional_network = ConvolutionalNetwork()

        return convolutional_network

    @staticmethod
    def train_network(network_object, training_data):

        network_object.compile_network()

        # network_object.train_network(training_data[0], training_data[1], training_data[2], training_data[3])