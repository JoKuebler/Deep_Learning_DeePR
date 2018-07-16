from src.Network.ConvolutionalNetwork import ConvolutionalNetwork


class Convolutional:

    @staticmethod
    def init_network():

        # Train network with Training and Test Set
        convolutional_network = ConvolutionalNetwork()

        return convolutional_network

    @staticmethod
    def train_network(network_object, training_data):

        network_object.compile_network()

        # network_object.train_network(training_data[0], training_data[1], training_data[2], training_data[3])