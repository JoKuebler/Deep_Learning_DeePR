from src.Data_Processing.file_reader import FileReader
from src.Data_Processing.encoding import Encoder
from src.Network.ConvolutionalNetwork import ConvolutionalNetwork

if __name__ == '__main__':

    # Init FileReader object
    file_read_train = FileReader('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/positive_data.txt',
                                 '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/negative_data.txt')

    # TODO BUILD VALIDATION SET OF
    # file_read_val = FileReader('')

    # Get Training Data as list
    pos_data, neg_data = file_read_train.read_training_data()

    # Init Encoder object
    encoder = Encoder()
    # Encode Training Data
    enc_data, target = encoder.encode(pos_data, neg_data)

    print(enc_data.shape)
    print(target.shape)

    # Init convolutinal network
    conv_net = ConvolutionalNetwork()
    # Train network
    conv_net.train_network(enc_data, target)



