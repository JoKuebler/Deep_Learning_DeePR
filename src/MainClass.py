from src.Data_Processing.file_reader import FileReader
from src.Data_Processing.encoding import Encoder

if __name__ == '__main__':

    # Create FileReader object
    file_reader = FileReader('/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/positive_data.txt',
                             '/ebio/abt1_share/update_tprpred/data/Convolutional/TrainingData/negative_data.txt')
    # Get Training Data as list
    pos_data, neg_data = file_reader.read_training_data()

    # Create Encoder object
    encoder = Encoder()
    # Encode Training Data
    enc_data, target = encoder.encode(pos_data, neg_data)




