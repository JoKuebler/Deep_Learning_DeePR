# TEXT ENCODED APPROACH
# Read in positive Samples
# with open(training_set) as raw_trps:
#     content = raw_trps.readlines()
#
# # Positive training samples
# train_samples = [x.strip() for x in content]
#
# # Read in negative samples
# with open(training_negative_set) as raw_neg:
#     content_neg = raw_neg.readlines()
#
# # Negative Training Samples
# train_samples_neg = [x.strip() for x in content_neg]
#
# # create labels
# pos_labels = [1] * 193
# neg_labels = [0] * 30
# labels = pos_labels + neg_labels
#
# # Combine positive and negative samples
# final_training_set = train_samples + train_samples_neg
#
# print(final_training_set)
#
# text_encoded_tprs = []
# # Preprocess Training Sample
# for sample in final_training_set:
#     text_encoded = PreprocessorObject.text_encode_fragment(sample)
#     text_encoded_tprs.append(text_encoded)
#
# # Convert to numpy Array for training
# np_samples = np.array(text_encoded_tprs)
# np_labels = np.array(labels)
#
# print(np_samples[0])
# print(np_labels)
#
# # create Neural Network Object
# NeuralNetworkObject = NeuralNetwork(np_samples, np_labels)
#
# # Train neural network
# NeuralNetworkObject.train_network()




# BIT Encoded Approach
# Read in positive Samples
# with open(training_set) as raw_trps:
#     content = raw_trps.readlines()
#
# # Read in negative Samples
# with open(training_negative_set) as neg_set:
#     content_neg = neg_set.readlines()
#
# encoded_fragments = []
#
# for fragment in content:
#     encoded = PreprocessorObject.encodeFragment(fragment.rstrip())
#     encoded_fragments.append(encoded)
#

# data_file = open('encoded.csv', 'w')
#
# for frag in encoded_fragments:
#     data_file.write(frag + ',1' + '\n')
#
# encoded_fragments = []
#
# for fragment in content_neg:
#     encoded = PreprocessorObject.encodeFragment(fragment.rstrip())
#     encoded_fragments.append(encoded)
#
# for frag in encoded_fragments:
#     data_file.write(frag + ',0' + '\n')
#
# data_file.close()
#
# dataset = numpy.loadtxt("encoded.csv", delimiter=",")
# training_data = dataset[:, 0:34]
# labels = dataset[:, 34]
#
# reshaped_labes = labels.reshape(223, 1)
#
# labels1 = np.random.randint(2, size=(223, 1))
# data = np.random.random((223, 34))
#
# print(training_data[0])
# print(data.shape)
# print(data[0])
# print(training_data.shape)
# print(training_data.reshape(223,34,1))

#
# NeuralNetworkObject = NeuralNetwork(training_data.reshape(223,34,1), reshaped_labes)
# #
# NeuralNetworkObject.train_network()





# PREDICTION
# Preprocessing Test Protein
# PreprocessorObject = Preprocessor(inputProtein)
#
# # Cut Input Protein into windows
# fragments = PreprocessorObject.cutInWindows(34)
# text_encoded_test_samples = []
#
#
# # To predict for probability for each fragment each one has to be encoded
# encoded_fragments = []
#
# for fragment in fragments:
#     encoded = PreprocessorObject.encodeFragment(fragment.rstrip())
#     encoded_fragments.append(encoded)
#
# data_file = open('pred_encoded.csv', 'w')
#
# for frag in encoded_fragments:
#     data_file.write(frag  + '\n')
#
# data_file.close()
#
# dataset = numpy.loadtxt("pred_encoded.csv", delimiter=",")
#
# print(len(dataset[0]))
# print(dataset[0])

# # Results contain probabilites for each class
#results = NeuralNetworkObject.predict(dataset[0])