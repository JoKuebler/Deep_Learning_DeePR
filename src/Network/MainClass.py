from src.Network.FeedForwardMain import FeedForward
from src.Network.ConvolutionalMain import Convolutional

if __name__ == '__main__':

    # ------------ Initializing Objects -------------
    # Initalize Main Runner Object
    # main_runner_ff = FeedForward()
    main_runner_conv = Convolutional()

    # Preprocessor with input data to predict
    # Preprocessor_ff = main_runner_ff.init_preprocessor(main_runner_ff.test_four)
    Preprocessor_conv = main_runner_conv.init_preprocessor()

    # Feed forward Network
    # ffNet = main_runner_ff.init_network()
    # Convolutinal Network
    convNet = main_runner_conv.init_network()

    # Feature Selector
    # FeatureSelector = FeatureSelection()
    # Result Evaluator
    # Result_Evaluator = Result_Evaluator()
    # PDB Parser
    # PDBParser = PDB_Parser()
    # -----------------------------------------------

    # -------------- Preprocessing ------------------
    # PDBParser.single_chains('/ebio/abt1_share/update_tprpred/data/TestPDB')
    # PDBParser.download_build('/ebio/abt1_share/update_tprpred/data/PDB_Approach/PDS/')

    # --------------- Training ----------------------
    # Initialize objects when training model
    # Training_data = main_runner_ff.init_training_data(Preprocessor_ff)
    Training_data_conv = main_runner_conv.init_training_data(Preprocessor_conv)
    # Traing network on data and apply cross validation
    # main_runner_ff.train_network(ffNet, Training_data)
    # main_runner_ff.cross_validate(Network, Training_data)
    # -----------------------------------------------

    # --------------- Safe & Load -------------------
    # If model has to be saved after training
    # main_runner_ff.save_network(ffNet)

    # Load network from saved json
    # main_runner_ff.load_network(ffNet)
    # -----------------------------------------------

    # --------------- Prediction --------------------
    # Predict raw training data (fragments)
    # main_runner_ff.predict_training_data(Network, Preprocessor_two, 0.8)

    # predict desired data with threshold (sequences)
    # main_runner_ff.single_predict(ffNet, Preprocessor_ff, 0.5)
    # -----------------------------------------------

    # --------------Evaluation ----------------------
    # Select Features
    # FeatureSelector.fit(Training_data[0], Training_data[1])

    # Predict Test Data and calculate Sensitivity and Specifity
    # results = Network.predict(Training_data[2])
    # print(Result_Evaluator.calculate_metrices(results, 7757))
    # ------------------------------------------------