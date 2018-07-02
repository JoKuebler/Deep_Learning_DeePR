from src.Network.PredictorMain import MainPredictor

if __name__ == '__main__':

    # Initalize Main Runner Object
    main_runner = MainPredictor()

    # Initialize necessary objects
    Preprocessor = main_runner.init_preprocessor(main_runner.test_two)
    Network = main_runner.init_network()

    # Initialize objects when training model
    Training_data = main_runner.init_training_data(Preprocessor)
    main_runner.train_network(Network, Training_data)

    # If model has to be saved after training
    main_runner.save_network(Network)

    # Load network from saved json
    # main_runner.load_network(Network)

    # Predict desired data with threshold
    # main_runner.predict_training_data(Network, Preprocessor, 0.8)

    main_runner.single_predict(Network, Preprocessor, 0.93)