from src.Network.PredictorMain import MainPredictor

if __name__ == '__main__':

    # Initalize Main Runner Object
    main_runner = MainPredictor()

    # Initialize necessary objects
    Preprocessor = main_runner.init_preprocessor(main_runner.test_three)
    Network = main_runner.init_network()

    # Initialize objects when training model
    # Training_data = init_training_data(Preprocessor)
    # train_network(Network, Training_data)

    # If model has to be saved after training
    # save_network(Network)

    # Load network from saved json
    main_runner.load_network(Network)

    # Predict desired data with threshold
    main_runner.predict_training_data(Network, Preprocessor, 0.0)