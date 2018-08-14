# Created by Jonas Kuebler on 07/02/2018
# Class to choose featues
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from src.Network.FeedForwardMain import MainPredictor
from sklearn import datasets
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE


class FeatureSelection:
    """Class to apply feature selection on the sets of features

       :param
     """

    def __init__(self):
        self.model = RFECV(RandomForestClassifier(), scoring='accuracy')

    def fit(self, training_data, training_labels):

        dataset_size = len(training_data)
        print(dataset_size)
        TwoDim_dataset = training_data.reshape(dataset_size, -1)

        # # create a base classifier used to evaluate a subset of attributes
        # model = LogisticRegression()
        # # # create the RFE model and select 3 attributes
        # rfe = RFE(model, 3)
        # rfe = rfe.fit(TwoDim_dataset, np.ravel(training_labels))
        # # # summarize the selection of the attributes
        # print(rfe.support_)
        # print(rfe.ranking_)

        # load the iris datasets
        dataset = datasets.load_iris()
        # create a base classifier used to evaluate a subset of attributes
        model = LogisticRegression()
        # create the RFE model and select 3 attributes
        print(dataset.data)
        rfe = RFE(model, 3)
        rfe = rfe.fit(dataset.data, dataset.target)
        # summarize the selection of the attributes
        print(rfe.support_)
        print(rfe.ranking_)



