# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

try:
    import os
except(Exception):
    raise("[ERROR] os was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import sys
except(Exception):
    raise("[ERROR] sys was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import csv
except(Exception):
    raise("[ERROR] csv was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import numpy as np
except(Exception):
    raise("[ERROR] numpy was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import math
except(Exception):
    raise("[ERROR] math was not imported properly. Exiting.")
    sys.exit(-1)

try:
    import itertools
except(Exception):
    raise("[ERROR] itertools was not imported properly")
    sys.exit(-1)

try:
    from sklearn.neighbors import KernelDensity
except(Exception):
    raise("[ERROR] sklearn was not imported properly")

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
except(Exception):
    raise "[ERROR] plt from matplotlib.pyplot was" + \
          "not imported properly. Exiting."
    sys.exit(-1)

try:
    import pandas as pd
except(Exception):
    raise("[ERROR] pandas was not imported properly. Exiting.")
    sys.exit(-1)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# CLASSES
# -----------------------------------------------------------------------------


class BayesClassifier(object):
    """
    Class that represents a Bayes classifier.
        :rtype: BayesClassifier object
    *Class members*
        *background_likelihood*
            Dictionary. Holds feature names as keys and the
                        corresponding distributions as values.
    """

    def __init__(self):
        # the classes
        self.class_col = None
        self.classes = list()
        self.class_dict = dict()

    def estimate_kernel(self, x):
        """
        Method...
        """

        # ---------------------------------------------------------------------
        # get the kernel
        # ---------------------------------------------------------------------
        # Function from scikit:
        # http://scikit-learn.org/stable/modules/density.html
        # bandwith = h
        # ---------------------------------------------------------------------

        kernel = KernelDensity(
            kernel='exponential',
            bandwidth=np.std(x, ddof=1) / 5 + np.finfo(float).eps).fit(
                np.array(x)[:, np.newaxis]
        )

        return(kernel)

    def fit(
        self,
        X_train,
        y_train
    ):
        """
        Method...
        """

        # initialize the variables we need
        self.class_col = y_train.name
        self.classes = y_train.unique().tolist()

        # join the data together
        feature_data_frame = X_train.join(y_train)

        # go over all classes
        for current_class in self.classes:

            # add the current class to the class_dict
            if current_class not in self.class_dict:
                self.class_dict[current_class] = dict()
            else:
                sys.stderr.write(
                    "ERROR: class {} \
                     already exists in the \
                     BayesClassifier.class_dict. \
                     {}".format(current_class,
                                os.linesep)
                )
                sys.exit(-1)

            # add for each feature the kernel density
            for current_feature in X_train.columns.values:
                # get the values that we observe for the current
                # class in our training data
                feature_values = \
                    feature_data_frame[feature_data_frame[
                        self.class_col] == current_class][
                            current_feature].tolist()
                # estimate the kernel
                kernel = self.estimate_kernel(feature_values)

                # store the kernel in the dictionary
                if current_feature not in self.class_dict[current_class]:
                    self.class_dict[current_class][current_feature] = kernel
                else:
                    sys.stderr.write(
                        """ERROR: feature {}
                         already exists in the
                         BayesClassifier.class_dict.
                         {}""".format(current_feature,
                                      os.linesep)
                    )
                    sys.exit(-1)

    def predict(self, X):
        """
        Method that predicts the class of each feature data
        set (row) in a given dataframe X.
        """

        # the list that will contain the predicted classes for each row
        predicted_label = []

        # do the prediction row by row (=dataset per dataset)
        for index, row in X.iterrows():

            # initialize the variables for the current dataset
            curr_max_probability = 0.0
            curr_max_class = None

            for current_class in self.class_dict:

                # calculate the probability to belong to the current class
                current_probability = 1.0
                # for current_feature, kernel in self.class_dict[current_class].iteritems():
                for current_feature, kernel in list(self.class_dict[current_class].items()):

                    # get the value that we have observed for the
                    # current feature in the current dataset (=row)
                    x = np.array([[row[current_feature]]])

                    # multiply it with the current probability
                    # HINT: the kernel.score_samples function gives back
                    #       the ln of the probability to belong to the
                    #       current kernel (=distribution).
                    current_probability *= np.exp(kernel.score_samples(x))

                if current_probability >= curr_max_probability:
                    curr_max_probability = current_probability
                    curr_max_class = current_class

            predicted_label.append(curr_max_class)

        return np.array(predicted_label)

    def predict_proba(self, row):
        """
        Method that predicts the class of each feature data set
        (row) in a given dataframe X.
        """

        # initialize the variables for the current dataset
        curr_max_probability = 0.0
        curr_max_class = None

        list_of_probabilities = []

        class_order = []

        for current_class in self.class_dict:

            class_order.append(current_class)

            # calculate the probability to belong to the current class
            current_probability = 1.0
            # for current_feature, kernel in self.class_dict[current_class].iteritems():
            for current_feature, kernel in list(self.class_dict[current_class].items()):

                # get the value that we have observed for the current feature
                # in the current dataset (=row)
                x = np.array([[row[current_feature]]])

                # multiply it with the current probability
                # HINT: the kernel.score_samples function gives back
                #       the ln of the probability to belong to the
                #       current kernel (=distribution).
                current_probability *= np.exp(kernel.score_samples(x))

            list_of_probabilities.append(current_probability)

            if current_probability >= curr_max_probability:
                curr_max_probability = current_probability
                curr_max_class = current_class

        normalized_probabilities = \
            [p / sum(list_of_probabilities) for p in list_of_probabilities]

        results = dict()

        for i in range(len(class_order)):
            results[class_order[i] + "_probability"] = \
                normalized_probabilities[i][0]
        results['classification'] = curr_max_class

        return(pd.Series(results))
